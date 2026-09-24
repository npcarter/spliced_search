import torch
import torch.nn as nn
"""
Contains loss functions for gene finding
"""

EPSILON = 1e-10  #used to prevent division by zero

# Maps each of the 15 FSM states to its reduced 5-state equivalent (index = 15-state id, value = 5-state id)
_REDUCED_STATE_LOOKUP = torch.tensor(
    [
        0,  # 0: intergenic
        2,  # 1: start -> Exon, frame 0
        4,  # 2: stop -> Exon, frame 2
        2,  # 3: Exon, frame 0
        3,  # 4: Exon, frame 1
        4,  # 5: Exon, frame 2
        1,  # 6: Intron, frame irrelevant
        1,  # 7: Intron, frame irrelevant
        1,  # 8: Intron, frame irrelevant
        2,  # 9: ASS0 -> Exon, frame 0
        3,  # 10: ASS1 -> Exon, frame 1
        4,  # 11: ASS2 -> Exon, frame 2
        2,  # 12: DSS0 -> Exon, frame 0
        3,  # 13: DSS1 -> Exon, frame 1
        4,  # 14: DSS2 -> Exon, frame 2
    ],
    dtype=torch.long,
)

class TiberiusLoss(nn.Module):
    """
    Reproduces the custom loss function used by Tiberius, as described in their work.  
    Loss = categorical cross entropy loss + f1_factor * (F1 score of states that represent codon positions + factor to penalize false positives in sequences with no codons)
    :param use_f1: if True, compute the F1 loss of the states that correspond to codon positions and include that in the loss, defaults to True
    :param f1_factor: Weight of the F1 loss and false-positive penalty relative to the CCE loss, defaults to 2.0
    """
    def __init__(self, use_f1=True, f1_factor=2.0):
        """
        Saves the provided values of use_f1 and f1_factor for later access
        """
        self.use_f1 = use_f1
        self.f1_factor = f1_factor

    @staticmethod
    def _reduce_states(fsm_states):
        """
        Converts a 2-D tensor of shape sequence_length x 15 into a one-D tensor of length sequence_length, where each element in the output 
        tensor contains the integer representation of the 5-state reduced FSM state corresponding to the state described by the corresponding 
        column of the input tensor.
        :param fsm_states: The 2-D tensor of one-hot FSM states
        """ 
        
        def _reduce_state(state_onehot):
            """
            Converts a single one-hot tensor into an integer representing the corresponding reduced state.
            Uses tensor indexing rather than Python-level branching (match/case), which is required
            for compatibility with torch.vmap since a batched tensor has no single concrete value to
            branch on.
            param: state_onehot: The input tensor. Must be of length 15
            """
            which_state = torch.argmax(state_onehot, dim=-1)
            return _REDUCED_STATE_LOOKUP.to(which_state.device)[which_state]
        return _reduce_state(fsm_states) # Just map the per_state reduction function across the states

    @staticmethod
    def _reduce_logits(logits_fifteen):
        """
        Converts a set of logits in 15-state format into 5-state reduced format.
        :param logits_fifteen: Tensor of shape (..., 15)
        """
        state_0 = logits_fifteen[..., 0]
        state_1 = torch.logsumexp(torch.stack([logits_fifteen[..., 6],logits_fifteen[..., 7], logits_fifteen[..., 8]]), dim=0)
        state_2 = torch.logsumexp(torch.stack([logits_fifteen[..., 1], logits_fifteen[..., 3],  logits_fifteen[..., 9], logits_fifteen[..., 12]]), dim = 0)
        state_3 = torch.logsumexp(torch.stack([logits_fifteen[..., 4], logits_fifteen[..., 10], logits_fifteen[..., 13]]), dim=0)
        state_4 = torch.logsumexp(torch.stack([logits_fifteen[..., 2], logits_fifteen[..., 5], logits_fifteen[..., 11], logits_fifteen[..., 14]]), dim=0)
        return torch.stack([state_0, state_1, state_2, state_3, state_4], dim=-1)

    @staticmethod
    def _predict_codon(logit_five):
        #return one if the sum of the logits for coding states (2, 3, 4) is greater than the sum of the logits for non-coding states
        # (0, 1)
        coding = torch.logsumexp(logit_five[..., 2:5], dim= -1)
        non_coding = torch.logsumexp(logit_five[..., 0:2], dim = -1)
        return(coding > non_coding)
    
    def forward(self, predictions_fifteenstate, targets_onehot):
        """
        Computes the loss for a set of predictions relative to the target categorization.  
        :param predictions_fifteenstate: The raw (un-normalized with Softmax) logits of the predicted probabilities that a given position is in each   
            FSM state.  Must be a 2-D tensor, where the first dimension has length equal to the number of nucleotides in the sequence and the second dimension is 15, and contains the logits for the probability that the corresponding nucleotide is in each of the model's FSM states
        :param targets_onehot: The FSM state that each nucleotide should be in if predicted properly. Should be a 2-D tensor of dimensions 
            sequence_length x 15, where each column is a one-hot vector with the position corresponding to the correct FSM state set to 1.0
        """
        print("starting forward")
        print(f"targets_onehot = {targets_onehot}")
        print(f"predictions_fifteenstate = {predictions_fifteenstate}")
        #first, compute the categorical cross-entropy of the predictions
        targets = self._reduce_states(targets_onehot)
        targets_flat = targets.reshape(-1)
        predictions = self._reduce_logits(predictions_fifteenstate)
        predictions_flat = predictions.reshape(-1, predictions.shape[-1])
        cce_loss = nn.functional.cross_entropy(predictions_flat, targets_flat, reduction='mean')
        print(f"targets_flat = {targets_flat}")
        print(f"predictions_flat = {predictions_flat}")
        print(f"cce_loss = {cce_loss}")
        if self.use_f1:
            # Find the number of sequence positions that held codons
            true_codons = torch.where(targets > 1, 1, 0) # states 2-4 in the five-state rep are coding 
            num_ground_truth_codons = torch.sum(true_codons)

            if num_ground_truth_codons >0: # sequence had at least one codon, compute F1


                # Find the sequence positions that were predicted to have codons
                codons_pred = torch.vmap(self._predict_codon, in_dims=0)(predictions)
                true_positives = true_codons * codons_pred  # will be 1 where a codon was both predicted and in ground truth


                num_true_positives = torch.sum(true_positives)
                num_predicted_positives = torch.sum(codons_pred)

                assert (num_predicted_positives >= num_true_positives), "Somehow got more true positives than predicted"
                assert(num_ground_truth_codons >= num_true_positives), "Somehow got more true positives than positives in ground truth"

                #extract the batch size from the input dimensions.  Handles the case where the last batch has different size than the rest
                batch_size = predictions_fifteenstate.shape[0] if predictions_fifteenstate.dim()==3 else 1
                precision = num_true_positives/(num_predicted_positives + EPSILON)
                recall = num_true_positives/(num_ground_truth_codons + EPSILON)
                f1_score = (2 * precision * recall)/(precision+ recall+ EPSILON)
                f1_loss = (1- f1_score) /batch_size # Subtract from one because a high f1 score is good, and we want a high loss 
                # when the model is doing a bad job.

                return cce_loss + (self.f1_factor * f1_loss) 
            else: #penalize for false positives.
                total_positions = targets_onehot.numel() // targets_onehot.shape[-1] #total number of sequence positions in batch
                 # Find the sequence positions that were predicted to have
                 # codons 
                codons_pred = torch.vmap(self._predict_codon, in_dims=0)(predictions)
                fpr = torch.sum(codons_pred)/total_positions
                return cce_loss + (self.f1_factor * fpr) #In Tiberius code, f1_loss is 0 if no ground truth positives, so don't need it in this code

        else: # not using f1, just return cce
            return cce_loss

        
  
  
  
  
  
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  

  

  
  
  
  
  
  
  
  
  
  

  

  

  
  
  
  
  
  
  

  
  
  
  
  
  
  

















