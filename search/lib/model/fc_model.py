import torch
import torch.nn as nn
"""
Contains a simple fully-connected baseline model for gene finding
"""

class TwoLayerFC(nn.Module):
    """
    Fully-connected network with one hidden layer (two weight layers) that maps an n-element sequence of one-hot
    nucleotides to an n-element sequence of per-position state logits.
    Every output position sees every input position, so the model is tied to a fixed sequence length.
    :param sequence_length: n, the number of elements in each input/output sequence, defaults to 128
    :param hidden_size: width of the hidden layer, defaults to 256
    :param num_inputs: size of the one-hot vector for each input element, defaults to 4 (A, C, G, T)
    :param num_outputs: number of output classes (FSM states) for each sequence element, defaults to 15
    :param dropout: dropout probability applied after the hidden activation, defaults to 0.0
    """
    def __init__(self, sequence_length=128, hidden_size=256, num_inputs=4, num_outputs=15, dropout=0.0):
        super().__init__()
        self.sequence_length = sequence_length
        self.num_inputs = num_inputs
        self.num_outputs = num_outputs
        self.net = nn.Sequential(
            nn.Flatten(start_dim=1),                                   # (batch, n, 4) -> (batch, n*4)
            nn.Linear(sequence_length * num_inputs, hidden_size),      # level 1
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, sequence_length * num_outputs),     # level 2
        )

    def forward(self, x):
        """
        Computes per-position logits for a batch of sequences.
        :param x: Tensor of shape (batch, sequence_length, num_inputs), or (sequence_length, num_inputs) for a single sample
        :return: Raw (un-normalized with Softmax) logits of shape (batch, sequence_length, num_outputs), or
            (sequence_length, num_outputs) if a single unbatched sample was passed in
        """
        unbatched = x.dim() == 2
        if unbatched:
            x = x.unsqueeze(0)
        if x.shape[1:] != (self.sequence_length, self.num_inputs):
            raise ValueError(f"Expected input of shape (batch, {self.sequence_length}, {self.num_inputs}), got {tuple(x.shape)}")
        logits = self.net(x).view(x.shape[0], self.sequence_length, self.num_outputs)
        return logits.squeeze(0) if unbatched else logits
