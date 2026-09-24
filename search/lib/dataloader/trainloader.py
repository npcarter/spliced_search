# trainloader.py -- functions for reading training data from disk and creating minibatches for training

import torch
from torch.utils.data import Dataset, DataLoader
from . import compress

def fsm_state_to_one_hot(fsmState):
    """
    Convert a character representing a state in the FSM to a one-hot training vector
    
    Expects fsmState to be a single character string, in the set "0-9,A-E". 
    Returns a one-hot tensor of float32 values, with 1.0 in the position corresponding to the state, and 0.0 elsewhere
    Returns "None" if the input is not a valid state character
 
    """
    if fsmState not in "0123456789ABCDE":
        return None
    oneHot = [0.0] * 15
    index = "0123456789ABCDE".index(fsmState)
    oneHot[index] = 1.0
    return torch.tensor(oneHot, dtype=torch.float32)


def nucleotide_to_one_hot(nucleotide):
    """
    Convert a character representing a nucleotide to a one-hot training vector
    
    Expects nucleotide to be a single character string, in the set "ACGT". 
    Returns a one-hot tensor of float32 values, with 1.0 in the position corresponding to the nucleotide, and 0.0 elsewhere
    Returns "None" if the input is not a valid nucleotide character
    """
    oneHot = [0.0] * 4
    match nucleotide:
        case 'A':
            oneHot[0] = 1.0
        case 'C':
            oneHot[1] = 1.0   
        case 'T':
            oneHot[2] = 1.0
        case 'G':
            oneHot[3] = 1.0
        case 'B':
            oneHot[1] = 0.333
            oneHot[2] = 0.333
            oneHot[3] = 0.333
        case 'D':
            oneHot[0] = 0.333
            oneHot[2] = 0.333
            oneHot[3] = 0.333
        case 'H':
            oneHot[0] = 0.333
            oneHot[1] = 0.333
            oneHot[2] = 0.333
        case 'V':
            oneHot[0] = 0.333
            oneHot[1] = 0.333
            oneHot[3] = 0.333
        case 'K':
            oneHot[2] = 0.5
            oneHot[3] = 0.5          
        case 'M':
            oneHot[0] = 0.5
            oneHot[1] = 0.5         
        case 'R':
            oneHot[0] = 0.5
            oneHot[3] = 0.5         
        case 'S':
            oneHot[1] = 0.5
            oneHot[3] = 0.5         
        case 'W':
            oneHot[0] = 0.5
            oneHot[2] = 0.5         
        case 'Y':
            oneHot[1] = 0.5
            oneHot[2] = 0.5
        case 'N' | 'X':
            oneHot[0] = 0.25
            oneHot[1] = 0.25
            oneHot[2] = 0.25
            oneHot[3] = 0.25
        case _:
            return None
    
    return torch.tensor(oneHot, dtype=torch.float32)

class TrainDataset(Dataset):
    """
    Dataset class for training data. Samples are stored on disk in a Zarr object, in compressed form (one byte/nucleotide).
    __getitem__() reads a sample from disk, converts it to one-hot tensor sample and label tensors, and returns it.
    """
    def __init__(self, zarr_path, sequence_length=16384, name_length=40):
        """
        Initialize the dataset with the path to the Zarr object containing the training data.
        """
        import zarr
        self.zarr_data = zarr.open(zarr_path, mode='r')
        (self.num_samples, _) = self.zarr_data.shape
        self.sequence_length = sequence_length
        self.name_length = name_length

    def __len__(self):
        """
        Return the number of samples in the dataset.
        """
        return self.num_samples

    def __getitem__(self, idx):
        """
        Read a sample from disk, convert it to one-hot tensor sample and label tensors, and return it.
        
        Args:
            idx (int): Index of the sample to retrieve.
        """

        # Read the sample from the Zarr object

        if (idx < 0 ) or (idx > self.num_samples):
            raise ValueError(f"Attempted to fetch training sample {idx}. Sample indices must be between 0 and {self.num_samples}.")
        compressed_item = self.zarr_data[idx]
        (name, sequence, fsm_states) = compress.decompress_chunk(compressed_item, self.sequence_length, self.name_length)
        sequence_onehot = []
        for nucleotide in sequence:
            sequence_onehot.append(nucleotide_to_one_hot(nucleotide))

        fsm_states_onehot = []
        for state in fsm_states:
            fsm_states_onehot.append(fsm_state_to_one_hot(state))

        return(name, sequence_onehot, fsm_states_onehot)

        
    
    
