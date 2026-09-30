# trainloader.py -- functions for reading training data from disk and creating minibatches for training

import torch
from torch.utils.data import Dataset
from . import compress

_FSM_STATE_TO_ONE_HOT_LOOKUP = torch.tensor([[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                                            [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
                                            [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
                                            [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
                                            [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
                                            [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
                                            [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
                                            [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
                                            [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
                                            [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
                                            [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
                                            [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
                                            [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
                                            [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
                                            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
                                            ], dtype=torch.float)

def fsm_state_to_one_hot(fsmState):
    """
    Convert a character representing a state in the FSM to a one-hot training vector
    
    Expects fsmState to be a single character string, in the set "0-9,A-E". 
    Returns a one-hot tensor of float32 values, with 1.0 in the position corresponding to the state, and 0.0 elsewhere
    Returns "None" if the input is not a valid state character
 
    """
    return _FSM_STATE_TO_ONE_HOT_LOOKUP[fsmState]

_NUCLEOTIDE_TO_ONE_HOT_LOOKUP = torch.tensor([[1.0, 0.0, 0.0, 0.0], #A
                                             [0.0, 1.0, 0.0, 0.0], #C
                                             [0.0, 0.0, 1.0, 0.0], #G
                                             [0.0, 0.0, 0.0, 1.0], #T
                                             [0.0, 0.333, 0.333, 0.333], #B
                                             [0.333, 0.0, 0.333, 0.333], #D
                                             [0.333, 0.333, 0.0, 0.333], #H
                                             [0.0, 0.0, 0.5, 0.5], #K
                                             [0.5, 0.5, 0.0, 0.0], #m
                                             [0.25, 0.25, 0.25, 0.25], #N/X        
                                             [0.5, 0.0, 0.5, 0.0], #R
                                             [0.0, 0.5, 0.5, 0.0], #S
                                             [0.333, 0.333, 0.333, 0.0], #V
                                             [0.5, 0.0, 0.0, 0.5], #W
                                             [0.0, 0.5, 0.0, 0.5], #Y
                                             [0.0, 0.0, 0.0, 0.0] #Z
                                            ], dtype = torch.float)
def nucleotide_to_one_hot(nucleotide):
    """
    Convert a character representing a nucleotide to a one-hot training vector
    
    Expects nucleotide to be a single character string, in the set "ACGT". 
    Returns a one-hot tensor of float32 values, with 1.0 in the position corresponding to the nucleotide, and 0.0 elsewhere
    Returns "None" if the input is not a valid nucleotide character
    """
    return _NUCLEOTIDE_TO_ONE_HOT_LOOKUP[nucleotide]

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

        sequence_onehot = torch.vmap(nucleotide_to_one_hot, in_dims = 0)(sequence)
        fsm_states_onehot = torch.vmap(fsm_state_to_one_hot, in_dims = 0)(fsm_states)

        return(name, sequence_onehot, fsm_states_onehot)

        
    
    
