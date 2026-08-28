# trainloader.py -- functions for reading training data from disk and creating minibatches for training

import torch

def fsmStateToOneHot(fsmState):
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


def nucleotideToOneHot(nucleotide):
    """
    Convert a character representing a nucleotide to a one-hot training vector
    
    Expects nucleotide to be a single character string, in the set "ACGT". 
    Returns a one-hot tensor of float32 values, with 1.0 in the position corresponding to the nucleotide, and 0.0 elsewhere
    Returns "None" if the input is not a valid nucleotide character
    """
    if nucleotide not in "ACGT":
        return None
    oneHot = [0.0] * 4
    index = "ACGT".index(nucleotide)
    oneHot[index] = 1.0
    return torch.tensor(oneHot, dtype=torch.float32)

