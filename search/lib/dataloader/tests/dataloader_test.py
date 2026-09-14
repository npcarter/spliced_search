#dataloader_tests.py -- tests for the dataloader module.

from .. import trainloader

def test_fsm_state_to_one_hot():
    for state in "0123456789ABCDE":
        oneHot = trainloader.fsm_state_to_one_hot(state)
        assert oneHot is not None, f"fsm_state_to_one_hot returned None for valid state {state}"
        assert oneHot.sum() == 1.0, f"fsm_state_to_one_hot returned a vector that does not sum to 1 for state {state}"
        assert oneHot["0123456789ABCDE".index(state)] == 1.0, f"fsm_state_to_one_hot returned a vector that does not have 1 in the correct position for state {state}"
    invalid = trainloader.fsm_state_to_one_hot("Z")
    assert invalid is None, "fsm_state_to_one_hot did not return None for invalid state 'Z'"
    print("All tests passed for fsm_state_to_one_hot")

def test_nucleotide_to_one_hot():
    for nucleotide in "ACGT":
        oneHot = trainloader.nucleotide_to_one_hot(nucleotide)
        assert oneHot is not None, f"nucleotide_to_one_hot returned None for valid nucleotide {nucleotide}"
        assert oneHot.sum() == 1.0, f"nucleotide_to_one_hot returned a vector that does not sum to 1 for nucleotide {nucleotide}"
        assert oneHot["ACGT".index(nucleotide)] == 1.0, f"nucleotide_to_one_hot returned a vector that does not have 1 in the correct position for nucleotide {nucleotide}"
    invalid = trainloader.nucleotide_to_one_hot("Z")
    assert invalid is None, "nucleotide_to_one_hot did not return None for invalid nucleotide 'Z'"
    print("All tests passed for nucleotide_to_one_hot")