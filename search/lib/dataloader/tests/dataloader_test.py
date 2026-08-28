#dataloader_tests.py -- tests for the dataloader module.

from .. import trainloader

def test_fsmStateToOneHot():
    for state in "0123456789ABCDE":
        oneHot = trainloader.fsmStateToOneHot(state)
        assert oneHot is not None, f"fsmStateToOneHot returned None for valid state {state}"
        assert oneHot.sum() == 1.0, f"fsmStateToOneHot returned a vector that does not sum to 1 for state {state}"
        assert oneHot["0123456789ABCDE".index(state)] == 1.0, f"fsmStateToOneHot returned a vector that does not have 1 in the correct position for state {state}"
    invalid = trainloader.fsmStateToOneHot("Z")
    assert invalid is None, "fsmStateToOneHot did not return None for invalid state 'Z'"
    print("All tests passed for fsmStateToOneHot")

def test_nucleotideToOneHot():
    for nucleotide in "ACGT":
        oneHot = trainloader.nucleotideToOneHot(nucleotide)
        assert oneHot is not None, f"nucleotideToOneHot returned None for valid nucleotide {nucleotide}"
        assert oneHot.sum() == 1.0, f"nucleotideToOneHot returned a vector that does not sum to 1 for nucleotide {nucleotide}"
        assert oneHot["ACGT".index(nucleotide)] == 1.0, f"nucleotideToOneHot returned a vector that does not have 1 in the correct position for nucleotide {nucleotide}"
    invalid = trainloader.nucleotideToOneHot("Z")
    assert invalid is None, "nucleotideToOneHot did not return None for invalid nucleotide 'Z'"
    print("All tests passed for nucleotideToOneHot")