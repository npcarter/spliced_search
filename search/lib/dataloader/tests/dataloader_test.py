#dataloader_tests.py -- tests for the dataloader module.

from .. import trainloader
from .. import compress
import torch 

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

    # Test the single-nucleotide codes
    for nucleotide in "ACTG":
        oneHot = trainloader.nucleotide_to_one_hot(nucleotide)
        assert oneHot is not None, f"nucleotide_to_one_hot returned None for valid nucleotide {nucleotide}"
        assert oneHot.sum() == 1.0, f"nucleotide_to_one_hot returned a vector that does not sum to 1 for nucleotide {nucleotide}"
        assert oneHot["ACTG".index(nucleotide)] == 1.0, f"nucleotide_to_one_hot returned a vector that does not have 1 in the correct position for nucleotide {nucleotide}"
    invalid = trainloader.nucleotide_to_one_hot("Z")
    assert invalid is None, "nucleotide_to_one_hot did not return None for invalid nucleotide 'Z'"

    # Now, the multi-nucleotide ones
    oneHot = trainloader.nucleotide_to_one_hot('B')
    assert torch.equal(oneHot, torch.tensor([0.0, 0.333, 0.333, 0.333], dtype=torch.float32)), f"nucleotide_to_one_hot returned incorrect values for 'B'"
    oneHot = trainloader.nucleotide_to_one_hot('D')
    assert torch.equal(oneHot, torch.tensor([0.333, 0.0, 0.333, 0.333], dtype=torch.float32)), f"nucleotide_to_one_hot returned incorrect values for 'D'"
    oneHot = trainloader.nucleotide_to_one_hot('H')
    assert torch.equal(oneHot, torch.tensor([0.333, 0.333, 0.333, 0.0], dtype=torch.float32)), f"nucleotide_to_one_hot returned incorrect values for 'H'"
    oneHot = trainloader.nucleotide_to_one_hot('V')
    assert torch.equal(oneHot, torch.tensor([0.333, 0.333, 0.0, 0.333], dtype=torch.float32)), f"nucleotide_to_one_hot returned incorrect values for 'V'"
    oneHot = trainloader.nucleotide_to_one_hot('N')
    assert torch.equal(oneHot, torch.tensor([0.25, 0.25, 0.25, 0.25], dtype=torch.float32)), f"nucleotide_to_one_hot returned incorrect values for 'N'"
    oneHot = trainloader.nucleotide_to_one_hot('X')
    assert torch.equal(oneHot, torch.tensor([0.25, 0.25, 0.25, 0.25], dtype=torch.float32)), f"nucleotide_to_one_hot returned incorrect values for 'X'"
    oneHot = trainloader.nucleotide_to_one_hot('K')
    assert torch.equal(oneHot, torch.tensor([0.0, 0.0, 0.5, 0.5], dtype=torch.float32)), f"nucleotide_to_one_hot returned incorrect values for 'K'"
    oneHot = trainloader.nucleotide_to_one_hot('M')
    assert torch.equal(oneHot, torch.tensor([0.5, 0.5, 0.0, 0.0], dtype=torch.float32)), f"nucleotide_to_one_hot returned incorrect values for 'M'"
    oneHot = trainloader.nucleotide_to_one_hot('R')
    assert torch.equal(oneHot, torch.tensor([0.5, 0.0, 0.0, 0.5], dtype=torch.float32)), f"nucleotide_to_one_hot returned incorrect values for 'R'"
    oneHot = trainloader.nucleotide_to_one_hot('S')
    assert torch.equal(oneHot, torch.tensor([0.0, 0.5, 0.0, 0.5], dtype=torch.float32)), f"nucleotide_to_one_hot returned incorrect values for 'S'"
    oneHot = trainloader.nucleotide_to_one_hot('W')
    assert torch.equal(oneHot, torch.tensor([0.5, 0.0, 0.5, 0.0], dtype=torch.float32)), f"nucleotide_to_one_hot returned incorrect values for 'W'"
    oneHot = trainloader.nucleotide_to_one_hot('Y')
    assert torch.equal(oneHot, torch.tensor([0.0, 0.5, 0.5, 0.0], dtype=torch.float32)), f"nucleotide_to_one_hot returned incorrect values for 'Y'"

    print("All tests passed for nucleotide_to_one_hot")


def test_TrainDataset():
    import tempfile
    import os
    import zarr
    import shutil


    # Build a small file of compressed training data
    tempdir = tempfile.TemporaryDirectory(delete=False).name

    input_file_path = os.path.join(tempdir, "input.txt")
     # Create a temporary input file
    with open(input_file_path, "wb") as temp_input:
        temp_input.write(b"Seq1\nACGTACGT\n01234567\n")
        temp_input.write(b"Seq2\nTGCATGCA\n76543210\n")
        temp_input.write(b"Seq3\nACGTACGT\n01234567\n")
        temp_input.write(b"Seq4\nTGCATGCA\n76543210\n")
        temp_input.write(b"Seq5\nACGTACGT\n01234567\n")
        temp_input.write(b"Seq6\nTGCATGCA\n76543210\n")


    # Create a temporary output file path
    output_file_path = os.path.join(tempdir, "output.zarr")

    try:
        # Compress the input file
        compress.compress_file(input_file_path, output_file_path, sequence_chunk_size=8, name_length=10, file_chunk_size=4)

        # Read the compressed data from the output file using zarr
        compressed_data = zarr.open(output_file_path, mode='r')

        test_dataset = trainloader.TrainDataset(output_file_path, sequence_length=8, name_length=10)

        # Read elements from the dataset in mixed order to test random access.
        (name, sequence, states) = test_dataset[3]
        assert name == "Seq4_0    ", "Got wrong name {name} for item 3 in test_TrainDataset"

        assert torch.equal(sequence[0], torch.tensor([0.0, 0.0, 1.0, 0.0])), "Wrong value returned for position 0 of item 3 sequence in test_TrainDataset"
        assert torch.equal(sequence[1], torch.tensor([0.0, 0.0, 0.0, 1.0])), "Wrong value returned for position 1 of item 3 sequence in test_TrainDataset"
        assert torch.equal(sequence[2], torch.tensor([0.0, 1.0, 0.0, 0.0])), "Wrong value returned for position 2 of item 3 sequence in test_TrainDataset"
        assert torch.equal(sequence[3], torch.tensor([1.0, 0.0, 0.0, 0.0])), "Wrong value returned for position 3 of item 3 sequence in test_TrainDataset"
        assert torch.equal(sequence[4], torch.tensor([0.0, 0.0, 1.0, 0.0])), "Wrong value returned for position 4 of item 3 sequence in test_TrainDataset"
        assert torch.equal(sequence[5], torch.tensor([0.0, 0.0, 0.0, 1.0])), "Wrong value returned for position 5 of item 3 sequence in test_TrainDataset"
        assert torch.equal(sequence[6], torch.tensor([0.0, 1.0, 0.0, 0.0])), "Wrong value returned for position 6 of item 3 sequence in test_TrainDataset"
        assert torch.equal(sequence[7], torch.tensor([1.0, 0.0, 0.0, 0.0])), "Wrong value returned for position 7 of item 3 sequence in test_TrainDataset"

        assert torch.equal(states[0], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 0 of item 3 in test_TrainDataset"
        assert torch.equal(states[1], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 1 of item 3 in test_TrainDataset"
        assert torch.equal(states[2], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 2 of item 3 in test_TrainDataset"
        assert torch.equal(states[3], torch.tensor([0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 3 of item 3 in test_TrainDataset"
        assert torch.equal(states[4], torch.tensor([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 4 of item 3 in test_TrainDataset"
        assert torch.equal(states[5], torch.tensor([0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 5 of item 3 in test_TrainDataset"
        assert torch.equal(states[6], torch.tensor([0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 6 of item 3 in test_TrainDataset"
        assert torch.equal(states[7], torch.tensor([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 7 of item 3 in test_TrainDataset"

        (name, sequence, states) = test_dataset[4]
        assert name == "Seq5_0    ", "Got wrong name {name} for item 4 in test_TrainDataset"
        assert torch.equal(sequence[7], torch.tensor([0.0, 0.0, 1.0, 0.0])), "Wrong value returned for position 7 of item 4 sequence in test_TrainDataset"
        assert torch.equal(sequence[6], torch.tensor([0.0, 0.0, 0.0, 1.0])), "Wrong value returned for position 6 of item 4 sequence in test_TrainDataset"
        assert torch.equal(sequence[5], torch.tensor([0.0, 1.0, 0.0, 0.0])), "Wrong value returned for position 5 of item 4 sequence in test_TrainDataset"
        assert torch.equal(sequence[4], torch.tensor([1.0, 0.0, 0.0, 0.0])), "Wrong value returned for position 4 of item 4 sequence in test_TrainDataset"
        assert torch.equal(sequence[3], torch.tensor([0.0, 0.0, 1.0, 0.0])), "Wrong value returned for position 3 of item 4 sequence in test_TrainDataset"
        assert torch.equal(sequence[2], torch.tensor([0.0, 0.0, 0.0, 1.0])), "Wrong value returned for position 2 of item 4 sequence in test_TrainDataset"
        assert torch.equal(sequence[1], torch.tensor([0.0, 1.0, 0.0, 0.0])), "Wrong value returned for position 1 of item 4 sequence in test_TrainDataset"
        assert torch.equal(sequence[0], torch.tensor([1.0, 0.0, 0.0, 0.0])), "Wrong value returned for position 0 of item 4 sequence in test_TrainDataset"
        assert torch.equal(states[7], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 7 of item 4 in test_TrainDataset"
        assert torch.equal(states[6], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 6 of item 4 in test_TrainDataset"
        assert torch.equal(states[5], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 5 of item 4 in test_TrainDataset"
        assert torch.equal(states[4], torch.tensor([0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 4 of item 4 in test_TrainDataset"
        assert torch.equal(states[3], torch.tensor([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 3 of item 4 in test_TrainDataset"
        assert torch.equal(states[2], torch.tensor([0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 2 of item 4 in test_TrainDataset"
        assert torch.equal(states[1], torch.tensor([0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 1 of item 4 in test_TrainDataset"
        assert torch.equal(states[0], torch.tensor([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 0 of item 4 in test_TrainDataset"

        (name, sequence, states) = test_dataset[5]
        assert name == "Seq6_0    ", "Got wrong name {name} for item 5 in test_TrainDataset"
        assert torch.equal(sequence[0], torch.tensor([0.0, 0.0, 1.0, 0.0])), "Wrong value returned for position 0 of item 5 sequence in test_TrainDataset"
        assert torch.equal(sequence[1], torch.tensor([0.0, 0.0, 0.0, 1.0])), "Wrong value returned for position 1 of item 5 sequence in test_TrainDataset"
        assert torch.equal(sequence[2], torch.tensor([0.0, 1.0, 0.0, 0.0])), "Wrong value returned for position 2 of item 5 sequence in test_TrainDataset"
        assert torch.equal(sequence[3], torch.tensor([1.0, 0.0, 0.0, 0.0])), "Wrong value returned for position 3 of item 5 sequence in test_TrainDataset"
        assert torch.equal(sequence[4], torch.tensor([0.0, 0.0, 1.0, 0.0])), "Wrong value returned for position 4 of item 5 sequence in test_TrainDataset"
        assert torch.equal(sequence[5], torch.tensor([0.0, 0.0, 0.0, 1.0])), "Wrong value returned for position 5 of item 5 sequence in test_TrainDataset"
        assert torch.equal(sequence[6], torch.tensor([0.0, 1.0, 0.0, 0.0])), "Wrong value returned for position 6 of item 5 sequence in test_TrainDataset"
        assert torch.equal(sequence[7], torch.tensor([1.0, 0.0, 0.0, 0.0])), "Wrong value returned for position 7 of item 5 sequence in test_TrainDataset"


        assert torch.equal(states[0], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 0 of item 5 in test_TrainDataset"
        assert torch.equal(states[1], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 1 of item 5 in test_TrainDataset"
        assert torch.equal(states[2], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 2 of item 5 in test_TrainDataset"
        assert torch.equal(states[3], torch.tensor([0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 3 of item 5 in test_TrainDataset"
        assert torch.equal(states[4], torch.tensor([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 4 of item 5 in test_TrainDataset"
        assert torch.equal(states[5], torch.tensor([0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 5 of item 5 in test_TrainDataset"
        assert torch.equal(states[6], torch.tensor([0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 6 of item 5 in test_TrainDataset"
        assert torch.equal(states[7], torch.tensor([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 7 of item 5 in test_TrainDataset"


        (name, sequence, states) = test_dataset[0]
        assert name == "Seq1_0    ", "Got wrong name {name} for item 0 in test_TrainDataset"
        assert torch.equal(sequence[7], torch.tensor([0.0, 0.0, 1.0, 0.0])), "Wrong value returned for position 7 of item 0 sequence in test_TrainDataset"
        assert torch.equal(sequence[6], torch.tensor([0.0, 0.0, 0.0, 1.0])), "Wrong value returned for position 6 of item 0 sequence in test_TrainDataset"
        assert torch.equal(sequence[5], torch.tensor([0.0, 1.0, 0.0, 0.0])), "Wrong value returned for position 5 of item 0 sequence in test_TrainDataset"
        assert torch.equal(sequence[4], torch.tensor([1.0, 0.0, 0.0, 0.0])), "Wrong value returned for position 4 of item 0 sequence in test_TrainDataset"
        assert torch.equal(sequence[3], torch.tensor([0.0, 0.0, 1.0, 0.0])), "Wrong value returned for position 3 of item 0 sequence in test_TrainDataset"
        assert torch.equal(sequence[2], torch.tensor([0.0, 0.0, 0.0, 1.0])), "Wrong value returned for position 2 of item 0 sequence in test_TrainDataset"
        assert torch.equal(sequence[1], torch.tensor([0.0, 1.0, 0.0, 0.0])), "Wrong value returned for position 1 of item 0 sequence in test_TrainDataset"
        assert torch.equal(sequence[0], torch.tensor([1.0, 0.0, 0.0, 0.0])), "Wrong value returned for position 0 of item 0 sequence in test_TrainDataset"
        assert torch.equal(states[7], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 7 of item 0 in test_TrainDataset"
        assert torch.equal(states[6], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 6 of item 0 in test_TrainDataset"
        assert torch.equal(states[5], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 5 of item 0 in test_TrainDataset"
        assert torch.equal(states[4], torch.tensor([0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 4 of item 0 in test_TrainDataset"
        assert torch.equal(states[3], torch.tensor([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 3 of item 0 in test_TrainDataset"
        assert torch.equal(states[2], torch.tensor([0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 2 of item 0 in test_TrainDataset"
        assert torch.equal(states[1], torch.tensor([0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 1 of item 0 in test_TrainDataset"
        assert torch.equal(states[0], torch.tensor([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 0 of item 0 in test_TrainDataset"

        (name, sequence, states) = test_dataset[1]
        assert name == "Seq2_0    ", "Got wrong name {name} for item 1 in test_TrainDataset"
        assert torch.equal(sequence[0], torch.tensor([0.0, 0.0, 1.0, 0.0])), "Wrong value returned for position 0 of item 1 sequence in test_TrainDataset"
        assert torch.equal(sequence[1], torch.tensor([0.0, 0.0, 0.0, 1.0])), "Wrong value returned for position 1 of item 1 sequence in test_TrainDataset"
        assert torch.equal(sequence[2], torch.tensor([0.0, 1.0, 0.0, 0.0])), "Wrong value returned for position 2 of item 1 sequence in test_TrainDataset"
        assert torch.equal(sequence[3], torch.tensor([1.0, 0.0, 0.0, 0.0])), "Wrong value returned for position 3 of item 1 sequence in test_TrainDataset"
        assert torch.equal(sequence[4], torch.tensor([0.0, 0.0, 1.0, 0.0])), "Wrong value returned for position 4 of item 1 sequence in test_TrainDataset"
        assert torch.equal(sequence[5], torch.tensor([0.0, 0.0, 0.0, 1.0])), "Wrong value returned for position 5 of item 1 sequence in test_TrainDataset"
        assert torch.equal(sequence[6], torch.tensor([0.0, 1.0, 0.0, 0.0])), "Wrong value returned for position 6 of item 1 sequence in test_TrainDataset"
        assert torch.equal(sequence[7], torch.tensor([1.0, 0.0, 0.0, 0.0])), "Wrong value returned for position 7 of item 1 sequence in test_TrainDataset"




        assert torch.equal(states[0], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 0 of item 1 in test_TrainDataset"
        assert torch.equal(states[1], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 1 of item 1 in test_TrainDataset"
        assert torch.equal(states[2], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 2 of item 1 in test_TrainDataset"
        assert torch.equal(states[3], torch.tensor([0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 3 of item 1 in test_TrainDataset"
        assert torch.equal(states[4], torch.tensor([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 4 of item 1 in test_TrainDataset"
        assert torch.equal(states[5], torch.tensor([0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 5 of item 1 in test_TrainDataset"
        assert torch.equal(states[6], torch.tensor([0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 6 of item 1 in test_TrainDataset"
        assert torch.equal(states[7], torch.tensor([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 7 of item 1 in test_TrainDataset"




        (name, sequence, states) = test_dataset[2]
        assert name == "Seq3_0    ", "Got wrong name {name} for item 0 in test_TrainDataset"
        assert torch.equal(sequence[7], torch.tensor([0.0, 0.0, 1.0, 0.0])), "Wrong value returned for position 7 of item 2 sequence in test_TrainDataset"
        assert torch.equal(sequence[6], torch.tensor([0.0, 0.0, 0.0, 1.0])), "Wrong value returned for position 6 of item 2 sequence in test_TrainDataset"
        assert torch.equal(sequence[5], torch.tensor([0.0, 1.0, 0.0, 0.0])), "Wrong value returned for position 5 of item 2 sequence in test_TrainDataset"
        assert torch.equal(sequence[4], torch.tensor([1.0, 0.0, 0.0, 0.0])), "Wrong value returned for position 4 of item 2 sequence in test_TrainDataset"
        assert torch.equal(sequence[3], torch.tensor([0.0, 0.0, 1.0, 0.0])), "Wrong value returned for position 3 of item 2 sequence in test_TrainDataset"
        assert torch.equal(sequence[2], torch.tensor([0.0, 0.0, 0.0, 1.0])), "Wrong value returned for position 2 of item 2 sequence in test_TrainDataset"
        assert torch.equal(sequence[1], torch.tensor([0.0, 1.0, 0.0, 0.0])), "Wrong value returned for position 1 of item 2 sequence in test_TrainDataset"
        assert torch.equal(sequence[0], torch.tensor([1.0, 0.0, 0.0, 0.0])), "Wrong value returned for position 0 of item 2 sequence in test_TrainDataset"
        assert torch.equal(states[7], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 7 of item 2 in test_TrainDataset"
        assert torch.equal(states[6], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 6 of item 2 in test_TrainDataset"
        assert torch.equal(states[5], torch.tensor([0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 5 of item 2 in test_TrainDataset"
        assert torch.equal(states[4], torch.tensor([0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 4 of item 2 in test_TrainDataset"
        assert torch.equal(states[3], torch.tensor([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 3 of item 2 in test_TrainDataset"
        assert torch.equal(states[2], torch.tensor([0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 2 of item 2 in test_TrainDataset"
        assert torch.equal(states[1], torch.tensor([0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 1 of item 2 in test_TrainDataset"
        assert torch.equal(states[0], torch.tensor([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])), "Wrong FSM states returned for position 0 of item 2 in test_TrainDataset"
        print("All tests passed for TrainDataset")


    finally:
        shutil.rmtree(tempdir)