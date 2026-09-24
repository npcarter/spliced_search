from .. import compress
import numpy as np
def test_fsm_char_to_nybble():
  for state in "0123456789ABCDE":
          nybble = compress.fsm_char_to_nybble(state)
          assert nybble is not None, f"fsm_char_to_nybble returned None for valid state {state}"
          assert nybble == "0123456789ABCDE".index(state), f"fsm_char_to_nybble returned incorrect value for state {state}"
  invalid = compress.fsm_char_to_nybble("Z")
  assert invalid is None, "fsm_char_to_nybble did not return None for invalid state 'Z'"
  print("All tests passed for fsm_char_to_nybble")

def test_fsm_nybble_to_char():
  for nybble in range(15):
          state = compress.fsm_nybble_to_char(nybble)
          assert state is not None, f"fsm_nybble_to_char returned None for valid nybble {nybble}"
          assert state == "0123456789ABCDE"[nybble], f"fsm_nybble_to_char returned incorrect value for nybble {nybble}"
  invalid = compress.fsm_nybble_to_char(15)
  assert invalid is None, "fsm_nybble_to_char did not return None for invalid nybble 15"
  print("All tests passed for fsm_nybble_to_char")

def test_fsm_char_to_nybble_and_fsm_nybble_to_char():
  for state in "0123456789ABCDE":
    nybble = compress.fsm_char_to_nybble(state)
    assert nybble is not None, f"fsm_char_to_nybble returned None for valid state {state}"
    assert nybble == "0123456789ABCDE".index(state), f"fsm_char_to_nybble returned incorrect value for state {state}"
        
    recovered_state = compress.fsm_nybble_to_char(nybble)
    assert recovered_state is not None, f"fsm_nybble_to_char returned None for valid nybble {nybble}"
    assert recovered_state == state, f"fsm_nybble_to_char did not recover the original state {state} from nybble {nybble}"
    
  print("All tests passed for fsm_char_to_nybble and fsm_nybble_to_char")

def test_dna_char_to_nybble():
  for state in "ACGTBDHKMNRSVWYXZ":
          nybble = compress.dna_char_to_nybble(state)
          assert nybble is not None, f"dna_char_to_nybble returned None for valid state {state}"
          assert nybble == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 9, 15]["ACGTBDHKMNRSVWYXZ".index(state)], f"dna_char_to_nybble returned incorrect value for state {state}"
  invalid = compress.dna_char_to_nybble("Q")
  assert invalid is None, "dna_char_to_nybble did not return None for invalid state 'Q'"
  print("All tests passed for dna_char_to_nybble")

def test_dna_nybble_to_char():
  for nybble in range(16):
          state = compress.dna_nybble_to_char(nybble)
          assert state is not None, f"dna_nybble_to_char returned None for valid nybble {nybble}"
          assert state == "ACGTBDHKMNRSVWYZ"[nybble], f"dna_nybble_to_char returned incorrect value for nybble {nybble}"
  invalid = compress.dna_nybble_to_char(16)
  assert invalid is None, "dna_nybble_to_char did not return None for invalid nybble 16"
  print("All tests passed for dna_nybble_to_char")

def test_dna_char_to_nybble_and_dna_nybble_to_char():
  for state in "ACGTBDHKMNRSVWYXZ":
    nybble = compress.dna_char_to_nybble(state)
    assert nybble is not None, f"dna_char_to_nybble returned None for valid state {state}"
    assert nybble == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 9, 15]["ACGTBDHKMNRSVWYXZ".index(state)], f"dna_char_to_nybble returned incorrect value for state {state}"
        
    recovered_state = compress.dna_nybble_to_char(nybble)
    assert recovered_state is not None, f"dna_nybble_to_char returned None for valid nybble {nybble}"
    assert (recovered_state == state or (state == 'X' and recovered_state == 'N')), f"dna_nybble_to_char did not recover the original state {state} from nybble {nybble}"
    
  print("All tests passed for dna_char_to_nybble and dna_nybble_to_char")

def test_fsm_nybble_to_char_and_fsm_char_to_nybble():
  for nybble in range(15):
    state = compress.fsm_nybble_to_char(nybble)
    assert state is not None, f"fsm_nybble_to_char returned None for valid nybble {nybble}"
    assert state == "0123456789ABCDE"[nybble], f"fsm_nybble_to_char returned incorrect value for nybble {nybble}"
        
    recovered_nybble = compress.fsm_char_to_nybble(state)
    assert recovered_nybble is not None, f"fsm_char_to_nybble returned None for valid state {state}"
    assert recovered_nybble == nybble, f"fsm_char_to_nybble did not recover the original nybble {nybble} from state {state}"
    
  print("All tests passed for fsm_nybble_to_char and fsm_char_to_nybble")

def test_compress_chunk():
  name = "Fnord"
  name_bytes =np.array([ord(c) for c in name])
  sequence = "ACGTACGT"
  fsmStates = "01234567"
  compressed = compress.compress_chunk(name, sequence, fsmStates, sequence_chunk_size=8, name_length=5)
  name_bytes = np.append(name_bytes, 0x00)  # hand-generate compressed values
  name_bytes = np.append(name_bytes, 0x11)
  name_bytes = np.append(name_bytes, 0x22)
  name_bytes = np.append(name_bytes, 0x33)
  name_bytes = np.append(name_bytes, 0x04)
  name_bytes = np.append(name_bytes, 0x15)
  name_bytes = np.append(name_bytes, 0x26)
  name_bytes = np.append(name_bytes, 0x37)
  assert np.array_equal(compressed,name_bytes), f"Unexpected compressed output: {compressed}"
  print("All tests passed for compress_chunk")

def test_decompress_chunk():
  name = "Fnord"
  sequence = "ACGTACGT"
  fsmStates = "01234567"
  compressed = np.array([ord(c) for c in name])  # Convert name to bytes
  compressed = np.append(compressed, 0x00)  # hand-generate compressed values
  compressed = np.append(compressed, 0x11)
  compressed = np.append(compressed, 0x22)
  compressed = np.append(compressed, 0x33)
  compressed = np.append(compressed, 0x04)
  compressed = np.append(compressed, 0x15)
  compressed = np.append(compressed, 0x26)
  compressed = np.append(compressed, 0x37)
  decompressed_name, decompressed_sequence, decompressed_fsmStates = compress.decompress_chunk(compressed, sequence_chunk_size=8, name_length=5)
  assert decompressed_name == name, f"Decompressed name '{decompressed_name}' does not match original '{name}'"
  assert decompressed_sequence == sequence, f"Decompressed sequence '{decompressed_sequence}' does not match original '{sequence}'"
  assert decompressed_fsmStates == fsmStates, f"Decompressed FSM states '{decompressed_fsmStates}' do not match original '{fsmStates}'"
  print("All tests passed for decompress_sequence")

def test_compress_chunk_and_decompress_chunk():
  name = "Fnord"
  sequence = "ACGTACGT"
  fsmStates = "01234567"

  # Compress the training data
  compressed = compress.compress_chunk(name, sequence, fsmStates, sequence_chunk_size=8, name_length=5)

  # Decompress the training data
  decompressed_name, decompressed_sequence, decompressed_fsmStates = compress.decompress_chunk(compressed, sequence_chunk_size=8, name_length=5)

  # Verify that the decompressed data matches the original data
  assert decompressed_name == name, f"Decompressed name '{decompressed_name}' does not match original '{name}'"
  assert decompressed_sequence == sequence, f"Decompressed sequence '{decompressed_sequence}' does not match original '{sequence}'"
  assert decompressed_fsmStates == fsmStates, f"Decompressed FSM states '{decompressed_fsmStates}' do not match original '{fsmStates}'"

  print("All tests passed for compress_sequence_and_decompress_sequence")

def test_split_sequence_into_chunks():
  #test base case (sequence longer than chunk size)
  name = "Fnord"
  sequence = "ACGTACGTACGTAC"
  fsmStates = "01234567012345"
  chunk_size = 8
  name_length = 10
  chunks = compress.split_sequence_into_chunks(name, sequence, fsmStates, chunk_size, name_length)
  expected_chunks = [
    ("Fnord_0   ", "ACGTACGT", "01234567"),
    ("Fnord_1   ", "GTACGTAC", "67012345") #Sequence is last chunk_size elements of sequence and FSM states
  ]
    
  assert len(chunks) == len(expected_chunks), f"Expected {len(expected_chunks)} chunks, got {len(chunks)}"
    
  for i, (chunk_name, seq_chunk, fsm_chunk) in enumerate(chunks):
    expected_name, expected_seq, expected_fsm = expected_chunks[i]
    assert len(chunk_name) == name_length, f"Chunk {i} name length {len(chunk_name)} does not match expected {name_length}"
    assert chunk_name == expected_name, f"Chunk {i} name '{chunk_name}' does not match expected '{expected_name}'"
    assert seq_chunk == expected_seq, f"Chunk {i} sequence '{seq_chunk}' does not match expected '{expected_seq}'"
    assert fsm_chunk == expected_fsm, f"Chunk {i} FSM states '{fsm_chunk}' do not match expected '{expected_fsm}'"

  #test edge case (sequence shorter than chunk size)
  name = "Fnord"
  sequence = "ACGT"
  fsmStates = "0123"
  chunk_size = 8
  name_length = 10
  chunk = compress.split_sequence_into_chunks(name, sequence, fsmStates, chunk_size, name_length)
  expected_chunk = ("Fnord_0   ", "ACGTZZZZ", "01230000") #Sequence and FSM states are padded to chunk_size 
  assert len(chunk) == 1, f"Expected 1 chunk, got {len(chunk)}"
  assert chunk[0] == expected_chunk, f"Chunk '{chunk}' does not match expected '{expected_chunk}'"
  print("All tests passed for split_sequence_into_chunks")

def test_compress_file():
   # This test will create a temporary input file with known content, compress it, and then check the output file for expected content.
  import tempfile
  import os
  import zarr
  import shutil

  tempdir = tempfile.TemporaryDirectory(delete=False).name
  input_file_path = os.path.join(tempdir, "input.txt")
  # Create a temporary input file
  with open(input_file_path, "wb") as temp_input:
    temp_input.write(b"Seq1\nACGTACGT\n01234567\n")
    temp_input.write(b"Seq2\nTGCATGCA\n76543210\n")
    temp_input.write(b"Seq1\nACGTACGT\n01234567\n")
    temp_input.write(b"Seq2\nTGCATGCA\n76543210\n")
    temp_input.write(b"Seq1\nACGTACGT\n01234567\n")
    temp_input.write(b"Seq2\nTGCATGCA\n76543210\n")

  # Create a temporary output file path
  output_file_path = os.path.join(tempdir, "output.zarr")

  # Compress the input file
  compress.compress_file(input_file_path, output_file_path, sequence_chunk_size=8, name_length=10, file_chunk_size=4)

  # Read the compressed data from the output file using zarr
  compressed_data = zarr.open(output_file_path, mode='r')

  # Check that the compressed data has the expected shape and dtype
  assert compressed_data.shape == (6, 18), f"Expected shape (6, 18), got {compressed_data.shape}"
  assert compressed_data.dtype == np.uint8, f"Expected dtype uint8, got {compressed_data.dtype}"

  print("All tests passed for compress_file")

  shutil.rmtree(tempdir)