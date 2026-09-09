from .. import compress
import numpy as np
def test_FsmCharToNybble():
  for state in "0123456789ABCDE":
          nybble = compress.FsmCharToNybble(state)
          assert nybble is not None, f"FsmCharToNybble returned None for valid state {state}"
          assert nybble == "0123456789ABCDE".index(state), f"FsmCharToNybble returned incorrect value for state {state}"
  invalid = compress.FsmCharToNybble("Z")
  assert invalid is None, "FsmCharToNybble did not return None for invalid state 'Z'"
  print("All tests passed for FsmCharToNybble")

def test_FsmNybbleToChar():
  for nybble in range(15):
          state = compress.FsmNybbleToChar(nybble)
          assert state is not None, f"FsmNybbleToChar returned None for valid nybble {nybble}"
          assert state == "0123456789ABCDE"[nybble], f"FsmNybbleToChar returned incorrect value for nybble {nybble}"
  invalid = compress.FsmNybbleToChar(15)
  assert invalid is None, "FsmNybbleToChar did not return None for invalid nybble 15"
  print("All tests passed for FsmNybbleToChar")

def test_FsmCharToNybble_and_FsmNybbleToChar():
  for state in "0123456789ABCDE":
    nybble = compress.FsmCharToNybble(state)
    assert nybble is not None, f"FsmCharToNybble returned None for valid state {state}"
    assert nybble == "0123456789ABCDE".index(state), f"FsmCharToNybble returned incorrect value for state {state}"
        
    recovered_state = compress.FsmNybbleToChar(nybble)
    assert recovered_state is not None, f"FsmNybbleToChar returned None for valid nybble {nybble}"
    assert recovered_state == state, f"FsmNybbleToChar did not recover the original state {state} from nybble {nybble}"
    
  print("All tests passed for FsmCharToNybble and FsmNybbleToChar")

def test_FsmNybbleToChar_and_FsmCharToNybble():
  for nybble in range(15):
    state = compress.FsmNybbleToChar(nybble)
    assert state is not None, f"FsmNybbleToChar returned None for valid nybble {nybble}"
    assert state == "0123456789ABCDE"[nybble], f"FsmNybbleToChar returned incorrect value for nybble {nybble}"
        
    recovered_nybble = compress.FsmCharToNybble(state)
    assert recovered_nybble is not None, f"FsmCharToNybble returned None for valid state {state}"
    assert recovered_nybble == nybble, f"FsmCharToNybble did not recover the original nybble {nybble} from state {state}"
    
  print("All tests passed for FsmNybbleToChar and FsmCharToNybble")

def test_compressTraining():
  name = "Fnord"
  name_bytes =np.array([ord(c) for c in name])
  sequence = "ACGTACGT"
  fsmStates = "01234567"
  compressed = compress.CompressTraining(name, sequence, fsmStates)
  name_bytes = np.append(name_bytes, 0xff)  # Append the terminator for name
  name_bytes = np.append(name_bytes, 0x00)  # hand-generate compressed values
  name_bytes = np.append(name_bytes, 0x11)
  name_bytes = np.append(name_bytes, 0x22)
  name_bytes = np.append(name_bytes, 0x33)
  name_bytes = np.append(name_bytes, 0x04)
  name_bytes = np.append(name_bytes, 0x15)
  name_bytes = np.append(name_bytes, 0x26)
  name_bytes = np.append(name_bytes, 0x37)
  print(name_bytes)
  assert np.array_equal(compressed,name_bytes), f"Unexpected compressed output: {compressed}"
  print("All tests passed for compressTraining")

def test_decompressTraining():
  name = "Fnord"
  sequence = "ACGTACGT"
  fsmStates = "01234567"
  compressed = np.array([ord(c) for c in name])  # Convert name to bytes
  compressed = np.append(compressed, 0xff)  # Append the terminator for name
  compressed = np.append(compressed, 0x00)  # hand-generate compressed values
  compressed = np.append(compressed, 0x11)
  compressed = np.append(compressed, 0x22)
  compressed = np.append(compressed, 0x33)
  compressed = np.append(compressed, 0x04)
  compressed = np.append(compressed, 0x15)
  compressed = np.append(compressed, 0x26)
  compressed = np.append(compressed, 0x37)
  decompressed_name, decompressed_sequence, decompressed_fsmStates = compress.DecompressTraining(compressed)
  assert decompressed_name == name, f"Decompressed name '{decompressed_name}' does not match original '{name}'"
  assert decompressed_sequence == sequence, f"Decompressed sequence '{decompressed_sequence}' does not match original '{sequence}'"
  assert decompressed_fsmStates == fsmStates, f"Decompressed FSM states '{decompressed_fsmStates}' do not match original '{fsmStates}'"
  print("All tests passed for decompressTraining")

  def test_compress_and_decompress_training():
    name = "Fnord"
    sequence = "ACGTACGT"
    fsmStates = "01234567"

    # Compress the training data
    compressed = compress.CompressTraining(name, sequence, fsmStates)

    # Decompress the training data
    decompressed_name, decompressed_sequence, decompressed_fsmStates = compress.DecompressTraining(compressed)

    # Verify that the decompressed data matches the original data
    assert decompressed_name == name, f"Decompressed name '{decompressed_name}' does not match original '{name}'"
    assert decompressed_sequence == sequence, f"Decompressed sequence '{decompressed_sequence}' does not match original '{sequence}'"
    assert decompressed_fsmStates == fsmStates, f"Decompressed FSM states '{decompressed_fsmStates}' do not match original '{fsmStates}'"

    print("All tests passed for compress_and_decompress_training")