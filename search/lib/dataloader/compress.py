import numpy as np
import zarr

def fsm_char_to_nybble(char):
  """
  Convert a character representing a state in the FSM to its corresponding nybble value.
  Returns None for invalid characters.
  """
  match char:
    case '0':
      return 0
    case '1':
      return 1
    case '2':
      return 2
    case '3':
      return 3
    case '4':
      return 4
    case '5':
      return 5
    case '6':
      return 6
    case '7':
      return 7
    case '8':
      return 8
    case '9':
      return 9
    case 'A':
      return 10
    case 'B':
      return 11
    case 'C':
      return 12
    case 'D':
      return 13
    case 'E':
      return 14
    case _:
      return None

def fsm_nybble_to_char(nybble):
  """
  Convert a nybble value back to its corresponding character representing a state in the FSM.
  Returns None for invalid nybble values.
  """
  match nybble:
    case 0:
      return '0'
    case 1:
      return '1'
    case 2:
      return '2'
    case 3:
      return '3'
    case 4:
      return '4'
    case 5:
      return '5'
    case 6:
      return '6'
    case 7:
      return '7'
    case 8:
      return '8'
    case 9:
      return '9'
    case 10:
      return 'A'
    case 11:
      return 'B'
    case 12:
      return 'C'
    case 13:
      return 'D'
    case 14:
      return 'E'
    case _:
      return None

def dna_char_to_nybble(char):

  match char:
    case 'A':
      return 0
    case 'C':
      return 1
    case 'G':
      return 2
    case 'T':
      return 3
    case 'B':
      return 4
    case 'D':
      return 5
    case 'H':
      return 6
    case 'K':
      return 7
    case 'M':
      return 8
    case 'N': # N and X are both "unknown"
      return 9
    case 'R':
      return 10
    case 'S':
      return 11
    case 'V':
      return 12
    case 'W':
      return 13
    case 'Y':
      return 14
    case 'X': # N and X are both "unknown"
      return 9
    case 'Z': # Z will be the padding character
      return 15
    case _:
      return None

def dna_nybble_to_char(nybble):
  """
  Convert a nybble value back to its corresponding character representing a nucleotide.
  Returns None for invalid nybble values.
  """
  match nybble:
    case 0:
      return 'A'
    case 1:
      return 'C'
    case 2:
      return 'G'
    case 3:
      return 'T'
    case 4:
      return 'B'
    case 5:
      return 'D'
    case 6:
      return 'H'
    case 7:
      return 'K'
    case 8:
      return 'M'
    case 9:
      return 'N'
    case 10:
      return 'R'
    case 11:
      return 'S'
    case 12:
      return 'V'
    case 13:
      return 'W'
    case 14:
      return 'Y'
    case 15:
      return 'Z'
    case _:
      return None

def compress_chunk(name, sequence, fsmStates, sequence_chunk_size=16384, name_length=40):
  """
  Converts the strings representing a chunk of sequence's name, DNA sequence, and FSM states 
  into a vector of bytes in compressed format.  
  The compressed data begins with the ASCII values of the characters in the name, followed by a 0xff terminator. This is followed by a compressed representation of the DNA sequence and FSM states, where each nucleotide and FSM state are represented by four-bit values that are packed into bytes. A byte of value 0xff indicates the boundary between the sequence name and the compressed sequence/FSM state data.  The function returns a numpy array of type uint8 containing the compressed data.  
  It is the responsibility of other code to ensure that the DNA sequence and FSM states are the same length and that they have been trimmed/padded to the correct length.  The function does not perform any validation of the input strings beyond checking that they are the same length.
  """
  if len(name) != name_length:
    raise ValueError(f"Name must be exactly {name_length} characters long.")
  if len(sequence) != sequence_chunk_size:
    raise ValueError(f"Sequence must be exactly {sequence_chunk_size} characters long.")
  if len(fsmStates) != sequence_chunk_size:
    raise ValueError(f"FSM states must be exactly {sequence_chunk_size} characters long.")

  compressed_data = np.zeros([sequence_chunk_size + name_length], dtype=np.uint8)

  #step one: copy the sequence name into the compressed data

  for index, ch in enumerate(name):
    by = ord(ch)
    if by > 255:
      raise ValueError(f"Character '{ch}' in name has ASCII value greater than 255.")
    compressed_data[index] = by

  #Now, create the compressed representation of the sequence and FSM states
  for i in range(0, len(sequence)):
    # Convert nucleotide and FSM state to nybble values
    nucleotide_nybble = dna_char_to_nybble(sequence[i])
    fsm_nybble = fsm_char_to_nybble(fsmStates[i])

    if nucleotide_nybble is None:
      raise ValueError(f"Invalid nucleotide character '{sequence[i]}' at position {i}.")
    if fsm_nybble is None:
      raise ValueError(f"Invalid FSM state character '{fsmStates[i]}' at position {i}.")

    # Pack two nybbles into one byte
    packed_byte = (nucleotide_nybble << 4) | fsm_nybble
    compressed_data[name_length+ i] = packed_byte

  assert compressed_data.shape[0] == name_length + sequence_chunk_size, f"Compressed data length {compressed_data.shape[0]} does not match expected length {name_length + sequence_chunk_size}"
  return compressed_data

def decompress_chunk(compressed_data, sequence_chunk_size=16384, name_length=40):
  """
  Decompresses the compressed training data back into its original components: name, sequence, and FSM states.  
  The function takes a numpy array of type uint8 containing the compressed data and returns a tuple containing the name (string), sequence (string), and FSM states (string).  
  """

  if len(compressed_data) != sequence_chunk_size + name_length:
    raise ValueError(f"Compressed data has incorrect length {len(compressed_data)}. Expected length is {sequence_chunk_size + name_length}.")


  # Extract the name from the compressed data
  name_bytes = compressed_data[:name_length]
  name = ''.join(chr(by) for by in name_bytes)

  # Extract the compressed sequence and FSM states
  compressed_sequence_fsm = compressed_data[name_length:name_length + sequence_chunk_size]

  sequence = []
  fsm_states = []

  for byte in compressed_sequence_fsm:
    nucleotide_nybble = (byte >> 4) & 0x0F
    fsm_nybble = byte & 0x0F

    nucleotide_char = dna_nybble_to_char(nucleotide_nybble)
    fsm_char = fsm_nybble_to_char(fsm_nybble)

    if nucleotide_char is None:
      raise ValueError(f"Invalid nucleotide nybble value '{nucleotide_nybble}' in compressed data.")
    if fsm_char is None:
      raise ValueError(f"Invalid FSM state nybble value '{fsm_nybble}' in compressed data.")

    sequence.append(nucleotide_char)
    fsm_states.append(fsm_char)

  return name, ''.join(sequence), ''.join(fsm_states)

def split_sequence_into_chunks(name, sequence, fsmStates, chunk_size=16384, name_length=40):
  """
  Splits the sequence and FSM states into chunks of specified size.  
  Returns a list of tuples, where each tuple contains the name of a chunk, a chunk of the sequence and the corresponding chunk of FSM states. 
  The name of the chunk is the original name with a suffix "_N" where N is the chunk number starting from 0, padded with space characters to reach the specified name_length.
  If the sequence length is shorter than the chunk size, pads the sequence with 'Z' and the FSM states with '0' (intergenic) to reach the chunk size.
  If the sequence length is not a multiple of the chunk size, the last chunk will be the last
  chunk_size characters of the sequence and FSM states, even though this creates overlap with the previous chunk. 
  """
  if len(sequence) != len(fsmStates):
    raise ValueError("Sequence and FSM states must be of the same length.")

  chunks = []
  if len(sequence) < chunk_size:
    new_seq = sequence.ljust(chunk_size, 'Z')
    new_fsm = fsmStates.ljust(chunk_size, '0')
    new_name = _make_chunk_name(name, 0, name_length)
    return([(new_name, new_seq, new_fsm)])


  for i in range(0, len(sequence), chunk_size):
    if i + chunk_size < len(sequence):
      seq_chunk = sequence[i:i + chunk_size]
      fsm_chunk = fsmStates[i:i + chunk_size]
    else:
      seq_chunk = sequence[len(sequence)-chunk_size:]
      fsm_chunk = fsmStates[len(fsmStates)-chunk_size:]
    new_name = _make_chunk_name(name, i // chunk_size, name_length)
    chunks.append((new_name, seq_chunk, fsm_chunk))


  return chunks

def _make_chunk_name(name, chunk_number, name_length):
  """
  Builds a chunk name in the form "{name}_{chunk_number}" padded with spaces to name_length.
  If the suffix alone doesn't fit within name_length, raises ValueError. Otherwise, name is
  truncated as needed to leave room for the suffix, since the suffix must always be present
  to keep chunk names unique.
  """
  suffix = f"_{chunk_number}"
  if len(suffix) > name_length:
    raise ValueError(f"Chunk suffix '{suffix}' alone exceeds the specified name length of {name_length}.")
  truncated_name = name[:name_length - len(suffix)]
  return f"{truncated_name}{suffix}".ljust(name_length)

def compress_file(input_file_path, output_file_path, sequence_chunk_size=16384, name_length=40, file_chunk_size=1024):
  """
  Compresses the contents of a file containing training data and writes the compressed data to an output file in Zarr format.  Each sequence in the input file is expected to be represented by three lines: the first line contains the name of the sequence, the second line contains the DNA sequence, and the third line contains the FSM states.  The function reads each sequence from the input file, calls SplitSequenceIntoChunks to split the sequence into chunks of the specified size, and then calls CompressSequence to compress each chunk.  The compressed data for each chunk is written to the output file in Zarr format.  The function does not perform any validation of the input file.
  """
  line_count = 0
  
  outfile = zarr.open(
    store = output_file_path,
    shape = (0, sequence_chunk_size + name_length),
    chunks = (file_chunk_size, sequence_chunk_size + name_length),
    dtype = np.uint8,
    mode = 'w'
  )

  index = 0
  with open(input_file_path, 'r') as infile:

    line_count = 0
    outchunk = np.zeros((file_chunk_size, sequence_chunk_size + name_length), dtype=np.uint8)
    for line in infile:
      if line_count == 0:
        name = line.rstrip('\n')
      if line_count == 1:
        sequence = line.rstrip('\n')
      if line_count == 2:
        fsmStates = line.rstrip('\n')
        chunks = split_sequence_into_chunks(name, sequence, fsmStates, chunk_size=sequence_chunk_size, name_length=name_length)
        for chunk in chunks:
          compressed_chunk = compress_chunk(chunk[0], chunk[1], chunk[2], sequence_chunk_size=sequence_chunk_size, name_length=name_length)
          outchunk[index, :] = compressed_chunk
          index += 1
          if index == file_chunk_size:
            outfile.append(outchunk)
            index = 0
      line_count = (line_count + 1) %3

  if index > 0: #There was a partial chunk at the end of the file, so write it to the output file
    outfile.append(outchunk[:index, :])



    
  
