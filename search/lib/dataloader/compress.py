import numpy as np

def FsmCharToNybble(char):
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

def FsmNybbleToChar(nybble):
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

def DnaCharToNybble(char):

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
    case 'N':
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
    case 'X':
      return 15
    case _:
      return None

def DnaNybbleToChar(nybble):
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
      return 'X'
    case _:
      return None

def CompressTraining(name, sequence, fsmStates):
  """
  Converts the strings representing a sequence's name, DNA sequence, and FSM states 
  into a vector of bytes in compressed format.  
  The compressed data begins with the ASCII values of the characters in the name, followed by a 0xff terminator. This is followed by a compressed representation of the DNA sequence and FSM states, where each nucleotide and FSM state are represented by four-bit values that are packed into bytes. A byte of value 0xff indicates the boundary between the sequence name and the compressed sequence/FSM state data.  The function returns a numpy array of type uint8 containing the compressed data.  
  """
  if len(sequence) != len(fsmStates):
    raise ValueError("Sequence and FSM states must be of the same length.")
  compressed_data = np.array([], dtype=np.uint8)

  #step one: copy the sequence name into the compressed data
  for ch in name:
    by = ord(ch)
    if by > 255:
      raise ValueError(f"Character '{ch}' in name has ASCII value greater than 255.")
    compressed_data = np.append(compressed_data, by)

  #Mark the transition to the compressed sequence and FSM states with a 0xff byte
  compressed_data = np.append(compressed_data, 0xff)  # Terminator for name

  #Now, create the compressed representation of the sequence and FSM states
  for i in range(0, len(sequence)):
    # Convert nucleotide and FSM state to nybble values
    nucleotide_nybble = DnaCharToNybble(sequence[i])
    fsm_nybble = FsmCharToNybble(fsmStates[i])

    if nucleotide_nybble is None:
      raise ValueError(f"Invalid nucleotide character '{sequence[i]}' at position {i}.")
    if fsm_nybble is None:
      raise ValueError(f"Invalid FSM state character '{fsmStates[i]}' at position {i}.")

    # Pack two nybbles into one byte
    packed_byte = (nucleotide_nybble << 4) | fsm_nybble
    compressed_data = np.append(compressed_data, packed_byte)
  return compressed_data

def DecompressTraining(compressed_data):
  """
  Decompresses the compressed training data back into its original components: name, sequence, and FSM states.  
  The function takes a numpy array of type uint8 containing the compressed data and returns a tuple containing the name (string), sequence (string), and FSM states (string).  
  """
  if len(compressed_data) == 0:
    raise ValueError("Compressed data is empty.")

  # Find the index of the 0xff terminator that separates the name from the compressed sequence/FSM states
  terminator_index = np.where(compressed_data == 0xff)[0]
  if len(terminator_index) == 0:
    raise ValueError("No terminator (0xff) found in compressed data.")
  
  terminator_index = terminator_index[0]

  # Extract the name from the compressed data
  name_bytes = compressed_data[:terminator_index]
  name = ''.join(chr(by) for by in name_bytes)

  # Extract the compressed sequence and FSM states
  compressed_sequence_fsm = compressed_data[terminator_index + 1:]

  sequence = []
  fsm_states = []

  for byte in compressed_sequence_fsm:
    nucleotide_nybble = (byte >> 4) & 0x0F
    fsm_nybble = byte & 0x0F

    nucleotide_char = DnaNybbleToChar(nucleotide_nybble)
    fsm_char = FsmNybbleToChar(fsm_nybble)

    if nucleotide_char is None:
      raise ValueError(f"Invalid nucleotide nybble value '{nucleotide_nybble}' in compressed data.")
    if fsm_char is None:
      raise ValueError(f"Invalid FSM state nybble value '{fsm_nybble}' in compressed data.")

    sequence.append(nucleotide_char)
    fsm_states.append(fsm_char)

  return name, ''.join(sequence), ''.join(fsm_states)

def CompressFile(input_file_path, output_file_path):
  """
  Compresses the contents of a file containing training data and writes the compressed data to an output file.  
  The input file is expected to contain three lines: the first line is the name, the second line is the DNA sequence, and the third line is the FSM states.  
  The function reads these lines, compresses them using the CompressTraining function, and writes the resulting compressed data to the specified output file in binary format.  
  """
  line_count = 0

  with open(input_file_path, 'r') as infile, open(output_file_path, 'wb') as outfile:
    for line in infile:
      if line_count == 0:
        name = line.strip()
      elif line_count == 1:
        sequence = line.strip()
      elif line_count == 2:
        fsm_states = line.strip()
        compressed_data = CompressTraining(name, sequence, fsm_states)
        outfile.write(compressed_data.tobytes())
      line_count = (line_count + 1)% 3
