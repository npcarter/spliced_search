import torch
from torch.utils.data import DataLoader
from lib.dataloader import trainloader
import argparse

parser = argparse.ArgumentParser(description='Top-level training loop')
parser.add_argument('traindb')
parser.add_argument('checkpointloc')
parser.add_argument('output_weightsfile')
parser.add_argument('--batch_size', '-b', default=16, type=int)
parser.add_argument('--epochs', '-e', default = 1, type=int)
parser.add_argument('--num_workers', '-n', default = 4, type=int)
parser.add_argument('--shuffle', '-s', default=True, type=bool)
parser.add_argument('--pin_memory', '-p', default=True, type=bool)
parser.add_argument('--drop_last', '-d', default = False, type = bool)
parser.add_argument('--prefetch_factor', '-f', default = 2, type=int)
parser.add_argument('--sequence_length', '-l', default = 16384, type=int)
parser.add_argument('--name_length', '-m', default = 40, type=int)


args = parser.parse_args()

def synchronize_device(device):
    if device.type == "cuda":
        torch.cuda.synchronize()
    elif device.type == "mps":
        torch.mps.synchronize()
    elif device.type == "cpu":
        pass # CPU operations are synchronous by default

def main(args):
  device = torch.device("cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu")
  if device != "cuda": #pinning memory isn't supported on mps
     pin_memory = False
  else:
     pin_memory = args.pin_memory
  epochs = 1
  the_dataset = trainloader.TrainDataset(args.traindb, sequence_length=args.sequence_length, name_length=args.name_length)
  train_loader = DataLoader(the_dataset, batch_size=args.batch_size, shuffle=args.shuffle, num_workers=args.num_workers,
                                pin_memory=pin_memory, drop_last=args.drop_last)#, prefetch_factor=args.prefetch_factor)

  found_names = {}
  for epoch in range(epochs):
    for batch_idx, (names, sequences, labels) in enumerate(train_loader):
      sequences.to(device)
      labels.to(device)
      for name in names:
        #print(f"{batch_idx} {name}")
        assert name not in found_names, "repeat name {name} found"
        found_names[name] = 1
      print(batch_idx)
      print(names[0])
  synchronize_device(device)

if __name__=='__main__':
  main(args)