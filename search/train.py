import torch
from torch.utils.data import DataLoader
import trainloader


device = torch.device("cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu")

epochs = 1

train_loader = torch.DataLoader()

for epoch in range(epochs):
