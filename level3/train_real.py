import os
import sys
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, random_split

# ---------------- PATH SETUP ----------------
ROOT = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(ROOT, "model"))
sys.path.append(os.path.join(ROOT, "data"))

from gnn_model import SimpleGNN
from dataset import ReactionDataset

# ---------------- DATASET ----------------
dataset = ReactionDataset(os.path.join(ROOT, "data", "uspto_sample.csv"))

train_size = int(0.8 * len(dataset))
val_size = len(dataset) - train_size
train_ds, val_ds = random_split(dataset, [train_size, val_size])

train_loader = DataLoader(train_ds, batch_size=1, shuffle=True)
val_loader = DataLoader(val_ds, batch_size=1)

# ---------------- MODEL ----------------
node_dim = 5
edge_dim = 6
model = SimpleGNN(node_dim, edge_dim)

criterion = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=1e-3)

# ---------------- TRAINING LOOP ----------------
for epoch in range(1, 31):
    model.train()
    train_loss = 0.0

    for g, y in train_loader:
        x = g["x"].squeeze(0)
        edge_index = g["edge_index"].squeeze(0)
        edge_attr = g["edge_attr"].squeeze(0)
        y = y.squeeze(0)  # FIX SHAPE

        optimizer.zero_grad()
        preds = model(x, edge_index, edge_attr)
        loss = criterion(preds, y)
        loss.backward()
        optimizer.step()

        train_loss += loss.item()

    # ---------------- VALIDATION ----------------
    model.eval()
    val_loss = 0.0
    with torch.no_grad():
        for g, y in val_loader:
            x = g["x"].squeeze(0)
            edge_index = g["edge_index"].squeeze(0)
            edge_attr = g["edge_attr"].squeeze(0)
            y = y.squeeze(0)  # FIX SHAPE

            preds = model(x, edge_index, edge_attr)
            val_loss += criterion(preds, y).item()

    print(
        f"Epoch {epoch:02d} | "
        f"Train Loss {train_loss:.4f} | "
        f"Val Loss {val_loss:.4f}"
    )

# ---------------- SAVE MODEL ----------------
os.makedirs(os.path.join(ROOT, "checkpoints"), exist_ok=True)
torch.save(model.state_dict(), os.path.join(ROOT, "checkpoints", "rc_gnn.pt"))
print("Model saved to checkpoints/rc_gnn.pt ✅")
