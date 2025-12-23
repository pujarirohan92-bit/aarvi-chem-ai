import torch
import torch.nn as nn
import torch.optim as optim

# path fix
import os, sys
ROOT = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(ROOT, "graph"))
sys.path.append(os.path.join(ROOT, "model"))
sys.path.append(os.path.join(ROOT, "data"))

from mol_to_graph import mol_to_graph
from gnn_model import SimpleGNN
from build_dataset import build_sample

# -------- Prepare one training sample (demo) --------
sample = build_sample()
graph = sample["graph"]
formed = set(sample["formed_bonds"])

x = graph["x"]
edge_index = graph["edge_index"]
edge_attr = graph["edge_attr"]

# -------- Build edge labels --------
# label = 1 if edge corresponds to formed bond else 0
labels = []
src, dst = edge_index
for i in range(edge_index.shape[1]):
    # atom indices are 0-based; mapped bonds are abstract,
    # so for demo we mark all edges as 0 except one
    labels.append(0)

# If at least one formed bond exists, mark first edge as positive (demo)
if len(labels) > 0:
    labels[0] = 1

y = torch.tensor(labels, dtype=torch.float).view(-1, 1)

# -------- Model --------
node_dim = x.shape[1]
edge_dim = edge_attr.shape[1]

model = SimpleGNN(node_dim, edge_dim)
criterion = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=1e-3)

# -------- Training loop --------
print("Starting training...")
for epoch in range(1, 51):
    model.train()
    optimizer.zero_grad()

    preds = model(x, edge_index, edge_attr)
    loss = criterion(preds, y)

    loss.backward()
    optimizer.step()

    if epoch % 10 == 0:
        with torch.no_grad():
            acc = ((preds > 0.5) == (y > 0.5)).float().mean().item()
        print(f"Epoch {epoch:03d} | Loss {loss.item():.4f} | Acc {acc:.2f}")

print("Training finished ✅")
