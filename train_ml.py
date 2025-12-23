import torch
import torch.nn as nn
import torch.optim as optim
from ml_model import ReactionFeasibilityModel

# 🔹 Dummy training dataset (STEP-20.3)
# features: [reaction_len, reactants, products, warning, condition_len]
X = torch.tensor([
    [20, 2, 1, 0, 15],
    [35, 3, 1, 1, 30],
    [15, 1, 1, 0, 10],
    [40, 3, 2, 1, 35],
    [25, 2, 1, 0, 20]
], dtype=torch.float32)

# Labels: feasible (1) / risky (0)
y = torch.tensor([[1], [0], [1], [0], [1]], dtype=torch.float32)

model = ReactionFeasibilityModel()
criterion = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=0.01)

print("🔄 Training ML model...")

for epoch in range(500):
    optimizer.zero_grad()
    output = model(X)
    loss = criterion(output, y)
    loss.backward()
    optimizer.step()

print("✅ Training complete")

torch.save(model.state_dict(), "model.pt")
print("💾 Model saved as model.pt")
