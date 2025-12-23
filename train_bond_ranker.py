import torch
import torch.nn as nn
import torch.optim as optim
import csv


class BondRankModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(3, 8),
            nn.ReLU(),
            nn.Linear(8, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        return self.net(x)


# Load dataset
X = []
y = []

with open("dataset/bond_training.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        X.append([
            float(row["priority"]),
            float(row["bond_type_len"]),
            float(row["reason_len"])
        ])
        y.append([float(row["label"])])

X = torch.tensor(X, dtype=torch.float32)
y = torch.tensor(y, dtype=torch.float32)

model = BondRankModel()
criterion = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=0.01)

print("🔄 Training bond ranker...")

for epoch in range(800):
    optimizer.zero_grad()
    output = model(X)
    loss = criterion(output, y)
    loss.backward()
    optimizer.step()

print("✅ Training complete")

torch.save(model.state_dict(), "bond_ranker.pt")
print("💾 Model saved as bond_ranker.pt")
