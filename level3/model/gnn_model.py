import torch
import torch.nn as nn
import torch.nn.functional as F


class SimpleGNN(nn.Module):
    def __init__(self, node_dim, edge_dim, hidden_dim=64):
        super().__init__()

        self.node_fc = nn.Linear(node_dim, hidden_dim)
        self.edge_fc = nn.Linear(edge_dim, hidden_dim)

        self.out = nn.Linear(hidden_dim * 2, 1)

    def forward(self, x, edge_index, edge_attr):
        x = F.relu(self.node_fc(x))
        e = F.relu(self.edge_fc(edge_attr))

        src, dst = edge_index
        edge_feat = torch.cat([x[src], e], dim=1)

        logits = self.out(edge_feat)
        return torch.sigmoid(logits)
