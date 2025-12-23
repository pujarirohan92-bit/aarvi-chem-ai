import torch
from gnn_model import SimpleGNN

node_dim = 5
edge_dim = 6

model = SimpleGNN(node_dim, edge_dim)

x = torch.rand(7, node_dim)
edge_index = torch.tensor([[0,1],[1,0]])
edge_attr = torch.rand(2, edge_dim)

out = model(x, edge_index, edge_attr)
print("Output shape:", out.shape)
