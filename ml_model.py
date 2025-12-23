import torch
import torch.nn as nn


class ReactionFeasibilityModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(5, 16),
            nn.ReLU(),
            nn.Linear(16, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        return self.net(x)


def featurize_step(step: dict):
    return torch.tensor([
        len(step.get("reaction", "")),
        len(step.get("reactants", [])),
        len(step.get("products", [])),
        1 if step.get("warning") else 0,
        len(step.get("conditions", ""))
    ], dtype=torch.float32)


def load_model(path="model.pt"):
    model = ReactionFeasibilityModel()
    try:
        model.load_state_dict(torch.load(path))
        model.eval()
    except FileNotFoundError:
        pass
    return model


_model = load_model()


def predict_feasibility(step: dict) -> float:
    with torch.no_grad():
        x = featurize_step(step)
        y = _model(x)
        return round(float(y.item()), 2)
