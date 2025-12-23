import csv
from pathlib import Path


DATASET_PATH = Path("dataset/reactions.csv")


def load_reaction_dataset():
    reactions = []
    with open(DATASET_PATH, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            reactions.append(row)
    return reactions


_DATASET = load_reaction_dataset()


def classify_reaction(step: dict) -> str:
    """
    Very simple matcher:
    Later → ML multiclass classifier
    """
    reaction_text = step.get("reaction", "").lower()

    for row in _DATASET:
        if row["reaction_class"].split("_")[0] in reaction_text:
            return row["reaction_class"]

    return "unknown_reaction"
