import random
from reaction_labeler import label_reaction
from data_cleaner import clean_reaction_record


def build_dataset(raw_records, seed=42):
    cleaned = []

    for r in raw_records:
        cr = clean_reaction_record(r)
        if not cr:
            continue

        rxn_class = label_reaction(cr["reactants"], cr["products"])
        cr["reaction_class"] = rxn_class
        cleaned.append(cr)

    random.seed(seed)
    random.shuffle(cleaned)

    split = int(0.8 * len(cleaned))
    train = cleaned[:split]
    val = cleaned[split:]

    return train, val
