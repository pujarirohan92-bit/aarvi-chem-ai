from train_model import train_model
from ml_engine import load_model


def initialize_ml(raw_records):
    model = train_model(raw_records)
    load_model(model)
    print("✅ ML model loaded into engine")
