from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from dataset_builder import build_dataset
from feature_builder import build_features


def train_model(raw_records):
    train, val = build_dataset(raw_records)

    X_train, y_train = build_features(train)
    X_val, y_val = build_features(val)

    model = RandomForestClassifier(
        n_estimators=100,
        random_state=42
    )

    model.fit(X_train, y_train)

    preds = model.predict(X_val)
    acc = accuracy_score(y_val, preds)

    print(f"Validation Accuracy: {acc:.2f}")

    return model
