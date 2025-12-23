from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error
from dataset_builder import build_dataset
from yield_feature_builder import build_yield_features


def train_yield_model(raw_records):
    train, val = build_dataset(raw_records)

    X_train, y_train = build_yield_features(train)
    X_val, y_val = build_yield_features(val)

    model = RandomForestRegressor(
        n_estimators=150,
        random_state=42
    )

    model.fit(X_train, y_train)

    preds = model.predict(X_val)
    mae = mean_absolute_error(y_val, preds)

    print(f"Yield MAE: {mae:.2f}")

    return model
