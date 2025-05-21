import argparse
import os
import pickle
import warnings
from typing import List

import joblib
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch import nn
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm

warnings.filterwarnings("ignore")

class M6ADataset(Dataset):
    def __init__(self, X: np.ndarray):
        self.X = torch.tensor(X, dtype=torch.float32)
    def __len__(self):
        return len(self.X)
    def __getitem__(self, idx):
        return self.X[idx]

class Attention(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int):
        super().__init__()
        self.attention = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, 1, bias=False),
        )
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        scores = self.attention(x)
        weights = torch.softmax(scores, dim=1)
        return x * weights

class M6APredictor(nn.Module):
    def __init__(
        self,
        hidden_layers: List[int],
        input_dim: int = 56,
        num_classes: int = 2,
        attention_hidden_dim: int = 32,
        dropout_rate: float = 0.5,
    ):
        super().__init__()
        self.attention = Attention(input_dim, attention_hidden_dim)
        layers: List[nn.Module] = []
        prev_dim = input_dim
        for h in hidden_layers:
            layers.extend([nn.Linear(prev_dim, h), nn.ReLU(), nn.Dropout(dropout_rate)])
            prev_dim = h
        self.hidden_layers = nn.Sequential(*layers)
        self.output_layer = nn.Linear(prev_dim, num_classes)
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.attention(x)
        x = self.hidden_layers(x)
        return self.output_layer(x)

def read_test_data(file_path: str) -> pd.DataFrame:
    print(f"Start reading data files: {file_path}")
    df = pd.read_csv(file_path, sep="\t", header=0, on_bad_lines="skip", engine="python")
    print(f"Read successfully: {df.shape[0]} rows, {df.shape[1]} columns")
    return df

def preprocess_data(df: pd.DataFrame, scaler, feature_names: List[str]):
    print("Start data preprocessing ...")
    contig = df["contig"].astype(str)
    position = df["position"]
    X = df[feature_names].copy()
    nan_mask = X.isna().any(axis=1)
    inf_mask = np.isinf(X.to_numpy()).any(axis=1)
    invalid_mask = nan_mask | inf_mask
    if invalid_mask.any():
        drop_n = int(invalid_mask.sum())
        print(f"[INFO] Remove {drop_n} lines containing NaN or inf/-inf.")
        X = X[~invalid_mask]
        contig = contig[~invalid_mask]
        position = position[~invalid_mask]
    if X.empty:
        raise ValueError("All samples are invalid and cannot be predicted")
    X_scaled = scaler.transform(X)
    print("Preprocessing completed")
    return X_scaled, contig.reset_index(drop=True), position.reset_index(drop=True)

def main():
    parser = argparse.ArgumentParser(description="M6A predictor")
    parser.add_argument("--test_file", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--preprocessed_train_data_file", default="preprocessed_train_data.pkl")
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--gpu_id", type=int, default=None)
    args = parser.parse_args()

    prep_path = os.path.join(args.output_dir, args.preprocessed_train_data_file)
    if not os.path.exists(prep_path):
        raise FileNotFoundError(f"Cannot find preprocessed training data file: {prep_path}")
    with open(prep_path, "rb") as f:
        train_info = pickle.load(f)
    scaler = train_info["scaler"]
    feature_names: List[str] = train_info["feature_names"]
    hidden_layers = train_info.get("hidden_layers", [256, 128])
    attention_hidden_dim = train_info.get("attention_hidden_dim", 32)
    dropout_rate = train_info.get("dropout_rate", 0.5)

    df_test = read_test_data(args.test_file)
    X_test, contigs, positions = preprocess_data(df_test, scaler, feature_names)

    test_loader = DataLoader(M6ADataset(X_test), batch_size=32, shuffle=False, num_workers=4)

    rf_path = os.path.join(args.output_dir, "random_forest.pkl")
    nn_path = os.path.join(args.output_dir, "best_model.pth")
    if not os.path.exists(rf_path) or not os.path.exists(nn_path):
        raise FileNotFoundError("Lacking random_forest.pkl or best_model.pth")
    rf_model = joblib.load(rf_path)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if args.gpu_id is not None and torch.cuda.is_available():
        if args.gpu_id < torch.cuda.device_count():
            device = torch.device(f"cuda:{args.gpu_id}")
        else:
            raise ValueError(f"GPU id {args.gpu_id} does not exist")

    nn_model = M6APredictor(
        hidden_layers=hidden_layers,
        input_dim=len(feature_names),
        num_classes=2,
        attention_hidden_dim=attention_hidden_dim,
        dropout_rate=dropout_rate,
    )
    nn_model.load_state_dict(torch.load(nn_path, map_location=device), strict=False)
    nn_model.to(device).eval()

    print("Random forest prediction ...")
    rf_probs = rf_model.predict_proba(X_test)[:, 1]

    print("Neural network prediction ...")
    nn_probs: List[float] = []
    with torch.no_grad():
        for batch in tqdm(test_loader, desc="NN"):
            batch = batch.to(device)
            nn_out = nn_model(batch)
            nn_probs.extend(torch.softmax(nn_out, 1)[:, 1].cpu().numpy())

    nn_probs = np.array(nn_probs)
    combined = (rf_probs + nn_probs) / 2

    out_df = pd.DataFrame(
        {
            "contig": contigs,
            "position": positions,
            "NeuralNetwork_Prob": nn_probs,
            "RandomForest_Prob": rf_probs,
            "Combined_Prob": combined,
        }
    )
    out_df.to_csv(args.output_file, sep="\t", index=False)
    print(f"The predicted results have been saved to {args.output_file}")

if __name__ == "__main__":
    main()
