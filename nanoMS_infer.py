import argparse
import os
import pickle
import warnings
from typing import List

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch import nn
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm
from sklearn.metrics import roc_auc_score, classification_report, confusion_matrix, accuracy_score, precision_recall_curve, average_precision_score

warnings.filterwarnings("ignore")

class M6ADataset(Dataset):
	def __init__(self, X: np.ndarray):
		self.X = torch.tensor(X, dtype=torch.float32)
	def __len__(self):
		return len(self.X)
	def __getitem__(self, idx):
		return self.X[idx]

class mlp(nn.Module):
	def __init__(self, 
				 hidden_layers, 
				 input_dim=56, 
				 num_classes=2,
				 dropout_rate=0.3):
		super(mlp, self).__init__()
		# 验证输入维度可被num_heads整除
		# # 全连接层保持不变
		layers = []
		prev_dim = input_dim
		for hidden_dim in hidden_layers:
			layers.append(nn.Linear(prev_dim, hidden_dim))
			layers.append(nn.ReLU())
			layers.append(nn.Dropout(dropout_rate))
			prev_dim = hidden_dim
		self.hidden_layers = nn.Sequential(*layers)
		self.output_layer = nn.Linear(prev_dim, num_classes)
		
	def forward(self, x):
		hidden = self.hidden_layers(x)  # (batch_size, hidden_layers[-1])
		logits = self.output_layer(hidden)  # (batch_size, num_classes)
		return logits

def read_test_data(file_path: str) -> pd.DataFrame:
	print(f"Start reading data files: {file_path}")
	df = pd.read_csv(file_path, sep="\t", header=0, on_bad_lines="skip", engine="python")
	print(f"Read successfully: {df.shape[0]} rows, {df.shape[1]} columns")
	return df

def get_device(args):
	if args.gpu_id and torch.cuda.is_available():
		device_ids = [int(x) for x in args.gpu_id.split(',')]
		device = torch.device(f'cuda:{device_ids[0]}')  # You can modify this logic based on how you want to use multiple GPUs
		return device, device_ids
	else:
		return torch.device('cpu'), []

def preprocess_data(df: pd.DataFrame, scaler, feature_names: List[str], label=None):
	print("Start data preprocessing ...")
	contig = df["contig"].astype(str)
	if label:
		add_label = df["label"].astype(str)
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
		add_label = add_label[~invalid_mask]
	if X.empty:
		raise ValueError("All samples are invalid and cannot be predicted")
	X_scaled = scaler.transform(X)
	print("Preprocessing completed")
	if not label:
		return X_scaled, contig.reset_index(drop=True), position.reset_index(drop=True)
	else:
		return X_scaled, contig.reset_index(drop=True), position.reset_index(drop=True), add_label.reset_index(drop=True)

def main():
	parser = argparse.ArgumentParser(description="nanoMS predictor")
	parser.add_argument("--test_file", type=str, required=True, help='The input file path')
	parser.add_argument("--model_dir", type=str, required=True, help='Output directory for nanoMS training')
	parser.add_argument("--preprocessed_train_data_file", type=str, default="preprocessed_train_data.pkl", help='Pre processed training data save/load file name')
	parser.add_argument("--output_file", type=str, required=True, help='Path of output results')
	parser.add_argument("--gpu_id", type=str, default=None, help='Comma-separated list of GPU IDs to use (e.g., "0,1,2"). If not specified, CPU will be used.')
	args = parser.parse_args()

	prep_path = os.path.join(args.model_dir, args.preprocessed_train_data_file)
	if not os.path.exists(prep_path):
		raise FileNotFoundError(f"Cannot find preprocessed training data file: {prep_path}")
	with open(prep_path, "rb") as f:
		train_info = pickle.load(f)
	scaler = train_info["scaler"]
	feature_names: List[str] = train_info["feature_names"]
	hidden_layers = train_info.get("hidden_layers", [256, 128])
	dropout_rate = train_info.get("dropout_rate", 0.5)

	df_test = read_test_data(args.test_file)
	X_test, contigs, positions = preprocess_data(df_test, scaler, feature_names)
	# X_test, contigs, positions, label = preprocess_data(df_test, scaler, feature_names, label=True)
	test_loader = DataLoader(M6ADataset(X_test), batch_size=1024, shuffle=False, num_workers=4)

	nn_path = os.path.join(args.model_dir, "best_model.pth")
	device, device_ids = get_device(args)
	
	nn_model = mlp(
		hidden_layers=hidden_layers,
		input_dim=len(feature_names),  # 56 numerical features
		num_classes=2,
		dropout_rate=dropout_rate
	).to(device)

	state_dict = torch.load(nn_path, map_location=device)
	
	if list(state_dict.keys())[0].startswith('module.'):
		state_dict = {k.replace('module.', ''): v for k, v in state_dict.items()}

	nn_model.load_state_dict(state_dict, strict=True) 
	
	if device_ids and len(device_ids) > 1:
		nn_model = torch.nn.DataParallel(nn_model, device_ids=device_ids)
		
	nn_model.eval()

	print("Neural network prediction ...")
	nn_probs: List[float] = []
	with torch.no_grad():
		for batch in tqdm(test_loader, desc="NN"):
			batch = batch.to(device)
			nn_out = nn_model(batch)
			nn_probs.extend(torch.softmax(nn_out, 1)[:, 1].cpu().numpy())

	nn_probs = np.array(nn_probs)
	
	
	if 'label' in locals() and len(list(set(label.values))) > 1:
		val_auc = roc_auc_score(label, nn_probs)
		print(f"AUC:{val_auc}")
		print(f"The predicted results have been saved to {args.output_file}")
		out_df = pd.DataFrame(
			{
				"contig": contigs,
				"position": positions,
				"label":label,
				"NeuralNetwork_Prob": nn_probs
			}
		)
	else:
		out_df = pd.DataFrame(
			{
				"contig": contigs,
				"position": positions,
				"NeuralNetwork_Prob": nn_probs
			}
		)
	out_df.to_csv(args.output_file, sep="\t", index=False)
	

if __name__ == "__main__":
	main()
