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
from sklearn.metrics import roc_auc_score, average_precision_score

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
		hidden = self.hidden_layers(x)
		logits = self.output_layer(hidden)
		return logits

def read_test_data(file_path: str) -> pd.DataFrame:
	print(f"Start reading data file: {file_path}")
	# 使用 \s+ 支持空格或Tab分隔，与训练读取逻辑保持一致
	df = pd.read_csv(file_path, sep=r'\s+', header=0, on_bad_lines="skip", engine="python")
	print(f"Read successfully: {df.shape[0]} rows, {df.shape[1]} columns")
	return df

def get_device(args):
	if args.gpu_id and torch.cuda.is_available():
		gpu_str = str(args.gpu_id)
		device_ids = [int(x.strip()) for x in gpu_str.split(',') if x.strip().isdigit()]
		if len(device_ids) > 0:
			device = torch.device(f'cuda:{device_ids[0]}')
			return device, device_ids
	elif torch.cuda.is_available():
		return torch.device('cuda:0'), [0]
	return torch.device('cpu'), []

def preprocess_data(df: pd.DataFrame, scaler, feature_names: List[str], has_label: bool = False):
	print("Start data preprocessing ...")
	
	# 检查必要的标识列是否存在
	if "contig" in df.columns:
		contig = df["contig"].astype(str)
	else:
		contig = df.iloc[:, 0].astype(str)

	if "position" in df.columns:
		position = df["position"]
	else:
		position = df.iloc[:, 1]

	y = None
	if has_label and "label" in df.columns:
		y = df["label"].values
	elif has_label and df.shape[1] == len(feature_names) + 3:
		# contig + position + features + label
		y = df.iloc[:, -1].values

	# 检查特征列是否都能获取
	missing_cols = [c for c in feature_names if c not in df.columns]
	if missing_cols:
		# 若测试集列名与训练集不一致但维度匹配，按位置取列
		print("[Warning] Column names do not completely match training features. Slicing by index positions.")
		X = df.iloc[:, 2:2 + len(feature_names)].copy()
		X.columns = feature_names
	else:
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
		if y is not None:
			y = y[~invalid_mask]

	if X.empty:
		raise ValueError("All samples are invalid and cannot be predicted.")

	X_scaled = scaler.transform(X)
	print("Preprocessing completed.")

	if y is not None:
		return X_scaled, contig.reset_index(drop=True), position.reset_index(drop=True), pd.Series(y).reset_index(drop=True)
	else:
		return X_scaled, contig.reset_index(drop=True), position.reset_index(drop=True), None

def main():
	parser = argparse.ArgumentParser(description="nanoMS predictor supporting Kit 002 (56 features) and Kit 004 (88 features)")
	parser.add_argument('-k', '--kits', type=str, choices=['002', '004'], default='002', 
						help='Kits version: 002 or 004. If not set, will load from metadata file.  (default: 002)')
	parser.add_argument("--test_file", type=str, required=True, help='The input test/inference file path.')
	parser.add_argument("--model_dir", type=str, required=True, help='Output directory containing trained model and metadata.')
	parser.add_argument("--preprocessed_train_data_file", type=str, default="preprocessed_train_data.pkl", help='Inference metadata save/load file name.')
	parser.add_argument("--output_file", type=str, required=True, help='Path of output prediction results.')
	parser.add_argument("--batch_size", type=int, default=1024, help='Batch size for prediction. default: 1024')
	parser.add_argument("--gpu_id", type=str, default=None, help='Comma-separated list of GPU IDs to use (e.g., "0" or "0,1,2").')
	args = parser.parse_args()

	prep_path = os.path.join(args.model_dir, args.preprocessed_train_data_file)
	if not os.path.exists(prep_path):
		raise FileNotFoundError(f"Cannot find preprocessed training metadata file: {prep_path}")
	
	with open(prep_path, "rb") as f:
		train_info = pickle.load(f)
		
	# 优先以命令行参数为准，其次读取元数据配置，缺省为 002
	kits = args.kits if args.kits is not None else train_info.get("kits", "002")
	expected_input_dim = 88 if kits == '004' else 56
	print(f"Running inference kit mode: {kits}, expected input feature dimension: {expected_input_dim}")

	scaler = train_info["scaler"]
	feature_names: List[str] = train_info["feature_names"]
	hidden_layers = train_info.get("hidden_layers", [256, 128])
	dropout_rate = train_info.get("dropout_rate", 0.5)

	if len(feature_names) != expected_input_dim:
		print(f"[Warning] Feature count in metadata ({len(feature_names)}) differs from kit expectation ({expected_input_dim}). Using metadata dimension.")
		expected_input_dim = len(feature_names)

	df_test = read_test_data(args.test_file)
	
	# 检查是否包含标签列
	has_label = ("label" in df_test.columns) or (df_test.shape[1] == expected_input_dim + 3)

	X_test, contigs, positions, labels = preprocess_data(
		df_test, 
		scaler=scaler, 
		feature_names=feature_names, 
		has_label=has_label
	)
	
	test_loader = DataLoader(M6ADataset(X_test), batch_size=args.batch_size, shuffle=False, num_workers=4)

	nn_path = os.path.join(args.model_dir, "best_model.pth")
	if not os.path.exists(nn_path):
		raise FileNotFoundError(f"Model weights file not found: {nn_path}")

	device, device_ids = get_device(args)
	print(f"Using device: {device}, device_ids: {device_ids}")

	# 动态传入 input_dim (适配 002/004 模型的对应维度)
	nn_model = mlp(
		hidden_layers=hidden_layers,
		input_dim=expected_input_dim,
		num_classes=2,
		dropout_rate=dropout_rate
	).to(device)

	state_dict = torch.load(nn_path, map_location=device)
	# 移除可能的 module. 前缀
	if any(k.startswith('module.') for k in state_dict.keys()):
		state_dict = {k.replace('module.', ''): v for k, v in state_dict.items()}

	nn_model.load_state_dict(state_dict, strict=True) 
	
	if device_ids and len(device_ids) > 1:
		nn_model = torch.nn.DataParallel(nn_model, device_ids=device_ids)
		
	nn_model.eval()

	print("Neural network predicting ...")
	nn_probs: List[float] = []
	with torch.no_grad():
		for batch in tqdm(test_loader, desc="NN Inference"):
			batch = batch.to(device)
			nn_out = nn_model(batch)
			probs = torch.softmax(nn_out, dim=1)[:, 1]
			nn_probs.extend(probs.cpu().numpy())

	nn_probs = np.array(nn_probs)
	
	# 输出结果整理
	out_dict = {
		"contig": contigs,
		"position": positions,
		"NeuralNetwork_Prob": nn_probs
	}

	if labels is not None:
		out_dict["label"] = labels
		# 如果标签有正负两类，计算评估指标
		unique_labels = np.unique(labels)
		if len(unique_labels) > 1:
			try:
				auc = roc_auc_score(labels, nn_probs)
				pr_auc = average_precision_score(labels, nn_probs)
				print(f"[Evaluation] Test Set AUC: {auc:.4f}, PR AUC: {pr_auc:.4f}")
			except Exception as e:
				print(f"[Warning] Failed to calculate metrics: {e}")

	out_df = pd.DataFrame(out_dict)
	
	# 保证输出目录存在
	out_dir = os.path.dirname(args.output_file)
	if out_dir:
		os.makedirs(out_dir, exist_ok=True)
		
	out_df.to_csv(args.output_file, sep="\t", index=False)
	print(f"Prediction results successfully saved to {args.output_file}")

if __name__ == "__main__":
	main()