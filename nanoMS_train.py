import argparse
import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, classification_report, confusion_matrix, accuracy_score, precision_recall_curve, average_precision_score
from sklearn.utils.class_weight import compute_class_weight
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader
import torch.nn.functional as F
from tqdm import tqdm
import warnings
import pickle

warnings.filterwarnings("ignore")

class FocalLoss(nn.Module):
	def __init__(self, alpha=1, gamma=2, reduction='mean'):
		"""
		:param alpha: Weight factor used to balance the loss of different classes
		:param gamma: Modulating factor focusing on hard-to-classify samples
		:param reduction: Loss aggregation method, 'none' | 'mean' | 'sum'
		"""
		super(FocalLoss, self).__init__()
		self.alpha = alpha
		self.gamma = gamma
		self.reduction = reduction

	def forward(self, inputs, targets):
		# Calculate cross-entropy loss
		ce_loss = F.cross_entropy(inputs, targets, reduction='none')
		# Calculate predicted probabilities
		pt = torch.exp(-ce_loss)
		# Calculate focal loss
		focal_loss = self.alpha * (1 - pt) ** self.gamma * ce_loss

		if self.reduction == 'mean':
			return focal_loss.mean()
		elif self.reduction == 'sum':
			return focal_loss.sum()
		else:
			return focal_loss

# 1. Data Preprocessing

def read_data(file_path, ncol):
	"""
	Read training or test data files and convert them into dataframes.
	Each row contains 56 numerical features and 1 label, totaling 57 columns.
	Automatically skip rows that do not contain exactly 57 columns.
	"""
	print(f"Start reading data file: {file_path}")
	
	try:
		# Assume data is space or tab separated, skip rows with inconsistent column counts
		data = pd.read_csv(file_path, sep='\s+', header=0, on_bad_lines='skip', engine='python')
	except Exception as e:
		print(f"Error reading file: {e}")
		return pd.DataFrame()
	
	expected_num_columns = ncol  # contig + position + 56 numerical features + 1 label
	actual_num_columns = data.shape[1]
	
	# Calculate skipped rows
	total_lines = sum(1 for _ in open(file_path)) - 1  # Subtract header row
	skipped_lines = total_lines - data.shape[0]
	
	if actual_num_columns != expected_num_columns:
		print(f"Warning: The number of columns in the data file does not match the expected. Expected {expected_num_columns} columns, but found {actual_num_columns} columns.")
	
	if skipped_lines > 0:
		print(f"Note: Skipped {skipped_lines} rows because their column count is not {expected_num_columns}.")
	
	if data.empty:
		print("Data is empty. Please check the data file format.")
	else:
		print(f"Successfully read data file. Number of samples: {data.shape[0]}, Number of columns: {data.shape[1]}")
	return data

def preprocess_data(df, scaler=None, fit=True):
	"""
	Standardize the data and drop rows containing NaN or inf/-inf.
	"""
	print("Starting data preprocessing...")

	# Extract numerical features and labels
	numerical_features = df.columns[:-1]  # Extract the first 58 feature columns
	# Select numerical features
	X_numerical = df[numerical_features].copy()

	# Drop rows containing NaN or inf/-inf
	print("Dropping rows with NaN or inf/-inf in numerical features...")
	initial_count = X_numerical.shape[0]
	# First drop NaN
	X_numerical = X_numerical.dropna(subset=numerical_features)
	# Then drop inf/-inf
	X_numerical = X_numerical[~np.isinf(X_numerical).any(axis=1)]
	final_count = X_numerical.shape[0]
	removed_count = initial_count - final_count
	if removed_count > 0:
		print(f"Dropped {removed_count} samples containing NaN or inf/-inf.")
	else:
		print("No samples found containing NaN or inf/-inf.")

	df = df.loc[X_numerical.index].reset_index(drop=True)
	X_numerical = X_numerical.reset_index(drop=True)

   
	if X_numerical.shape[0] == 0:
		raise ValueError("All samples contain NaN or inf/-inf values, unable to train or test.")

	# Standardize numerical features
	if fit:
		print("Fitting StandardScaler...")
		scaler = StandardScaler().fit(X_numerical)
		print("StandardScaler fitting completed.")
	print("Standardizing numerical features...")
	X_numerical = scaler.transform(X_numerical)
	print("Numerical features standardization completed.")

	
	X = X_numerical

   
	y = df['label'].values

	print("Data preprocessing completed.")

	if fit:
		return {
			'X': X,
			'y': y,
			'scaler': scaler,
			'feature_names': numerical_features
		}
	else:
		return {
			'X': X,
			'y': y,
			'feature_names': numerical_features
		}

class M6ADataset(Dataset):
	def __init__(self, X, y):
		self.X = torch.tensor(X, dtype=torch.float)
		self.y = torch.tensor(y, dtype=torch.long)

	def __len__(self):
		return len(self.y)

	def __getitem__(self, idx):
		return self.X[idx], self.y[idx]

class mlp(nn.Module):
	def __init__(self, 
				 hidden_layers, 
				 input_dim=56, 
				 num_classes=2,
				 dropout_rate=0.3):
		super(mlp, self).__init__()
		# Validate input dimensions
		# Fully connected layers remain unchanged
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

def get_device(args):
	if args.gpu_id and torch.cuda.is_available():
		device_ids = [int(x) for x in args.gpu_id.split(',')]
		device = torch.device(f'cuda:{device_ids[0]}')  # You can modify this logic based on how you want to use multiple GPUs
		return device, device_ids
	else:
		return torch.device('cpu'), []

def train_and_evaluate(args, data_train, data_test):
	print("Creating datasets and dataloaders...")
	train_X, train_y = data_train['X'], data_train['y']
	test_X, test_y = data_test['X'], data_test['y']

	# Print the number of training and testing samples
	print(f"Number of training samples: {train_X.shape[0]}")
	print(f"Number of testing samples: {test_X.shape[0]}")

	train_dataset = M6ADataset(train_X, train_y)
	test_dataset = M6ADataset(test_X, test_y)

	print("Creating DataLoaders...")
	train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True, num_workers=4)
	test_loader = DataLoader(test_dataset, batch_size=args.batch_size, shuffle=False, num_workers=4)

	nn_model = mlp(
		hidden_layers=args.hidden_layers,
		input_dim=train_X.shape[1],  # 56 numerical features
		num_classes=2,
		dropout_rate=args.dropout_rate
	)
	print("Neural network model instantiated.")

	device, device_ids = get_device(args)
	print(f"Using device: {device}, device_ids: {device_ids}")

	nn_model.to(device)

	if len(device_ids) > 1:
		print(f"Enabling DataParallel using GPUs: {device_ids}")
		nn_model = nn.DataParallel(nn_model, device_ids=device_ids)


	class_weights = compute_class_weight(
		class_weight='balanced',
		classes=np.unique(train_y),
		y=train_y
	)
	class_weights = torch.tensor(class_weights, dtype=torch.float).to(device)
	
	# Use Focal Loss
	criterion = FocalLoss(alpha=args.focal_alpha, gamma=args.focal_gamma)
	optimizer = torch.optim.AdamW(
		nn_model.parameters(),
		lr=args.learning_rate,
		weight_decay=args.l2  # Use args.l2 as weight_decay
	)
	print("Loss function and optimizer defined.")

	scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
		optimizer,
		mode='min',
		factor=0.5,
		patience=5,
		verbose=True
	)
	print("Learning rate scheduler set up.")

	best_val_loss = float('inf')
	best_epoch = 0
	patience_counter = 0
	best_val_auc = 0
	
	history = {
		'train_loss': [],
		'train_auc': [],
		'train_pr_auc': [],  # PR AUC for training
		'val_loss': [],
		'val_auc': [],
		'val_pr_auc': []     # PR AUC for validation
	}
	
	train_res_file = open(os.path.join(args.output_dir, "train_res.txt"), "w")
	print("Starting model training...")
	for epoch in range(1, args.epochs + 1):
		print(f"\n=== Epoch {epoch}/{args.epochs} ===")
		nn_model.train()
		running_loss = 0.0
		all_preds = []
		all_labels = []
		all_probs = []

		print("Starting training phase...")
		for batch_idx, (X_batch, y_batch) in enumerate(tqdm(train_loader, desc="Training", leave=False)):
			X_batch = X_batch.to(device)
			y_batch = y_batch.to(device)

			optimizer.zero_grad()
			outputs = nn_model(X_batch)
			loss = criterion(outputs, y_batch)
			loss.backward()
			optimizer.step()

			running_loss += loss.item() * X_batch.size(0)
			probs = F.softmax(outputs, dim=1)[:, 1]
			all_probs.extend(probs.detach().cpu().numpy())
			_, predicted = torch.max(outputs, 1)
			all_preds.extend(predicted.detach().cpu().numpy())
			all_labels.extend(y_batch.detach().cpu().numpy())

		train_loss = running_loss / len(train_dataset)

		# Calculate PR AUC
		train_precision, train_recall, _ = precision_recall_curve(all_labels, all_probs)
		train_pr_auc = average_precision_score(all_labels, all_probs)

		try:
			train_auc = roc_auc_score(all_labels, all_probs)
		except ValueError as ve:
			print(f"Error calculating training AUC: {ve}")
			train_auc = 0.0

		history['train_loss'].append(train_loss)
		history['train_auc'].append(train_auc)
		history['train_pr_auc'].append(train_pr_auc)

		print(f"Training set - Loss: {train_loss:.4f}, AUC: {train_auc:.4f}, PR AUC: {train_pr_auc:.4f}")

		train_res_file.write(f"Epoch {epoch} Training set - Loss: {train_loss:.4f}, AUC: {train_auc:.4f}, PR AUC: {train_pr_auc:.4f}\n")
		print("Starting validation phase...")
		nn_model.eval()
		val_loss = 0.0
		val_preds = []
		val_labels_list = []
		val_probs = []
		with torch.no_grad():
			for X_batch, y_batch in tqdm(test_loader, desc="Validation", leave=False):
				X_batch = X_batch.to(device)
				y_batch = y_batch.to(device)

				outputs = nn_model(X_batch)
				loss = criterion(outputs, y_batch)
				val_loss += loss.item() * X_batch.size(0)

				probs = F.softmax(outputs, dim=1)[:, 1]
				val_probs.extend(probs.detach().cpu().numpy())
				_, predicted = torch.max(outputs, 1)
				val_preds.extend(predicted.detach().cpu().numpy())
				val_labels_list.extend(y_batch.detach().cpu().numpy())

		val_loss /= len(test_dataset)

		
		val_precision, val_recall, _ = precision_recall_curve(val_labels_list, val_probs)
		val_pr_auc = average_precision_score(val_labels_list, val_probs)

		try:
			val_auc = roc_auc_score(val_labels_list, val_probs)
		except ValueError as ve:
			print(f"Error calculating validation AUC: {ve}")
			val_auc = 0.0

		history['val_loss'].append(val_loss)
		history['val_auc'].append(val_auc)
		history['val_pr_auc'].append(val_pr_auc)

		print(f"Validation set - Loss: {val_loss:.4f}, AUC: {val_auc:.4f}, PR AUC: {val_pr_auc:.4f}")
		train_res_file.write(f"Epoch {epoch} Validation set - Loss: {val_loss:.4f}, AUC: {val_auc:.4f}, PR AUC: {val_pr_auc:.4f}\n")

		scheduler.step(val_loss)
		print("Learning rate scheduler updated.")

		if val_auc > best_val_auc:
			best_val_loss = val_loss
			best_val_prc = val_pr_auc
			best_val_auc = val_auc
			best_epoch = epoch
			patience_counter = 0
			torch.save(nn_model.state_dict(), os.path.join(args.output_dir, 'best_model.pth'))
			print(f"Save the current best model at epoch {epoch}.")
			test_res = {"val_pred":val_preds, "val_label":val_labels_list, "val_probs":val_probs}
			test_res_df = pd.DataFrame(test_res)
			test_res_df.to_csv(os.path.join(args.output_dir, "val_pred_res.txt"), sep="\t",index=False)
		else:
			patience_counter += 1
			print(f"The verification loss has not improved. Patience counter: {patience_counter}/{args.patience}")
			if patience_counter >= args.patience:
				print("Trigger early stop. Stop training.")
				break
		if args.save_epoch_model:
			checkpoint_path = os.path.join(args.output_dir, f'model_epoch_{epoch}.pth')
			torch.save(nn_model.state_dict(), checkpoint_path)
			print(f"Saved model checkpoint for epoch {epoch}.")

		print(f"\nTraining completed. The best model is in the {best_epoch} epoch, the validation loss: {best_val_loss:.4f}, AUC: {best_val_auc:.4f}, PR AUC: {best_val_prc:.4f}")
	
	# history_df = pd.DataFrame(history)
	train_res_file.close()
	

def parse_arguments():
	parser = argparse.ArgumentParser(description="Training a Neural Network to Predict RNA Structural Modifications and m6A Sites from Nanopore Sequencing Data")
	parser.add_argument('--train_file', type=str, required=True, help='Training data file path.')
	parser.add_argument('--valid_ratio', type=float, default=0.2, help='Validation set ratio. default (0.2)')
	parser.add_argument('--valid_file', type=str, default=None, help='Validation data file path. If no file is provided, it will be split from the training file according to the ratio set by --valid_ratio.')
	parser.add_argument('--output_dir', type=str, required=True, help='Directory to save model and results.')
	parser.add_argument('--ncol', type=int, default=59, help='Number of data columns')

	parser.add_argument('--epochs', type=int, default=200, help='Number of training epochs.')
	parser.add_argument('--batch_size', type=int, default=32, help='Training batch size.')
	parser.add_argument('--learning_rate', type=float, default=0.0001, help='Optimizer learning rate.')
	parser.add_argument('--weight_decay', type=float, default=1e-4, help='Weight decay (L2 regularization).')  # Deprecated, use --l2
	parser.add_argument('--patience', type=int, default=20, help='Patience counter for early stopping.')
	parser.add_argument('--seed', type=int, default=42, help='Random seed.')
	parser.add_argument('--dropout_rate', type=float, default=0.5, help='Dropout rate for regularization.')
	parser.add_argument('--focal_alpha', type=float, default=1.0, help='Alpha parameter for focal loss.')
	parser.add_argument('--focal_gamma', type=float, default=2.0, help='Gamma parameter for focal loss.')
	parser.add_argument("--gpu_id", type=str, default=None, help='Comma-separated list of GPU IDs to use (e.g., "0,1,2"). If not specified, CPU will be used.')
	parser.add_argument('--preprocessed_train_data_file', type=str, default='preprocessed_train_data.pkl', help='Path to save/load preprocessed training data.')
	parser.add_argument('--preprocessed_test_data_file', type=str, default='preprocessed_test_data.pkl', help='Path to save/load preprocessed testing data.')
	parser.add_argument('--save_epoch_model', action='store_true', help='Save the model for each epoch')

 
	parser.add_argument('--l2', type=float, default=1e-4, help='L2 regularization (weight_decay).')
	parser.add_argument('--hidden_layers', type=int, nargs='+', default=[256, 128], help='List of neural network hidden layer sizes. e.g., --hidden_layers 256 128')
	parser.add_argument('--use_preprocessed', action='store_true', help='Whether to use preprocessed data files.')
	return parser.parse_args()

def main():
	args = parse_arguments()
	os.makedirs(args.output_dir, exist_ok=True)
	print(f"Output directory created or already exists: {args.output_dir}")

	torch.manual_seed(args.seed)
	np.random.seed(args.seed)
	if torch.cuda.is_available():
		torch.cuda.manual_seed_all(args.seed)
	print(f"Random seed set to: {args.seed}")

	preprocessed_train_data_path = os.path.join(args.output_dir, args.preprocessed_train_data_file)
	preprocessed_test_data_path = os.path.join(args.output_dir, args.preprocessed_test_data_file)

	if args.use_preprocessed:
		if os.path.exists(preprocessed_train_data_path) and os.path.exists(preprocessed_test_data_path):
			print(f"Loading preprocessed training data from {preprocessed_train_data_path}...")
			with open(preprocessed_train_data_path, 'rb') as f:
				preprocessed_train = pickle.load(f)
			print("Preprocessed training data loaded successfully.")

			print(f"Loading preprocessed testing data from {preprocessed_test_data_path}...")
			with open(preprocessed_test_data_path, 'rb') as f:
				preprocessed_test = pickle.load(f)
			print("Preprocessed testing data loaded successfully.")
		else:
			print("Preprocessed data files not found. Please ensure they exist, or do not use the --use_preprocessed option.")
			return
	else:
		
		print("Reading training data...")
		train_df = read_data(args.train_file, args.ncol)
		if train_df.shape[1] != args.ncol:
			print(f"Training data column count mismatch. Expected {args.ncol} columns (contig + position + 56 numerical features + 1 label), but found {train_df.shape[1]} columns.")
			return
		if train_df.empty:
			print("No valid data read from training file. Please check data file format.")
			return

		if args.valid_file:
			print("Reading testing data...")
			test_df = read_data(args.valid_file, args.ncol)
			if test_df.shape[1] != args.ncol:
				print(f"Testing data column count mismatch. Expected {args.ncol} columns (contig + position + 56 numerical features + 1 label), but found {test_df.shape[1]} columns.")
				return
			if test_df.empty:
				print("No valid data read from testing file. Please check data file format.")
				return
			if train_df.shape[1] != 57:
				test_df = test_df.iloc[:,2:]
				train_df = train_df.iloc[:,2:]
		else:
			if train_df.shape[1] == 57:
				test_df = train_df.iloc[:,:].sample(frac=args.valid_ratio)
				train_df = train_df.iloc[:,:].drop(test_df.index)
			else:
				test_df = train_df.iloc[:,2:].sample(frac=args.valid_ratio)
				train_df = train_df.iloc[:,2:].drop(test_df.index)

		preprocessed_train = preprocess_data(train_df, fit=True)
		
		preprocessed_test = preprocess_data(
			test_df,
			scaler=preprocessed_train['scaler'],
			fit=False
		)

		
		with tqdm(total=2, desc="Saving preprocessed data") as pbar:
			with open(preprocessed_train_data_path, 'wb') as f:
				pickle.dump(preprocessed_train, f)
			pbar.update(1)
			with open(preprocessed_test_data_path, 'wb') as f:
				pickle.dump(preprocessed_test, f)
			pbar.update(1)
		print(f"Preprocessed training data saved to {preprocessed_train_data_path}.")
		print(f"Preprocessed testing data saved to {preprocessed_test_data_path}.")

	train_and_evaluate(args, preprocessed_train, preprocessed_test)

if __name__ == "__main__":
	main()