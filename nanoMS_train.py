
import argparse
import os
import pandas as pd
import numpy as np
import re
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, classification_report, confusion_matrix, accuracy_score, precision_recall_curve, average_precision_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils.class_weight import compute_class_weight
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader
import torch.nn.functional as F
from tqdm import tqdm
import warnings
import pickle
import joblib
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


class FocalLoss(nn.Module):
    def __init__(self, alpha=1, gamma=2, reduction='mean'):
        """
        :param alpha: 权重因子，用于平衡不同类别的损失
        :param gamma: 调制因子，聚焦于难分类样本
        :param reduction: 损失的聚合方式，'none' | 'mean' | 'sum'
        """
        super(FocalLoss, self).__init__()
        self.alpha = alpha
        self.gamma = gamma
        self.reduction = reduction

    def forward(self, inputs, targets):
        # 计算交叉熵损失
        ce_loss = F.cross_entropy(inputs, targets, reduction='none')
        # 计算预测概率
        pt = torch.exp(-ce_loss)
        # 计算焦点损失
        focal_loss = self.alpha * (1 - pt) ** self.gamma * ce_loss

        if self.reduction == 'mean':
            return focal_loss.mean()
        elif self.reduction == 'sum':
            return focal_loss.sum()
        else:
            return focal_loss

# 1. 数据预处理

def read_data(file_path):
    """
    读取训练或测试数据文件，并将其转换为数据框。
    每行包含contig,position, 56个数值特征和1个标签，总计57列。
    自动跳过不包含正好57列的行。
    """
    print(f"开始读取数据文件: {file_path}")
    
    try:
        # 假设数据以空格或制表符分隔，跳过不符合列数的行
        data = pd.read_csv(file_path, sep='\s+', header=0, on_bad_lines='skip', engine='python')
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return pd.DataFrame()
    
    expected_num_columns = 59  # contig + position + 56features + 1label
    actual_num_columns = data.shape[1]
    
    # 计算跳过的行数
    total_lines = sum(1 for _ in open(file_path)) - 1  # 减去标题行
    skipped_lines = total_lines - data.shape[0]
    
    if actual_num_columns != expected_num_columns:
        print(f"Warning: The number of columns in the data file does not match the expected number. Expected {expected_num_columns} column, but found {actual_num_columns} columns.")
    
    if skipped_lines > 0:
        print(f"Note: {skipped_lines} rows has been skipped because the number of columns in these rows is not {expected_num_columns}.")
    
    if data.empty:
        print("The data is empty. Please check the format of the data file.")
    else:
        print(f"Successfully read data file. Sample size: {data.shape[0]}, Columns:{data.shape[1]}")
    return data

def preprocess_data(df, scaler=None, fit=True):
    """
    对数据进行标准化，删除包含 NaN 或 inf/-inf 的行。
    """
    print("Start data preprocessing...")

    # 提取数值特征和标签
    numerical_features = df.columns[:-1]  # 提取前56列特征
    # 选择数值特征
    X_numerical = df[numerical_features].copy()

    # 删除包含 NaN 或 inf/-inf 的行
    print("Delete rows containing NaN or inf/- nf in numerical features...")
    initial_count = X_numerical.shape[0]
    # 先删除 NaN
    X_numerical = X_numerical.dropna(subset=numerical_features)
    # 再删除 inf/-inf
    X_numerical = X_numerical[~np.isinf(X_numerical).any(axis=1)]
    final_count = X_numerical.shape[0]
    removed_count = initial_count - final_count
    if removed_count > 0:
        print(f"Remov {removed_count} lines containing NaN or inf/-inf.")
    else:
        print("No samples containing NaN or nf/- nf were found.")

    df = df.loc[X_numerical.index].reset_index(drop=True)
    X_numerical = X_numerical.reset_index(drop=True)

   
    if X_numerical.shape[0] == 0:
        raise ValueError("All samples contain NaN or nf/- nf values and cannot be trained or tested.")

    # 标准化数值特征
    if fit:
        print("Start fitting StandardScaler...")
        scaler = StandardScaler().fit(X_numerical)
        print("StandardScaler fitting completed.")
    print("Start standardizing numerical features...")
    X_numerical = scaler.transform(X_numerical)
    print("Standardization of numerical features completed.")

    
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

class Attention(nn.Module):
    def __init__(self, input_dim, hidden_dim):
        super(Attention, self).__init__()
        self.attention = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, 1, bias=False)
        )

    def forward(self, x):
        # x: (batch_size, input_dim)
        scores = self.attention(x)  # (batch_size, 1)
        weights = F.softmax(scores, dim=1)  # (batch_size, 1)
        output = x * weights  # (batch_size, input_dim)
        return output

class M6APredictor(nn.Module):
    def __init__(self,
                 hidden_layers,  # List of hidden layer sizes
                 input_dim=56,  # 56 numerical features
                 num_classes=2,
                 attention_hidden_dim=32,
                 dropout_rate=0.5):
        super(M6APredictor, self).__init__()

        # 注意力机制
        self.attention = Attention(input_dim, attention_hidden_dim)

        # 全连接层
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
        # x: (batch_size, input_dim)
        x = self.attention(x)  # (batch_size, input_dim)
        hidden = self.hidden_layers(x)  # (batch_size, hidden_layers[-1])
        logits = self.output_layer(hidden)  # (batch_size, num_classes)
        return logits

def train_and_evaluate(args, data_train, data_test):
    print("Create datasets and data loaders...")
    train_X, train_y = data_train['X'], data_train['y']
    test_X, test_y = data_test['X'], data_test['y']

    # 打印训练和测试样本数量
    print(f"Number of training samples: {train_X.shape[0]}")
    print(f"Number of valid samples: {test_X.shape[0]}")

    train_dataset = M6ADataset(train_X, train_y)
    test_dataset = M6ADataset(test_X, test_y)

    print("Create DataLoaders...")
    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True, num_workers=4)
    test_loader = DataLoader(test_dataset, batch_size=args.batch_size, shuffle=False, num_workers=4)

    print("Training a Random Forest Model...")
    rf_model_path = os.path.join(args.output_dir, 'random_forest.pkl')
    if os.path.exists(rf_model_path):
        print(f"Load the existing random forest model from{rf_model_path}...")
        rf_model = joblib.load(rf_model_path)
        print("The random forest model has been loaded.")
    else:
        rf_model = RandomForestClassifier(
            n_estimators=args.rf_n_estimators,
            max_depth=args.rf_max_depth,
            min_samples_split=args.rf_min_samples_split,
            min_samples_leaf=args.rf_min_samples_leaf,
            class_weight=args.rf_class_weight,
            random_state=args.seed,
            n_jobs=-1
        )
        with tqdm(total=1, desc="Training a Random Forest Model") as pbar:
            rf_model.fit(train_X, train_y)
            pbar.update(1)
        joblib.dump(rf_model, rf_model_path)
        print(f"The random forest model has been trained and saved to {rf_model_path}.")

    nn_model = M6APredictor(
        hidden_layers=args.hidden_layers,
        input_dim=train_X.shape[1],  # 56 numerical features
        num_classes=2,
        attention_hidden_dim=args.attention_hidden_dim,
        dropout_rate=args.dropout_rate
    )
    print("The neural network model has been instantiated.")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if args.gpu_id is not None:
        if torch.cuda.is_available():
            if args.gpu_id < torch.cuda.device_count():
                device = torch.device(f"cuda:{args.gpu_id}")
                print(f"Useing GPU: {device}")
            else:
                raise ValueError(f"The specified {args.gpu_id}  exceeds the number of available GPUs ({torch.cuda.device_count()})。")
        else:
            raise ValueError("CUDA is not available.")
    else:
        if torch.cuda.is_available():
            print(f"Use default GPU: {device}")
        else:
            print("CUDA is not available, Using CPU.")
    nn_model.to(device)


    class_weights = compute_class_weight(
        class_weight='balanced',
        classes=np.unique(train_y),
        y=train_y
    )
    class_weights = torch.tensor(class_weights, dtype=torch.float).to(device)
    
    # 使用焦点损失
    criterion = FocalLoss(alpha=args.focal_alpha, gamma=args.focal_gamma)
    optimizer = torch.optim.AdamW(
        nn_model.parameters(),
        lr=args.learning_rate,
        weight_decay=args.l2  # 使用 args.l2 作为 weight_decay
    )
    print("The loss function and optimizer have been defined.")

    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode='min',
        factor=0.5,
        patience=5,
        verbose=True
    )
    print("The learning rate scheduler has been set.")

    best_val_loss = float('inf')
    best_epoch = 0
    patience_counter = 0

    history = {
        'train_loss': [],
        'train_auc': [],
        'train_pr_auc': [],  # PR AUC for training
        'val_loss': [],
        'val_auc': [],
        'val_pr_auc': []     # PR AUC for validation
    }

    print("Start model training ..")
    for epoch in range(1, args.epochs + 1):
        print(f"\n=== Epoch {epoch}/{args.epochs} ===")
        nn_model.train()
        running_loss = 0.0
        all_preds = []
        all_labels = []
        all_probs = []

        print("Start the training phase ..")
        for batch_idx, (X_batch, y_batch) in enumerate(tqdm(train_loader, desc="training", leave=False)):
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

        # 计算 PR AUC
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

        print("Start validation phase...")
        nn_model.eval()
        val_loss = 0.0
        val_preds = []
        val_labels_list = []
        val_probs = []
        with torch.no_grad():
            for X_batch, y_batch in tqdm(test_loader, desc="validation", leave=False):
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

        scheduler.step(val_loss)
        print("The learning rate scheduler has been updated.")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_epoch = epoch
            patience_counter = 0
            torch.save(nn_model.state_dict(), os.path.join(args.output_dir, 'best_model.pth'))
            print(f"Save the current best model at epoch {epoch}.")
        else:
            patience_counter += 1
            print(f"The verification loss has not improved. Patience counter:: {patience_counter}/{args.patience}")
            if patience_counter >= args.patience:
                print("Trigger early stop. Stop training.")
                break

        if args.save_epoch_model:
            checkpoint_path = os.path.join(args.output_dir, f'model_epoch_{epoch}.pth')
            torch.save(nn_model.state_dict(), checkpoint_path)
            print(f"The {epoch} epoch of model checkpoints has been saved")

    print(f"\nTraining completed. The best model is in the {best_epoch} epoch, the validation loss is {best_val_loss:.4f}.")
    np.save(os.path.join(args.output_dir, "history.npy"), history)

def parse_arguments():
    parser = argparse.ArgumentParser(description="RNA m6A Binding Site Predictor with Neural Network and Random Forest")
    parser.add_argument('--train_file', type=str, required=True, help='Train data file path')
    parser.add_argument('--valid_ratio', type=float, default=0.2, help='Validation set ratio')
    parser.add_argument('--output_dir', type=str, required=True, help='Save the directory of models and results.')
    parser.add_argument('--epochs', type=int, default=200, help='Number of training epochs')
    parser.add_argument('--batch_size', type=int, default=32, help='Training batch size')
    parser.add_argument('--learning_rate', type=float, default=0.0001, help='learning_rate')
    # parser.add_argument('--weight_decay', type=float, default=1e-4, help='权重衰减 (L2 正则化)。')  # Deprecated, use --l2
    parser.add_argument('--patience', type=int, default=20, help='Patient early stop counter')
    parser.add_argument('--seed', type=int, default=42, help='seed')
    parser.add_argument('--attention_hidden_dim', type=int, default=32, help='The hidden dimensions of attention mechanism')
    parser.add_argument('--dropout_rate', type=float, default=0.5, help='dropout ratio')
    parser.add_argument('--focal_alpha', type=float, default=1.0, help='Focal alpha')
    parser.add_argument('--focal_gamma', type=float, default=2.0, help='Focal gamma')
    parser.add_argument('--gpu_id', type=int, default=None, help='GPU ID')
    parser.add_argument('--preprocessed_train_data_file', type=str, default='preprocessed_train_data.pkl', help='Pre processed training data save/load file name')
    parser.add_argument('--preprocessed_test_data_file', type=str, default='preprocessed_test_data.pkl', help='Pre processed test data save/load file name')
    parser.add_argument('--save_epoch_model', action='store_true', help='Save the model for each epoch')

    parser.add_argument('--l2', type=float, default=1e-4, help='L2 (weight_decay)')
    parser.add_argument('--hidden_layers', type=int, nargs='+', default=[256, 128], help='Size of hidden layers in neural networks. --hidden_layers 256 128')
    parser.add_argument('--rf_n_estimators', type=int, default=200, help='The number of estimators for random forests.')
    parser.add_argument('--rf_max_depth', type=int, default=10, help='The maximum depth of a random forest.')
    parser.add_argument('--rf_class_weight', type=str, default='balanced', help='The category weights of random forests.')
    parser.add_argument('--rf_min_samples_split', type=int, default=5, help='The minimum number of sample partitions for a random forest.')
    parser.add_argument('--rf_min_samples_leaf', type=int, default=2, help='The minimum number of leaf node samples in a random forest.')
    parser.add_argument('--use_preprocessed', action='store_true', help='Whether to use preprocessed data files.')
    return parser.parse_args()

def main():
    args = parse_arguments()
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Output directory has been created or already exists: {args.output_dir}")

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)
    print(f"The random seed has been set to: {args.seed}")

    preprocessed_train_data_path = os.path.join(args.output_dir, args.preprocessed_train_data_file)
    preprocessed_test_data_path = os.path.join(args.output_dir, args.preprocessed_test_data_file)

    if args.use_preprocessed:
        if os.path.exists(preprocessed_train_data_path) and os.path.exists(preprocessed_test_data_path):
            print(f"Load preprocessed training data from {preprocessed_train_data_path} ...")
            with open(preprocessed_train_data_path, 'rb') as f:
                preprocessed_train = pickle.load(f)
            print("The preprocessed training data has been successfully loaded.")

            print(f"Load preprocessed valid data from {preprocessed_test_data_path} ...")
            with open(preprocessed_test_data_path, 'rb') as f:
                preprocessed_test = pickle.load(f)
            print("The preprocessed valid data has been successfully loaded.")
        else:
            print("Preprocessed data file not found. Please ensure that the preprocessed data file exists or do not use the -- use_preprocessed option.")
            return
    else:
        
        print("Read training data...")
        train_df = read_data(args.train_file)
        if train_df.empty:
            print("No valid data was read from the training file. Please check the format of the data file.")
            return
        
        if train_df.shape[1] != 59:
            print(f"The number of training data columns is incorrect. Expected 59 columns (contig+position+56features + 1label), but found {train_df.shape[1]} columns.")
            return

        test_df = train_df.iloc[:,2:].sample(frac=args.valid_ratio)
        train_df = train_df.iloc[:,2:].drop(test_df.index)
        
        preprocessed_train = preprocess_data(train_df, fit=True)
        
        preprocessed_test = preprocess_data(
            test_df,
            scaler=preprocessed_train['scaler'],
            fit=False
        )

        
        with tqdm(total=2, desc="Save preprocessed data") as pbar:
            with open(preprocessed_train_data_path, 'wb') as f:
                pickle.dump(preprocessed_train, f)
            pbar.update(1)
            with open(preprocessed_test_data_path, 'wb') as f:
                pickle.dump(preprocessed_test, f)
            pbar.update(1)
        print(f"The preprocessed training data has been saved to {preprocessed_train_data_path}")
        print(f"The preprocessed validvalid data has been saved to{preprocessed_test_data_path}。")

    train_and_evaluate(args, preprocessed_train, preprocessed_test)

if __name__ == "__main__":
    main()
