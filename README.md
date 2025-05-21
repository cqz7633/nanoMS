<p align="center"><img src="./images/nanoMS_logo.png" width="500px" height="250px" style="float: left;" ></p>

******************

# nanoMS
Simultaneous detection of RNA m6A and structure from direct RNA-seq data.

## Create Environment with Conda
First, download the repository and create the environment.

```
git clone https://github.com/cqz7633/nanoMS.git
cd ./nanoMS
conda env create -f environment.yml
```

Then, activate the `nanoMS` environment.

```
conda activate nanoMS
```

## Data processing
We have provided demo data in the `./Data`.

### 1. Basecalling
Guppy performs data trimming, filtering and basecalling, using FAST5 format files as input.

Example usage is as follows:
```
guppy_basecaller -i /path/to/FAST5 -s /path/to/output --config /path/to/configuration
```
### 2. Align nanopore current signals to reference k-mers
Guppy performs data trimming, filtering and basecalling, using FAST5 format files as input.

Example usage is as follows:
```
nanopolish eventalign --reads /path/to/reads/FASTA --bam /path/to/aligned/bam --genome /path/to/genome/FASTA --scale-events > /path/to/events/align
```
*NOTE:* nanopolish has been integrated into the nanoMS environment.

### 3. Clean events
Clean the current information file obtained from nanopolis using the `clean_event.py` script

The parameters of the `clean_event.py` script is provided as below:
```
usage: clean_event.py [-h] --input INPUT --output OUTPUT
                      [--processes PROCESSES]

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT         Input data file path.
  --output OUTPUT       Output result file path.
  --processes PROCESSES
                        Number of processes used.
```
An example of running a command is provided as below:
```
python ./scripts/clean_event.py --input ./data/Demo_H9_nanopolish_events.tsv --output /PATH/to/clean_events.txt --processes 12
```

### 4.1 Generate training data for m6A sites

First, prepare a file containing known m6A sites as training labels, where the first column is the transcript name and the second column is the relative position of the m6A site within the transcript. 

The format is as follows:
| Trans| Position |
|--------|---------|
| ENST00000416718.2 | 82 |
| ENST00000416718.2 | 145 |
| ENST00000416718.2 | 157 |
| ENST00000327044.7 | 1335 |
| ENST00000477976.5 | 2766 |
| ENST00000379370.7 | 6241 |

*NOTE:* This file contains a header and is divided by tabs.

Generate training data from cleaned event data using the `generate_m6a_train.py` script.

The parameters of the `generate_m6a_train.py` script is provided as below:
```
usage: generate_m6a_train.py [-h] --input_file INPUT_FILE --output_file
                             OUTPUT_FILE --ref_pos_file REF_POS_FILE

optional arguments:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE
                        Input file after clean_event.py process
  --output_file OUTPUT_FILE
                        Output dir path
  --ref_pos_file REF_POS_FILE
                        Reference position file
```
An example of running a command is provided as below:
```
python ./scripts/generate_m6a_train.py --input_file /PATH/to/clean_events.txt --output_file /PATH/to/m6A_train_data.tsv --ref_pos_file ./data/Demo_H9_ref_position.tsv
```
### 4.2 Generate training data for secondary structure

First, prepare a file containing known RNA secondary structures as training labels. Here, we use icSHAPE data, where the first column is the transcript name, the second column is the transcript length, the third column is ‘*’, and each subsequent column represents the pairing score for each position in the transcript. The format is as follows:
| Trans| Length | * | Pos1 | Pos2 | Pos3 | Pos4 | ... |
|--------|---------|---------|---------|---------|---------|---------|---------|
| ENST00000574232.5 | 2129 | * | 0.189 | 0.529 | 0.043 | 0.087 | ... |
| ENST00000576646.7 | 629 | * | 0.176 | 0.347 | 0.347 | 0.871 | ... |
| ENST00000592202.5 | 1469 | * | 0.445 | 0.012 | 0.312 | 0.546 | ... |

*NOTE:* This file does not contains a header and is divided by tabs.

Generate training data from cleaned event data using the `generate_struct_train.py` script.

The parameters of the `generate_struct_train.py` script is provided as below:
```
usage: generate_struct_train.py [-h] --input_file INPUT_FILE --output_file
                                OUTPUT_FILE --shape_file SHAPE_FILE

optional arguments:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE
                        Input file after clean_event.py process
  --output_file OUTPUT_FILE
                        Output file path
  --shape_file SHAPE_FILE
                        Reference icshape file
```
An example of running a command is provided as below:
```
python ./scripts/generate_m6a_train.py --input_file /PATH/to/clean_events.txt --output_file /PATH/to/struct_train_data.tsv --shape_file ./data/Demo_H9_shape.tsv
```
### 4.3 Generate m6A sites and secondary structure data for inference

Simultaneously generate m6A sites and secondary structure data for inference using the `generate_infer_data.py` script.

The parameters of the `generate_infer_data.py` script is provided as below:
```
usage: generate_infer_data.py [-h] --input_file INPUT_FILE --output_dir
                              OUTPUT_DIR --file_prefix FILE_PREFIX

optional arguments:
  -h, --help            show this help message and exit
  --input_file INPUT_FILE
                        Input file after clean_event.py process
  --output_dir OUTPUT_DIR
                        Output dir path
  --file_prefix FILE_PREFIX
                        Output file prefix
```
An example of running a command is provided as below:
```
python ./scripts/generate_infer_data.py --input_file /PATH/to/clean_events.txt --output_dir /PATH/to/output/dir --file_prefix Demo_prefix
```

## nanoMS training

Train the nanoMS model using the nanoMS_train.py script, which accepts either m6A site data or secondary structure data as input. 

The parameters for nanoMS_train.py are as follows:
```
usage: nanoMS_train.py [-h] --train_file TRAIN_FILE [--valid_ratio VALID_RATIO] --output_dir OUTPUT_DIR
                       [--epochs EPOCHS] [--batch_size BATCH_SIZE] [--learning_rate LEARNING_RATE]
                       [--patience PATIENCE] [--seed SEED] [--attention_hidden_dim ATTENTION_HIDDEN_DIM]
                       [--dropout_rate DROPOUT_RATE] [--focal_alpha FOCAL_ALPHA] [--focal_gamma FOCAL_GAMMA]
                       [--gpu_id GPU_ID] [--preprocessed_train_data_file PREPROCESSED_TRAIN_DATA_FILE]
                       [--preprocessed_test_data_file PREPROCESSED_TEST_DATA_FILE] [--save_epoch_model] [--l2 L2]
                       [--hidden_layers HIDDEN_LAYERS [HIDDEN_LAYERS ...]] [--rf_n_estimators RF_N_ESTIMATORS]
                       [--rf_max_depth RF_MAX_DEPTH] [--rf_class_weight RF_CLASS_WEIGHT]
                       [--rf_min_samples_split RF_MIN_SAMPLES_SPLIT] [--rf_min_samples_leaf RF_MIN_SAMPLES_LEAF]
                       [--use_preprocessed]

optional arguments:
  -h, --help            show this help message and exit
  --train_file TRAIN_FILE
                        Train data file path
  --valid_ratio VALID_RATIO
                        Validation set ratio
  --output_dir OUTPUT_DIR
                        Save the directory of models and results.
  --epochs EPOCHS       Number of training epochs
  --batch_size BATCH_SIZE
                        Training batch size
  --learning_rate LEARNING_RATE
                        learning_rate
  --patience PATIENCE   Patient early stop counter
  --seed SEED           seed
  --attention_hidden_dim ATTENTION_HIDDEN_DIM
                        The hidden dimensions of attention mechanism
  --dropout_rate DROPOUT_RATE
                        dropout ratio
  --focal_alpha FOCAL_ALPHA
                        Focal alpha
  --focal_gamma FOCAL_GAMMA
                        Focal gamma
  --gpu_id GPU_ID       GPU ID
  --preprocessed_train_data_file PREPROCESSED_TRAIN_DATA_FILE
                        Pre processed training data save/load file name
  --preprocessed_test_data_file PREPROCESSED_TEST_DATA_FILE
                        Pre processed test data save/load file name
  --save_epoch_model    Save the model for each epoch
  --l2 L2               L2 (weight_decay)
  --hidden_layers HIDDEN_LAYERS [HIDDEN_LAYERS ...]
                        Size of hidden layers in neural networks. --hidden_layers 256 128
  --rf_n_estimators RF_N_ESTIMATORS
                        The number of estimators for random forests.
  --rf_max_depth RF_MAX_DEPTH
                        The maximum depth of a random forest.
  --rf_class_weight RF_CLASS_WEIGHT
                        The category weights of random forests.
  --rf_min_samples_split RF_MIN_SAMPLES_SPLIT
                        The minimum number of sample partitions for a random forest.
  --rf_min_samples_leaf RF_MIN_SAMPLES_LEAF
                        The minimum number of leaf node samples in a random forest.
  --use_preprocessed    Whether to use preprocessed data files.
```
An example of running a command is provided as below:
```
python ./nanoMS_train.py \
	--batch_size 32 \
	--train_file /PATH/to/m6a/or/structure/train/data \
	--epochs 150 \
	--patience 15 \
	--gpu_id 0 \
	--output_dir /PATH/to/output/dir \
```

## nanoMS inference
Perform inference using the trained nanoMS model via the nanoMS_infer.py script. 

The parameters of the nanoMS_infer.py script are as follows:
```
usage: nanoMS_infer.py [-h] --test_file TEST_FILE --model_dir MODEL_DIR
                       [--preprocessed_train_data_file PREPROCESSED_TRAIN_DATA_FILE] --output_file OUTPUT_FILE
                       [--gpu_id GPU_ID]

optional arguments:
  -h, --help            show this help message and exit
  --test_file TEST_FILE
                        The input file path
  --model_dir MODEL_DIR
                        Output directory for nanoMS training
  --preprocessed_train_data_file PREPROCESSED_TRAIN_DATA_FILE
                        Pre processed training data save/load file name
  --output_file OUTPUT_FILE
                        Path of output results
  --gpu_id GPU_ID       GPU ID
```
An example of running a command is provided as below:
```
python ./nanoMS_infer.py \
	--model_dir /PATH/to/nanoMS/trained/output/dir \
	--test_file /PATH/to/m6a/or/structure/inference/data \
	--output_file /PATH/to/output/file \
	--gpu_id 0 \
```
