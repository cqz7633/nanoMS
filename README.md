<p align="center"><img src="./images/nanoMS_logo.png" width="500px" height="250px" style="float: left;" ></p>

******************

# nanoMS
Simultaneous detection of RNA m6A and structure from direct RNA-seq data

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
### 1. Clean events
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

### 2.1 Generate training data for m6A sites

First, prepare a file containing known m6A sites as training labels, where the first column is the transcript name and the second column is the relative position of the m6A site within the transcript. The format is as follows:
| Trans| Position |
|--------|---------|
| ENST00000416718.2 | 82 |
| ENST00000416718.2 | 145 |
| ENST00000416718.2 | 157 |
| ENST00000327044.7 | 1335 |
| ENST00000477976.5 | 2766 |
| ENST00000379370.7 | 6241 |

*NOTE:* This file contains a header.

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
### 2.2 Generate training data for secondary structure

First, prepare a file containing known m6A sites as training labels, where the first column is the transcript name and the second column is the relative position of the m6A site within the transcript. The format is as follows:
| Trans| Position |
|--------|---------|
| ENST00000416718.2 | 82 |
| ENST00000416718.2 | 145 |
| ENST00000416718.2 | 157 |
| ENST00000327044.7 | 1335 |
| ENST00000477976.5 | 2766 |
| ENST00000379370.7 | 6241 |

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

