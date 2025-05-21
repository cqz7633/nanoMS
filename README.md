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

The parameters of the `FIAAU_integrate.R` script is provided as below:
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

### 2. Samples information files

The pipeline requires `FastQ` and `BAM` files as inputs, so we need to create two separate files to provide their locations and corresponding sample information. Please note that the sample information file does not require a header.

#### FastQ information file
The sample file of FastQ, the first column is the absolute path of the FastQ file. If it is `pair-end` data, the `R1` and `R2` files of each sample are arranged together. The second column is sample information, used to provide control and treatment group information in the FastQ files.

The FastQ sample file is provided as below:

| column1| column2 |
|--------|---------|
| /PATH/control_rep1_R1.fq.gz | control |
| /PATH/control_rep1_R2.fq.gz | control |
| /PATH/control_rep2_R1.fq.gz | control |
| /PATH/control_rep2_R2.fq.gz | control |
| /PATH/treatment_rep1_R1.fq.gz | treatment |
| /PATH/treatment_rep1_R2.fq.gz | treatment |
| /PATH/treatment_rep2_R1.fq.gz | treatment |
| /PATH/treatment_rep2_R2.fq.gz | treatment |

#### BAM information file
The sample file of FastQ, the first column is the absolute path of the BAM file. The second column is sample information, used to provide control and treatment group information in the BAM files.

The BAM sample file is provided as below:

| column1| column2 |
|--------|---------|
| /PATH/control_rep1.bam| control |
| /PATH/control_rep2.bam | control |
| /PATH/treatment_rep1.bam | treatment |
| /PATH/treatment_rep2.bam | treatment |

## Run FIAAU pipline

FIAAU consists of two parts, one for running six methods and the other for integrating the results of the six methods. Please do not change the structure of the original FIAAU directory.

### 6 methods process

The operation of six methods is controlled by a Python script, with each method generating a Bash script and running simultaneously.  

The parameters of the `FIAAU_process.py` script is provided as below:

```
usage: FIAAU_process.py [-h] [-f F] [-b B] [-c C] [-t T] [-p P] [-m M] [-r R] [-o O] [-l L] [-bs BS] [-ct CT]

optional arguments:
  -h, --help  show this help message and exit
  -f F        absolute paths of fastq information file
  -b B        absolute paths of bam information file
  -c C        control
  -t T        treatment
  -p P        paired or single ('y' or 'n') default: y
  -m M        core numbers default: 4
  -r R        reads length  default: 150
  -o O        out put dir  default: ./FIAAU_Y-m-d_H-M-S
  -l L        run the step 'makeTFfasta' of LABRAT or not ('y' or 'n') default: n
  -bs BS      bin size for calculating big wig file default: 10
  -ct CT      coverage cut off default: 0.5
```
An example of running a command is provided as below:

```
python FIAAU_process.py -f /PATH/FastQ_info.txt -b /PATH/Bam_info.txt -c control -t treatment -o /PATH/output/
```

### Integration results process

The integration process is controlled by an R script, which will generate the FIAAU_integate directory under the output parameter directory of FIAAU_process.py and save the results here. Please be careful not to change the structure of the six method directories and the names of the result files in the FIAAU process output directory. The integration section will read them in a relatively fixed path. 

The parameters of the `FIAAU_integrate.R` script is provided as below:
```
Usage: FIAAU_integrate.R [-[-help|h]] [-[-fiaau_dir|f] <character>] [-[-control|c] <character>] [-[-treatment|t] <character>] [-[-over_num|n] [<integer>]] [-[-qapa_cut|q] [<double>]] [-[-p_met|pm] [<double>]] [-[-p_int|pi] [<double>]]
    -h|--help         help
    -f|--fiaau_dir    FIAAU dir, the output dir for FIAAU_process.py
    -c|--control      control name
    -t|--treatment    treatment name
    -n|--over_num     method overlap number, default: 2
    -q|--qapa_cut     QAPA diff cutoff, default: 20
    -pm|--p_met        p value cutoff of each method, default: 0.05
    -pi|--p_int        p_integrate cutoff, default: 0.05
```
An example of running a command is provided as below:

```
Rscript FIAAU_integrate.R -f /PATH/FIAAU_process_output_dir -c control -t treatment
```
