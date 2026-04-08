import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
import os
import multiprocessing as mp
from multiprocessing import Pool, cpu_count
import warnings
import gc

warnings.filterwarnings('ignore')

def process_group(args):
    """Function to process a single data group"""
    group_key, group_data = args
    read_name, contig = group_key
    
    # Sort by position
    group_data = group_data.sort_values('position')
    
    # Extract core data as numpy arrays to improve lookup speed
    positions = group_data['position'].values
    kmers = group_data['reference_kmer'].values
    
    feature_columns = [
        'event_level_mean', 'event_stdv', 'event_length', 'standardized_level',
        'event_level_median', 'event_stdv_median', 'event_length_median', 'diff_stdv'
    ]
    feature_matrix = group_data[feature_columns].values
    
    # Build fast index 
    pos_to_idx = {pos: idx for idx, pos in enumerate(positions)}
    pos_to_kmer = {pos: kmers[idx] for idx, pos in enumerate(positions)}
    
    result_rows = []
    
    # [Modification 1]: Iterate using sorted(pos_to_idx.keys())
    # This not only strictly ensures ascending order processing but also inherently deduplicates 
    # redundant positions, perfectly aligning with the code below
    all_positions = sorted(pos_to_idx.keys())
    
    for center_pos in all_positions:
        # 1. Check for continuous upstream and downstream data (-3 to +3)
        required_positions = [center_pos + i for i in range(-3, 4)]
        if not all(p in pos_to_idx for p in required_positions):
            continue
            
        # 2. Check if position-2 exists
        pos_minus_2 = center_pos - 2
        if pos_minus_2 not in pos_to_kmer:
            continue
            
        # 3. Extract features of 7 windows
        features = []
        valid = True
        for offset in range(-3, 4):
            pos = center_pos + offset
            idx = pos_to_idx[pos]
            row_features = feature_matrix[idx]
            
            # Check for invalid values
            if np.any(np.isnan(row_features)) or np.any(np.isinf(row_features)):
                valid = False
                break
            features.extend(row_features)
            
        if not valid:
            continue
            
        # Keep output column order: contig, position, read_name, features...
        result_rows.append([contig, center_pos, read_name] + features)
        
    return result_rows

def get_output_columns():
    """Construct header for the output file"""
    feature_columns = [
        'event_level_mean', 'event_stdv', 'event_length', 'standardized_level',
        'event_level_median', 'event_stdv_median', 'event_length_median', 'diff_stdv'
    ]
    output_columns = ['contig', 'position', 'read_name']
    for offset in ['-3', '-2', '-1', '0', '+1', '+2', '+3']:
        for feat in feature_columns:
            output_columns.append(f"{offset}_{feat}")
    return output_columns

def process_data_parallel(input_file, output_file, num_processes=None):
    if num_processes is None:
        num_processes = max(1, cpu_count() - 1)
    
    print(f"Using {num_processes} processes for parallel processing")
    
    dtype_dict = {
        'contig': 'category',
        'position': 'int32',
        'reference_kmer': 'category',
        'read_name': 'category',
        'event_level_mean': 'float32',
        'event_stdv': 'float32',
        'event_length': 'float32',
        'standardized_level': 'float32',
        'event_level_median': 'float32',
        'event_stdv_median': 'float32',
        'event_length_median': 'float32',
        'diff_stdv': 'float32'
    }

    print("Reading and loading data into memory...")
    df = pd.read_csv(input_file, sep='\t', dtype=dtype_dict)
    
    print("Creating group iterator...")
    group_iter = df.groupby(['read_name', 'contig'])
    total_groups = group_iter.ngroups
    
    # 3. Prepare output file
    header = get_output_columns()
    
    # [Modification 2]: First use pandas to generate an empty file with a header to ensure header consistency
    pd.DataFrame(columns=header).to_csv(output_file, sep='\t', index=False)
    
    print(f"Starting parallel processing of {total_groups} data groups...")
    
    # Open file in append mode
    with open(output_file, 'a', newline='') as f:
        with Pool(processes=num_processes) as pool:
            pbar = tqdm(total=total_groups, desc="Processing progress")
            for result in pool.imap_unordered(process_group, group_iter, chunksize=200):
                if result:
                    # [Modification 3]: Use Pandas to_csv instead of native csv.writer for streaming write.
                    # This retains the original "compute and release memory" optimization strategy,
                    # and ensures that the underlying np.float32 format matches exactly with Pandas unified saving representation.
                    chunk_df = pd.DataFrame(result)
                    chunk_df.to_csv(f, sep='\t', header=False, index=False)
                pbar.update(1)
            pbar.close()

    # Release memory
    del df
    gc.collect()
    print(f"Processing complete, results saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Extract context features of the structural modification at the single-molecule level for subsequent inference.")
    parser.add_argument("-i", "--input", required=True, help="Input file path")
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    parser.add_argument("-p", "--processes", type=int, default=4, help="Number of parallel processes")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist")
        return
    
    print(f"Input file: {args.input}")
    process_data_parallel(args.input, args.output, args.processes)

if __name__ == "__main__":
    main()