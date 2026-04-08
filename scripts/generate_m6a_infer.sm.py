import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse
import os
from multiprocessing import Pool, cpu_count
import warnings
import gc

warnings.filterwarnings('ignore')

def is_drach_motif(kmer):
    """Check if it is a DRACH motif (D=A/G/T, R=A/G, A=A, C=C, H=A/T/C)"""
    if len(kmer) != 5:
        return False
    
    # Use sets for fast lookup
    d_positions = {'A', 'G', 'T'}
    r_positions = {'A', 'G'}
    a_position = {'A'}
    c_position = {'C'}
    h_positions = {'A', 'T', 'C'}
    
    return (kmer[0] in d_positions and 
            kmer[1] in r_positions and 
            kmer[2] in a_position and 
            kmer[3] in c_position and 
            kmer[4] in h_positions)

def process_group(args):
    """Function to process a single data group, used for multiprocessing"""
    group_key, group_data = args
    read_idx, contig = group_key
    
    # Sort each group by position
    group_data = group_data.sort_values('position')
    
    # Use numpy arrays to improve performance
    positions = group_data['position'].values
    kmers = group_data['reference_kmer'].values
    
    # Feature columns
    feature_columns = [
        'event_level_mean', 'event_stdv', 'event_length', 'standardized_level',
        'event_level_median', 'event_stdv_median', 'event_length_median', 'diff_stdv'
    ]
    
    # Extract feature matrix
    feature_matrix = group_data[feature_columns].values
    
    # Create mapping dictionary from position to index
    pos_to_idx = {pos: idx for idx, pos in enumerate(positions)}
    pos_to_kmer = {pos: kmers[idx] for idx, pos in enumerate(positions)}
    
    result_rows = []
    
    # Pre-calculate all possible center positions
    all_positions = sorted(pos_to_idx.keys())
    
    for center_pos in all_positions:
        # Check for continuous upstream and downstream data (-3 to +3)
        required_positions = [center_pos + i for i in range(-3, 4)]
        if not all(pos in pos_to_idx for pos in required_positions):
            continue
        
        # Check if position-2 is a DRACH motif
        pos_minus_2 = center_pos - 2
        if pos_minus_2 not in pos_to_kmer:
            continue
            
        kmer_minus_2 = pos_to_kmer[pos_minus_2]
        if not is_drach_motif(kmer_minus_2):
            continue
        
        # Extract features - use numpy array operations to improve performance
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
        
        # Build result row
        result_row = [contig, center_pos, read_idx] + features
        result_rows.append(result_row)
    
    return result_rows

def get_output_columns():
    """Get output column names (consistent with the original code)"""
    feature_columns = [
        'event_level_mean', 'event_stdv', 'event_length', 'standardized_level',
        'event_level_median', 'event_stdv_median', 'event_length_median', 'diff_stdv'
    ]
    output_columns = ['contig', 'position', 'read_name']
    for offset in ['-3', '-2', '-1', '0', '+1', '+2', '+3']:
        for feature in feature_columns:
            output_columns.append(f"{offset}_{feature}")
    return output_columns

def process_data_parallel(input_file, output_file, num_processes=None):
    """Memory optimization + precision alignment parallel processing version"""
    
    if num_processes is None:
        num_processes = max(1, cpu_count() - 1)
    
    print(f"Using {num_processes} processes for parallel processing")
    
    # Force type conversion to retain the original code's float32 precision
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
    
    print("Reading data via Pandas...")
    # One-time read: heavily compressed relying on dtype_dict, 10 million lines take up only ~300MB of memory
    df = pd.read_csv(input_file, sep='\t', dtype=dtype_dict)
    print(f"Original data row count: {len(df)}")
    
    print("Building group index...")
    # Core optimization 1: Absolutely do not use list(df.groupby()), only generate an iterator to prevent memory explosion
    group_iter = df.groupby(['read_name', 'contig'])
    total_groups = group_iter.ngroups
    
    output_columns = get_output_columns()
    
    print("Starting parallel processing and streaming writes...")
    valid_rows_count = 0
    
    # Open file ready for writing while computing
    with open(output_file, 'w', encoding='utf-8') as f_out:
        # Write header
        f_out.write('\t'.join(output_columns) + '\n')
        
        with Pool(processes=num_processes) as pool:
            # Core optimization 2: imap_unordered returns results to be processed immediately, no memory backlog
            pbar = tqdm(pool.imap_unordered(process_group, group_iter, chunksize=100), 
                        total=total_groups, desc="Parallel processing progress")
            
            for result_rows in pbar:
                if result_rows:
                    # Core optimization 3: Convert each small batch of results to a DataFrame, write using Pandas' native to_csv.
                    # This ensures the process of converting floats to strings is exactly identical to the original code, absolutely no precision differences.
                    temp_df = pd.DataFrame(result_rows, columns=output_columns)
                    temp_df.to_csv(f_out, sep='\t', index=False, header=False)
                    valid_rows_count += len(result_rows)
    
    # Clean memory
    del df
    gc.collect()
    
    print(f"\nProcessing complete. Valid extracted rows: {valid_rows_count}")
    print(f"Results saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Extract context features of the m6A-associated DRACH motif at the single-molecule level for subsequent inference.")
    parser.add_argument("-i", "--input", required=True, help="Input file path")
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    parser.add_argument("-p", "--processes", type=int, default=4, 
                       help="Number of parallel processes (default: 4)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist")
        return
    
    process_data_parallel(args.input, args.output, args.processes)

if __name__ == "__main__":
    main()