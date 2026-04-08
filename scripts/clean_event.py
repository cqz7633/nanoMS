import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import os
from time import perf_counter
import psutil
from collections import defaultdict
from statistics import mode
import warnings
warnings.filterwarnings('ignore')

def parse_arguments():
    parser = argparse.ArgumentParser(description='Merge entries at the same location and perform data cleaning.')
    parser.add_argument('-i','--input', type=str, required=True, help='Input data file path.')
    parser.add_argument('-o','--output', type=str, required=True, help='Output result file path.')
    parser.add_argument('-p','--processes', type=int, default=4, help='Number of processes used.')
    parser.add_argument('-c','--chunk_size', type=int, default=1000000, help='chunk size')
    return parser.parse_args()

def filter_and_average_vectorized(values, std_multiplier=3):
    if len(values) == 1:
        return values[0]
    
    values = np.asarray(values, dtype=np.float64)
    median = np.median(values)
    std = np.std(values)
    
    if std == 0 or np.isnan(std):
        return median
    
    # Vectorized operations
    mask = np.abs(values - median) <= std_multiplier * std
    filtered = values[mask]
    
    return np.mean(filtered) if len(filtered) > 0 else median

def process_chunk_vectorized(chunk_df):
    """Process chunk using vectorized operations"""
    # Use groupby for efficient grouping
    grouped = chunk_df.groupby(['contig', 'position'])
    
    results = []
    for (contig, position), group in grouped:
        # Extract numpy arrays to speed up computation
        event_level_mean = group['event_level_mean'].values
        event_stdv = group['event_stdv'].values
        event_length = group['event_length'].values
        standardized_level = group['standardized_level'].values
        
        # Calculate the mode of reference_kmer
        ref_kmers = group['reference_kmer'].values
        if len(ref_kmers) == 1:
            ref_kmer = ref_kmers[0]
        else:
            # Use numpy's unique function to calculate the mode
            unique, counts = np.unique(ref_kmers, return_counts=True)
            ref_kmer = unique[counts.argmax()]
        
        result = {
            'contig': contig,
            'position': position,
            'reference_kmer': ref_kmer,
            'event_level_mean': filter_and_average_vectorized(event_level_mean),
            'event_stdv': filter_and_average_vectorized(event_stdv),
            'event_length': filter_and_average_vectorized(event_length),
            'model_stdv': group['model_stdv'].iloc[0],
            'standardized_level': filter_and_average_vectorized(standardized_level),
            'event_level_median': np.median(event_level_mean),
            'event_stdv_median': np.median(event_stdv),
            'event_length_median': np.median(event_length)
        }
        results.append(result)

    return results

def process_file_in_chunks(file_path, chunk_size, num_processes):
    """Process file in chunks with progress printing"""
    processed_chunks = 0
    
    def update_progress(result):
        """Callback function to update progress"""
        nonlocal processed_chunks
        processed_chunks += 1
        print(f"processing chunk {processed_chunks}, processed {(processed_chunks)*chunk_size} rows...", flush=True)
        return result

    with Pool(processes=num_processes) as pool:
        chunks = pd.read_csv(file_path, sep='\t', chunksize=chunk_size)
        results = []
        
        # Use callback function
        for chunk_result in pool.imap_unordered(process_chunk_vectorized, chunks):
            results.extend(update_progress(chunk_result))
        
        return merge_duplicate_positions(results)

def merge_duplicate_positions(results):
    """Merge identical positions that may appear in different chunks"""
    position_dict = defaultdict(list)
    
    # Group by position
    for result in results:
        key = (result['contig'], result['position'])
        position_dict[key].append(result)
    
    # Merge data at the same position
    final_results = []
    for key, group_results in tqdm(position_dict.items(), desc="Merging duplicate positions"):
        if len(group_results) == 1:
            final_results.append(group_results[0])
        else:
            # Need to merge multiple results
            merged = merge_multiple_results(group_results)
            final_results.append(merged)
    
    return final_results

def merge_multiple_results(results):
    """Merge multiple results at the same position"""
    # Collect all values
    all_values = defaultdict(list)
    for result in results:
        for key, value in result.items():
            if key not in ['contig', 'position']:
                if key == 'reference_kmer':
                    all_values[key].append(value)
                elif key in ['event_level_mean', 'event_stdv', 'event_length', 'standardized_level']:
                    all_values[key].append(value)
                elif key == 'model_stdv':
                    all_values[key] = value  # Keep only one value
                elif key.endswith('_median'):
                    all_values[key].append(value)
    
    # Calculate final values
    merged = {
        'contig': results[0]['contig'],
        'position': results[0]['position']
    }
    
    # Mode of reference_kmer
    ref_kmers = all_values['reference_kmer']
    if len(ref_kmers) == 1:
        merged['reference_kmer'] = ref_kmers[0]
    else:
        unique, counts = np.unique(ref_kmers, return_counts=True)
        merged['reference_kmer'] = unique[counts.argmax()]
    
    # Numerical fields
    for key in ['event_level_mean', 'event_stdv', 'event_length', 'standardized_level']:
        merged[key] = filter_and_average_vectorized(all_values[key])
    
    merged['model_stdv'] = all_values['model_stdv']
    
    # Median fields
    for key in ['event_level_median', 'event_stdv_median', 'event_length_median']:
        merged[key] = np.median(all_values[key])
    
    return merged

def main():
    args = parse_arguments()
    # Check if the file exists
    if not os.path.exists(args.input):
        print(f"Input file does not exist: {args.input}")
        return
    
    print(f"Start processing file: {args.input}")
    print(f"Using {args.processes} processes")
    
    # Process file
    results = process_file_in_chunks(args.input, args.chunk_size, args.processes)
    
    print(f"Processing completed, generating output file...")
    
    # Convert to DataFrame
    result_df = pd.DataFrame(results)
    
    # Data cleaning
    result_df["diff_stdv"] = (result_df["event_stdv"] - result_df["model_stdv"]) / result_df["model_stdv"]
    
    # Drop model_stdv column
    if 'model_stdv' in result_df.columns:
        result_df = result_df.drop(columns=['model_stdv'])
    
    # Sort by contig and position
    result_df = result_df.sort_values(['contig', 'position'])
    
    # Write to output file
    result_df.to_csv(args.output, sep='\t', index=False)
    
    print(f"Output file saved to: {args.output}")

if __name__ == "__main__":
    start = perf_counter()
    main()
    print(f"Total elapsed time: {perf_counter() - start:.6f} seconds")