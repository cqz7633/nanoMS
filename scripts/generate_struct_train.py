#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
from functools import partial
import multiprocessing as mp
import os
import shutil
import sys
import tempfile
import numpy as np
import pandas as pd
from tqdm import tqdm

# Numerical feature columns to extract
DATA_COLS = [
    "event_level_mean", "event_stdv", "event_length", "standardized_level", 
    "event_level_median", "event_stdv_median", "event_length_median", "diff_stdv"
]

# Declare global variable for child processes to share data
global_struct_dict = {}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract context features of the structural modification at the site level for training.")
    parser.add_argument("-i", "--input_file", required=True, type=str, help="Input file after clean_event.py process")
    parser.add_argument("-o", "--output_file", required=True, type=str, help="Output file path")
    parser.add_argument("-r", "--ref_shape_file", required=True, type=str, help="Reference icshape file")
    parser.add_argument("-k", "--kits", choices=["002", "004"], default="002", help="002 for 7-mer, 004 for 11-mer (default: 002)")
    parser.add_argument("-p", "--processes", type=int, default=4, help="Number of parallel processes (default: 4)")
    return parser.parse_args()

def init_worker(shared_dict):
    """Initialize the global dictionary for child processes"""
    global global_struct_dict
    global_struct_dict = shared_dict

def read_icshape_optimized(shape_file):
    """Optimized version of icshape reading, using numpy float32 to save a large amount of memory"""
    struct_dict = {}
    with open(shape_file, "r") as sf:
        for line in tqdm(sf, desc="Loading icshape"):
            line = line.strip()
            if not line:
                continue
            cols = line.split('\t')
            contig = cols[0].split('|')[0].strip()
            
            if len(cols) < 4:
                struct_dict[contig] = np.array([], dtype=np.float32)
                continue
                
            scores = []
            for val in cols[3:]:
                if val == 'NULL':
                    scores.append(np.nan)
                else:
                    try:
                        scores.append(float(val))
                    except ValueError:
                        scores.append(np.nan)
            struct_dict[contig] = np.array(scores, dtype=np.float32)
    return struct_dict

def get_file_info(filepath, offset):
    """Retrieve column index mappings and construct the dynamic output header."""
    with open(filepath, 'r') as f:
        header = f.readline().rstrip('\r\n').split('\t')
        
    contig_idx = header.index("contig")
    pos_idx = header.index("position")
    
    # Get column indices for DATA_COLS
    data_col_indices = [header.index(col) for col in DATA_COLS]
    
    new_header = ["contig", "position"]
    for i in range(-offset, offset + 1):
        prefix = f"+{i}" if i > 0 else str(i)
        for col_name in DATA_COLS:
            new_header.append(f"{prefix}_{col_name}")
    new_header.append("label")
    
    return contig_idx, pos_idx, data_col_indices, new_header

def chunk_generator(filepath, temp_dir, chunk_size, offset, contig_idx, pos_idx, data_col_indices):
    """Stream-read files in chunks and package tasks."""
    overlap = 2 * offset
    with open(filepath, 'r') as f:
        next(f)  # Skip header
        chunk = []
        batch_id = 0
        for line in f:
            chunk.append(line)
            if len(chunk) == chunk_size:
                tmp_out = os.path.join(temp_dir, f"part_{batch_id}.tsv")
                yield (chunk, tmp_out, offset, contig_idx, pos_idx, data_col_indices)
                chunk = chunk[-overlap:]  # Keep overlap to prevent window truncation
                batch_id += 1
        if len(chunk) > overlap:
            tmp_out = os.path.join(temp_dir, f"part_{batch_id}.tsv")
            yield (chunk, tmp_out, offset, contig_idx, pos_idx, data_col_indices)

def process_chunk(task):
    """
    Subprocess task: parse windows, filter/validate, label, and write results to a temporary file.
    """
    global global_struct_dict
    chunk_lines, tmp_out, offset, contig_idx, pos_idx, data_col_indices = task
    
    window_size = 2 * offset + 1
    # Use fast string splitting instead of DictReader for performance
    rows = [line.rstrip('\r\n').split('\t') for line in chunk_lines]
    n_rows = len(rows)
    
    results = []
    
    for i in range(n_rows - 2 * offset):
        window = rows[i : i + window_size]
        
        # 1. Check consistency of the contig
        contig0 = window[0][contig_idx].split('|')[0].strip()
        if any(row[contig_idx].split('|')[0].strip() != contig0 for row in window[1:]):
            continue
            
        # 2. Check for continuity in positions
        try:
            pos0 = int(window[0][pos_idx])
            if any(int(window[j][pos_idx]) != pos0 + j for j in range(1, window_size)):
                continue
        except ValueError:
            continue
            
        # 3. Label extraction & filtering (center index is the offset)
        pos_center = pos0 + offset
        contig_center = contig0
        
        scores = global_struct_dict.get(contig_center)
        if scores is None or pos_center < 0 or pos_center >= len(scores):
            continue
            
        label_val = scores[pos_center]
        if np.isnan(label_val):
            continue
            
        if label_val <= 0.2:
            final_label = 0
        elif label_val >= 0.8:
            final_label = 1
        else:
            continue  # Discard samples between 0.2 and 0.8
            
        # 4. Extract features & filter NaN/Inf
        is_valid = True
        features = []
        for row in window:
            for idx in data_col_indices:
                try:
                    val = float(row[idx])
                    if np.isnan(val) or np.isinf(val):
                        is_valid = False
                        break
                    features.append(val)
                except ValueError:
                    is_valid = False
                    break
            if not is_valid:
                break
                
        if is_valid:
            results.append([contig_center, pos_center] + features + [final_label])
            
    # Write to temporary file if data is valid to save RAM
    if results:
        with open(tmp_out, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(results)
        return tmp_out
    return None

def get_total_lines(filepath):
    """Estimate total lines in the file quickly."""
    with open(filepath, 'rb') as f:
        return sum(1 for _ in f)

def main():
    args = parse_arguments()
    input_file = args.input_file
    output_file = args.output_file
    shape_file = args.ref_shape_file
    threads = args.processes
    offset = 3 if args.kits == "002" else 5

    if not os.path.exists(input_file):
        print(f"Error: The file {input_file} does not exist.")
        return

    # 1. Memory optimization: Load struct dict
    print("\n[1/4] Loading structure dict...")
    struct_dict = read_icshape_optimized(shape_file)

    # 2. Prepare parallel processing and temporary directory
    print(f"\n[2/4] Processing modules with {threads} processes | Kit: {args.kits} (Offset: {offset})...")
    contig_idx, pos_idx, data_col_indices, new_header = get_file_info(input_file, offset)
    
    total_lines = get_total_lines(input_file) - 1
    chunk_size = 50000
    overlap = 2 * offset
    total_chunks = (total_lines // (chunk_size - overlap)) + 1
    
    temp_dir = tempfile.mkdtemp(prefix="struct_proc_")
    
    pool = mp.Pool(processes=threads, initializer=init_worker, initargs=(struct_dict,))
    generator = chunk_generator(input_file, temp_dir, chunk_size, offset, contig_idx, pos_idx, data_col_indices)
    
    temp_files = []
    with tqdm(total=total_chunks, desc="Processing chunks") as pbar:
        for res_tmp_file in pool.imap_unordered(process_chunk, generator):
            if res_tmp_file:
                temp_files.append(res_tmp_file)
            pbar.update(1)
            
    pool.close()
    pool.join()

    # Free memory used by the large dictionary
    del struct_dict 

    if not temp_files:
        print("Error: The filtered and merged module is empty. Exiting.")
        shutil.rmtree(temp_dir)
        sys.exit()

    # 3. Merge results and remove duplicates
    print(f"\n[3/4] Merging {len(temp_files)} part files & Dropping duplicates...")
    df_list = [pd.read_csv(tf, sep="\t", names=new_header) for tf in temp_files]
    full_df = pd.concat(df_list, ignore_index=True)
    
    # Clean up temporary directory
    shutil.rmtree(temp_dir)
    
    print(f"Extracted {len(full_df)} valid samples before duplicate removal.")
    
    # Deduplicate based on all feature columns (excluding contig, position, and label)
    feature_cols = full_df.columns[2:-1]
    duplicates_mask = full_df.duplicated(subset=feature_cols, keep=False)
    df_flt = full_df[~duplicates_mask]
    
    # 4. Save output
    print("\n[4/4] Saving output...")
    df_flt.to_csv(output_file, sep="\t", index=False)
    print(f"Success! Filtered dataframe shape: {df_flt.shape}")
    print(f"Output saved as {output_file}.")

if __name__ == "__main__":
    main()