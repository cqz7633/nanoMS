import sys
import csv
import argparse
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import multiprocessing as mp
from functools import partial

# Numerical feature columns to extract
DATA_COLS = [
    "event_level_mean", "event_stdv", "event_length", "standardized_level", 
    "event_level_median", "event_stdv_median", "event_length_median", "diff_stdv"
]

# Declare a global variable for child processes to share data (avoids huge memory consumption from passing large dictionaries in multiprocessing)
global_struct_dict = {}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract context features of the structural modification at the site level for training.")
    parser.add_argument("-i", "--input_file", required=True, type=str, help="Input file after clean_event.py process")
    parser.add_argument("-o", "--output_file", required=True, type=str, help="Output file path")
    parser.add_argument("-r", "--ref_shape_file", required=True, type=str, help="Reference icshape file")
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

def chunk_generator(filepath, chunk_size=100000, overlap=6):
    """
    Streaming generator for large files.
    Keep overlap (6 rows) to ensure 7-mer modules are not truncated at chunk boundaries
    """
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        chunk = []
        for row in reader:
            chunk.append(row)
            if len(chunk) == chunk_size:
                yield chunk
                chunk = chunk[-overlap:]  # Keep the last 6 rows as the beginning of the next chunk
        if len(chunk) > overlap:
            yield chunk

def process_chunk(chunk_rows):
    """Multiprocessing worker function: Completes continuity check, flatten, scoring, filtering, NaN removal, etc. within a window"""
    global global_struct_dict
    results = []
    
    for i in range(len(chunk_rows) - 6):
        window = chunk_rows[i:i+7]
        
        # 1. Check contig consistency
        contig0 = window[0]['contig'].split('|')[0].strip()
        if any(row['contig'].split('|')[0].strip() != contig0 for row in window[1:]):
            continue
            
        # 2. Check continuity
        try:
            pos0 = int(window[0]['position'])
            if any(int(window[j]['position']) != pos0 + j for j in range(1, 7)):
                continue
        except ValueError:
            continue
            
        # 3. Label value extraction & filtering (7-mer center is the 4th, i.e., index 3)
        pos4 = pos0 + 3
        contig4 = contig0
        
        scores = global_struct_dict.get(contig4)
        if scores is None or pos4 < 0 or pos4 >= len(scores):
            continue
            
        label_val = scores[pos4]
        if np.isnan(label_val):
            continue
            
        # Filter label based on original logic
        if label_val <= 0.2:
            final_label = 0
        elif label_val >= 0.8:
            final_label = 1
        else:
            continue  # Discard data between 0.2 and 0.8, and invalid data
            
        # 4. Feature flattening extraction & real-time NaN/Inf filtering
        is_valid = True
        features = []
        for row in window:
            for col in DATA_COLS:
                try:
                    val = float(row[col])
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
            # Assemble results: [contig, position] + features + [label]
            results.append([contig4, pos4] + features + [final_label])
            
    return results

def get_total_lines(filepath):
    """Quickly estimate the number of file lines to configure the progress bar"""
    with open(filepath, 'rb') as f:
        return sum(1 for _ in f)

def main():
    args = parse_arguments()
    input_file = args.input_file
    output_file = args.output_file
    shape_file = args.ref_shape_file
    threads = args.processes

    if not os.path.exists(input_file):
        print(f"Error: The file {input_file} does not exist.")
        return

    # 1. Memory optimization: Load scoring dictionary
    print("\n[1/4] Loading structure dict...")
    struct_dict = read_icshape_optimized(shape_file)

    # 2. Prepare for parallel processing
    print(f"\n[2/4] Processing modules with {threads} processes...")
    total_lines = get_total_lines(input_file)
    chunk_size = 50000
    total_chunks = (total_lines // chunk_size) + 1
    
    # Build header (strictly consistent with the original code in format and order)
    new_header = ["contig", "position"]
    offsets = range(-3, 4)
    for i in offsets:
        prefix = f"+{i}" if i > 0 else str(i)
        for col_name in DATA_COLS:
            new_header.append(f"{prefix}_{col_name}")
    new_header.append("label")

    # Create multiprocessing pool
    pool = mp.Pool(processes=threads, initializer=init_worker, initargs=(struct_dict,))
    generator = chunk_generator(input_file, chunk_size=chunk_size, overlap=6)
    
    # Parallel processing and real-time result collection
    all_valid_results = []
    with tqdm(total=total_chunks, desc="Processing chunks") as pbar:
        for res_chunk in pool.imap_unordered(process_chunk, generator):
            all_valid_results.extend(res_chunk)
            pbar.update(1)
            
    pool.close()
    pool.join()

    # Clean up the huge shared dictionary to release memory
    del struct_dict 

    if not all_valid_results:
        print("Error: The filtered and merged module is empty. Exiting.")
        sys.exit()

    print(f"\nExtracted {len(all_valid_results)} valid samples before duplicate removal.")

    # 3. Result merging and deduplication
    print("\n[3/4] Converting to DataFrame & Dropping duplicates...")
    df = pd.DataFrame(all_valid_results, columns=new_header)
    
    # Determine feature column range for deduplication (strictly aligned with original code: from the third column to the second to last)
    feature_cols = df.columns[2:-1]
    
    # Mark all duplicate rows (including all repeated instances) and discard them
    duplicates_mask = df.duplicated(subset=feature_cols, keep=False)
    df_flt = df[~duplicates_mask]
    
    # 4. Save file
    print("\n[4/4] Saving output...")
    df_flt.to_csv(output_file, sep="\t", index=False)
    print(f"Success! Filtered dataframe shape: {df_flt.shape}")
    print(f"Output saved as {output_file}.")

if __name__ == "__main__":
    main()