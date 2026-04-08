import sys
import csv
import argparse
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import multiprocessing as mp

# Define feature column names to extract
DATA_COLS = [
    "event_level_mean", "event_stdv", "event_length", "standardized_level", 
    "event_level_median", "event_stdv_median", "event_length_median", "diff_stdv"
]

# DRACH rules dictionary
VALID_FIRST_BASES = {
    'D': {'A', 'G', 'T'},  # D -> A, G, T (U)
    'R': {'A', 'G'},       # R -> A, G
    'A': {'A'},            # A -> A
    'C': {'C'},            # C -> C
    'H': {'A', 'C', 'T'}   # H -> A, C, T (U)
}

global_binding_sites = set()

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract context features of the m6A-associated DRACH motif at the site level for training.")
    parser.add_argument("-i", "--input_file", required=True, type=str, help="Input clean_event.py file")
    parser.add_argument("-o", "--output_file", required=True, type=str, help="Output file path")
    parser.add_argument("-r", "--ref_pos_file", required=True, type=str, help="Reference position file")
    parser.add_argument("-b", "--base", type=int, choices=[0, 1], default=0, help="Label base type: 0 for 0-base, 1 for 1-base")
    parser.add_argument("-p", "--processes", type=int, default=4, help="Number of parallel processes (default: 4)")
    return parser.parse_args()

def init_worker(shared_sites):
    global global_binding_sites
    global_binding_sites = shared_sites

def is_valid_drach_sequence(kmer):
    if len(kmer) != 5:
        return False
    return (kmer[0] in VALID_FIRST_BASES['D'] and 
            kmer[1] in VALID_FIRST_BASES['R'] and 
            kmer[2] == 'A' and 
            kmer[3] == 'C' and 
            kmer[4] in VALID_FIRST_BASES['H'])

def read_binding_sites(binding_file, base_type):
    binding_sites = set()
    with open(binding_file, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        next(reader, None)
        for row in reader:
            if not row: continue
            transcript_id = row[0].split('|')[0].strip()
            transcript_pos = int(row[1])
            if base_type == 1:
                transcript_pos -= 1
            binding_sites.add((transcript_id, transcript_pos))
    return binding_sites

def chunk_generator(filepath, chunk_size=100000, overlap=6):
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        chunk = []
        for row in reader:
            chunk.append(row)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = chunk[-overlap:]  
        if len(chunk) > overlap:
            yield chunk

def process_chunk(chunk_rows):
    global global_binding_sites
    results = []
    
    for i in range(len(chunk_rows) - 6):
        window = chunk_rows[i:i+7]
        
        # 1. Check contig consistency
        contig_list = [row['contig'].split('|')[0].strip() for row in window]
        if len(set(contig_list)) != 1:
            continue
        contig_id = contig_list[0]
            
        # 2. Check position continuity
        try:
            positions = [int(row['position']) for row in window]
            if any(positions[j] != positions[0] + j for j in range(1, 7)):
                continue
        except ValueError:
            continue
            
        # 3. DRACH validation: Strictly corresponds to module[1][2] in the original code
        # The original code checks the reference_kmer of the 2nd row (index 1) within the window
        kmer_to_check = window[1].get('reference_kmer', '')
        if not is_valid_drach_sequence(kmer_to_check):
            continue
            
        # 4. Extract output site and Label: Strictly corresponds to module[3] in the original code
        # The original code uses the 4th row (index 3) as the representative position of the module for matching and output
        center_row = window[3]
        output_pos = int(center_row['position'])
        
        if (contig_id, output_pos) in global_binding_sites:
            label = 1
        else:
            label = 0
            
        # 5. Feature flattening and NaN filtering
        is_valid_features = True
        features = []
        for row in window:
            for col in DATA_COLS:
                try:
                    val = float(row[col])
                    if np.isnan(val) or np.isinf(val):
                        is_valid_features = False
                        break
                    features.append(val)
                except ValueError:
                    is_valid_features = False
                    break
            if not is_valid_features:
                break
                
        if is_valid_features:
            results.append([contig_id, output_pos] + features + [label])
            
    return results

def get_total_lines(filepath):
    with open(filepath, 'rb') as f:
        return sum(1 for _ in f)

def main():
    args = parse_arguments()

    if not os.path.exists(args.input_file):
        print(f"Error: {args.input_file} not found.")
        return

    print(f"\n[1/4] Loading reference sites (Base-{args.base})...")
    binding_sites = read_binding_sites(args.ref_pos_file, args.base)

    print(f"\n[2/4] Processing with {args.processes} processes...")
    total_lines = get_total_lines(args.input_file)
    chunk_size = 50000
    total_chunks = (total_lines // chunk_size) + 1
    
    new_header = ["contig", "position"]
    for i in range(-3, 4):
        prefix = f"+{i}" if i > 0 else str(i)
        for col in DATA_COLS:
            new_header.append(f"{prefix}_{col}")
    new_header.append("label")

    pool = mp.Pool(processes=args.processes, initializer=init_worker, initargs=(binding_sites,))
    generator = chunk_generator(args.input_file, chunk_size=chunk_size, overlap=6)
    
    all_valid_results = []
    with tqdm(total=total_chunks, desc="Stream Processing") as pbar:
        for res_chunk in pool.imap_unordered(process_chunk, generator):
            all_valid_results.extend(res_chunk)
            pbar.update(1)
            
    pool.close()
    pool.join()

    if not all_valid_results:
        print("Error: No data matched.")
        return

    print("\n[3/4] Converting to DataFrame & Deduplicating...")
    df = pd.DataFrame(all_valid_results, columns=new_header)
    feature_cols = df.columns[2:-1]
    df_flt = df[~df.duplicated(subset=feature_cols, keep=False)]
    
    print("\n[4/4] Shuffling and Saving...")
    df_shuffled = df_flt.sample(frac=1, random_state=42).reset_index(drop=True)
    df_shuffled.to_csv(args.output_file, sep="\t", index=False)
    print(f"Success! Output: {args.output_file}")

if __name__ == "__main__":
    main()