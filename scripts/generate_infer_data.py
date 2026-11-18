#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import csv
import argparse
import os
import math
import threading
from tqdm import tqdm
import multiprocessing as mp

# Global DRACH rule dictionary
valid_first_bases = {
    'D': {'A', 'G', 'T'},
    'R': {'A', 'G'},
    'A': {'A'},
    'C': {'C'},
    'H': {'A', 'C', 'T'}
}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate m6A sites and secondary structure data for inference at the site level.")
    parser.add_argument("-i", "--input_file", type=str, required=True, help="Input file after clean_event.py process")
    parser.add_argument("-o", "--output_dir", type=str, required=True, help="Output dir path")
    parser.add_argument("-f", "--file_prefix", type=str, required=True, help="Output file prefix")
    parser.add_argument("-p", "--processes", type=int, default=4, help="Number of processes to use (default: 4)")
    parser.add_argument("-k", "--kits", choices=['002', '004'], default='002', help="002 for 7-mer, 004 for 11-mer (default: 002)")
    parser.add_argument("-m", "--mode", type=str, choices=['both', 'm6a', 'structure'], default='both',
                        help="Choose outputs: 'both' (default), 'm6a' (only m6a sites), 'structure' (only structure)")
    return parser.parse_args()

def is_valid_drach_sequence(kmer_sequence):
    if len(kmer_sequence) != 5: return False
    if kmer_sequence[0] not in valid_first_bases['D']: return False
    if kmer_sequence[1] not in valid_first_bases['R']: return False
    if kmer_sequence[2] != 'A': return False
    if kmer_sequence[3] != 'C': return False
    if kmer_sequence[4] not in valid_first_bases['H']: return False
    return True

def get_file_info(filepath, offset):
    """Extract header information and generate new header based on offset"""
    with open(filepath, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        fields = next(reader)
        
    contig_idx = fields.index('contig')
    pos_idx = fields.index('position')
    kmer_idx = fields.index('reference_kmer')
    
    data_indices = [i for i, f in enumerate(fields) if f not in ('contig', 'position', 'reference_kmer')]
    data_cols = [fields[i] for i in data_indices]
    
    new_header = ["contig", "position"]
    # Dynamically generate header
    for i in range(-offset, offset + 1):
        prefix = f"+{i}" if i > 0 else str(i)
        for col_name in data_cols:
            new_header.append(f"{prefix}_{col_name}")
    new_header.append("label")
    
    return contig_idx, pos_idx, kmer_idx, data_indices, new_header

def reader_thread_func(filepath, chunk_size, overlap, task_q, num_workers):
    """Streaming file reader with overlap support."""
    with open(filepath, 'r') as f:
        next(f)  # Skip header
        chunk = []
        for line in f:
            chunk.append(line)
            if len(chunk) == chunk_size:
                task_q.put(chunk)
                # Keep overlap rows as the start of the next chunk
                chunk = chunk[-overlap:]
        if len(chunk) > overlap:
            task_q.put(chunk)
            
    for _ in range(num_workers):
        task_q.put(None)

def worker_process_func(task_q, res_q, mode, contig_idx, pos_idx, kmer_idx, data_indices, offset):
    """Core processor with dynamic sliding window based on offset."""
    window_size = 2 * offset + 1
    # Verify the DRACH kmer offset in the window.
    # Regardless of window size, the DRACH 5-mer always starts 2 positions before the center.
    # Center index of window is 'offset'. DRACH start index is 'offset - 2'.
    drach_check_idx = offset - 2 

    while True:
        chunk_lines = task_q.get()
        if chunk_lines is None:
            res_q.put(None)
            break
            
        chunk_rows = [line.rstrip('\n').split('\t') for line in chunk_lines]
        
        stru_res = []
        m6a_res = []
        
        # Iterate through chunk with a sliding window of size window_size
        for i in range(len(chunk_rows) - 2 * offset):
            window = chunk_rows[i : i + window_size]
            
            # 1. Check contig consistency
            contig0 = window[0][contig_idx]
            if any(row[contig_idx] != contig0 for row in window[1:]): continue
            
            # 2. Check position continuity
            try:
                pos0 = int(window[0][pos_idx])
                if any(int(window[j][pos_idx]) != pos0 + j for j in range(1, window_size)): continue
            except ValueError:
                continue
                
            # 3. Extract features and validate NaNs
            is_valid = True
            features = []
            for row in window:
                for idx in data_indices:
                    val_str = row[idx]
                    # Quick validation for NaN/inf
                    if val_str in ('nan', 'NaN', 'inf', '-inf', ''):
                        is_valid = False; break
                    features.append(val_str)
                if not is_valid: break
            if not is_valid: continue
                
            # Extract info from center row
            center_row = window[offset]
            out_str = f"{center_row[contig_idx]}\t{center_row[pos_idx]}\t" + '\t'.join(features) + '\t-1\n'
            
            if mode in ['both', 'structure']:
                stru_res.append(out_str)
                
            if mode in ['both', 'm6a']:
                # Get kmer for DRACH validation (corresponding position in center window)
                kmer_sequence = window[drach_check_idx][kmer_idx]
                if is_valid_drach_sequence(kmer_sequence):
                    m6a_res.append(out_str)
                    
        res_q.put((stru_res, m6a_res))

def main():
    args = parse_arguments()
    offset = 3 if args.kits == '002' else 5
    
    if not os.path.exists(args.input_file):
        print(f"Error: The file {args.input_file} does not exist.")
        sys.exit(1)

    site_output = os.path.join(args.output_dir, f"{args.file_prefix}.m6a.tsv") 
    stru_output = os.path.join(args.output_dir, f"{args.file_prefix}.structure.tsv") 

    contig_idx, pos_idx, kmer_idx, data_indices, new_header = get_file_info(args.input_file, offset)
    
    # Estimate progress bar
    total_lines = sum(1 for _ in open(args.input_file, 'rb')) - 1
    chunk_size = 50000
    overlap = 2 * offset
    total_chunks = (total_lines // (chunk_size - overlap)) + 1

    f_m6a = open(site_output, 'w') if args.mode in ['both', 'm6a'] else None
    f_stru = open(stru_output, 'w') if args.mode in ['both', 'structure'] else None
    
    header_str = '\t'.join(new_header) + '\n'
    if f_m6a: f_m6a.write(header_str)
    if f_stru: f_stru.write(header_str)

    print(f"Processing with {args.processes} processes | Kit: {args.kits} (Offset: {offset}) | Mode: {args.mode.upper()}")
    
    task_queue = mp.Queue(maxsize=args.processes * 2)
    result_queue = mp.Queue(maxsize=args.processes * 2)

    reader_thread = threading.Thread(
        target=reader_thread_func, 
        args=(args.input_file, chunk_size, overlap, task_queue, args.processes)
    )
    reader_thread.start()

    workers = []
    for _ in range(args.processes):
        p = mp.Process(
            target=worker_process_func,
            args=(task_queue, result_queue, args.mode, contig_idx, pos_idx, kmer_idx, data_indices, offset)
        )
        p.start()
        workers.append(p)

    m6a_count = 0
    stru_count = 0
    active_workers = args.processes

    with tqdm(total=total_chunks, desc="Processing") as pbar:
        while active_workers > 0:
            res = result_queue.get()
            if res is None:
                active_workers -= 1
            else:
                stru_res, m6a_res = res
                if f_stru and stru_res:
                    f_stru.writelines(stru_res)
                    stru_count += len(stru_res)
                if f_m6a and m6a_res:
                    f_m6a.writelines(m6a_res)
                    m6a_count += len(m6a_res)
                pbar.update(1)

    reader_thread.join()
    for p in workers:
        p.join()

    if f_m6a: f_m6a.close()
    if f_stru: f_stru.close()

    print("\nProcessing Complete!")
    if args.mode in ['both', 'structure']: print(f"Saved structure to {stru_output}")
    if args.mode in ['both', 'm6a']: print(f"Saved m6A sites to {site_output}")

if __name__ == "__main__":
    main()