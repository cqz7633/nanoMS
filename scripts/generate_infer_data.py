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

def get_file_info(filepath):
    """Extract header information and generate new header"""
    with open(filepath, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        fields = next(reader)
        
    contig_idx = fields.index('contig')
    pos_idx = fields.index('position')
    kmer_idx = fields.index('reference_kmer')
    
    data_indices = [i for i, f in enumerate(fields) if f not in ('contig', 'position', 'reference_kmer')]
    data_cols = [fields[i] for i in data_indices]
    
    new_header = ["contig", "position"]
    for i in range(-3, 4):
        prefix = f"+{i}" if i > 0 else str(i)
        for col_name in data_cols:
            new_header.append(f"{prefix}_{col_name}")
    new_header.append("label")
    
    return contig_idx, pos_idx, kmer_idx, data_indices, new_header

def reader_thread_func(filepath, chunk_size, overlap, task_q, num_workers):
    """Independent thread: responsible for streaming the file. Blocks automatically when the queue is full, perfectly controlling memory."""
    with open(filepath, 'r') as f:
        next(f)  # Skip header
        chunk = []
        for line in f:
            chunk.append(line)
            if len(chunk) == chunk_size:
                task_q.put(chunk)  # If the queue is full, it will pause and wait here
                chunk = chunk[-overlap:]  # Keep the last 6 lines as overlap
        if len(chunk) > overlap:
            task_q.put(chunk)
            
    # Send termination signal to all worker processes
    for _ in range(num_workers):
        task_q.put(None)

def worker_process_func(task_q, res_q, mode, contig_idx, pos_idx, kmer_idx, data_indices):
    """Core of multiprocessing: only process assigned task chunks"""
    while True:
        chunk_lines = task_q.get()
        if chunk_lines is None:
            res_q.put(None)  # Tell the main process that this process has finished
            break
            
        # Split strings in the child process to reduce main process overhead
        chunk_rows = [line.rstrip('\n').split('\t') for line in chunk_lines]
        
        stru_res = []
        m6a_res = []
        
        for i in range(len(chunk_rows) - 6):
            window = chunk_rows[i:i+7]
            
            contig0 = window[0][contig_idx]
            if any(row[contig_idx] != contig0 for row in window[1:]): continue
                
            try:
                pos0 = int(window[0][pos_idx])
                if any(int(window[j][pos_idx]) != pos0 + j for j in range(1, 7)): continue
            except ValueError:
                continue
                
            is_valid = True
            features = []
            for row in window:
                for idx in data_indices:
                    try:
                        val = float(row[idx])
                        if math.isnan(val) or math.isinf(val):
                            is_valid = False
                            break
                        features.append(str(val))
                    except ValueError:
                        is_valid = False
                        break
                if not is_valid: break
            if not is_valid: continue
                
            center_row = window[3]
            target_contig = center_row[contig_idx]
            target_pos = center_row[pos_idx]  # Directly take the position string from original data
            
            # Concatenate directly into the final string for writing, saving massive memory
            out_str = target_contig + '\t' + target_pos + '\t' + '\t'.join(features) + '\t-1\n'
            
            if mode in ['both', 'structure']:
                stru_res.append(out_str)
                
            if mode in ['both', 'm6a']:
                kmer_sequence = window[1][kmer_idx]
                if is_valid_drach_sequence(kmer_sequence):
                    m6a_res.append(out_str)
                    
        # Push results into the output queue
        res_q.put((stru_res, m6a_res))

def get_total_lines(filepath):
    with open(filepath, 'rb') as f:
        return sum(1 for _ in f) - 1

def main():
    args = parse_arguments()

    if not os.path.exists(args.input_file):
        print(f"Error: The file {args.input_file} does not exist.")
        sys.exit(1)

    site_output = os.path.join(args.output_dir, args.file_prefix + ".m6a.tsv") 
    stru_output = os.path.join(args.output_dir, args.file_prefix + ".structure.tsv") 

    contig_idx, pos_idx, kmer_idx, data_indices, new_header = get_file_info(args.input_file)
    
    print(f"Counting lines in {args.input_file}...")
    total_lines = get_total_lines(args.input_file)
    chunk_size = 50000
    overlap = 6
    total_chunks = (total_lines // (chunk_size - overlap)) + 1

    # Prepare for file writing
    f_m6a = open(site_output, 'w') if args.mode in ['both', 'm6a'] else None
    f_stru = open(stru_output, 'w') if args.mode in ['both', 'structure'] else None
    
    header_str = '\t'.join(new_header) + '\n'
    if f_m6a: f_m6a.write(header_str)
    if f_stru: f_stru.write(header_str)

    print(f"Processing with {args.processes} processes (Mode: {args.mode.upper()})")
    
    # === Core optimization: Use bounded queues to strictly limit memory ===
    task_queue = mp.Queue(maxsize=args.processes * 2)
    result_queue = mp.Queue(maxsize=args.processes * 2)

    # 1. Start an independent reader thread
    reader_thread = threading.Thread(
        target=reader_thread_func, 
        args=(args.input_file, chunk_size, overlap, task_queue, args.processes)
    )
    reader_thread.start()

    # 2. Start multiple consumer processes
    workers = []
    for _ in range(args.processes):
        p = mp.Process(
            target=worker_process_func,
            args=(task_queue, result_queue, args.mode, contig_idx, pos_idx, kmer_idx, data_indices)
        )
        p.start()
        workers.append(p)

    m6a_count = 0
    stru_count = 0
    active_workers = args.processes

    # 3. The main process is responsible for listening to results and safely writing to disk in a single thread
    with tqdm(total=total_chunks, desc="Processing Modules") as pbar:
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

    # Cleanup
    reader_thread.join()
    for p in workers:
        p.join()

    if f_m6a: f_m6a.close()
    if f_stru: f_stru.close()

    print("\nProcessing Complete!")
    if args.mode in ['both', 'structure']:
        print(f"Extracted {stru_count} Structure modules saved to {stru_output}")
    if args.mode in ['both', 'm6a']:
        print(f"Extracted {m6a_count} Filtered DRACH sites saved to {site_output}")

if __name__ == "__main__":
    main()