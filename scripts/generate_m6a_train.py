#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import multiprocessing as mp
import os
import shutil
import sys
import tempfile
import numpy as np
import pandas as pd
from tqdm import tqdm

# Feature columns
DATA_COLS = [
    "event_level_mean", "event_stdv", "event_length", "standardized_level", 
    "event_level_median", "event_stdv_median", "event_length_median", "diff_stdv"
]

# DRACH validation sets
DRACH_D = {"A", "G", "T"}
DRACH_R = {"A", "G"}
DRACH_H = {"A", "C", "T"}

def is_valid_drach_sequence(kmer_sequence):
    """Verify if the middle 5 bases match the DRACH motif."""
    if len(kmer_sequence) != 5:
        return False
    return (
        kmer_sequence[0] in DRACH_D
        and kmer_sequence[1] in DRACH_R
        and kmer_sequence[2] == "A"
        and kmer_sequence[3] == "C"
        and kmer_sequence[4] in DRACH_H
    )

def process_contig_records(records, binding_sites, offset):
    """Perform sliding window matching and feature extraction within a contig."""
    n = len(records)
    window_size = 2 * offset + 1
    if n < window_size:
        return []

    positions = [r["pos"] for r in records]
    first_bases = [r["base"] for r in records]
    results = []

    # Iterate through valid windows
    for i in range(offset, n - offset):
        # 1. Continuity check
        is_continuous = True
        for j in range(1, window_size):
            if positions[i - offset + j] != positions[i - offset] + j:
                is_continuous = False
                break
        if not is_continuous:
            continue

        # 2. DRACH verification: take the 5 positions at the window center
        # The window center index is i; the central 5 positions are [i-2, i-1, i, i+1, i+2]
        drach_kmer = "".join(first_bases[i - 2 : i + 3])
        if not is_valid_drach_sequence(drach_kmer):
            continue

        # 3. Label determination
        center_record = records[i]
        contig_center = center_record["raw_contig"]
        pos_center = positions[i]
        clean_contig = contig_center.split("|")[0]
        label = "1" if (clean_contig, pos_center) in binding_sites else "0"

        # 4. Extract flattened features
        merged_row = [contig_center, pos_center]
        for idx in range(i - offset, i + offset + 1):
            merged_row.extend(records[idx]["features"])
        merged_row.append(label)

        results.append(merged_row)
    return results

def worker_task(args):
    """Single process to handle a batch of data."""
    chunk_lines, binding_sites, tmp_out, offset = args
    records = []
    current_contig = None
    extracted_rows = []

    for line in chunk_lines:
        parts = line.rstrip("\r\n").split("\t")
        contig = parts[0]
        pos = int(parts[1])
        base = parts[2][0] # Assuming reference kmer is in the 3rd column
        # Get features, skipping the first three columns (contig, pos, ref_kmer)
        features = parts[3:] 

        if current_contig is not None and contig != current_contig:
            extracted_rows.extend(process_contig_records(records, binding_sites, offset))
            records = []

        current_contig = contig
        records.append({"raw_contig": contig, "pos": pos, "base": base, "features": features})

    if records:
        extracted_rows.extend(process_contig_records(records, binding_sites, offset))

    if extracted_rows:
        with open(tmp_out, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(extracted_rows)
    return tmp_out

def build_new_header(offset):
    new_header = ["contig", "position"]
    for i in range(-offset, offset + 1):
        prefix = f"+{i}" if i > 0 else str(i)
        for col_name in DATA_COLS:
            new_header.append(f"{prefix}_{col_name}")
    new_header.append("label")
    return new_header

def parse_arguments():
    parser = argparse.ArgumentParser(description="Efficient m6A site data extraction.")
    parser.add_argument("-i", "--input_file", required=True)
    parser.add_argument("-o", "--output_file", required=True)
    parser.add_argument("-r", "--ref_pos_file", required=True)
    parser.add_argument("-k", "--kits", choices=["002", "004"], default="002", help="002 for 7-mer, 004 for 11-mer (default: 002)")
    parser.add_argument("-p", "--processes", type=int, default=4)
    parser.add_argument("-b", "--base", type=int, choices=[0, 1], default=0, help="The ref_pos_file uses a 0-based or 1-based coordinate (default: 0)")
    return parser.parse_args()

def main():
    args = parse_arguments()
    offset = 3 if args.kits == "002" else 5
    
    print(f"Loading binding sites (Base-{args.base})...")
    binding_sites = set()
    with open(args.ref_pos_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        for row in reader:
            pos = int(row[1]) - args.base
            binding_sites.add((row[0].split('|')[0], pos))

    temp_dir = tempfile.mkdtemp(prefix="m6a_proc_")
    
    # Generate batch tasks
    tasks = []
    chunk_lines = []
    current_contig = None
    CHUNK_SIZE = 500 # Number of contigs processed per task
    batch_id = 0

    with open(args.input_file, "r") as f:
        next(f) # Skip header
        for line in f:
            contig = line.split("\t", 1)[0]
            if current_contig != contig:
                if len(chunk_lines) > 0 and len(tasks) % CHUNK_SIZE == 0:
                    tmp_out = os.path.join(temp_dir, f"part_{batch_id}.tsv")
                    tasks.append((chunk_lines, binding_sites, tmp_out, offset))
                    chunk_lines = []
                    batch_id += 1
                current_contig = contig
            chunk_lines.append(line)
        if chunk_lines:
            tmp_out = os.path.join(temp_dir, f"part_{batch_id}.tsv")
            tasks.append((chunk_lines, binding_sites, tmp_out, offset))

    print(f"Processing in {len(tasks)} chunks...")
    with mp.Pool(processes=args.processes) as pool:
        temp_files = list(tqdm(pool.imap(worker_task, tasks), total=len(tasks)))

    # Merge results
    print("Merging results...")
    new_header = build_new_header(offset)
    df_list = [pd.read_csv(tf, sep="\t", names=new_header) for tf in temp_files if os.path.exists(tf)]
    
    if not df_list:
        print("No data extracted.")
        return
        
    full_df = pd.concat(df_list, ignore_index=True)
    
    # Filter and save
    print("Deduplicating and filtering...")
    feature_cols = full_df.columns[2:-1]
    # Remove duplicates
    full_df = full_df[~full_df.duplicated(subset=feature_cols, keep=False)]
    full_df = full_df.sample(frac=1, random_state=42)
    full_df.to_csv(args.output_file, sep="\t", index=False)
    
    shutil.rmtree(temp_dir)
    print(f"Done. Saved to {args.output_file}")

if __name__ == "__main__":
    main()