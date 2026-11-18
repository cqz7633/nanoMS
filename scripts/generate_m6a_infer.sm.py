import argparse
import os
import warnings
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
from tqdm import tqdm

warnings.filterwarnings("ignore")

FEATURE_COLUMNS = [
    "event_level_mean",
    "event_stdv",
    "event_length",
    "standardized_level",
    "event_level_median",
    "event_stdv_median",
    "event_length_median",
    "diff_stdv",
]

DRACH_D = {"A", "G", "T"}
DRACH_R = {"A", "G"}
DRACH_H = {"A", "T", "C"}


def is_drach_motif(kmer):
    """Check whether a 5-base sequence matches DRACH."""
    if not isinstance(kmer, str) or len(kmer) != 5:
        return False

    return (
        kmer[0] in DRACH_D
        and kmer[1] in DRACH_R
        and kmer[2] == "A"
        and kmer[3] == "C"
        and kmer[4] in DRACH_H
    )


def process_group(args):
    """Process one read/contig group, preserving single-molecule results."""
    group_key, group_data, offset = args
    read_idx, contig = group_key

    group_data = group_data.sort_values("position", kind="mergesort")

    positions = group_data["position"].to_numpy(copy=False)
    kmers = group_data["reference_kmer"].to_numpy(copy=False)
    feature_matrix = group_data[FEATURE_COLUMNS].to_numpy(copy=False)

    window_size = 2 * offset + 1
    if len(positions) < window_size:
        return []

    # Candidate centers must have a complete, consecutive position window.
    centers = np.arange(offset, len(positions) - offset)
    window_positions = positions[
        centers[:, None] + np.arange(-offset, offset + 1)
    ]
    continuous = np.all(
        window_positions == positions[centers, None] + np.arange(-offset, offset + 1),
        axis=1,
    )
    centers = centers[continuous]

    if len(centers) == 0:
        return []

    # Preserve the original motif check at center position -2.
    motif_indices = centers - 2
    motif_valid = np.fromiter(
        (is_drach_motif(kmers[idx]) for idx in motif_indices),
        dtype=bool,
        count=len(motif_indices),
    )
    centers = centers[motif_valid]

    if len(centers) == 0:
        return []

    feature_offsets = np.arange(-offset, offset + 1)
    selected_features = feature_matrix[
        centers[:, None] + feature_offsets
    ]
    finite_rows = np.isfinite(selected_features).all(axis=(1, 2))
    centers = centers[finite_rows]
    selected_features = selected_features[finite_rows]

    if len(centers) == 0:
        return []

    rows = []
    for center_idx, features in zip(centers, selected_features):
        rows.append(
            [contig, positions[center_idx], read_idx, *features.reshape(-1)]
        )

    return rows


def get_output_columns(offset):
    output_columns = ["contig", "position", "read_name"]
    for window_offset in range(-offset, offset + 1):
        prefix = f"+{window_offset}" if window_offset > 0 else str(window_offset)
        output_columns.extend(
            f"{prefix}_{feature}" for feature in FEATURE_COLUMNS
        )
    return output_columns


def process_data_parallel(input_file, output_file, num_processes=None, kits="002"):
    if num_processes is None:
        num_processes = max(1, cpu_count() - 1)

    offset = 3 if kits == "002" else 5
    print(f"Using {num_processes} processes for parallel processing")
    print(f"Kit {kits}: extracting {2 * offset + 1} positions per site")

    dtype_dict = {
        "contig": "category",
        "position": "int32",
        "reference_kmer": "category",
        "read_name": "category",
        **{column: "float32" for column in FEATURE_COLUMNS},
    }

    print("Reading data via Pandas...")
    df = pd.read_csv(input_file, sep="\t", dtype=dtype_dict)
    print(f"Original data row count: {len(df)}")

    print("Building group index...")
    group_iter = df.groupby(["read_name", "contig"], sort=False, observed=True)
    total_groups = group_iter.ngroups
    output_columns = get_output_columns(offset)

    print("Starting parallel processing and streaming writes...")
    valid_rows_count = 0

    with open(output_file, "w", encoding="utf-8") as f_out:
        f_out.write("\t".join(output_columns) + "\n")

        with Pool(processes=num_processes) as pool:
            tasks = (
                (group_key, group_data, offset)
                for group_key, group_data in group_iter
            )
            progress = tqdm(
                pool.imap_unordered(process_group, tasks, chunksize=100),
                total=total_groups,
                desc="Parallel processing progress",
            )

            for result_rows in progress:
                if result_rows:
                    temp_df = pd.DataFrame(result_rows, columns=output_columns)
                    temp_df.to_csv(
                        f_out, sep="\t", index=False, header=False
                    )
                    valid_rows_count += len(result_rows)

    del df
    print(f"\nProcessing complete. Valid extracted rows: {valid_rows_count}")
    print(f"Results saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract context features of the m6A-associated DRACH motif "
            "at the single-molecule level."
        )
    )
    parser.add_argument("-i", "--input", required=True, help="Input file path")
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=4,
        help="Number of parallel processes (default: 4)",
    )
    parser.add_argument(
        "-k",
        "--kits",
        choices=["002", "004"],
        default="002",
        help="002 extracts 7 positions; 004 extracts 11 positions (default: 002)",
    )

    args = parser.parse_args()

    if not os.path.exists(args.input):
        parser.error(f"Input file does not exist: {args.input}")

    process_data_parallel(
        args.input,
        args.output,
        args.processes,
        args.kits,
    )


if __name__ == "__main__":
    main()