import argparse
import gc
import os
import warnings
from itertools import repeat
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


def process_group(args):
    """Process one read/contig group and extract its context features."""
    (read_name, contig), group_data, offset = args

    group_data = group_data.sort_values("position")
    positions = group_data["position"].to_numpy()
    feature_matrix = group_data[FEATURE_COLUMNS].to_numpy()
    row_count = len(positions)

    if row_count < 2 * offset + 1:
        return []

    # Preserve the original behavior for duplicate positions: the last row
    # for a position is used by the position-to-index mapping.
    pos_to_idx = {pos: idx for idx, pos in enumerate(positions)}
    unique_positions = np.fromiter(
        sorted(pos_to_idx), dtype=positions.dtype, count=len(pos_to_idx)
    )
    index_by_position = np.fromiter(
        (pos_to_idx[pos] for pos in unique_positions),
        dtype=np.intp,
        count=len(unique_positions),
    )

    # Only centers with a full positional range can have a valid window.
    if len(unique_positions) < 2 * offset + 1:
        return []

    center_start = np.searchsorted(unique_positions, unique_positions + offset)
    # Locate candidate center indices by their position in the sorted unique array.
    candidate_indices = np.arange(offset, len(unique_positions) - offset)
    center_positions = unique_positions[candidate_indices]

    # Require every position in the window, including the center, to exist.
    expected_positions = (
        center_positions[:, None]
        + np.arange(-offset, offset + 1, dtype=positions.dtype)
    )
    window_indices = np.searchsorted(unique_positions, expected_positions)
    in_bounds = window_indices < len(unique_positions)
    continuous = in_bounds.all(axis=1)

    if continuous.any():
        rows = np.flatnonzero(continuous)
        continuous[rows] = np.all(
            unique_positions[window_indices[rows]] == expected_positions[rows],
            axis=1,
        )

    # Original logic requires position -2 to be present; it does not test a motif.
    minus_two_indices = np.searchsorted(unique_positions, center_positions - 2)
    minus_two_valid = minus_two_indices < len(unique_positions)
    valid_rows = np.flatnonzero(minus_two_valid)
    minus_two_valid[valid_rows] = (
        unique_positions[minus_two_indices[valid_rows]]
        == center_positions[valid_rows] - 2
    )

    keep = continuous & minus_two_valid
    candidate_indices = candidate_indices[keep]
    center_positions = center_positions[keep]
    window_indices = window_indices[keep]

    if len(candidate_indices) == 0:
        return []

    # Map unique-position indices back to the original sorted group rows.
    feature_indices = index_by_position[window_indices]
    window_features = feature_matrix[feature_indices]
    finite = np.isfinite(window_features).all(axis=(1, 2))

    result_rows = []
    for center_pos, features in zip(center_positions[finite], window_features[finite]):
        result_rows.append(
            [contig, center_pos, read_name, *features.reshape(-1)]
        )

    return result_rows


def get_output_columns(offset):
    columns = ["contig", "position", "read_name"]
    for window_offset in range(-offset, offset + 1):
        prefix = f"+{window_offset}" if window_offset > 0 else str(window_offset)
        columns.extend(f"{prefix}_{feature}" for feature in FEATURE_COLUMNS)
    return columns


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

    print("Reading and loading data into memory...")
    df = pd.read_csv(input_file, sep="\t", dtype=dtype_dict)
    print(f"Input rows: {len(df)}")

    print("Creating group iterator...")
    group_iter = df.groupby(
        ["read_name", "contig"],
        sort=False,
        observed=True,
    )
    total_groups = group_iter.ngroups
    header = get_output_columns(offset)

    print(f"Starting parallel processing of {total_groups} data groups...")
    with open(output_file, "w", encoding="utf-8", newline="") as output:
        output.write("\t".join(header) + "\n")

        with Pool(processes=num_processes) as pool:
            tasks = (
                (group, group_data, offset)
                for group, group_data in group_iter
            )
            results = pool.imap_unordered(process_group, tasks, chunksize=100)

            for result in tqdm(results, total=total_groups, desc="Processing progress"):
                if result:
                    chunk_df = pd.DataFrame(result, columns=header)
                    chunk_df.to_csv(
                        output,
                        sep="\t",
                        header=False,
                        index=False,
                    )

    del df
    gc.collect()
    print(f"Processing complete, results saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract context features at the single-molecule level "
            "for subsequent inference."
        )
    )
    parser.add_argument("-i", "--input", required=True, help="Input file path")
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=4,
        help="Number of parallel processes",
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

    print(f"Input file: {args.input}")
    process_data_parallel(
        args.input,
        args.output,
        args.processes,
        args.kits,
    )


if __name__ == "__main__":
    main()