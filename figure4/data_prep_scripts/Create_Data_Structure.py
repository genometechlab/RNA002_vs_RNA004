import pysam
import numpy as np
import os
import time
from numba import njit
import pickle  
from tqdm import tqdm
import argparse

@njit
def get_modification_positions(numerical_sequence, target_number, deltas):
    """
    Convert 'deltas' into absolute positions for a particular 'target_number'.
    """
    deltas[1:] += 1
    target_indices = np.where(numerical_sequence == target_number)[0]
    return target_indices[np.cumsum(deltas)]

def mm_ml(mm_tag, ml_tag):
    """
    Yields tuples of (mod_label, deltas_array, ml_scores).
    """
    mm_parts = mm_tag.split(';')[:-1]  # remove last empty split
    mm_labels = [part.split(',')[0] for part in mm_parts]  # e.g. 'A+a'
    mm_deltas = [np.array(list(map(int, part.split(',')[1:]))) for part in mm_parts]

    ml_index = 0
    ml_tag_values = []
    for deltas in mm_deltas:
        length = len(deltas)
        ml_tag_values.append(ml_tag[ml_index : ml_index + length])
        ml_index += length

    for mod_label, deltas_array, ml_array in zip(mm_labels, mm_deltas, ml_tag_values):
        yield mod_label, deltas_array, ml_array

def parse_read_data(read):
    """
    Extract key info from a pysam.AlignedSegment.
    """
    mm_tag = read.get_tag('MM') if read.has_tag('MM') else None
    ml_tag = read.get_tag('ML') if read.has_tag('ML') else None
    if ml_tag is not None and not isinstance(ml_tag, list):
        ml_tag = list(ml_tag)

    return (
        read.query_name,
        read.query_sequence,
        read.get_aligned_pairs(with_seq=True, matches_only=True),
        mm_tag,
        ml_tag,
        read.reference_name
    )

def read_to_numpy_combined(read_tuple, valid_mod_labels, clamp_scores=True):
    """
    Return a single (num_bases, 3 + len(valid_mod_labels)) array for this read.
    Columns:
        0 = read position
        1 = reference position
        2 = base identity (A=0, C=1, G=2, T=3, N=4)
        3.. = ML scores (one column per mod_label), or -1 if no data for that mod
    """
    read_id, read_sequence, aligned_pairs, mm_tag, ml_tag, chromosome = read_tuple
    base_to_number = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    num_bases = len(read_sequence)

    data = np.full((num_bases, 3 + len(valid_mod_labels)), -1, dtype=np.int16)

    # Column 2 = base identity
    numeric_sequence = np.array([base_to_number.get(b, 4) for b in read_sequence], dtype=np.uint8)

    # Fill read_pos and ref_pos
    raw_read_pos, raw_ref_pos, ref_base = zip(*aligned_pairs) if aligned_pairs else ([], [], [])
    raw_read_pos = np.array([rp if rp is not None else -1 for rp in raw_read_pos], dtype=np.int16)
    raw_ref_pos  = np.array([rp if rp is not None else -1 for rp in raw_ref_pos], dtype=np.int16)

    # Only place them if in_range
    in_range = (raw_read_pos >= 0) & (raw_read_pos < num_bases)
    valid_indices = np.where(in_range)[0]

    data[raw_read_pos[valid_indices], 0] = raw_read_pos[valid_indices]  # read_pos
    data[raw_read_pos[valid_indices], 1] = raw_ref_pos[valid_indices]   # ref_pos
    data[:, 2] = numeric_sequence

    # If no MM/ML tags, done.
    if not mm_tag or not ml_tag:
        return data

    # For each mod_label present in the read, place its ML scores in the correct column.
    for mod_label, deltas_array, ml_scores in mm_ml(mm_tag, ml_tag):
        # If this mod_label is not one of our valid_mod_labels, skip it
        if mod_label not in valid_mod_labels:
            continue

        col_index = 3 + valid_mod_labels.index(mod_label)

        # Figure out which numeric base we are targeting
        if mod_label.startswith("A+"):
            target_number = 0
        elif mod_label.startswith("C+"):
            target_number = 1
        elif mod_label.startswith("G+"):
            target_number = 2
        elif mod_label.startswith("T+"):
            target_number = 3
        else:
            target_number = 4  # 'N'

        mod_positions = get_modification_positions(numeric_sequence, target_number, deltas_array)

        # Insert ML scores
        for pos_idx, score in zip(mod_positions, ml_scores):
            if clamp_scores:
                if score > 32767:
                    score = 32767
                if score < -32768:
                    score = -32768
            data[pos_idx, col_index] = score

    return data

def parse_args():
    parser = argparse.ArgumentParser(
        description="Process all user-specified mods in a BAM, store them by read_id -> single array."
    )
    parser.add_argument(
        "-b",
        "--bam_path",
        type=str,
        required=True,
        help="Path to the input BAM file."
    )
    parser.add_argument(
        "-t",
        "--tags",
        nargs="+",
        default=["A+a.", "A+17596.", "C+m.", "T+17802."],
        help="List of modification labels to look for."
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        default=".",
        help="Directory to save the output data."
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=True,
        help="Filename for output data structure"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="If set, only process the first 50 reads for debugging."
    )
    parser.add_argument(
        "--file_type",
        type=str,
        default="pickle",
        choices=["pickle", "npz", "hdf5"],
        help="Output file format."
    )
    parser.add_argument(
        "--strict_length",
        type=int,
        default=0,
        help="If > 0, only keep reads of this exact length."
    )
    
    return parser.parse_args()

def main():
    args = parse_args()
    bam_path = args.bam_path
    valid_mod_labels = args.tags

    print(f"\nProcessing BAM: {bam_path}")
    print(f"Looking for these modifications: {valid_mod_labels}")
    print("Debug mode is:", args.debug)

    results = {}

    bam = pysam.AlignmentFile(bam_path, "rb")

    read_iter = enumerate(bam)
    for idx, read in tqdm(read_iter):
        if args.debug and idx >= 50:
            break

        # Skip unwanted reads
        if (read.is_secondary or read.is_unmapped or
            read.is_supplementary or read.has_tag('pi') or read.is_reverse):
            continue

        # If user wants to strictly keep reads of certain length
        if args.strict_length > 0 and len(read.query_sequence) != args.strict_length:
            continue

        # Get the read info
        read_tuple = parse_read_data(read)
        read_id, read_sequence, aligned_pairs, mm_tag, ml_tag, chromosome = read_tuple

        # If no aligned_pairs or any None in them, skip
        if not aligned_pairs:
            continue
        if any(rp is None for _, rp, base in aligned_pairs):
            continue

        # Build the NumPy array
        combined_array = read_to_numpy_combined(read_tuple, valid_mod_labels, clamp_scores=True)
        results[read_id] = combined_array

    bam.close()

    # Save the dictionary
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Save a placehold for creating different filetype outputs
    base_filename = args.prefix
    if args.file_type == "pickle":
        outfile_path = os.path.join(args.output_dir, base_filename + ".pkl")
        print(f"\nSaving final dictionary to {outfile_path} (pickle)...")
        with open(outfile_path, "wb") as out_f:
            pickle.dump(results, out_f)
            
    print("Done.")

if __name__ == "__main__":
    main()
