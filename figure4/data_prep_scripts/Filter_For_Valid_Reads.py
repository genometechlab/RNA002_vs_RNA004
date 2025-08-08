import pysam
import pickle
import numba
import numpy as np
from collections import defaultdict
from tqdm import tqdm
import argparse
import os
import subprocess


# Collect reference names from the reference fasta file
def get_reference_names_from_fasta(fasta_file):
    """
    Extracts reference names from a FASTA file.
    """
    with pysam.FastaFile(fasta_file) as fasta:
        return list(fasta.references)


# Count and filter read_id's within a bam
def count_and_filter(
    bam_file,
    reference_names,
    fasta_file,
    output_pkl_path,
    reverse_lookup_pkl_path,
    exact_match=False,
    debug=False
):
    """
    Counts reads that fully span a region and optionally enforces exact matches.

    - If exact_match=False, collects all reads that cover [start..stop).
    - If exact_match=True, collects only reads that cover [start..stop) AND whose
      aligned sequence in that region exactly matches the reference.

    Returns a tuple: (counts, read_ids_dict, reverse_lookup_dict)
    """
    counts = {}
    read_ids_dict = defaultdict(list)
    reverse_lookup_dict = {}

    # Open fasta file and loop through for every reference in fasta file
    with pysam.FastaFile(fasta_file) as fasta, pysam.AlignmentFile(bam_file, "rb") as bam:
        for reference in tqdm(reference_names, desc="Processing References"):
            if debug:
                print(f"Processing reference: {reference}")

            # Filter IVT and prefixified oligos separately
            if reference.endswith("_IVT"):
                start, stop = 49, 124  
            else:
                start, stop = 46, 121  
                
            if debug:
                print(f"Using interval: {reference}:{start}-{stop} (pysam half-open)")

            try:
                reference_sequence = fasta.fetch(reference, start, stop)
                if debug:
                    print(f"Reference sequence length: {len(reference_sequence)}")
                    print(f"Reference sequence: {reference_sequence}")

                total_reads_in_region = bam.count(reference=reference, start=start, end=stop)
                if debug:
                    print(f"bam.count() => {total_reads_in_region} reads overlap any part of {reference}:{start}-{stop}")

            except ValueError:
                counts[reference] = "Reference not found in BAM file"
                if debug:
                    print(f"Reference '{reference}' not found in BAM file.")
                continue

            read_positions = []
            coverage_filter_count = 0
            exact_match_count = 0

            # Loop through each read, ensure the reads are primary aligned and they are aligned within the specific regions
            for read_idx, read in enumerate(bam.fetch(reference=reference, start=start, stop=stop)):
                if (
                    not read.is_secondary
                    and not read.is_unmapped
                    and read.reference_start <= start
                    and read.reference_end >= stop
                ):
                    coverage_filter_count += 1

                    # If an exact match is wanted (to know essentially "perfect" reads) -- Not likely since 5' end is always cut 10-15 bases
                    if exact_match:
                        aligned_sequence = ''.join(
                            read.query_alignment_sequence[q_i]
                            for r_i, q_i in read.get_aligned_pairs(matches_only=True)
                            if (
                                r_i is not None
                                and q_i is not None
                                and 0 <= q_i < len(read.query_alignment_sequence)
                                and start <= r_i < stop
                            )
                        )
                        if aligned_sequence == reference_sequence:
                            exact_match_count += 1
                            read_positions.append(read.reference_start)
                            read_ids_dict[reference].append(read.query_name)
                            reverse_lookup_dict[read.query_name] = reference
                            if debug and read_idx < 10:
                                print(f"Read {read.query_name} matched exactly in {reference}:{start}-{stop}")
                    else:
                        read_positions.append(read.reference_start)
                        read_ids_dict[reference].append(read.query_name)
                        reverse_lookup_dict[read.query_name] = reference
                        if debug and read_idx < 10:
                            print(f"Read {read.query_name} covers region {reference}:{start}-{stop}")

            if exact_match:
                counts[reference] = exact_match_count
            else:
                counts[reference] = len(read_positions)

            if debug:
                if exact_match:
                    print(
                        f"For {reference}, {coverage_filter_count} reads fully covered the region. "
                        f"{exact_match_count} of those matched the reference exactly."
                    )
                else:
                    print(
                        f"For {reference}, {coverage_filter_count} reads fully covered the region. "
                        f"Collected {len(read_positions)} reads (no exact match enforced)."
                    )

    if debug:
        print("Saving read IDs to PKL files...")

    # Dump contents to pickle files
    with open(output_pkl_path, "wb") as pkl_file:
        pickle.dump(read_ids_dict, pkl_file)
    with open(reverse_lookup_pkl_path, "wb") as reverse_pkl_file:
        pickle.dump(reverse_lookup_dict, reverse_pkl_file)

    if debug:
        print("Done saving pickle files.")
        print("Summary counts:", counts)

    return counts, read_ids_dict, reverse_lookup_dict


def filter_and_index(num_threads, read_id_list, input_bam, output_bam):
    """
    Filters a BAM file using a given list of read IDs and then indexes the filtered BAM file.

    Parameters:
        num_threads (int): Number of threads to use.
        read_id_list (str): Path to the Read_ID_List.txt file.
        input_bam (str): Path to the BAM file needing filtering.
        output_bam (str): Path to the output filtered BAM file.
    """
    # Build the samtools view command.
    view_command = [
        "samtools", "view",
        "-@", str(num_threads),
        "-N", read_id_list,
        "-b", input_bam
    ]
    
    # Run samtools view and write its output to the output BAM file.
    try:
        with open(output_bam, "wb") as out_file:
            print(f"Running: {' '.join(view_command)} > {output_bam}")
            subprocess.run(view_command, stdout=out_file, check=True)
    except subprocess.CalledProcessError as e:
        print("Error during samtools view:")
        raise e

    # Build the samtools index command.
    index_command = [
        "samtools", "index",
        "-@", str(num_threads),
        output_bam
    ]
    
    # Run samtools index to index the newly created filtered BAM file.
    try:
        print(f"Running: {' '.join(index_command)}")
        subprocess.run(index_command, check=True)
    except subprocess.CalledProcessError as e:
        print("Error during samtools index:")
        raise e


def parse_args():
    parser = argparse.ArgumentParser(
        description="Read through a bamfile, and filter based on specific filtering criteria."
    )
    parser.add_argument(
        "-b",
        "--bam_path",
        type=str,
        required=True,
        help="Path to the input BAM file."
    )
    parser.add_argument(
        "-r",
        "--ref_path",
        type=str,
        required=True,
        help="Path to the fasta reference"
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=str,
        required=True,
        help="Path to output directory for all outfiles"
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=True,
        help="Name for the output file"        
    )
    parser.add_argument(
        "--exact_match",
        action="store_true",
        help="If an exact match to the reference is wanted, input as True"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Input True if you'd want to see the verbose process."
    )
    parser.add_argument(
        "-@",
        "--num_threads",
        type=int,
        default=1,
        help="The number of threads for samtools processing"
    )
    
    return parser.parse_args()


def main():
    args = parse_args()
    
    # Check that the output directory exists
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # Parse Arguments
    fasta_file_path = args.ref_path
    bam_file_path = args.bam_path
    output_pkl_path = os.path.join(args.outdir, f"{args.prefix}_read_ids_by_reference.pkl")
    reverse_lookup_pkl_path = os.path.join(args.outdir, f"{args.prefix}_read_id_to_reference.pkl")
    read_ids_list_path = os.path.join(args.outdir, f"{args.prefix}_read_ids_list.txt")
    num_threads = args.num_threads
    output_bam = os.path.join(args.outdir, f"{args.prefix}_bam_filtered_for_valid_reads.bam")

    # Change exact_match to True if you want to collect only reads that match exactly.
    counts, read_ids_dict, reverse_lookup_dict = count_and_filter(
        bam_file=bam_file_path,
        reference_names=get_reference_names_from_fasta(fasta_file_path),
        fasta_file=fasta_file_path,
        output_pkl_path=output_pkl_path,
        reverse_lookup_pkl_path=reverse_lookup_pkl_path,
        exact_match=args.exact_match,
        debug=args.debug
    )

    # Find all valid read IDs
    all_read_ids = []
    for ref_name, rids in read_ids_dict.items():
        all_read_ids.extend(rids)

    with open(read_ids_list_path, "w") as txt_file:
        for rid in all_read_ids:
            txt_file.write(rid + "\n")

    print("Preview of collected Read IDs:")
    for rid in all_read_ids[:10]:
        print(rid)
        
    print(f"Total unique read IDs: {len(set(all_read_ids))}")
    print(f"List of read IDs saved to: {read_ids_list_path}")

    # Run Samtools
    print("Running samtools view to filter and index reads")
    filter_and_index(num_threads, read_ids_list_path, bam_file_path, output_bam)
    print("Process Completed")

if __name__ == "__main__":
    main()
    