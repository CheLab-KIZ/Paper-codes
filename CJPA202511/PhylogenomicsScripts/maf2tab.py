#!/usr/bin/env python3
# -*- coding:utf-8 -*-

# Func: Convert MAF file to a tabular format with specified individuals.
# Usage: maf2tab.py [-h] -i INPUT -I INDV [-r REFERENCE] [-o OUTPUT] [-c]
#                   [--out-freq OUT_FREQ] [--ref-gap] [--var]
#                   [--max-mis MAX_MIS] [--g2n | --n2g]
# By Wenjie DONG [2024/8/11]


import argparse
import os
import gzip
import sys


def parse_species_list(species_file):
    with open(species_file, 'r') as f:
        species_list = [line.strip() for line in f if line.strip()]
    return species_list


def open_maf_file(maf_file):
    if maf_file.endswith('.gz'):
        return gzip.open(maf_file, 'rt')
    else:
        return open(maf_file, 'r')


def open_output(outputfile, stdout):
    if stdout:
        return sys.stdout
    else:
        return open(outputfile, 'w')


def process_maf(maf_file, species_list, reference, output_file, stdout, out_freq, ref_gap, snp_only, max_mis, g2n, n2g):
    with open_maf_file(maf_file) as maf, open_output(output_file, stdout) as output:
        block = []
        chromosome = None
        ref_seq = None
        ref_start = None
        species_sequences = {species: None for species in species_list}

        # Buffer list to store lines before writing to file
        buffer = []
        count = 0

        # Write the header of the output table
        output.write("chr\tpos\tref\t" + "\t".join(species_list) + "\n")

        for line in maf:
            line = line.strip()

            # Process the alignment block when an empty line is encountered
            if line == "" and block:
                if ref_seq is not None:
                    ref_pos = ref_start + 1
                    gap_counter = 0  # Used to track the fractional positions for gaps

                    for i, base in enumerate(ref_seq):
                        # Count the base types excluding gaps and missing data
                        bases = [seq[i] for seq in species_sequences.values() if seq and seq[i] not in ['-', 'N']]
                        base_types = set(bases)
                        missing_count = len(species_list) - len(bases)
                        missing_ratio = missing_count / len(species_list)

                        # Filter out positions with missing data ratio higher than the threshold
                        if missing_ratio > max_mis:
                            continue

                        # Check if it's a variant position
                        if base != '-' and (not snp_only or len(base_types) > 1):
                            row = [chromosome, ref_pos, ref_seq[i]]

                            for species in species_list:
                                seq = species_sequences[species]
                                base = seq[i] if seq and seq[i] != '-' else 'N'
                                if g2n and base == '-':
                                    base = 'N'
                                elif n2g and base == 'N':
                                    base = '-'
                                row.append(base)

                            buffer.append("\t".join(map(str, row)))
                            count += 1
                            ref_pos += 1
                            gap_counter = 0  # Reset gap counter
                        elif base == '-' and ref_gap:
                            # Handle gaps in the reference with fractional positions
                            if ref_gap and (not snp_only or len(base_types) > 1):
                                gap_counter += 1
                                gap_position = f"{ref_pos - 1}.{gap_counter}"
                                row = [chromosome, gap_position, ref_seq[i]]

                                for species in species_list:
                                    seq = species_sequences[species]
                                    base = seq[i] if seq and seq[i] != '-' else 'N'
                                    if g2n and base == '-':
                                        base = 'N'
                                    elif n2g and base == 'N':
                                        base = '-'
                                    row.append(base)

                                buffer.append("\t".join(map(str, row)))
                                count += 1

                        # Write to file when the buffer reaches the specified frequency
                        if count >= out_freq:
                            output.write("\n".join(buffer) + "\n")
                            buffer = []
                            count = 0

                block = []
                ref_seq = None
                species_sequences = {species: None for species in species_list}
                continue

            # Add the line to the current block
            block.append(line)

            # Process each line in the alignment block
            if line.startswith("s "):
                tokens = line.split()
                seq_info = tokens[1].split('.')
                species_name = seq_info[0]
                seq_name = seq_info[1]
                seq_start = int(tokens[2])
                seq = tokens[6].upper()  # Convert bases to uppercase

                # Determine the reference sequence
                if ref_seq is None:
                    if reference and species_name == reference:
                        ref_seq = seq
                        ref_start = seq_start
                        chromosome = seq_name
                    elif not reference:
                        ref_seq = seq
                        ref_start = seq_start
                        chromosome = seq_name

                if species_name in species_sequences:
                    species_sequences[species_name] = seq

        # Ensure any remaining buffered content is written to the file
        if buffer:
            output.write("\n".join(buffer) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="maf2tab: Convert MAF file to a tabular format with specified individuals.")
    parser.add_argument("-i", "--input", required=True, help="input MAF file")
    parser.add_argument("-I", "--indv", required=True, help="file with list of individual names (one per line)")
    parser.add_argument("-r", "--reference", help="reference species name (default: first sequence in the block)")
    parser.add_argument("-o", "--output", help="output file name (default: input file prefix + '.tab')")
    parser.add_argument("-c", "--stdout", action="store_true", help="output to standard output instead of a file")
    parser.add_argument("--out-freq", type=int, default=1, help="frequency of writing output to file (default: 1)")
    parser.add_argument("--ref-gap", action="store_true", help="include gaps in reference with fractional positions")
    parser.add_argument("--var", action="store_true", help="only output variant sites")
    parser.add_argument("--max-mis", type=float, default=1.0, help="maximum missing proportion (0~1) to include a position (default: 1)")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--g2n", action="store_true", help="convert all gaps to 'N'")
    group.add_argument("--n2g", action="store_true", help="convert all missing characters (N) to gaps")

    args = parser.parse_args()

    # If the output file name is not specified, use the input file prefix with ".tab" suffix
    if args.output is None:
        input_prefix = os.path.splitext(os.path.basename(args.input))[0]
        args.output = f"{input_prefix}.tab"

    species_list = parse_species_list(args.indv)
    process_maf(args.input, species_list, args.reference, args.output, args.stdout, args.out_freq, args.ref_gap, args.var, args.max_mis, args.g2n, args.n2g)
