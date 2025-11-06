#!/usr/bin/env python3
# -*- coding:utf-8 -*-

# Func: HaploExtrMRNA: Extract mRNA sequences from a Haploid file.
# Usage: HaploExtrMRNA.py [-h] -i INPUT -g GFF [-o OUTDIR] [-p PREFIX] [--all]
#                         [--fasta] [--no-ref] [--min-len MIN_LEN]
#                         [--min-indv MIN_INDV] [--max-mis MAX_MIS]
#                         [--min-var MIN_VAR] [--min-pi MIN_PI]
#                         [--seqmask SEQMASK]
# By Wenjie DONG [2024/8/20]


import os
import sys
import argparse
import gzip
from datetime import datetime
from decimal import Decimal


def open_file(filename, mode):
    """Open uncompressed or 'gz' compressed file"""
    if filename.endswith(".gz"):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)


def parse_gff(gff_file):
    """Parse the GFF file to extract sorted mRNA and CDS information by chromosomes."""
    mRNA_info = {}
    with open_file(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            if fields[2] == "mRNA":
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                attributes = fields[8]
                mRNA_id = attributes.split("ID=")[-1].split(";")[0]
                if chrom not in mRNA_info:
                    mRNA_info[chrom] = []
                mRNA_info[chrom].append({"id": mRNA_id, "start": start, "end": end, "cds": []})
            elif fields[2] == "CDS":
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                attributes = fields[8]
                parent_id = attributes.split("Parent=")[-1].split(";")[0]
                for mRNA in mRNA_info.get(chrom, []):
                    if mRNA["id"] == parent_id:
                        mRNA["cds"].append((start, end))
                        break

    # Sort mRNA by start position within each chromosome
    for chrom in mRNA_info:
        mRNA_info[chrom].sort(key=lambda x: x["start"])

    return mRNA_info


def do_seqmask(mRNA_block, min_basenum, no_ref):
    """Mask the sequences whose base number is lower than the threshold value."""
    basenum = [0] * len(mRNA_block[0][2:]) if no_ref else [0] * len(mRNA_block[0][3:])
    ambiguous = ["R", "Y", "M", "K", "S", "W", "H", "B", "V", "D", "N", "-", "?"]

    # Count base number of each sequence
    for pos in mRNA_block:
        bases = pos[2:] if no_ref else pos[3:]
        for i, base in enumerate(bases):
            if base not in ambiguous:
                basenum[i] += 1

    # Get the indexes of sequences to be masked
    mask_index = [i + 2 for i, num in enumerate(basenum) if num < min_basenum] if no_ref else [i + 3 for i, num in enumerate(basenum) if num < min_basenum]

    # Do hard mask
    for i in range(len(mRNA_block)):
        for j in mask_index:
            mRNA_block[i][j] = "N"

    return len(mask_index)


def check_metrics(mRNA_block, num_individuals, min_indv, max_mis, min_var, min_pi, no_ref):
    """Calculate metrics for the extracted data to determine whether to export the mRNA."""
    sequence_length = len(mRNA_block)

    # Initialize variables
    non_missing_counts = [0] * num_individuals
    variant_counts = 0
    informative_counts = 0
    total_bases = sequence_length * num_individuals
    total_missing = 0

    ambiguous = ["R", "Y", "M", "K", "S", "W", "H", "B", "V", "D", "N", "-", "?"]

    for pos in mRNA_block:
        bases = pos[2:] if no_ref else pos[3:]
        alleles = [base for base in bases if base not in ambiguous]
        # Determine whether the site is a variable site or parsimony informative site
        if len(set(alleles)) > 1:
            variant_counts += 1
            if sum(1 for a in set(alleles) if alleles.count(a) > 1) > 1:  # is parsimony informative site
                informative_counts += 1
        # Count non-missing bases and missing data
        for i, base in enumerate(bases):
            if base not in ["N", "-", "?"]:
                non_missing_counts[i] += 1
            else:
                total_missing += 1

    indv_num = sum(1 for count in non_missing_counts if count > 0)
    mis_prop = total_missing / total_bases

    if min_indv and indv_num < min_indv:
        return f"sequence number: {indv_num} < {min_indv}", False

    if max_mis and mis_prop > max_mis:
        return f"total missing proportion: {mis_prop} > {max_mis}", False

    if min_var and variant_counts < min_var:
        return f"variant number: {variant_counts} < {min_var}", False

    if min_pi and informative_counts < min_pi:
        return f"parsimony informative site number: {informative_counts} < {min_pi}", False

    return "", True


def extract_mRNA(haploid_file, mRNA_info, file_prefix, output_dir, include_all, no_ref, min_indv, min_len, max_mis, min_var, min_pi, seqmask, sitefilt_mis, fasta_output):
    """Extract mRNA sequences from the Haploid file and do the filter and output."""
    outputnum, filtnum, mask_num, site_filtnum = 0, 0, 0, 0
    with open_file(haploid_file, 'rt') as f:
        header = f.readline().strip().split("\t")
        individuals = header[2:] if no_ref else header[3:]

        mrna_pointer = {chrom: 0 for chrom in mRNA_info}
        current_mRNA = {chrom: mRNA_info[chrom][0] if mRNA_info[chrom] else None for chrom in mRNA_info}
        mRNA_block = []

        for line in f:
            fields = line.strip().split("\t")
            chrom = fields[0]
            decimal_pos = Decimal(fields[1])
            position = int(decimal_pos)
            ref_base = None if no_ref else fields[2]
            bases = [b.upper() for b in fields[2:]] if no_ref else [b.upper() for b in fields[3:]]

            if chrom not in mRNA_info:
                continue

            mRNA = current_mRNA[chrom]

            if mRNA is None:
                continue

            # Finished extracting one mRNA and determine whether to output
            if position > mRNA["end"]:
                if mRNA_block:
                    sitefilt_msg = f"Filter out {site_filtnum} sites -> " if site_filtnum else ""
                    if len(mRNA_block) >= min_len:
                        if seqmask:
                            mask_num = do_seqmask(mRNA_block, seqmask, no_ref)  # mask sequences
                        msg, keep = check_metrics(mRNA_block, len(individuals), min_indv, max_mis, min_var, min_pi, no_ref)
                        mask_msg = f"Masked {mask_num} sequences -> " if mask_num else ""
                        if keep:  # don't filter the mRNA out
                            filename = f"{output_dir}/{file_prefix}.mRNA_{mRNA['id']}.{chrom}-{mRNA['start']}-{mRNA['end']}{'' if include_all else '.CDS'}.{'fasta' if fasta_output else 'haplo'}"
                            print(f"{sitefilt_msg}{mask_msg}Output mRNA {mRNA['id']} to {filename} ... ", end='', flush=True)
                            write_output(mRNA_block, individuals, header, filename, fasta_output, no_ref)
                            print("Done.")
                            outputnum += 1
                        else:
                            sys.stderr.write(f"{sitefilt_msg}{mask_msg}Filter mRNA {mRNA['id']} ({chrom}:{mRNA['start']}-{mRNA['end']}); [REASON] {msg}.\n")
                            filtnum += 1
                    else:
                        sys.stderr.write(f"{sitefilt_msg}Filter mRNA {mRNA['id']} ({chrom}:{mRNA['start']}-{mRNA['end']}); [REASON] sequence length: {len(mRNA_block)} < {min_len}.\n")
                        filtnum += 1
                elif site_filtnum:  # all site was filtered out
                    sitefilt_msg = f"Filter out {site_filtnum} sites -> " if site_filtnum else ""
                    sys.stderr.write(f"{sitefilt_msg}Filter mRNA {mRNA['id']} ({chrom}:{mRNA['start']}-{mRNA['end']}); [REASON] sequence length: {len(mRNA_block)} < {min_len}.\n")
                    filtnum += 1

                # clear temporary storage and move the pointer to the next mRNA
                mRNA_block = []
                mask_num = 0
                site_filtnum = 0
                while mRNA and position > mRNA["end"]:
                    mrna_pointer[chrom] += 1
                    if mrna_pointer[chrom] < len(mRNA_info[chrom]):
                        mRNA = mRNA_info[chrom][mrna_pointer[chrom]]
                    else:
                        mRNA = None
                    current_mRNA[chrom] = mRNA

            # Collect sites of current mRNA and store the data to the temporary storage
            if mRNA and mRNA["start"] <= position <= mRNA["end"]:
                if include_all or any(start <= position <= end for start, end in mRNA["cds"]):
                    # filter out the sites with much missing data
                    if sitefilt_mis and sum(1 for b in bases if b in ["N", "-", "?"]) / len(individuals) > sitefilt_mis:  # calculate missing proportion
                        site_filtnum += 1
                    else:
                        if no_ref:
                            mRNA_block.append([chrom, str(decimal_pos)] + bases)
                        else:
                            mRNA_block.append([chrom, str(decimal_pos), ref_base] + bases)

        # Processing the final mRNA block if it still in temporary storage
        if mRNA_block:
            sitefilt_msg = f"Filter out {site_filtnum} sites -> " if site_filtnum else ""
            if len(mRNA_block) >= min_len:
                if seqmask:
                    mask_num = do_seqmask(mRNA_block, seqmask, no_ref)
                msg, keep = check_metrics(mRNA_block, len(individuals), min_indv, max_mis, min_var, min_pi, no_ref)
                mask_msg = f"Masked {mask_num} sequences -> " if mask_num else ""
                if keep:
                    filename = f"{output_dir}/{file_prefix}.mRNA_{mRNA['id']}.{chrom}-{mRNA['start']}-{mRNA['end']}{'' if include_all else '.CDS'}.{'fasta' if fasta_output else 'haplo'}"
                    print(f"{sitefilt_msg}{mask_msg}Output mRNA {mRNA['id']} to {filename} ... ", end='', flush=True)
                    write_output(mRNA_block, individuals, header, filename, fasta_output, no_ref)
                    print("Done.")
                    outputnum += 1
                else:
                    sys.stderr.write(f"{sitefilt_msg}{mask_msg}Filter mRNA {mRNA['id']} ({chrom}:{mRNA['start']}-{mRNA['end']}); [REASON] {msg}.\n")
                    filtnum += 1
            else:
                sys.stderr.write(f"{sitefilt_msg}Filter mRNA {mRNA['id']} ({chrom}:{mRNA['start']}-{mRNA['end']}); [REASON] sequence length: {len(mRNA_block)} < {min_len}.\n")
                filtnum += 1
        elif site_filtnum:
            sitefilt_msg = f"Filter out {site_filtnum} sites -> " if site_filtnum else ""
            sys.stderr.write(f"{sitefilt_msg}Filter mRNA {mRNA['id']} ({chrom}:{mRNA['start']}-{mRNA['end']}); [REASON] sequence length: {len(mRNA_block)} < {min_len}.\n")
            filtnum += 1

    return outputnum, filtnum


def write_output(mRNA_data, individuals, header, filename, fasta_output, no_ref):
    """Write mRNA block to file in Haploid or FASTA format."""
    with open(filename, 'w') as output_file:
        if fasta_output:
            for i, individual in enumerate(individuals):
                sequence = ''.join(pos[i + 2] for pos in mRNA_data) if no_ref else ''.join(pos[i + 3] for pos in mRNA_data)
                if sequence.count("N") + sequence.count("-") + sequence.count("?") < len(sequence):  # Exclude fully missing sequences
                    output_file.write(f">{individual}\n{sequence}\n")
        else:
            # write header line
            output_file.write("\t".join(header) + "\n")
            # write sites
            for pos in mRNA_data:
                output_file.write("\t".join(pos) + "\n")


def main():
    parser = argparse.ArgumentParser(description="HaploExtrMRNA v1.0: Extract mRNA sequences from a Haploid file. Provided by Wenjie DONG [2024/8/20].", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", required=True, help="input Haploid file (uncompressed or 'gz' compressed)")
    parser.add_argument("-g", "--gff", required=True, help="input GFF annotation file (uncompressed or 'gz' compressed)")
    parser.add_argument("-o", "--outdir", default=".", help="output directory (default: current directory)")
    parser.add_argument("-p", "--prefix", help="prefixes of output files (default: the prefix of input Haploid file)")
    parser.add_argument("--all", action="store_true", help="keep all positions within mRNA, including noncoding regions")
    parser.add_argument("--fasta", action="store_true", help="output sequences in FASTA format instead of Haploid format")
    parser.add_argument("--no-ref", action="store_true", help="input Haploid file does not include reference (or 'major') column")

    filter_group = parser.add_argument_group("mRNA filter")
    filter_group.add_argument("--min-len", type=int, default=1, help="minimum number of positions required")
    filter_group.add_argument("--min-indv", type=int, default=1, help="minimum number of non-blank individuals required")
    filter_group.add_argument("--max-mis", type=float, help="maximum allowed total missing data proportion")
    filter_group.add_argument("--min-var", type=int, help="minimum number of variant positions required")
    filter_group.add_argument("--min-pi", type=int, help="minimum number of parsimony-informative sites required")

    sfilter_group = parser.add_argument_group("sequence/site filter")
    sfilter_group.add_argument("--seqmask", type=int, help="do hard mask to the sequences whose unambiguous base number is fewer than this value")
    sfilter_group.add_argument("--sitefilt", type=float, help="filter out bases whose missing proportion is higher than this value")

    args = parser.parse_args()

    # Set prefixes of output files
    prefix = args.prefix if args.prefix else os.path.splitext(os.path.basename(args.input))[0]

    start_time = datetime.now()
    print(f"## HaploExtrMRNA v1.0: Extract mRNA sequences from a Haploid file.\n## START RUN at [{start_time.strftime('%Y-%m-%d %H:%M:%S')}]")
    print(f"## Input Haploid file: {args.input}")
    print(f"## Input GFF annotation file: {args.gff}")
    print(f"## Output prefix: {prefix}")
    print(f"## Output directory: {args.outdir}")
    print(f"## Output format: {'FASTA' if args.fasta else 'Haploid'}")
    print(f"## Extract mode: {'all (CDS + non-coding)' if args.all else 'CDS'}")
    print("## Output log -> standard out; Filter log -> standard error")
    print("##### START PROCESSING #####")

    print("--> reading gff file ... ", end='', flush=True)
    mRNA_info = parse_gff(args.gff)
    total_mrna_num = sum(len(mrnas) for mrnas in mRNA_info.values())
    print("Done.")

    # Create the output directory if it doesn't exist
    print("--> creating output directory ... ", end='', flush=True)
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
        print("Done.")
    else:
        print("Already exists!")

    # Extract CDS sequences from the Haploid file
    print("--> start extracting mRNAs\n### OUTPUT LOG ###\n------------------")
    output_num, filtered_num = extract_mRNA(args.input, mRNA_info, prefix, args.outdir, args.all, args.no_ref,
                          args.min_indv, args.min_len, args.max_mis, args.min_var, args.min_pi, args.seqmask, args.sitefilt, args.fasta)
    print("------------------\n###  END LOG  ###")

    # Calculate running time
    end_time = datetime.now()
    timedelta = end_time - start_time
    hours, remainder = divmod(timedelta.seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"##### END PROCESSING #####\n## END at [{end_time.strftime('%Y-%m-%d %H:%M:%S')}]")
    print(f"## Running time: {hours}:{minutes}:{seconds}")

    # Make summary
    print("\n############### SUMMARY ###############")
    print(f"* Total mRNA number in GFF file: {total_mrna_num}")
    print(f"* Covered mRNA number: {output_num + filtered_num}")
    print(f"* Output mRNA number: {output_num}")
    print(f"* Filtered out mRNA number: {filtered_num}\n#######################################\n")


if __name__ == "__main__":
    main()
