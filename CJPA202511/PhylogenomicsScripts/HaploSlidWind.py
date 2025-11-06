#!/usr/bin/env python3
# -*- coding:utf-8 -*-

# Func: Split haploid file into non-overlapping sliding windows
# Usage: HaploSlidWind-v2.3.py [-h] -i INPUT_FILE -o PREFIX (-w SIZE | -s NUM) [-c CHR_LENGTH_FILE] [-H CHROM_HEADER] [--no-ref] [--non-int]
#                              [-f GFF_FILE | -I FILE] [--extract [FEATURE ...] | --sfilter [FEATURE ...] | --wfilter [FEATURE ...]]
#                              [--seqfilt BASE_NUM] [--sitefilt MIS_PROP] [--no-dec] [--ambig]
#                              [--min-len LENGTH] [--min-var NUM] [--min-pi NUM] [--max-mis MIS_PROP] [--min-seqn NUM] [--dist DISTANCE]
#                              [--keep-indv INDV_FILE | --remove-indv INDV_FILE]
#                              [--fasta] [--keep-ref] [--filt-blank]
#                              [--g2n | --n2g | --mis2quest | --mis2star]
# The first version is made by Wenjie DONG at 2024/3/21
# Author contact: dongwenjie@mail.kiz.ac.cn
# Institute: Kunming Institute of Zoology, CAS
## Update log:
##   v2.3 [2024/11/12]:
##         Added functions to keep specified individuals or remove some individuals when reading the input.
##         Optimized some parts.
##   v2.2 [2024/9/7]:
##         Added function to read feature intervals from a simple interval file instead of GFF file.
##         Added function to generate site windows instead of physical genomic windows.
##         Added option to keep the reference (or major) sequence in FASTA alignment.
##         Optimized some parts.
##   v2.1 [2024/8/11]:
##         Added a method for generating windows.
##         Fixed an issue where the number of output sequences was lower than the set minimum value.
##         Position is allowed to be a decimal (for the output file of maf2tab.py with option "--ref-gap").
##         User is allowed to set the starting column of samples.
##   v2.0 [2024/7/24]:
##         Added function to filter/extract data according to a GFF annotation file.
##         Added function to output FASTA format.
##         Added function to filter sequences and sites which have high missing proportion.
##         Reorganized arguments.


import sys
import argparse
import gzip
import time
from collections import defaultdict
from decimal import Decimal


def replace_char(txt, mode=None):
    if mode == "g2n":
        return txt.replace('-', 'N')
    elif mode == "n2g":
        return txt.replace('N', '-').replace('n', '-')
    elif mode == "mis2quest":
        return txt.replace('N', '?').replace('n', '?').replace('-', '?')
    elif mode == "mis2star":
        return txt.replace('N', '*').replace('n', '*').replace('-', '*')
    else:
        return txt


def read_individuals(indv_file):
    with open(indv_file) as indfile:
        indv_list = [i for i in indfile.read().splitlines() if i]

    return indv_list


def read_chromosome_lengths(chromosome_lengths_file):
    chromosome_lengths = {}
    print("## Read chromosome lengths ...", end='', flush=True)
    with open(chromosome_lengths_file, 'r') as f:
        for line in f:
            chromosome, length = line.strip().split('\t')
            chromosome_lengths[chromosome] = int(length)
    print(f" Done.\n## -- Chromosome number: {len(chromosome_lengths)}")

    return chromosome_lengths


def read_intervals(intv_file):
    intervals_dict = {}

    print(f"## Reading intervals from file '{intv_file}' ...", end='', flush=True)
    # Determine if the file is gzipped by checking its extension
    open_func = gzip.open if intv_file.endswith('.gz') else open

    with open_func(intv_file, 'rt') as file:
        line_num = 1
        for line in file:
            # Skip comments and empty lines
            if line.startswith('#') or not line.strip():
                continue

            parts = line.split()

            # Check file format
            if len(parts) != 3:
                print(f" Failed. Error info: the number of columns is not 3 (at line {line_num})\n -- Program terminated.")
                exit(1)
            else:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])

                if chrom not in intervals_dict:
                    intervals_dict[chrom] = [[start, end]]
                else:
                    intervals_dict[chrom].append([start, end])

                line_num += 1

    print(" Done.")

    # Sort intervals
    print("## Sorting intervals ...", end='', flush=True)
    for chrom in intervals_dict.keys():
        intervals_dict[chrom].sort(key=lambda x: x[0])
    print(" Done.")

    return intervals_dict


def read_feature_intervals(gff_file):
    intervals_dict = {}

    print(f"## Reading features from file '{gff_file}' ...", end='', flush=True)
    # Determine if the file is gzipped by checking its extension
    open_func = gzip.open if gff_file.endswith('.gz') else open

    with open_func(gff_file, 'rt') as file:
        for line in file:
            # Skip comments and empty lines
            if line.startswith('#') or not line.strip():
                continue

            parts = line.strip().split('\t')

            # Ensure the line has the correct number of columns
            if len(parts) != 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attributes = parts

            # Convert start and end positions to integers
            start = int(start)
            end = int(end)

            # Add the interval to the dictionary
            if feature not in intervals_dict:
                intervals_dict[feature] = {seqid: []}
            if seqid not in intervals_dict[feature]:
                intervals_dict[feature][seqid] = []
            intervals_dict[feature][seqid].append([start, end])
    print(" Done.")

    # Sort intervals
    print("## Sorting intervals of each features ...", end='', flush=True)
    for feature in intervals_dict.keys():
        for seqid in intervals_dict[feature].keys():
            intervals_dict[feature][seqid].sort(key=lambda x: x[0])
    print(" Done.")

    return intervals_dict


def merge_intervals(interval_dict_list):  # merge intervals by chromosomes
    allintv = defaultdict(list)
    intv_merged = defaultdict(list)

    print("## Merge intervals ...", end='', flush=True)
    # concatenate all interval lists for each chromosome and sort intervals by starting position
    for intv_dict in interval_dict_list:
        for seqid, intv in intv_dict.items():
            allintv[seqid].extend(intv)
    for seqid in allintv.keys():
        allintv[seqid].sort(key=lambda x: x[0])

    # merge overlapped intervals
    for seqid in allintv.keys():
        if len(allintv[seqid]) > 1:
            current_intv = [i for i in allintv[seqid][0]]
            i = 1
            while i < len(allintv[seqid]):
                while i < len(allintv[seqid]) and allintv[seqid][i][0] <= current_intv[1]:
                    if allintv[seqid][i][1] > current_intv[1]:
                        current_intv[1] = allintv[seqid][i][1]
                    i += 1
                intv_merged[seqid].append(current_intv)
                if i < len(allintv[seqid]):
                    current_intv = [i for i in allintv[seqid][i]]
                i += 1
            if current_intv not in intv_merged[seqid]:
                intv_merged[seqid].append(current_intv)
        else:
            intv_merged[seqid] = [i for i in allintv[seqid]]

    print(" Done.")

    return intv_merged


def feature_filter_site(pos, intervals, start_idx):  # start_idx is a list to pass parameters by reference
    if start_idx[0] >= len(intervals):
        return False
    else:
        if pos < intervals[start_idx[0]][0]:
            return False
        elif intervals[start_idx[0]][0] <= pos <= intervals[start_idx[0]][1]:
            return True
        else:
            while start_idx[0] < len(intervals) and pos > intervals[start_idx[0]][1]:
                start_idx[0] += 1
            if start_idx[0] >= len(intervals):
                return False
            else:
                if intervals[start_idx[0]][0] <= pos <= intervals[start_idx[0]][1]:
                    return True
                else:
                    return False


def feature_extract_site(pos, intervals, start_idx):
    if start_idx[0] >= len(intervals):
        return False
    else:
        if pos < intervals[start_idx[0]][0]:
            return True  # filter out sites
        elif intervals[start_idx[0]][0] <= pos <= intervals[start_idx[0]][1]:
            return False  # keep sites
        else:
            while start_idx[0] < len(intervals) and pos > intervals[start_idx[0]][1]:
                start_idx[0] += 1
            if start_idx[0] >= len(intervals):
                return False
            else:
                if intervals[start_idx[0]][0] <= pos <= intervals[start_idx[0]][1]:
                    return False
                else:
                    return True


def feature_filter_window(wind_start, wind_end, intervals, start_idx):
    if start_idx[0] >= len(intervals):
        return False
    else:
        if wind_end < intervals[start_idx[0]][0]:
            return False  # keep the window
        elif wind_start <= intervals[start_idx[0]][0] <= wind_end:
            return True  # filter out
        elif wind_start > intervals[start_idx[0]][1]:
            while start_idx[0] < len(intervals) and wind_start > intervals[start_idx[0]][1]:
                start_idx[0] += 1
            if start_idx[0] >= len(intervals):
                return False
            else:
                if wind_start <= intervals[start_idx[0]][0] <= wind_end:
                    return True
                else:
                    return False


def filter_ss(lines, ambig_as_mis, seqfilt, sitefilt, start_col):  # Filter sequence and sites based on degree of missing
    mischar = ['-', 'N', 'n', '?', '*', 'R', 'r', 'Y', 'y', 'M', 'm', 'K', 'k', 'S', 's', 'W', 'w', 'H', 'h', 'B', 'b', 'V', 'v', 'D', 'd'] if ambig_as_mis else ['-', 'N', 'n', '?', '*']
    filt_seqnum, filt_sitenum = 0, 0
    sqlen = len(lines)
    matrix = [i.strip('\n').split('\t') for i in lines]

    seqnum = len(matrix[0]) - start_col

    # Filter sequences by base number
    if seqfilt:
        # count base number for each sequence
        seq_basecount = [0 for i in matrix[0][start_col:]]
        for row in matrix:
            for i_col in range(seqnum):
                if row[i_col + start_col] not in mischar:
                    seq_basecount[i_col] += 1

        # change bases of filtered sequences to N
        seqidx_to_filt = [i for i in range(seqnum) if seq_basecount[i] < seqfilt]
        filt_seqnum = len(seqidx_to_filt)
        for i in range(sqlen):
            for idx in seqidx_to_filt:
                matrix[i][idx + start_col] = "N"

    # Filter sites by missing proportion
    if sitefilt is not None:
        # function to count missing number of a site
        def miscount(base_list, mischar):
            misnum = 0
            for mc in mischar:
                misnum += base_list.count(mc)
            return misnum

        # filter site whose missing proportion is greater than threshold value
        ori_sitenum = len(matrix)
        matrix = [line for line in matrix if miscount(line[start_col:], mischar) / seqnum <= sitefilt]
        filt_sitenum = ori_sitenum - len(matrix)

    # Convert matrix to data lines
    datalines = ['\t'.join(line) + '\n' for line in matrix]

    return datalines, filt_seqnum, filt_sitenum


def filter_window(lines, ambig_as_mis, minlen, minvar, minpi, maxmis, minseqnum, start_col):
    mischar = ['-', 'N', '?', '*', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D'] if ambig_as_mis else ['-', 'N', '?', '*']
    sqlen = len(lines)

    # Filter sequence length
    if sqlen < minlen:
        return True, f'len = {sqlen} bp < {minlen} bp'
    var_num, pi_num = 0, 0

    matrix = [i.upper().strip('\n').split('\t')[start_col:] for i in lines]

    # Calculate variable sites and parsimony sites numbers
    if minvar > 0 or minpi > 0:
        for char in matrix:
            all_bases = [i for i in char if i not in mischar]
            bases_catgry = set(all_bases)
            if len(bases_catgry) > 1:
                var_num += 1
            if len([b for b in bases_catgry if all_bases.count(b) > 1]) > 1:
                pi_num += 1

        # filter var. and pi site numbers
        if var_num < minvar:
            return True, f'var. num = {var_num} < {minvar}'
        if pi_num < minpi:
            return True, f'pi num = {pi_num} < {minpi}'

    # Count the sequence number with missing less than {maxgap}%
    if minseqnum > 0 and maxmis < 1:
        seq_num = len(matrix[0])
        # count missing number of each sequence
        mis_num = [0 for i in range(seq_num)]
        for char in matrix:
            for i in range(seq_num):
                if char[i] in mischar:
                    mis_num[i] += 1
        # count the number of sequences that meet the requirement
        seqnum_req = 0
        for n in mis_num:
            if n / sqlen < maxmis:
                seqnum_req += 1
        if seqnum_req < minseqnum:
            return True, f'number of sequences whose proportion of missing characters is less than {maxmis}: {seqnum_req} < {minseqnum}'

    return False, ''


def output_window(header, datalines, file_prefix, hap2fasta, filt_blank, minsqnum, msg, start_col, keep_ref):
    if not hap2fasta:  # output haploid format
        print(msg, flush=True)
        with open(f'{file_prefix}.haplo', 'w') as window_file:
            window_file.write(header)
            window_file.write(''.join(datalines))
    else:  # output fasta format
        # Convert haploid format to FASTA format
        if keep_ref:
            seqnames = header.strip().split('\t')[2:]
        else:
            seqnames = header.strip().split('\t')[start_col:]
        sequences = [[] for i in seqnames]
        for dtl in datalines:
            if keep_ref:
                bases = dtl.strip().split('\t')[2:]
            else:
                bases = dtl.strip().split('\t')[start_col:]
            for i in range(len(bases)):
                sequences[i].append(bases[i])
        fasta_content = []
        for i in range(len(sequences)):
            if filt_blank:
                # filter blank sequences
                filt = True
                for char in set(sequences[i]):
                    if char in ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']:
                        filt = False
                        break
            else:
                filt = False
            if not filt:
                fasta_content.append(f">{seqnames[i]}\n{''.join(sequences[i])}\n")

        if len(fasta_content) < minsqnum:
            sys.stderr.write(f'Filter out the window "{file_prefix.split(".")[1].replace("-", ":", 1)}"; REASON: sequence number {len(fasta_content)} < {minsqnum}.\n')
        else:
            print(msg, flush=True)
            with open(f'{file_prefix}.fasta', 'w') as window_file:
                window_file.write(''.join(fasta_content))


def sitenum_sliding_window(haplo_file, output_prefix, site_window_size, chromosome_lengths, start_col, individuals, header_chr, seqfilt,
                           sitefilt, no_dec, ambig_as_mis, minlen, minvar, minpi, maxmis, minsqnum, distance, repchar_mode,
                           output_fasta, keep_ref, filt_blank, featr_filt, intvs):
    """ Generate sliding windows based on site numbers """
    # Initialize variables
    current_chromosome = None
    current_window_start = 1
    current_window_end = 0
    last_output_window_end = -distance
    site_num = 0
    start_idx = [0]
    window_total_num, window_filtered_num, window_output_num = 0, 0, 0

    with gzip.open(haplo_file, 'rt') if haplo_file.endswith('.gz') else open(haplo_file, 'r') as f:
        header = ''
        indexes = []
        indv_mode = None
        tmp_store = []  # temporarily store the data lines
        for line in f:
            cols = line.strip().split('\t')
            chrom, pos, *bases = cols

            if line.startswith(f'{header_chr}\t'):  # get header
                if not header:
                    indv_names = bases[start_col - 2:]
                    if individuals:
                        indv_list, indv_mode = individuals
                        indexes = [indv_names.index(i) for i in indv_list]
                        indexes.sort()
                        if indv_mode == 0:  # remove mode
                            bases = bases[:start_col - 2] + [indv_names[i] for i in range(len(bases)) if i not in indexes]
                        elif indv_mode == 1:  # keep mode
                            bases = bases[:start_col - 2] + [indv_names[i] for i in indexes]
                    header = f'{chrom}\t{pos}\t' + '\t'.join(bases) + '\n'
                continue

            if no_dec and '.' in pos:  # skip decimcal positions (the position of a site in the deletion region in the reference is assigned a decimal value)
                continue

            pos_dec = Decimal(pos)
            pos = int(pos_dec)

            if indexes:  # keep or remove some individuals
                indv_bases = bases[start_col - 2:]
                if indv_mode == 0:  # remove mode
                    bases = bases[:start_col - 2] + [indv_bases[i] for i in range(len(indv_bases)) if i not in indexes]
                elif indv_mode == 1:  # keep mode
                    bases = bases[:start_col - 2] + [indv_bases[i] for i in indexes]

            # Write the data of the last window to file if window changed
            if tmp_store and (chrom != current_chromosome or site_num == site_window_size):
                window_total_num += 1
                filt = False
                if featr_filt[2] and current_chromosome in intvs:  # "feature filter window" mode
                    filt = feature_filter_window(current_window_start, current_window_end, intvs[current_chromosome], start_idx)
                if filt:
                    sys.stderr.write(f'Filter out the window "{current_chromosome}:{current_window_start}-{current_window_end}"; REASON: feature filter.\n')
                else:
                    tmp_store, filt_seqnum, filt_mis_sitenum = filter_ss(tmp_store, ambig_as_mis, seqfilt, sitefilt, start_col)
                    filt, msg = filter_window(tmp_store, ambig_as_mis, minlen, minvar, minpi, maxmis, minsqnum, start_col)
                    if filt:
                        sys.stderr.write(f'Filter out the window "{current_chromosome}:{current_window_start}-{current_window_end}"; REASON: {msg}.\n')
                    else:
                        window_filtered_num += 1
                        if current_window_start - last_output_window_end > distance:
                            # print messages
                            filt_seq_msg = f'Filter out {filt_seqnum} sequences -> ' if filt_seqnum > 0 else ''
                            filt_mis_site_msg = f'Filter out {filt_mis_sitenum} sites by missing -> ' if filt_mis_sitenum > 0 else ''
                            outputmsg = f'{filt_seq_msg}{filt_mis_site_msg}Output window {current_chromosome}:{current_window_start}-{current_window_end}.'

                            # write window to file
                            filepre = f'{output_prefix}.{current_chromosome}-{current_window_start}-{current_window_end}'
                            output_window(header, tmp_store, filepre, output_fasta, filt_blank, minsqnum, outputmsg, start_col, keep_ref)
                            window_output_num += 1
                            last_output_window_end = current_window_end
                tmp_store.clear()
                site_num = 0

            # Update current chromosome if changed
            if chrom != current_chromosome:
                current_chromosome = chrom
                current_window_start = pos
                current_window_end = 0
                last_output_window_end = -distance
                site_num = 0
                start_idx = [0]

            # Filter or extract sites by feature
            if featr_filt[0] and chrom in intvs:  # "extract site" mode
                filt = feature_extract_site(pos, intvs[chrom], start_idx)
            elif featr_filt[1] and chrom in intvs:  # "filter site" mode
                filt = feature_filter_site(pos, intvs[chrom], start_idx)
            else:
                filt = False

            if not filt:
                if chromosome_lengths and current_chromosome not in chromosome_lengths:
                    continue
                else:
                    # Set the start position of the new window as the position of the first site
                    if not tmp_store:
                        current_window_start = pos

                    # Append the data line to the temporary list
                    if site_num <= site_window_size:
                        # Change the gap or missing character if necessary
                        bases_txt = replace_char('\t'.join(bases), repchar_mode)
                        tmp_store.append(f"{chrom}\t{pos_dec}\t{bases_txt}\n")
                        site_num += 1
                        current_window_end = pos  # set the end position

        # Write the last window to file
        window_total_num += 1
        filt = False
        if featr_filt[2] and current_chromosome in intvs:  # "feature filter window" mode
            filt = feature_filter_window(current_window_start, current_window_end, intvs[current_chromosome], start_idx)
        if filt:
            sys.stderr.write(f'Filter out the window "{current_chromosome}:{current_window_start}-{current_window_end}"; REASON: feature filter.\n')
        else:
            tmp_store, filt_seqnum, filt_mis_sitenum = filter_ss(tmp_store, ambig_as_mis, seqfilt, sitefilt, start_col)
            filt, msg = filter_window(tmp_store, ambig_as_mis, minlen, minvar, minpi, maxmis, minsqnum, start_col)
            if filt:
                sys.stderr.write(f'Filter out the window "{current_chromosome}:{current_window_start}-{current_window_end}"; REASON: {msg}.\n')
            else:
                window_filtered_num += 1
                if current_window_start - last_output_window_end > distance:
                    # print messages
                    filt_seq_msg = f'Filter out {filt_seqnum} sequences -> ' if filt_seqnum > 0 else ''
                    filt_mis_site_msg = f'Filter out {filt_mis_sitenum} sites by missing -> ' if filt_mis_sitenum > 0 else ''
                    outputmsg = f'{filt_seq_msg}{filt_mis_site_msg}Output window {current_chromosome}:{current_window_start}-{current_window_end}.'

                    # write window to file
                    filepre = f'{output_prefix}.{current_chromosome}-{current_window_start}-{current_window_end}'
                    output_window(header, tmp_store, filepre, output_fasta, filt_blank, minsqnum, outputmsg, start_col, keep_ref)
                    window_output_num += 1

    print(f'\n### Summary ###\nTotal number of windows generated: {window_total_num}\nNumber of windows after filtering: {window_filtered_num}\nNumber of output windows: {window_output_num}\n###############', flush=True)


def genomic_sliding_window(haplo_file, output_prefix, window_size, chromosome_lengths, start_col, individuals, header_chr, non_int, seqfilt,
                   sitefilt, no_dec, ambig_as_mis, minlen, minvar, minpi, maxmis, minsqnum, distance, repchar_mode,
                   output_fasta, keep_ref, filt_blank, featr_filt, intvs):
    """ Generate sliding windows using physical window size """
    # Initialize variables
    current_chromosome = None
    current_window_start = 1
    current_window_end = window_size
    last_output_window_end = -distance
    start_idx = [0]
    window_total_num, window_filtered_num, window_output_num = 0, 0, 0

    with gzip.open(haplo_file, 'rt') if haplo_file.endswith('.gz') else open(haplo_file, 'r') as f:
        header = ''
        tmp_store = []  # temporarily store the data lines
        indexes = []
        for line in f:
            cols = line.strip().split('\t')
            chrom, pos, *bases = cols

            if line.startswith(f'{header_chr}\t'):  # get header
                if not header:
                    indv_names = bases[start_col - 2:]
                    if individuals:
                        indv_list, indv_mode = individuals
                        indexes = [indv_names.index(i) for i in indv_list]
                        indexes.sort()
                        if indv_mode == 0:  # remove mode
                            bases = bases[:start_col - 2] + [indv_names[i] for i in range(len(bases)) if i not in indexes]
                        elif indv_mode == 1:  # keep mode
                            bases = bases[:start_col - 2] + [indv_names[i] for i in indexes]
                    header = f'{chrom}\t{pos}\t' + '\t'.join(bases) + '\n'
                continue

            if no_dec and '.' in pos:  # skip decimcal positions (the position of a site in the deletion region in the reference is assigned a decimal value)
                continue
            
            pos_dec = Decimal(pos)
            pos = int(pos_dec)

            if indexes:  # keep or remove some individuals
                indv_bases = bases[start_col - 2:]
                if indv_mode == 0:  # remove mode
                    bases = bases[:start_col - 2] + [indv_bases[i] for i in range(len(indv_bases)) if i not in indexes]
                elif indv_mode == 1:  # keep mode
                    bases = bases[:start_col - 2] + [indv_bases[i] for i in indexes]

            # Write the data of the last window to file if window changed
            if tmp_store and (chrom != current_chromosome or pos > current_window_end):
                window_total_num += 1
                filt = False
                if featr_filt[2] and current_chromosome in intvs:  # "feature filter window" mode
                    filt = feature_filter_window(current_window_start, current_window_end, intvs[current_chromosome], start_idx)
                if filt:
                    sys.stderr.write(f'Filter out the window "{current_chromosome}:{current_window_start}-{current_window_end}"; REASON: feature filter.\n')
                else:
                    tmp_store, filt_seqnum, filt_mis_sitenum = filter_ss(tmp_store, ambig_as_mis, seqfilt, sitefilt, start_col)
                    filt, msg = filter_window(tmp_store, ambig_as_mis, minlen, minvar, minpi, maxmis, minsqnum, start_col)
                    if filt:
                        sys.stderr.write(f'Filter out the window "{current_chromosome}:{current_window_start}-{current_window_end}"; REASON: {msg}.\n')
                    else:
                        window_filtered_num += 1
                        if current_window_start - last_output_window_end > distance:
                            # print messages
                            filt_seq_msg = f'Filter out {filt_seqnum} sequences -> ' if filt_seqnum > 0 else ''
                            filt_mis_site_msg = f'Filter out {filt_mis_sitenum} sites by missing -> ' if filt_mis_sitenum > 0 else ''
                            outputmsg = f'{filt_seq_msg}{filt_mis_site_msg}Output window {current_chromosome}:{current_window_start}-{current_window_end}.'

                            # write window to file
                            filepre = f'{output_prefix}.{current_chromosome}-{current_window_start}-{current_window_end}'
                            output_window(header, tmp_store, filepre, output_fasta, filt_blank, minsqnum, outputmsg, start_col, keep_ref)
                            window_output_num += 1
                            last_output_window_end = current_window_end
                tmp_store.clear()

            # Update current chromosome if changed
            if chrom != current_chromosome:
                current_chromosome = chrom
                current_window_start = 1
                current_window_end = window_size
                last_output_window_end = -distance
                start_idx = [0]

            # Filter or extract sites by feature
            if featr_filt[0] and chrom in intvs:  # "extract site" mode
                filt = feature_extract_site(pos, intvs[chrom], start_idx)
            elif featr_filt[1] and chrom in intvs:  # "filter site" mode
                filt = feature_filter_site(pos, intvs[chrom], start_idx)
            else:
                filt = False

            if not filt:
                if current_chromosome in chromosome_lengths:
                    # Determine if the position is in the current window
                    if pos > current_window_end:
                        if non_int:  # free start position based on the first site in window
                            current_window_start = pos
                            current_window_end = pos + window_size - 1
                        else:  # fixed start position (start-1 is an integer multiple of window size)
                            n = int((pos + window_size - 1) / window_size)
                            current_window_start = (n - 1) * window_size + 1
                            current_window_end = n * window_size

                    # Check the validity of the window and append the data line to the temporary list
                    if current_window_start < chromosome_lengths[current_chromosome]:
                        if current_window_end > chromosome_lengths[current_chromosome]:
                            current_window_end = chromosome_lengths[current_chromosome]
                        # Change the gap or missing character if necessary
                        bases_txt = replace_char('\t'.join(bases), repchar_mode)
                        tmp_store.append(f"{chrom}\t{pos_dec}\t{bases_txt}\n")

        # Write the last window to file
        window_total_num += 1
        filt = False
        if featr_filt[2] and current_chromosome in intvs:  # "feature filter window" mode
            filt = feature_filter_window(current_window_start, current_window_end, intvs[current_chromosome], start_idx)
        if filt:
            sys.stderr.write(f'Filter out the window "{current_chromosome}:{current_window_start}-{current_window_end}"; REASON: feature filter.\n')
        else:
            tmp_store, filt_seqnum, filt_mis_sitenum = filter_ss(tmp_store, ambig_as_mis, seqfilt, sitefilt, start_col)
            filt, msg = filter_window(tmp_store, ambig_as_mis, minlen, minvar, minpi, maxmis, minsqnum, start_col)
            if filt:
                sys.stderr.write(f'Filter out the window "{current_chromosome}:{current_window_start}-{current_window_end}"; REASON: {msg}.\n')
            else:
                window_filtered_num += 1
                if current_window_start - last_output_window_end > distance:
                    # print messages
                    filt_seq_msg = f'Filter out {filt_seqnum} sequences -> ' if filt_seqnum > 0 else ''
                    filt_mis_site_msg = f'Filter out {filt_mis_sitenum} sites by missing -> ' if filt_mis_sitenum > 0 else ''
                    outputmsg = f'{filt_seq_msg}{filt_mis_site_msg}Output window {current_chromosome}:{current_window_start}-{current_window_end}.'

                    # write window to file
                    filepre = f'{output_prefix}.{current_chromosome}-{current_window_start}-{current_window_end}'
                    output_window(header, tmp_store, filepre, output_fasta, filt_blank, minsqnum, outputmsg, start_col, keep_ref)
                    window_output_num += 1

    print(f'\n### Summary ###\nTotal number of windows generated: {window_total_num}\nNumber of windows after filtering: {window_filtered_num}\nNumber of output windows: {window_output_num}\n###############', flush=True)


def main():
    examples = '''######################################################################################################\n\nExamples:\n
    1.Generate continuous 100kb genomic sliding windows and output windows in Haploid format.\n
        HaploSlidWind-v2.3.py -i data.haplo -o ./windows/data_slidwind -c chrom_len.txt -w 100000 > output.log\n
    2.Generate 10kb physical genomic windows and ensure that the window spacing is at least 50kb and output the window in FASTA format.
      The start positions of the windows are not integer multiples of the window size.
      Input file is gzip/bgzip compressed.\n
        HaploSlidWind-v2.3.py -i mydata.haplo.gz -o ./my_windows/data_10kw_sp50k -c chrlen.txt -w 10000 --dist 50000 --non-int --fasta > output.log\n
    3.Generate 5kb physical genomic windows and ensure that the window spacing is at least 25kb and output windows in FASTA format.
      Remove the sequences which have fewer than 1000 bases, and filter out blank windows.\n
        HaploSlidWind-v2.3.py -i data.haplo -o ./windows/data_sw5k_25k -c chrom_len.txt -w 5000 --dist 25000 --fasta --seqfilt 1000 --filt-blank > output.log 2>filter.log\n
    4.Generate continuous 50kb genomic windows and only keep the sites in non-coding region, and output windows in FASTA format.\n
        HaploSlidWind-v2.3.py -i data.haplo -o ./windows/data_50k_noncoding -c chrlen.txt -w 50000 -f ref_annot.gff --sfilter CDS --fasta > output.log 2> filter.log\n
    5.Generate 10kb physical genomic windows whose start positions are not integer multiples of the window size and output the window in Haploid format.
      Filter out the sites with the proportion of missing characters larger than 50%.
      Filter out the windows with sequence length less than 5000bp after site filtering.
      Filter out the sites where over 50% of sequences have more than 20% missing data (total sequence number is 50 for instance).\n
        HaploSlidWind-v2.3.py -i data.haplo -o ./windows/data_10k_filt -c chrlen.txt -w 10000 --non-int --sitefilt 0.5 --min-len 5000 --max-mis 0.2 --min-seqn 25 >output.log 2>filter.log\n
    6. Generate 5kb site windows and don't keep the final windows of each chromosome if their lengths are shorter than 5kb.\n
        HpaloSlidWind-v2.3.py -i data.haplo -o ./windows/data_5k_site -s 5000 --min-len 5000\n
    7. Generate 1kb physical genomic windows and only keep individuals specified in file indv.txt and write windows in FASTA format.\n
        HpaloSlidWind-v2.3.py -i mydata.haplo.gz -o ./1kb_windows_filter_indv/mydata_1kw_filt -c chrlen.txt -w 1000 --keep-indv indv.txt --fasta 1> out.log 2> filt.log\n
######################################################################################################
[ Updated in 2024/11/12 ]
****************************************************
~~
~~~~ Author: Wenjie Dong
~~~~~ Institute: Kunming Institute of Zoology, CAS
~~~~ Contact: dongwenjie@mail.kiz.ac.cn
~~
****************************************************'''
    parser = argparse.ArgumentParser(description="HaploSlidWind v2.3: Split haploid file into non-overlapping sliding windows.",
                                     epilog=examples, formatter_class=argparse.RawTextHelpFormatter)
    maingrp = parser.add_argument_group("main options")
    maingrp.add_argument("-i", "--input", type=str, required=True, metavar="INPUT_FILE", help="input haploid file (text file or gzip compressed (*.gz))")
    maingrp.add_argument("-o", "--output", type=str, required=True, metavar="PREFIX", help="output file prefix")
    windgrp = maingrp.add_mutually_exclusive_group(required=True)
    windgrp.add_argument("-w", "--window-size", type=int, metavar="SIZE", help="constant physical genomic window size")
    windgrp.add_argument("-s", "--site-number", type=int, metavar="NUM", help="generate windows with constant site number, so the physical window size is variable, as is the starting points of these windows ('--non-int' is set as default)")
    maingrp.add_argument("-c", "--chrom-len", type=str, metavar="CHR_LENGTH_FILE", help="chromosome lengths file (two columns: <chromosome> <length>; no header)")
    maingrp.add_argument("-H", "--chr-header", default="chr", metavar="CHROM_HEADER", help="the header of the first column (chromosome) [chr]")
    maingrp.add_argument("--no-ref", action="store_true", help="the reference (or major) column (column three) doesn't exist")
    maingrp.add_argument("--non-int", action="store_true", help="the window start position is the position of the first site in that window, otherwise, the start position minus 1 is an integer multiple of the window size")

    featrgrp = parser.add_argument_group("feature extract/filter")
    ffile_megrp = featrgrp.add_mutually_exclusive_group()
    ffile_megrp.add_argument("-f", "--feature", metavar="GFF_FILE", help="feature file (GFF format) for structure annotations")
    ffile_megrp.add_argument("-I", "--intervals", metavar="FILE", help="file contains intervals for data extraction or filtration (three columns: <chromosome> <start> <end>; no header)")
    fmegrp = featrgrp.add_mutually_exclusive_group()
    fmegrp.add_argument("--extract", nargs="*", dest="extr_featr", metavar="FEATURE", help="keep the sites belong to specified features (e.g. CDS) while creating windows")
    fmegrp.add_argument("--sfilter", nargs="*", dest="sfilter", metavar="FEATURE", help="filter out the sites belong to specified features while creating windows")
    fmegrp.add_argument("--wfilter", nargs="*", dest="wfilter", metavar="FEATURE", help="only keep the windows outside the regions containing specified features")

    sfiltgrp = parser.add_argument_group("sequence/site filter for missing")
    sfiltgrp.add_argument("--seqfilt", type=int, metavar="BASE_NUM", help="filter out the sequences whose base number (not missing) is fewer than this value")
    sfiltgrp.add_argument("--sitefilt", type=float, metavar="MIS_PROP", help="filter out the sites whose missing proportion is greater than this value")
    sfiltgrp.add_argument("--no-dec", action="store_true", help="filter out all decimal positions (for a site in the deletion region in the reference, the position is assigned a decimal value)")
    sfiltgrp.add_argument("--ambig", action="store_true", help="treat all ambiguous sites as missing")

    filtgrp = parser.add_argument_group("window filter")
    filtgrp.add_argument("--min-len", default=1, type=int, metavar="LENGTH", help="filter out the windows with matrix (alignment) length lower than this value [1]")
    filtgrp.add_argument("--min-var", default=0, type=int, metavar="NUM", help="filter out the windows with variable sites number lower than this value")
    filtgrp.add_argument("--min-pi", default=0, type=int, metavar="NUM", help="filter out the windows with parsimony informative sites number lower than this value")
    filtgrp.add_argument("--max-mis", default=0.99999, type=float, dest="maxmis", metavar="MIS_PROP", help="the maximum proportion of missing characters of a sequence [0.99999]")
    filtgrp.add_argument("--min-seqn", default=1, type=int, dest="minsqnum", metavar="NUM", help="the minimum number of sequence whose proportion of missing characters is less than the value provided by '--max-mis' [1]")
    filtgrp.add_argument("--dist", default=0, type=int, metavar="DISTANCE", help="the minimum distance between two adjacent output windows")

    idvfiltgrp = parser.add_argument_group("individual filter")
    idvmegrp = idvfiltgrp.add_mutually_exclusive_group()
    idvmegrp.add_argument("--keep-indv", metavar="INDV_FILE", help="only the individuals specified in a text file are retained when reading the input")
    idvmegrp.add_argument("--remove-indv", metavar="INDV_FILE", help="remove the individuals specified in a text file when reading the input")

    fastagrp = parser.add_argument_group("FASTA output settings")
    fastagrp.add_argument("--fasta", action="store_true", help="output FASTA format instead of haploid format")
    fastagrp.add_argument("--keep-ref", action="store_true", help="keep the reference (or major) sequence in alignment")
    fastagrp.add_argument("--filt-blank", action="store_true", help="don't output blank sequences")

    misrepgrp = parser.add_argument_group("missing character replacement")
    megrp = misrepgrp.add_mutually_exclusive_group()
    megrp.add_argument("--g2n", action="store_true", help="replace gaps with 'N'")
    megrp.add_argument("--n2g", action="store_true", help="replace 'N' with gaps")
    megrp.add_argument("--mis2quest", action="store_true", help="replace missing characters ('N','-') with '?'")
    megrp.add_argument("--mis2star", action="store_true", help="replace missing characters ('N','-') with '*'")

    args = parser.parse_args()

    featr_filt = [False, False, False]  # switch off: extr_featr, sfilter, wfilter

    # Check if chromosome length file is provided if genomic window size is used
    if args.window_size and not args.chrom_len:
        parser.error(f"chromosome length file ('--chrom-len') must be provided if using genomic windows ('--window_size')")

    # Check argument conflict for site-window mode
    if args.site_number and args.min_len > args.site_number:
        parser.error(f"'--min-len {args.min_len}': minimum alignment length should not be larger than the window size of sites ('--site-number {args.site_number}')")

    # Set start column of data matrix
    if args.no_ref:
        start_col = 2
        if args.keep_ref:
            parser.error("arguments conflict: '--no-ref' and '--keep-ref'")
    else:
        start_col = 3

    # Print start time and command line
    start_time, start_struct_time = time.time(), time.localtime()
    print(f"## Start time: {time.strftime('%Y-%m-%d %H:%M:%S', start_struct_time)}")
    print("## CMD: " + " ".join(sys.argv))

    # Set the mode for replacing missing characters
    if args.g2n:
        repchar_mode = "g2n"
    elif args.n2g:
        repchar_mode = "n2g"
    elif args.mis2quest:
        repchar_mode = "mis2quest"
    elif args.mis2star:
        repchar_mode = "mis2star"
    else:
        repchar_mode = None

    # Read chromosome lengths
    if args.chrom_len:
        chrom_len = read_chromosome_lengths(args.chrom_len)
    else:
        chrom_len = None

    # Read the individual file if one is specified
    if args.keep_indv:
        indvs = [read_individuals(args.keep_indv), 1]  # "1": keep mode
    elif args.remove_indv:
        indvs = [read_individuals(args.remove_indv), 0]  # "0": remove mode
    else:
        indvs = None

    # If filter by features, get interval lists of all chromosomes
    if args.feature:  # read intervals from a GFF file
        featr_intv = read_feature_intervals(args.feature)
        if not featr_intv:
            sys.stderr.write("<ERROR> No features were obtained. Please check your GFF file.\n")
            exit(1)
    elif args.intervals:  # read intervals from a simple interval file
        featr_intv = read_intervals(args.intervals)
        if not featr_intv:
            sys.stderr.write("<ERROR> No intervals were obtained. Please check your interval file.\n")
            exit(1)
    else:
        featr_intv = None

    # Merge intervals if multiple features were specified
    merged_intvs = None
    if featr_intv:
        if args.feature:  # GFF file
            featr_filt = [args.extr_featr, args.sfilter, args.wfilter]
            if args.extr_featr or args.sfilter or args.wfilter:
                featr_names = [i for i in featr_filt if i][0]
                print(f"## Specified features: {', '.join(featr_names)}")
                if len(featr_names) > 1:
                    specified_featr_intv = [featr_intv[name] for name in featr_names]
                    merged_intvs = merge_intervals(specified_featr_intv)
                else:
                    merged_intvs = featr_intv[featr_names[0]]
            else:
                sys.stderr.write("<WARNING> No features were specified. Nothing to do with the given GFF file.\n")
        elif args.intervals:  # interval file
            featr_filt = ["--extract" in sys.argv, "--sfilter" in sys.argv, "--wfilter" in sys.argv]
            if featr_filt[0] or featr_filt[1] or featr_filt[2]:
                print("## Specified features: NONE -> all intervals")
                merged_intvs = featr_intv
            else:
                sys.stderr.write("<WARNING> No operation was specified (--extract/--sfilter/--wfilter). Nothing to do with the given intervals.\n")

    # Print window settings
    if args.window_size:
        print(f"## Generate physical genomic windows with constant size: {args.window_size} bp")
    else:
        print(f"## Generate site windows with site quantity : {args.site_number} bp")
    if args.dist:
        print(f"## The interval between windows: {args.dist} bp")

    # Print individual settings
    if args.keep_indv:
        print(f"## Only keep these individuals [num: {len(indvs[0])}]: {', '.join(indvs[0])}")
    elif args.remove_indv:
        print(f"## Remove these individuals [num: {len(indvs[0])}]: {', '.join(indvs[0])}")

    # Print filter settings
    if args.seqfilt:
        print(f"## Sequence filter: keep sequences with length >= {args.seqfilt}")
    if args.sitefilt is not None:
        if args.sitefilt > 1:
            args.sitefilt /= 100
        print(f"## Site filter: keep sites with missing proportion <= {args.sitefilt * 100}%")

    # Print output settings
    if args.fasta:
        print("## Output format: FASTA")
    else:
        print("## Output format: Haploid")

    # Start generating sliding windows
    print("#### Start making sliding windows ####")
    if args.window_size:
        genomic_sliding_window(args.input, args.output, args.window_size, chrom_len, start_col, indvs, args.chr_header, args.non_int, args.seqfilt,
                       args.sitefilt, args.no_dec, args.ambig, args.min_len, args.min_var, args.min_pi, args.maxmis, args.minsqnum,
                       args.dist, repchar_mode, args.fasta, args.keep_ref, args.filt_blank, featr_filt, merged_intvs)
    else:
        sitenum_sliding_window(args.input, args.output, args.site_number, chrom_len, start_col, indvs, args.chr_header, args.seqfilt,
                       args.sitefilt, args.no_dec, args.ambig, args.min_len, args.min_var, args.min_pi, args.maxmis, args.minsqnum,
                       args.dist, repchar_mode, args.fasta, args.keep_ref, args.filt_blank, featr_filt, merged_intvs)

    # Print end time
    end_time, end_strut_time = time.time(), time.localtime()
    time_intv = end_time - start_time
    m, s = divmod(time_intv, 60)
    h, m = divmod(m, 60)
    print(f"## Program finished at {time.strftime('%Y-%m-%d %H:%M:%S', end_strut_time)}")
    print(f"## Running time: {h:.0f}:{m:.0f}:{s:.2f}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("-- Program terminated by user.")
