#!/usr/bin/env python3
# -*- coding:utf-8 -*-

# Func: Process DFOIL output files and generate introgression summary
# Usage: python3 dfoil_summary.py [-h] (-f FILES [FILES ...] | -l FILE_LIST) [-o OUTPUT] [--output-tab] [--path PATH]
# By Wenjie DONG [2024/10/9]

import argparse
import os
import subprocess
import sys
import shutil


def run_dfoil_analyze(input_file, dfoil_path=None):
    if dfoil_path:
        if not os.path.isfile(dfoil_path):
            raise FileNotFoundError(f"The specified path '{dfoil_path}' is not a valid file.")
        cmd = [dfoil_path, input_file]
    else:
        dfoil_script = 'dfoil_analyze.py'
        # Check if the script is available in the environment path or current working directory
        script_path = shutil.which(dfoil_script) or os.path.abspath(dfoil_script)
        if not os.path.isfile(script_path):
            raise FileNotFoundError(f"{dfoil_script} not found in PATH or current directory. Please specify the --path.")
        sys.stderr.write(f"{input_file} > {script_path}")
        cmd = [script_path, input_file]
    
    # Capture the output of dfoil_analyze.py
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout


def extract_introgression_calls(dfoil_output):
    introgression_start = dfoil_output.find("# Introgression Calls:")
    if introgression_start == -1:
        raise ValueError("Introgression Calls section not found.")
    
    introgression_lines = dfoil_output[introgression_start:].strip().split("\n")[1:]
    introgression_calls = {}
    
    for line in introgression_lines:
        if line.strip():
            call_type, count = line.split()
            introgression_calls[call_type] = int(count)
    
    sys.stderr.write(" -> Extract introgression numbers")
    
    return introgression_calls


def extract_summary_statistics(dfoil_output):
    summary_start = dfoil_output.find("# DFOIL component summary statistics:")
    if summary_start == -1:
        raise ValueError("DFOIL component summary statistics section not found.")
    
    summary_end = dfoil_output.find("# Introgression Calls:")

    summary_lines = [i for i in dfoil_output[summary_start:summary_end].strip().split("\n")[1:] if i]
    
    return "\n".join(summary_lines)


def write_introgression_table(file_data, output_file):
    introgression_types = ['none', '13', '14', '23', '24', '31', '41', '32', '42', '123', '124']
    
    with open(output_file, 'w') if output_file else sys.stdout as f_out:
        # Write header
        f_out.write("file\t" + "\t".join(introgression_types) + "\n")
        
        # Write data for each file
        for file_name, calls in file_data.items():
            row = [file_name] + [str(calls.get(t, 0)) for t in introgression_types]
            f_out.write("\t".join(row) + "\n")


def write_summary_statistics(file_name, summary):
    output_path = os.path.splitext(file_name)[0] + ".dfoil_summary_stat.tsv"
    with open(output_path, 'w') as f_out:
        f_out.write(summary)
    sys.stderr.write(f" -> Write summary statistics to {output_path}\n")


def main():
    parser = argparse.ArgumentParser(description="Process DFOIL output files and generate introgression summary.")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-f", "--files", nargs='+', help="list of input files")
    group.add_argument("-l", "--file-list", help="file containing a list of input files")
    
    parser.add_argument("-o", "--output", help="output file for introgression summary table (default: stdout)")
    parser.add_argument("--output-tab", action='store_true', help="output the DFOIL component summary statistics")
    parser.add_argument("--path", help="path to dfoil_analyze.py")
    
    args = parser.parse_args()

    # Get the list of input files
    if args.files:
        input_files = args.files
    elif args.file_list:
        with open(args.file_list, 'r') as f:
            input_files = [line.strip() for line in f if line.strip()]
    
    # Store introgression data for all files
    introgression_data = {}
    
    for input_file in input_files:
        try:
            # Run dfoil_analyze and capture output
            dfoil_output = run_dfoil_analyze(input_file, args.path)
            
            # Extract Introgression Calls
            introgression_calls = extract_introgression_calls(dfoil_output)
            introgression_data[input_file] = introgression_calls
            
            # If --output-tab is set, also extract and write summary statistics
            if args.output_tab:
                summary_statistics = extract_summary_statistics(dfoil_output)
                write_summary_statistics(input_file, summary_statistics)
        except FileNotFoundError as e:
            print(e)
            sys.exit(1)
        except Exception as e:
            print(f"Error processing {input_file}: {e}", file=sys.stderr)
    
    # Write introgression summary table
    write_introgression_table(introgression_data, args.output)
    print(f"Write introgression counts table to file {args.output}.", file=sys.stderr)


if __name__ == "__main__":
    main()
