#!/bin/bash

## Func: call fasta2dfoil.py in parallel for multiple sample groups and fasta files
## Usage: bash fasta2dfoil_parallel.sh output_dir sample_groups_file fasta_list_file outgroup max_parallel_tasks
## By Wenjie DONG [2024/9/12]

# Show usage if the number of arguments doesn't match the requirement
[[ $# -ne 5 ]] && echo "usage: bash $0 <output_dir> <sample_groups_file> <fasta_list_file> <outgroup> <max_parallel_tasks>" && exit

outdir=$1            # First argument: output directory
sample_grps=$2       # Second argument: list file with four-sample groups at each line
fasta_list=$3        # Third argument: list file with paths to all FASTA files
outgroup=$4          # Fourth argument: outgrou sequence name
max_cpu=$5           # Fifth argument: maximum number of parallel tasks

mkdir -p $outdir

# Initialize the number of currently running tasks to 0
running_tasks=0

# Read each line from the sample groups file
cat $sample_grps | while read -r line; do
    # Create a task and run it in the background
    (
        # Initialize the index for the FASTA file sequence
        fasta_index=0

        # Create the merged file for the current line and write the header
        merged_file="${outdir}/$(echo "$line" | tr " " ".").merged_counts.tsv"
        header_written=false
        > $merged_file   # Clear or create the merged file

        # Read each FASTA file from the fasta_list
        cat $fasta_list | while read -r fasta; do
            fasta_prefix=$(basename "$fasta" .fasta)   # Get the prefix of the FASTA file (without extension)

            # Create a unique temporary file for this process
            tmp_counts=$(mktemp)

            # Generate the counts file
            fasta2dfoil.py <(selectSeqs.pl -f <(echo "$line" | tr " " $"\n"; echo "$outgroup") "$fasta") \
            --out ${outdir}/$(echo "$line" | tr " " ".").${fasta_index}.counts \
            --names $(echo "$line" | tr " " ","),"$outgroup"

            # Get the path to the generated counts file
            counts_file="${outdir}/$(echo "$line" | tr " " ".").${fasta_index}.counts"
            
            # Modify the counts file: replace the chromosome column with the FASTA prefix and position with the sequence index
            awk -v prefix="$fasta_prefix" -v pos="$fasta_index" 'NR==1 {print $0} NR>1 {$1=prefix; $2=pos; print}' $counts_file > "$tmp_counts"

            # Merge into the final file, keeping only one header
            if [ "$header_written" = false ]; then
                cat "$tmp_counts" >> $merged_file
                header_written=true
            else
                grep -v "^#" "$tmp_counts" >> $merged_file   # Skip the header
            fi

            # Delete the processed counts file and temporary file
            rm -f "$counts_file" "$tmp_counts"
            
            # Increment the FASTA file index
            (( fasta_index++ ))

        done
    ) &
    
    # Increment the number of background tasks
    (( running_tasks++ ))

    # If the number of running tasks reaches the maximum allowed, wait for all background tasks to finish
    if [[ "$running_tasks" -ge "$max_cpu" ]]; then
        wait
        running_tasks=0  # Reset task counter
    fi
done

# Ensure all tasks are completed
wait
