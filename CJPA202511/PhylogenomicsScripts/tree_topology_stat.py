#!/usr/bin/env python3

# Func: Count and extract topologies from gene trees.
# Usage: tree_topology_stat.py [-h] -i INPUT -o OUTPRE [-n TOP_N] [-s {count,weight}]
# By Wenjie DONG [2024/12/16]


import argparse
from collections import defaultdict
from Bio import Phylo
from io import StringIO
import re
import warnings


def normalize_topology(tree):
    """
    Removing branch lengths and node labels (e.g., support values), sorting tip labels for consistent topology.
    """
    # Remove branch lengths and node labels
    for clade in tree.find_clades():
        clade.branch_length = None
        clade.confidence = None  # Remove node labels (e.g., support values)

    # Convert tree back to a Newick string
    newick_str = StringIO()
    Phylo.write(tree, newick_str, "newick")
    newick_str = newick_str.getvalue().strip()

    # Remove branch lengths and node labels
    clean_newick = re.sub(r':\d+(\.\d+)?([eE][-+]?\d+)?', '', newick_str)  # Remove branch lengths
    clean_newick = re.sub(r'\)\d+', ')', clean_newick)  # Remove node labels

    return clean_newick


def calculate_support_weight(newick_str, tree_index):
    """
    Calculate the weight of a tree based on average node support values.
    """
    # Extract node support values using regex
    support_values = re.findall(r'\)(\d+)', newick_str)
    support_values = [float(value) for value in support_values]

    if not support_values:
        warnings.warn(f"Tree {tree_index + 1}: No support values found. Default weight set to 1.")
        return 1.0  # Default weight if no support values found

    # Calculate the average support and normalize to [0, 1]
    avg_support = sum(support_values) / len(support_values)
    weight = avg_support / 100.0  # Normalize to a fraction
    return weight


def process_trees(input_file, sort_by="count"):
    """
    Process input trees to count topologies and calculate weighted sums.
    """
    topology_data = []

    with open(input_file, "r") as file:
        for index, line in enumerate(file):
            newick_str = line.strip()
            if not newick_str:
                continue

            # Parse tree and normalize topology
            tree = Phylo.read(StringIO(newick_str), "newick")
            normalized_topology = normalize_topology(tree)

            # Calculate weight for the current tree
            weight = calculate_support_weight(newick_str, index)

            # Store topology and weight
            topology_data.append((normalized_topology, weight))

    # Aggregate data by topology
    topology_stats = defaultdict(lambda: {"count": 0, "weight_sum": 0.0})
    for topology, weight in topology_data:
        topology_stats[topology]["count"] += 1
        topology_stats[topology]["weight_sum"] += weight

    # Convert to list and sort
    sorted_topologies = sorted(
        topology_stats.items(),
        key=lambda x: x[1]["weight_sum"] if sort_by == "weight" else x[1]["count"],
        reverse=True
    )

    return sorted_topologies


def write_nexus(sorted_topologies, output_prefix):
    """
    Write the topologies to a Nexus file.
    """
    nexus_file = f"{output_prefix}.nexus"
    with open(nexus_file, "w") as nexus_out:
        nexus_out.write("#NEXUS\nbegin trees;\n")
        for idx, (topo, _) in enumerate(sorted_topologies, 1):
            tree_name = f"top{idx}"
            nexus_out.write(f"  tree {tree_name} = {topo}\n")
        nexus_out.write("end;\n")
    print(f"Topologies saved to: {nexus_file}")


def write_counts_table(sorted_topologies, output_prefix):
    """
    Write the topology counts and weight sums to a tab-delimited text file.
    """
    table_file = f"{output_prefix}_count.tsv"
    with open(table_file, "w") as table_out:
        table_out.write("Tree_name\tCount\tWeight_Sum\n")
        for idx, (topo, stats) in enumerate(sorted_topologies, 1):
            table_out.write(f"top{idx}\t{stats['count']}\t{stats['weight_sum']:.4f}\n")
    print(f"Counts table saved to: {table_file}")


def main():
    parser = argparse.ArgumentParser(description="Count and extract topologies from gene trees.")
    parser.add_argument("-i", "--input", required=True, help="Input gene tree file in Newick format")
    parser.add_argument("-o", "--outpre", required=True, help="Output file prefix for Nexus and table files")
    parser.add_argument("-n", "--top-n", type=int, default=None, help="Number of top topologies to output (default: all)")
    parser.add_argument("-s", "--sort-by", choices=["count", "weight"], default="count",
                        help="Sort by 'count' or 'weight' (default: count)")

    args = parser.parse_args()

    # Process input file and count topologies
    sorted_topologies = process_trees(args.input, sort_by=args.sort_by)

    # Limit to the top N topologies if specified
    if args.top_n is not None:
        sorted_topologies = sorted_topologies[:args.top_n]

    # Write results to output files
    write_nexus(sorted_topologies, args.outpre)
    write_counts_table(sorted_topologies, args.outpre)


if __name__ == "__main__":
    main()
