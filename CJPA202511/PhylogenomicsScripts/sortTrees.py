#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Func: Sort trees based on terminal node count, max root-to-tip distance, basal branch length, name character sum, and terminal names.
# Usage: python3 sortTrees.py [-h] -i INPUT [-o OUTPUT] [-f {nexus,newick}] -s {i,d,ai,ad,r} [-b]
# By Wenjie DONG [2025/8/16]
# Email: dwjdaniel@foxmail.com

import argparse
import os
import sys
import functools
import random
import copy
from Bio import Phylo
from Bio.Phylo import BaseTree
from io import StringIO


def detect_file_format(filename):
    """Detect file format based on content (first line)"""
    try:
        with open(filename, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#NEXUS'):
                return 'nexus'
        return 'newick'
    except Exception as e:
        sys.exit(f"Error reading file: {str(e)}")


def compute_subtree_features(clade, use_brlen):
    """Compute and cache features for each node"""
    if not clade.clades:  # Terminal node
        name_str = clade.name or ""
        name_sum = sum(ord(c) for c in name_str)
        clade.features = {
            'terminal_count': 1,
            'max_distance': 0 if use_brlen else None,
            'name_sum': name_sum,
            'terminal_names': [name_str] if name_str else []
        }
    else:
        terminal_count = 0
        max_distance = 0 if use_brlen else None
        name_sum = 0
        terminal_names = []
        for child in clade.clades:
            compute_subtree_features(child, use_brlen)
            child_feat = child.features
            terminal_count += child_feat['terminal_count']
            
            if use_brlen:
                child_branch = child.branch_length or 0.0
                current_distance = child_branch + child_feat['max_distance']
                if current_distance > max_distance:
                    max_distance = current_distance
            
            name_sum += child_feat['name_sum']
            terminal_names.extend(child_feat['terminal_names'])
        
        terminal_names.sort()
        clade.features = {
            'terminal_count': terminal_count,
            'max_distance': max_distance,
            'name_sum': name_sum,
            'terminal_names': terminal_names
        }
    return clade


def deepcopy_clade(clade):
    """Create a deep copy of a clade preserving all attributes and comments"""
    new_clade = BaseTree.Clade()
    
    # Copy all attributes including comments, confidence, etc.
    for attr, value in clade.__dict__.items():
        if attr != 'clades':  # We'll handle clades separately
            setattr(new_clade, attr, copy.deepcopy(value))
    
    # Recursively copy child clades
    new_clade.clades = [deepcopy_clade(child) for child in clade.clades]
    
    return new_clade


def deepcopy_tree(tree):
    """Create a deep copy of a tree preserving all attributes"""
    # Create a new tree object with the copied root clade
    new_tree = BaseTree.Tree(root=deepcopy_clade(tree.clade), rooted=tree.rooted)
    
    # Copy all tree attributes
    for attr, value in tree.__dict__.items():
        if attr not in ['clade', 'rooted']:  # Skip attributes we've already handled
            setattr(new_tree, attr, copy.deepcopy(value))
    
    return new_tree


def sort_tree(tree, mode, use_brlen):
    """Recursively sort tree nodes while preserving all attributes and comments"""
    # Create a deep copy of the tree
    tree_copy = deepcopy_tree(tree)
    compute_subtree_features(tree_copy.clade, use_brlen)
    
    def compare_clades(a, b, current_mode):
        """Custom comparison function for two clades"""
        a_feat = a.features
        b_feat = b.features
        
        # 1. Compare terminal node count
        if a_feat['terminal_count'] != b_feat['terminal_count']:
            # REVERSED: 'i' now means larger values first (decrease), 'd' means smaller values first (increase)
            return (a_feat['terminal_count'] - b_feat['terminal_count']) if current_mode == 'd' else (b_feat['terminal_count'] - a_feat['terminal_count'])
        
        # 2. Compare max root-to-tip distance (only if use_brlen is True)
        if use_brlen and a_feat['max_distance'] is not None and b_feat['max_distance'] is not None:
            if abs(a_feat['max_distance'] - b_feat['max_distance']) > 1e-6:
                # REVERSED: 'i' now means larger values first (decrease), 'd' means smaller values first (increase)
                return (a_feat['max_distance'] - b_feat['max_distance']) if current_mode == 'd' else (b_feat['max_distance'] - a_feat['max_distance'])
        
        # 3. Compare basal branch length (only if use_brlen is True)
        if use_brlen:
            a_branch = a.branch_length or 0.0
            b_branch = b.branch_length or 0.0
            if abs(a_branch - b_branch) > 1e-6:
                # REVERSED: 'i' now means larger values first (decrease), 'd' means smaller values first (increase)
                return (a_branch - b_branch) if current_mode == 'd' else (b_branch - a_branch)
        
        # 4. Compare name character sum
        if a.clades and b.clades and a_feat['name_sum'] != b_feat['name_sum']:
            # REVERSED: 'i' now means larger values first (decrease), 'd' means smaller values first (increase)
            return (a_feat['name_sum'] - b_feat['name_sum']) if current_mode == 'd' else (b_feat['name_sum'] - a_feat['name_sum'])
        
        # 5. Compare terminal names (REVERSED DIRECTION)
        if not a.clades and not b.clades:
            a_name = a.name or ""
            b_name = b.name or ""
            # REVERSED: Now 'i' means A-Z (alphabetical), 'd' means Z-A (reverse alphabetical)
            if current_mode == 'd':
                # Z to A (reverse alphabetical)
                return (a_name < b_name) - (a_name > b_name)
            else:
                # A to Z (alphabetical)
                return (a_name > b_name) - (a_name < b_name)
        
        # 6. Preserve original order
        return 0
    
    def recursive_fixed_sort(clade, current_mode):
        """Fixed sorting mode"""
        for child in clade.clades:
            recursive_fixed_sort(child, current_mode)
        if clade.clades:
            clade.clades.sort(key=functools.cmp_to_key(lambda x,y: compare_clades(x,y,current_mode)))
    
    def recursive_alternate_sort(clade, current_mode, depth=0):
        """Alternating sorting mode"""
        # REVERSED: next_mode logic remains the same, but meanings are reversed
        next_mode = 'd' if depth % 2 == 0 else 'i'
        if mode == 'ad':
            next_mode = 'i' if depth % 2 == 0 else 'd'
        for child in clade.clades:
            recursive_alternate_sort(child, next_mode, depth+1)
        if clade.clades:
            clade.clades.sort(key=functools.cmp_to_key(lambda x,y: compare_clades(x,y,current_mode)))
    
    def recursive_random_sort(clade):
        """Random sorting mode"""
        current_mode = random.choice(['i', 'd'])
        for child in clade.clades:
            recursive_random_sort(child)
        if clade.clades:
            clade.clades.sort(key=functools.cmp_to_key(lambda x,y: compare_clades(x,y,current_mode)))
    
    # Apply sorting based on mode
    if mode in ['i', 'd']:
        recursive_fixed_sort(tree_copy.clade, mode)
    elif mode in ['ai', 'ad']:
        # REVERSED: root_mode logic remains the same, but meanings are reversed
        root_mode = 'd' if mode == 'ai' else 'i'
        recursive_alternate_sort(tree_copy.clade, root_mode)
    elif mode == 'r':
        recursive_random_sort(tree_copy.clade)
    else:
        sys.exit(f"Error: Unsupported sort mode '{mode}'")
    
    return tree_copy


def main():
    parser = argparse.ArgumentParser(description='sortTrees: Sort phylogenetic tree by tip numbers, branch lengths and tip names. By Wenjie DONG [2025/8/14].', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='input tree file (NEXUS or Newick)')
    parser.add_argument('-o', '--output', help='output file (default: stdout)')
    parser.add_argument('-f', '--format', choices=['nexus', 'newick'], default='nexus', help='output format (NEXUS or Newick)')
    parser.add_argument('-s', '--sort', choices=['i', 'd', 'ai', 'ad', 'r'], required=True, 
                        help='sort mode: i(increase), d(decrease), ai(alternate increase), ad(alternate decrease), r(random)')
    parser.add_argument('-b', '--brlen', action='store_true', help='consider branch lengths in sorting')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        sys.exit(f"Error: Input file '{args.input}' not found")
    
    # Detect file format based on content
    file_format = detect_file_format(args.input)
    
    try:
        trees = list(Phylo.parse(args.input, format=file_format))
    except Exception as e:
        sys.exit(f"Error parsing input: {str(e)}")
    
    # Process each tree while preserving all attributes and comments
    processed_trees = [sort_tree(tree, args.sort, args.brlen) for tree in trees]
    
    # Prepare output
    output = StringIO()
    Phylo.write(processed_trees, output, args.format)
    result = output.getvalue()
    output.close()
    
    # Write output
    if args.output:
        with open(args.output, 'w') as f:
            f.write(result)
        print(f"Sorted tree written to {args.output} in {args.format} format")
    else:
        print(result)


if __name__ == "__main__":
    main()
