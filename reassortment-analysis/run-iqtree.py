#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 17:07:10 2021

@author: Alexey Markin
"""
from Bio import SeqIO
import subprocess
from dendropy import Tree
import sys
import re


def collapse_zero_branches(tree_path: str, threshold=1e-6) -> str:
    tree = Tree.get(path=tree_path, schema='newick')
    print("Original # of internal nodes: %d" % len(tree.internal_nodes()))
    tree.collapse_unweighted_edges(threshold)
    print("After cleaning: %d" % len(tree.internal_nodes()))
    updated_tree_path = tree_path
    tree_file = open(updated_tree_path, 'w+')
    tree_file.write(str(tree) + ';\n')
    tree_file.close()
    return updated_tree_path


def replace_names(tree_path: str, subbed_to_original: dict):
    tree = Tree.get(path=tree_path, schema='newick')
    tree_str = str(tree)
    for subbed in subbed_to_original.keys():
        original = subbed_to_original[subbed]
        tree_str = tree_str.replace(subbed, original)
    tree_file = open(tree_path, 'w+')
    tree_file.write(tree_str + ';\n')
    tree_file.close()
    print('Replaced names.')


def run_iqtree(alignment_path: str, iqtree_args=None):
    seqs = SeqIO.parse(alignment_path, 'fasta')
    subbed_to_original = {re.subn(r'(\||/|;|\+|\(|\)|:)', '_', seq.id)[0]: seq.id
                         for seq in seqs}
    
    if iqtree_args and iqtree_args.count('-pre') > 0:
        prefix = iqtree_args[iqtree_args.index('-pre') + 1]
    else:
        prefix = alignment_path
    tree_file = prefix + '.treefile'
    command = ['iqtree', '-s', alignment_path]
    if iqtree_args:
        command = command + iqtree_args
    print(' '.join(command))
    subprocess.call(command, stderr=subprocess.STDOUT)
    replace_names(tree_file, subbed_to_original)
    # tree_file = collapse_zero_branches(tree_file)


if __name__ == '__main__':
    args = sys.argv[1:]
    if len(args) == 1:
        run_iqtree(args[0])
    else:
        run_iqtree(args[0], args[1:])
