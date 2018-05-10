'''
    Python functions to run the pipeline of analyses for topological data
    analysis on genome sequence alignments for recombination detection.
    To be used with TREE.py.

    Authors: Devon Humphreys, Melissa McGuirl, Michael Miyagi
    Corresponding Author: Melissa McGuirl (melissa_mcguirl@brown.edu)
    Date: 05/10/2018

'''
from __future__ import print_function
import sys, re
import os
import numpy as np
import subprocess, fileinput

# @profile
def Hamming_distance(string1, string2):
    '''Computes the Hamming distances between two sequences for a given alignment
    of genome sequences.'''
    if len(string1) != len(string2):
        raise ValueError("Sequences must be of the same length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(string1, string2) if ch1 != "N" and ch2 != "N")

# @profile
def format_data(input_file, filetype):
    '''Formats input data for Hamming distance matrix computations'''
    if filetype != "fasta":
        with open(input_file, 'r') as f:
            lines = f.readlines()
        for i in xrange(len(lines)):
            lines[i] = lines[i][1:]
    else:
        with open(input_file, 'r') as f:
            lines = f.readlines()[1:]
    return lines

def empty_matrix(lines):
    '''Generate an empty matrix to fill with Hamming distances'''
    return np.ndarray(shape = (len(lines), len(lines)), dtype = float)

# @profile
def populate_matrix(matrix, lines):
    '''Populates matrix elements with pairwise Hamming distances'''
    for i in range(len(lines)):
        for j in range(len(lines)):
            matrix[i][j] = Hamming_distance(lines[i], lines[j])
    return matrix

def print_to_file(data, filename):
    '''Prints Hamming distance matrix to a file in usable format'''
    outfile = open(filename, "w")
    data = str(data.tolist())
    outfile.write(data)
    outfile.close()
    return None

# @profile
def reformat_Hamming(Ham_file):
    '''Returns reformatted output Hamming file
    for use with Ripser'''
    for line in fileinput.input(Ham_file, inplace = True):
        line = re.sub(']', '\n', line, flags = re.M)#.rstrip()
        line = re.sub('(\[|^,)', '', line, flags = re.M)
        line = re.sub('^\s+', '', line, flags = re.M)
      #  print(line)
    return None

