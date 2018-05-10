'''
    TREE.py

    Input: fasta formatted sequence alignments or hamming distance matrix.

    Output: Hamming distance matrices, Ripser persistent homology intervals,
    barcode statistics for each alignment input.

    Usage: python TREE.py  <path/to/sequence/files>

    Authors: Devon Humphreys, Melissa McGuirl, Michael Miyagi
    Corresponding Author: Melissa McGuirl (melissa_mcguirl@brown.edu)
    Updated: 05/10/18
'''
import glob, time
from data_processing import *
from barcode_stats import getBars, barStats
from predict_rho import *
import numpy as np
import argparse

def main():

    descriptor = '''A Topological Recombination Rate Efficient Estimator. This
    software takes as input either a collection of sequence alignments in FASTA
    format (*.fasta) or a distance matrix.'''
    
    parser = argparse.ArgumentParser(description = descriptor)

    parser.add_argument('-i', '--indir',
                        action = 'store',
                        required = True,
                        help = '''provide path to directory containing sequence
                        alignment file or distance matrix.''')
    parser.add_argument('-t', '--InType', action = 'store', required = False,
                        default = 'FASTA', help = '''specify FASTA or DIST
                        input.''')

    args = parser.parse_args()
    glob_start = time.time()

    inFile = args.indir
    inType = args.InType

    
    if inType == 'FASTA':

        HammingFile = "HammingMat" 
        lines = format_data(inFile, 'fasta')
        matrix = empty_matrix(lines)
        hamm_matrix = populate_matrix(matrix, lines)
        print("Hamming distance matrix computed.")
        print_to_file(hamm_matrix, HammingFile)
        reformat_Hamming(HammingFile)

    else:
        HammingFile = inFile

    RipserFile = "RipserFile"
    run_Ripser(HammingFile, RipserFile)
    print("Persistent Homology Computations complete.")
    reformat_Ripser(RipserFile)
    end_point, dim0, dim1 = getBars(RipserFile)
    StatsFile = "BCStats"
    barStats(end_point,dim0, dim1, StatsFile)
    avg0, var0, b1 = get_bstats(StatsFile)
    logrho_pred = log_rho(avg0, var0, b1)
    pred_rho = np.exp(logrho_pred)
    glob_end = time.time()
    
    print("TREE Prediction: rho={0}".format(pred_rho))
    print("Total time: {0}".format(glob_end - glob_start))


if __name__=="__main__":

    main()
