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
from swindow import *
import subprocess, os, time
from data_processing import *
from barcode_stats import getBars, barStats
from predict_rho import *
import numpy as np
import argparse

def main():

    descriptor = '''A Topological Recombination Rate Efficient Estimator. This
    software takes as input either a collection of sequence alignments in FASTA
    format or a distance matrix (specify if inputting distance matrix.)'''

    parser = argparse.ArgumentParser(description = descriptor)

    parser.add_argument('-i', '--indir',
                        action = 'store',
                        required = True,
                        help = '''provide path to directory containing sequence
                        alignment file or distance matrix.''')
    parser.add_argument('-t', '--InType', action = 'store', required = False,
                        default = 'FASTA', help = '''specify FASTA or DIST
                        input.''')
    parser.add_argument('-s', '--swindow', action = 'store_true', required = False,
                        help = '''run a sliding window analysis''')
    parser.add_argument('-w', '--WindowSize', action = 'store', required = False,
                        help = '''provide size of sliding window''')
    parser.add_argument('-o', '--output', action = 'store', required = False,
                        help = '''provide name for output file, if desired''')
    parser.add_argument('-n', '--name', action = 'store', required = False,
                        help =''' provide unique identifier for a batch of sliding
                                window runs''')

    args = parser.parse_args()
    glob_start = time.time()

    inFile = args.indir
    inType = args.InType


    if inType == 'FASTA':
        # compute Hamming distance matrix and reformat it for Ripser inputs.
        HammingFile = "HammingMat"
        print("Reading FASTA file...")
        lines = format_data(inFile, 'fasta')
        print("FASTA file read.")
        matrix = empty_matrix(lines)
        print("Computing Hamming distance matrix...")
        hamm_matrix = populate_matrix(matrix, lines)
        print("Hamming distance matrix computed.")
        print_to_file(hamm_matrix, HammingFile)
        reformat_Hamming(HammingFile)

    else:
        HammingFile = inFile

    # run ripser.
    RipserFile = "RipserFile"
    cmd = "ripser %s | cut -f 2 -d: | awk 'NR > 1 {print}' | cut -c 3- | sed 's/.$//' > %s" % (HammingFile,  RipserFile)
    os.system(cmd)
    print("Persistent Homology computations complete.")
    end_point, dim0, dim1 = getBars(RipserFile)
    # compute barcode statistics
    StatsFile = "BCStats"
    barStats(end_point,dim0, dim1, StatsFile)
    avg0, var0, b1 = get_bstats(StatsFile)
    # compute log rho
    logrho_pred = log_rho(avg0, var0, b1)
    # convert to TREE predictions
    pred_rho = np.exp(logrho_pred)
    glob_end = time.time()
    # print outputs
    print("TREE Prediction: rho^={0}".format(pred_rho))
    print("Total time: {0}".format(glob_end - glob_start))


    if args.swindow:
        inFile  = args.indir
        OUT = args.output
        N = args.windowSize
        name = args.name
        print("Processing data...")
        data = process_data(inFile)
        print("Data processing done.")
        # get sliding window data
        print("Begin sliding window...")
        s_data = segWindow(data, int(N))
        print("Sliding window analysis complete.")
        #get barcode stats
        print("Processing barcode statistics...")
        getBarCodeStats(s_data, name, OUT)
        print("Barcodes processed.")
        # plot results
        print("Plotting results...")
        Stats_files = glob.glob(OUT + '/BStats_' + name + '_window_*.txt')
        #plotSplitStats(Stats_files, name)
        plotPsi(Stats_files)
        print("TREE analysis complete.")



if __name__=="__main__":

    main()
