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
from ripser import ripser 


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
                        help = '''provide directory in which output 
                        files should be stored.''')
    parser.add_argument('-N', '--normFactor', action = 'store', required =False, 
                        default = 1000, help =''' provide normalization factor for scaling
                            rho. Default =1000.''')
    parser.add_argument('-n', '--name', action = 'store', required = False,
                        help =''' provide unique identifier for a batch of sliding
                                window runs''')
    parser.add_argument('-b', '--base', action = 'store_true', required = False,
                        help = '''if doing a sliding window analysis, this flag
                            sets the window to move across raw base pairs rather than
                            SNPs. The default behavior is to move over only
                            segregating sites and ignore constant sites.''')
    parser.add_argument('-f', '--offset', action = 'store', required = False,
                        help = '''supply the offset for a sliding window using raw
                            base pairs (how much overlap there will be).''')
    parser.add_argument('-g', '--graph', action = 'store', required = False,
                        help='''supply name of graph to write to file''')

    args = parser.parse_args()
    glob_start = time.time()

    inFile = args.indir
    inType = args.InType
    normFactor = float(args.normFactor)

    if not args.swindow:
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
        else:
            HammingFile = inFile
            hamm_matrix = np.loadtxt(HammingFile, dtype='float', delimiter = ',')
        # run ripser.
        print("Running Ripser...")
        bars = ripser(hamm_matrix,  distance_matrix=True, maxdim=1)
        print("Ripser analysis complete.")
        end_point = np.max(hamm_matrix)
        dim0 = bars['dgms'][0]
        dim1 = bars['dgms'][1]
        print("Persistent Homology computations complete.")
        # compute barcode statistics
        StatsFile = "BCStats"
        print("Computing barcode statistics...")
        barStats(end_point, dim0, dim1, StatsFile)
        avg0, var0, b1 = get_bstats(StatsFile)
        print("Barcode statistics complete.\n")
        print("Psi={0}\nVariance={1}\nBetti1={2}\n".format(avg0, var0, b1))
        # compute log rho
        print("Computing Rho^...")
        logrho_pred = log_rho(avg0, var0, b1)
        # convert to TREE predictions
        pred_rho = np.exp(logrho_pred)/normFactor
        print("Rho estimation complete.\n")
        glob_end = time.time()
        # print outputs
        print("TREE Prediction: rho^={0}".format(pred_rho))
        print("Total time: {0}".format(glob_end - glob_start))

    elif args.swindow:
        inFile  = args.indir
        OUT = args.output
        N = args.WindowSize
        name = args.name
        graph = args.graph
        print("Processing data...")
        data = process_data(inFile)
        print("Data processing done.")
        # get sliding window data
        if args.base:
            s_data = baseWindow(data, int(N), int(args.offset))
            print("Begin sliding window...")
        else:
            s_data = segWindow(data, int(N))
        print("Begin sliding window...")
        print("Sliding window analysis complete.")
        #get barcode stats
        print("Processing barcode statistics...")
        getBarCodeStats(s_data, name, OUT)
        print("Barcodes processed.")
        # plot results
        print("Plotting results...")
        Stats_files = glob.glob(OUT + '/BStats_' + name + '_window_*.txt')
        #plotSplitStats(Stats_files, name)
        plot_TREE(Stats_files, graph)
        print("TREE analysis complete.")



if __name__=="__main__":

    main()
