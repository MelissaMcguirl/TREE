# TREE
Topological REcombination Estimator

Authors: Devon P. Humphreys, Melissa R. McGuirl, Michael Miyagi, Andrew J. Blumberg

For questions/comments please contact Melissa R. McGuirl at melissa_mcguirl@brown.edu.

This software takes as input either (1) a collection of genomes in FASTA format  or (2) a distance matrix and predicts the underlying recombination rate from topological summary statistics of the data. The supported distance matrix formats are those formats that are currently supported by Ripser:

      -comma-separated values lower triangular distance matrix (preferred)
      -comma-separated values upper triangular distance matrix (MATLAB output from the function pdist)
      -comma-separated values full distance matrix

The user has the option getting recombination rate estimates over a sliding window analysis and can specify a normalization factor (default = 1/1000). 

This software is based upon the work presented in Humphreys, D.P., McGuirl, M.R., Miyagi, M., and Blumberg, A.J. Fast Estimation of Recombination Rates Using Topological Data Analysis (2018). Preprint: https://www.biorxiv.org/content/early/2018/08/20/395210

Dependencies: Ripser (https://github.com/Ripser/ripser), Python 2.7, matplotlib, numpy

Usage: python TREE.py -i INPUT_FILE [-t INPUT_TYPE] [-s SLIDING_WINDOW_FLAG] [-w WINDOW_SIZE] [-o OUTPUT_DIRECTORY] [-N NORMALIZATION_FACTOR] [-n FILENAME_IDENTIFIER] [-b BASE_FLAG] [-f OFFSET_VALUE] [-g PLOT_NAME]

Note, the default input type is a FASTA file and the default normalization factor is 1/1000. 

Examples:    

      1) python TREE.py -i ../examples/seq_example.fasta (input = FASTA file)   
      2) python TREE.py -i ../examples/hamming_example -t DIST (input = distance matrix)
      3) python TREE.py -i ../examples/seq_example.fasta -s -w 20 -o ../examples/outputs -n test -g ../examples/outputs/outputPlt (sliding window analysis)
      
Outputs:
      Sample expected output for the sliding window analysis is provided in examples/outputs/
      
Pipeline:

      0) Compute Hamming distance matrix
      1) Feed Hamming distance matrix into Ripser to compute dimension 0 and dimension 1 persistent homology barcodes
      2) Extract topological summary statistics (psi, b1, phi)
      3) Predict recombination rate using TREE
      4) Plot results (for sliding windows)

Help:
      0) python TREE.py -h
