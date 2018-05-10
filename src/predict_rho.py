import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt

# a function to extract barcode stats from file
def get_bstats(file):
    stats = open(file, 'r')
    stats = stats.readlines()

    avg0 = float(stats[0][:-1])
    var0 = float(stats[1][:-1])
    b1 = float(stats[2][:-1])
    return avg0, var0, b1

# a function to compute log(TREE)
def log_rho(avg0, var0, b1):
    A =  5.53046629e-02
    B =  -3.74380246e-04
    C = 1.81333434e-02
    D = -1.79713403e-04
    E =  -5.93368387e-05
    y_int = 2.24756003254
    logrho = A*avg0 + B*var0 + C*b1 + D*(avg0**2) + E*(b1**2) +  y_int
    return logrho 


def main():
    descriptor = "A python function that predict rho using psi, var0, b1."

    parser = argparse.ArgumentParser(description = descriptor)
    parser.add_argument('-i', '--indir', required = True, help = 'provide path to folder containing barcode stats, labeled as BCStats*')
    args = parser.parse_args()
    IN = args.indir
    
    files = glob.glob(IN + '/BCStats_*')
    rho_pred = []
    window = []
    psi = []
    # loop through files in directory and predict recombination 
    for file in files:
        avg0, var0, b1 = get_bstats(file)
        logrho_pred = log_rho(avg0, var0, b1)
        rho_pred.append(np.exp(logrho_pred))
        window_start = file.index('window_') + 7
        window_end = file.index('.txt')
        window.append(float(file[window_start:window_end]))
        psi.append(avg0)
        j = j + 1

    # plot results
    fig, ax = plt.subplots()
    ax.scatter(window, rho_pred, edgecolors=(0, 0, 0))
    ax.set_ylabel('rho')
    ax.set_xlabel('window')
    ax.set_title('TREE Estimates')
    plt.show()


if __name__ == "__main__":
    main()
