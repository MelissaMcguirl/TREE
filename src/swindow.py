#Michael Miyagi
#Sliding Window implementation on genome strings.
#It lets us count by segregating sites and by raw bp.

from numpy import *
import matplotlib.pyplot as plt
from data_processing import *
from barcode_stats import *
from predict_rho import *
#Restrictions:
#NumSites should be even. For reasonable behavior, make offset=length/2

def segWindow(inpMat,numsites):
	outputMat=[]
#For this, I assume the offset is 1/2 the length.
	tempMat=[''.join(seq) for seq in zip(*inpMat)]
	outputMat=list(i for i in tempMat if i!=len(i)*i[0])
	return	baseWindow([''.join(seq) for seq in zip(*outputMat)],numsites,numsites/2)

def baseWindow(inpMat,length,offset):
	windowmat=[]
	for inp in inpMat:
		windowmat.append(vstack(([inp[i:i+length] for i in range(0,len(inp)-length+1,length)],[inp[i:i+length]
		for i in range(offset,len(inp)-length+offset,length)])).reshape((-1,),order='F'))
	return windowmat

def process_data(inFile):
    infile = open(inFile, 'r')
    data = infile.readlines()[1::2]
    for i in data:
        i = i[3021:]
    for i in range(len(data)):
        data[i] =  ''.join([x for x in data[i] if not x.isdigit()])
        data[i] = data[i].replace(' ', '')
        data[i] = data[i][:-1]
    print("Data Input Done")
    return data

# get barcode stats per window, print to file
def getBarCodeStats(s_data, name, OUT):
    for i in xrange(len(s_data[0])):
    # Comput hamming distance and run ripser on each subpopulation
        HammingFile = OUT +  '/Hamm_' + name + '_window_' + str(i) + '.txt'
        RipserFile = OUT + '/Rip_' + name + '_window_' + str(i) + '.txt'
        StatsFile = OUT + '/BStats_' + name + '_window_' + str(i) + '.txt'
        subPop = [s_data[j][i] for j in xrange(len(s_data))]
        matrix = empty_matrix(subPop)
        hamm_matrix = populate_matrix(matrix, subPop)
        print_to_file(hamm_matrix, HammingFile)
        reformat_Hamming(HammingFile)
        run_Ripser(HammingFile, RipserFile)
        reformat_Ripser(RipserFile)
    # Get barcode stats
        end_point, dim0, dim1 = getBars(RipserFile)
        barStats(end_point, dim0, dim1, StatsFile)
        print("Window {0} complete".format(i))

# plot avg0, avg1, b0, and b1 as window slides across genomes
def plotSplitStats(files, name):
    b0 = []
    avg0 = []
    avg1 = []
    b1 = []
    window = []
    # get stats
    for file in files:
        stats = open(file, 'r')
        stats = stats.readlines()
        b0.append(float(stats[0]))
        avg0.append(float(stats[1]))
        b1.append(float(stats[5]))
        avg1.append(float(stats[6]))
        start = file.index('_window_') + 8
        end = file.index('.txt')
        window.append(float(file[start:end]))
    # plot
    f, axarr = plt.subplots(2, 2)
    axarr[0, 0].scatter(window, b0)
    axarr[0, 0].set_title('Betti0')
    axarr[0, 1].scatter(window, avg0)
    axarr[0, 1].set_title('Average dim 0 bar length')
    axarr[1, 0].scatter(window,b1 )
    axarr[1, 0].set_title('Betti1')
    axarr[1, 1].scatter(window, avg1)
    axarr[1, 1].set_title('Average dim 1 bar length')
    plt.suptitle('Barcode stats across a sliding window, ' + name)
    plt.savefig(name)

def read_stats(file):
	stats = open(file, 'r')
	stats = stats.readlines()
	avg0 = float(stats[0])
	var0 = float(stats[1])
	b1 = float(stats[2])
	return avg0, var0, b1

def plot_TREE(files, name):
    avg0 = []
    window = []
    rho_pred = []
    # get stats
    for file in files:
        avg0, var0, b1 = read_stats(file)
        logrho_pred = log_rho(avg0, var0, b1)
        rho_pred.append(np.exp(logrho_pred)/1000)

        window_start = file.index('window_') + 7
        window_end = file.index('.txt')
        window.append(float(file[window_start:window_end]))

    # plot
    #f, axarr = plt.subplots(2, 2)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    '''Create scatterplot with regression of betti_n on r'''
    axes = plt.gca()
    #axes.set_ylim([0, max(np.log(x)) + 1])
    #ax1.set_title(args.title)
    ax1.set_xlabel("Window")
    ax1.set_ylabel(r"$\rho$")
    ax1.plot(window, rho_pred, 'b.')
    plt.savefig(name)

# compute mean b1
def mean_b1(files):
    b1 = []
    for file in files:
        stats = open(file, 'r')
        stats = stats.readlines()
        b1.append(float(stats[5]))
    return np.mean(b1)

def rho_est(barstats, coeffs):
    rho_est = coeffs[0] + coeffs[1]*barstats[0] + coeffs[2]*barstats[1] + \
              coeffs[3]*barstats[0]*barstats[1] + coeffs[4]*(barstats[0]**2) + \
              coeffs[5]*(barstats[1]**2)
    return rho_est


#To Demonstrate Functionality
#for x in baseWindow(['stringaaa','stringsss'],4,2):
#	print x

#for x in segWindow(['ACTGGGAAAA','ACTCCCGGGG','ACTCCCATTT'],4):
#	print x
