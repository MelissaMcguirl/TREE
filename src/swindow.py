#Michael Miyagi
#Sliding Window implementation on genome strings.
#It lets us count by segregating sites and by raw bp.

from numpy import *

e#Restrictions:
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
		windowmat.append(vstack(([inp[i:i+length] for i in range(0,len(inp)-length+1,length)],[inp[i:i+length] for i in range(offset,len(inp)-length+offset,length)])).reshape((-1,),order='F'))
	return windowmat

def LDWindow(inpMat, LDFile, N):
	LDFile = open(LDFile, 'r')
	LDdata = LDFile.readlines()[3::]
	numSNPs = len(LDdata)
	numWindows = numSNPs/N
	windows = []
	for j in range(numWindows):
		windows.append((int(LDdata[j*N].split(' ')[0]),int(LDdata[(j+1)*N -1].split(' ')[1])))

	windowmat = []
	for inp in inpMat:
		windowmat.append(vstack(inp[windows[i][0]:windows[i][1]] for i in range(numWindows)).reshape(-1,))
	return windowmat


#To Demonstrate Functionality
#for x in baseWindow(['stringaaa','stringsss'],4,2):
#	print x

#for x in segWindow(['ACTGGGAAAA','ACTCCCGGGG','ACTCCCATTT'],4):
#	print x
