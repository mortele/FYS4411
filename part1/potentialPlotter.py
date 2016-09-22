import matplotlib.pyplot 	as plt
import numpy 				as np


R = []
V = []

fileName = '../build-part1-Desktop_Qt_5_7_0_clang_64bit-Release/VMCData.dat'

with open(fileName, 'r') as inFile :
	for line in inFile :
		line = line.split()
		R.append(line[0])
		V.append(line[1])


plt.figure()
plt.plot(R,V,'ro-')
plt.xlabel('R')
plt.ylabel('V')
plt.title('Potential as a function of intermolecular distance for H2 molecule.')
plt.show()