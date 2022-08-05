import sys, os
from os.path import dirname, abspath, exists
from numpy import histogram, array

# Setup import path

thisDir = dirname(abspath(__file__))
sys.path.remove(thisDir)

nucDir = dirname(thisDir)
sys.path.append(nucDir)

from NucApi import Nucleus, StrucCalcParams
from matplotlib import pyplot

nucFile = '../data/examples/NXT-33_50k.nuc'
nuc = Nucleus(nucFile)

groupName = 'singleCell'

distsExlude, distsInclude = nuc.calcBootstrapCrossValid(groupName)

binMax = 5.0

data, edges = histogram(distsExlude, 50, (0.0, binMax))
data = array(data, float)
data /= data.sum()
pyplot.plot(edges[:-1], data, color='red')

data, edges = histogram(distsInclude, 50, (0.0, binMax))
data = array(data, float)
data /= data.sum()
pyplot.plot(edges[:-1], data, color='green')

pyplot.show()

rmsdsBundle, rmsdsPartition, rmsdsResample = nuc.calcBootstrapRmsds(groupName)

binMax = 0.75

data, edges = histogram(rmsdsBundle, 50, (0.0, binMax))
data = array(data, float)
data /= data.sum()
pyplot.plot(edges[:-1], data, color='red')

data, edges = histogram(rmsdsPartition, 50, (0.0, binMax))
data = array(data, float)
data /= data.sum()
pyplot.plot(edges[:-1], data, color='green')

data, edges = histogram(rmsdsResample, 50, (0.0, binMax))
data = array(data, float)
data /= data.sum()
pyplot.plot(edges[:-1], data, color='blue')

pyplot.show()
