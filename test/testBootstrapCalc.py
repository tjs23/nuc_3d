import sys, os
from os.path import dirname, abspath, exists

# Setup import path

thisDir = dirname(abspath(__file__))
sys.path.remove(thisDir)

nucDir = dirname(thisDir)
sys.path.append(nucDir)

from NucApi import Nucleus, StrucCalcParams

nucFile = '../data/examples/NXT-33_50k.nuc'
nuc = Nucleus(nucFile)

print('Setup calc params')
calcParams = StrucCalcParams(distPowerLaw=-0.33, distLower=0.8, distUpper=1.2,
                             bboneLower=0.1, bboneUpper=1.1, seqScale=10000,
                             randSeed=7, randStart=True,  randWalk=False,
                             randRad=100.0)

calcParams.addAnnealStage(domainDict={}, bboneSep=[int(5e5), int(5e5)],
                          tempStart=5000.0, tempEnd=10.0, tempSteps=100,
                          dynSteps=1000, timeStep=0.001, useTrans=True)

print('Calculate structures')
groupName = 'singleCell'
chromosomes=['1']

strucCode = nuc.getNextStructureCode()
strucGroup = nuc.getStructureGroup(strucCode, 'Test1') 
 
nuc.annealStructure(groupName, chromosomes, numModels=10,
                    calcParams=calcParams, numCpus=10,
                    structure=strucCode)                         

nuc.annealStructureBootstrap(groupName, chromosomes, 
                             numModels=10, numPartitions=4, numResamples=2,
                             calcParams=calcParams, numCpus=10)

nuc.setBootstrapOverviewStructure(groupName)

nuc.save()
