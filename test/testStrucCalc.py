import sys, os
from os.path import dirname, abspath, exists

# Setup import path

thisDir = dirname(abspath(__file__))
sys.path.remove(thisDir)

nucDir = dirname(thisDir)
sys.path.append(nucDir)

# 

bamFile = "data/SiCUP/sample_1268/sample_1268_L008_TAAGGCGA_ACTGCATA_SiCUPPED.bam"
nucFile = 'StrucCalcTest_temp.nuc'
groupName='singleCell'

if exists(nucFile):
  os.remove(nucFile)

#

from NucApi import Nucleus, StrucCalcParams

print('New Nuc')
nuc = Nucleus(nucFile)

print('Import BAM contacts')
nuc.importContacts(bamFile, 'sam', groupName)

print('Remove noise')
nuc.removeIsolatedContacts(groupName, threshold=int(1e6))


print('Setup calc params')
calcParams = StrucCalcParams(distPowerLaw=-0.33, distLower=0.8, distUpper=1.2,
                             bboneLower=0.1, bboneUpper=1.1, seqScale=10000,
                             randSeed=7, randStart=True,  randWalk=False, randRad=100.0)

calcParams.addAnnealStage(domLstDict={}, bboneSep=[int(5e4), int(5e4)], tempStart=5000.0, tempEnd=10.0,
                          tempSteps=100, dynSteps=100, timeStep=0.001, useTrans=True)

                    
print('Calculate structure')
nuc.annealStructure(groupName, chromosomes=['1'], numModels=1,
                    calcParams=calcParams, numCpus=1)                         

