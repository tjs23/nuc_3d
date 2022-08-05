import sys

sys.path.append('/home/tjs23/nuc/')

from core.NucApi import readContacts

petFile = '/home/tjs23/nuc/data/raw/mouseTh1/TN_SC_1.fend_pairs'

nuc = readContacts(petFile, 'Test.nuc')

nuc.calcStructure(chromosomes=['X'], numModels=1, tMax=5000.0, tMin=10.0,
                  tempSteps=100, dynSteps=1000, randomStart=True, randomWalk=False,
                  trans=True, callback=None, randomSeed=None, randomWidth=400.0,
                  backboneSep=1.0/5e5)
