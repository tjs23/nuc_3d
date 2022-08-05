import sys, glob, colorsys, os
from os.path import abspath, dirname, join, splitext, split
sys.path.append(abspath(dirname(dirname(__file__))))

from NucApi import Nucleus, StrucCalcParams

def calcDomainStructures(combNucPath, filePaths, chromo, start, end):
  """
  Calculate topological domain structure from multiple single-cell datasets.
  Makes a combined Nuc file
  """
  
  start = int(start)
  end = int(end)
  
  os.unlink(combNucPath)
  
  nuc = Nucleus(combNucPath)
  nuc.save()
  nuc.addChromosomes({chromo:range(start, end, 1000)})
  
  nFiles = float(len(filePaths))
  groupNames = []
  
  for i, filePath in enumerate(filePaths):
    
    tempNuc = Nucleus(filePath)

    # get the contacts just for one domain
    # setup chromosomes for just one domain
    # add to combined Nuc
    
    groupName, group = tempNuc.getContactGroups()[0]
    
    cacheDict = tempNuc.getCachedContacts(groupName)
    
    
    if chromo not in cacheDict:
      continue
    
    if chromo not in cacheDict[chromo]:
      continue
    
    contacts = cacheDict[chromo][chromo] # (posA, posB, count) * n
    
    contacts = [(a, b, c) for a, b, c in contacts.T if (start <= a <= end) and (start <= b <= end)]
    
    print '%3d %s %s %4d' % (i, filePath, groupName, len(contacts))
     
    contactDict = {(chromo, chromo):contacts}
        
    newName = 'singleCell%03d' % (i+1)
    groupNames.append(newName)
    selected = True # i == 0
    color = colorsys.hsv_to_rgb(0.75 * i/nFiles, 1.0, 1.0)
    
    nuc.setContacts(newName, contactDict, replace=False,
                    selected=selected, color=color)

  
  nuc.save()
  
  # Calculate coordinates
  
  calcParams = StrucCalcParams(distPowerLaw=-0.33, distLower=0.8, distUpper=1.2,
                               bboneLower=0.1, bboneUpper=1.1, seqScale=10000,
                               randSeed=22, randStart=True,  randWalk=False, randRad=1000.0)
 
  
  calcParams.addAnnealStage(domainDict=None, bboneSep=None, # Native
                            tempStart=1000.0, tempEnd=10.0,
                            tempSteps=1000, dynSteps=100, timeStep=0.001, useTrans=True)
  
  for i, groupName in enumerate(groupNames):
     nuc.annealStructure(groupName, chromosomes=[chromo], numModels=10,
                         calcParams=calcParams, bgCalc=False, numCpus=10, structure=None)
 
  nuc.save()
  

if __name__ == '__main__':

  
  filePaths = glob.glob('/home/tjs23/Desktop/Orig_Nextera_Data/GoodNuc/S1112*.nuc')
 
  combNucFile = '/home/tjs23/nucleus/DomainsCombined.nuc'
  
  calcDomainStructures(combNucFile, filePaths, '1', start=69.2e6, end=78.55e6)
