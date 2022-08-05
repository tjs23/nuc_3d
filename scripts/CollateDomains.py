import sys, glob, colorsys, os
from os.path import abspath, dirname, join, splitext, split
sys.path.append(abspath(dirname(dirname(__file__))))

from numpy import array, float32
from NucApi import Nucleus, StrucCalcParams

def collateDomains(combNucPath, filePaths, chromo, start, end):
  """
  Combine structures and contacts for a given chromosome region from
  from multiple single-cell datasets. Makes a single, combined Nuc file.
  """
  
  start = int(start)
  end = int(end)
  
  if os.path.exists(combNucPath):
    os.unlink(combNucPath)
  
  nuc = Nucleus(combNucPath)
  nuc.removeStructure(0)
  nuc.save()
  
  nFiles = float(len(filePaths))
  groupNames = []
  
  for i, filePath in enumerate(filePaths):
    
    tempNuc = Nucleus(filePath)
    fileName = os.path.split(filePath)[1]
    
    newName = os.path.splitext(fileName)[0].split('_')[1]
    
    # get the contacts and coords just for one domain
    # setup chromosomes for just one domain
    # add to combined Nuc
    
    groupName, group = tempNuc.getContactGroups()[0]
    
    cacheDict = tempNuc.getCachedContacts(groupName)
    
    if chromo not in cacheDict:
      continue
    
    if chromo not in cacheDict[chromo]:
      continue
    
    # New structure slot
      
    nuc.getStructureGroup(i, newName)
    
    # Fetch cis contacts
    
    contacts = cacheDict[chromo][chromo] # (posA, posB, count) * n
    
    # Get region contacts
    
    contacts = [(a, b, c) for a, b, c in contacts.T if (start <= a <= end) and (start <= b <= end)]
    
    # Set region contacts
     
    contactDict = {(chromo, chromo):contacts}
        
    #newName = 'singleCell%03d' % (i+1)
    groupNames.append(newName)
    
    selected = True # i == 0
    color = colorsys.hsv_to_rgb(0.75 * i/nFiles, 1.0, 1.0)
    
    nuc.setContacts(newName, contactDict, replace=False,
                    selected=selected, color=color)
    
    # Fetch all coords
    
    coordsGroup = tempNuc._getCoordsGroup()
    
    modelCoords = coordsGroup[chromo]
    
    # Fetch all particle positions
    
    particGroup = tempNuc._getParticleGroup()
    
    positions = array(particGroup[chromo]['positions'])
    
    idxToPos = {}
    for j, p in enumerate(positions):
      idxToPos[j] = p 
        
    # Filter postion and coord regions
    
    indices = (start <= positions).nonzero()

    positions = positions[indices]
    
    modelCoords = modelCoords[:,indices[0]]

    indices = (positions <= end).nonzero()

    positions = positions[indices]
    
    modelCoords = modelCoords[:,indices[0]]

    # Set region positions
    
    nuc.addChromosomes({chromo:positions}, structure=i)
    
    posToIdx = {}
    for j, p in enumerate(positions):
      posToIdx[p] = j
    
    # Set region coords
    
    nuc.setModelChromosomeCoords(modelCoords, chromo, structure=i)
    nuc.getStructureGroup(i).attrs['displayModels'] = [0] # array(range(len(modelCoords)))
    
    # Set restraints
    
    restGroup = tempNuc._getRestraintsGroup()
    
    restSubGroup = tempNuc._getGroup(chromo, restGroup)
    
    restArray = array(restSubGroup[chromo]).T

    restraints = []
    
    for rest in restArray:
      a = int(rest[0])
      b = int(rest[1])
      
      posA = idxToPos[a]
      posB = idxToPos[b]
      
      if (start <= posA <= end) and (start <= posB <= end):
        rest[0] = float(posToIdx[posA])
        rest[1] = float(posToIdx[posB])
        restraints.append(rest)
        
    restraints = array(restraints, float32).T
    
    restGroup = nuc._getRestraintsGroup(structure=i)
    
    restSubGroup = nuc._getGroup(chromo, restGroup)
   
    nuc._setData(chromo, restSubGroup, float32, restraints)
    
    print 'Struc:%3d From:%s Group:%s Particles:%4d Contacts:%4d Restraints:%4d' % (i, split(filePath)[1], newName, len(positions), len(contacts), len(restraints[0]))

    nuc.save()
  
  
  nuc.structureAlign()
  nuc.mergeStructures(nModels=1)
  nuc.getStructureGroup(i+1).attrs['displayModels'] = array(range(nuc.getNumModels(i+1)))
  
  nuc.save()
  

if __name__ == '__main__':

  
  #filePaths = glob.glob('/home/tjs23/nucleus/Nextera01/S1112*400kb.nuc')
  filePaths = glob.glob('/home/tjs23/nucleus/Nextera_01/best_100kb/*_10x_100kb.nuc')
 
  combNucFile = '/home/tjs23/nucleus/DomainsCombined_Chr12_58-70Mb.nuc'
  
  collateDomains(combNucFile, filePaths, '12', start=58e6, end=70e6)
