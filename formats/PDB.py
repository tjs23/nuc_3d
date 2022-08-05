from formats.Util import getFileObj
from numpy import array, uint32, float32

def exportCoords(filePath, posDict, coordsDict, nuc, structure, scale=1.0, extended=True):

  alc = ' '
  ins = ' '
  prefix = 'HETATM'
  lFormat = '%-80.80s\n'
  
  if extended:
    pdbFormat = '%-6.6s%5.1d %4.4s%s%3.3s %s%4.1d%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  %10d\n'
    terFormat = '%-6.6s%5.1d      %s %s%4.1d%s                                                     %10d\n'
  else:
    pdbFormat = '%-6.6s%5.1d %4.4s%s%3.3s %s%4.1d%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  \n'
    terFormat = '%-6.6s%5.1d      %s %s%4.1d%s                                                     \n'

  fileObj = getFileObj(filePath, 'w')
  write = fileObj.write

  sample = nuc.sample
  title = 'Nuc3D genome structure export'
  infoDict = nuc.getSampleInfo()
  particGroup = nuc._getParticleGroup(structure)
  
  chromosomes = list(posDict.keys())
  sort_chromos = []
  
  for chromo in chromosomes:
    if chromo[:3] == 'chr':
      key = chromo[3:]
    
    else:
      key = chromo
    
    if key.isdigit():
      key = '%03d' % int(key)
    
    sort_chromos.append((key, chromo))
  
  sort_chromos.sort()
  sort_chromos = [x[1] for x in sort_chromos]
    
  chrA = sort_chromos[0]
  pos = posDict[chrA]
  seqSep = nuc.getBackboneSpacing()
  coords = coordsDict[chrA]
  numModels = len(coords)
  
  line = 'TITLE     %s' % title
  write(lFormat % line)
  write(lFormat % 'REMARK 210')
  
  #for key in sorted(infoDict):
  #  info = '%s' % infoDict[key]
  #  write(lFormat % 'REMARK 210 %s:%s' % (key,info.rstrip())) 
 
  write(lFormat % 'REMARK 210 Atom type C are backbone nodes')
  write(lFormat % 'REMARK 210 Atom type N are restrained nodes')
  write(lFormat % 'REMARK 210 Atom number increases every %s bases' % seqSep)
  write(lFormat % 'REMARK 210 Residue code indicates chromosome')
  write(lFormat % 'REMARK 210 Residue number represents which sequence Mb the atom is in')
  write(lFormat % 'REMARK 210 Chain letter is different every chromosome, where Chr1=a, Chr2=b etc.')
  #write(lFormat % 'REMARK 210 B-factor field is sequence mapability')
  #write(lFormat % 'REMARK 210 Coordinate scale is 1 unit to %.2f microns' % scale)
  
  if extended:
    fileObj.write(lFormat % 'REMARK 210 Extended PDB format with particle seq. pos. in last column')
  
  fileObj.write(lFormat % 'REMARK 210')
  
  pos_chromo = {}
  
  for m in range(numModels):
    line = 'MODEL     %4d' % (m+1)
    write(lFormat  % line)
    
    c = 0
    j = 1
    seqPrev = None
    for k, chromo in enumerate(sort_chromos):
      chain_code = chr(ord('a')+k)
      
      tlc = chromo
      while len(tlc) < 2:
        tlc = '_' + tlc
      
      if len(tlc) == 2:
        tlc = 'C' + tlc
      
      if len(tlc) > 3:
        tlc = tlc[:3]
      
      chromoModelCoords = coordsDict[chromo][m]
      
      if not len(coords):
        continue
      
      pos = posDict[chromo]
      backbone = array(particGroup[chromo[3:]]['backbone'])
      mapability = array(particGroup[chromo[3:]]['mapability'])
      
      for i, seqPos in enumerate(pos):
        c += 1
 
        seqMb = int(seqPos//1e6) + 1
        
        if seqMb == seqPrev:
          j += 1
        else:
          j = 1
        
        if backbone[i]:
          el = 'C'
          a = 'C%d' % j
        else:
          el = 'N'
          a = 'N%d' % j
          
        aName = '%-3s' % a
        b = mapability[i]
        x, y, z = chromoModelCoords[i] #XYZ coordinates
         
        seqPrev = seqMb
        pos_chromo[c] = chromo
        
        if extended:
          line  = pdbFormat % (prefix,c,aName,alc,tlc,chain_code,seqMb,ins,x,y,z,0.0,b,el,seqPos)
        else:
          line  = pdbFormat % (prefix,c,aName,alc,tlc,chain_code,seqMb,ins,x,y,z,0.0,b,el)
          
        write(line)
 
    write(lFormat  % 'ENDMDL')
 
  for i in range(c-2):
     if pos_chromo[i+1] == pos_chromo[i+2]:
       line = 'CONECT%5.1d%5.1d' % (i+1, i+2)
       write(lFormat  % line)
 
  write(lFormat  % 'END')
  fileObj.close()

