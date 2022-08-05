import gzip, sys
from numpy import array, zeros, savez, log2, int32, arange, dstack, uint32, float32
from glob import glob
from collections import defaultdict

from os.path import abspath, dirname
sys.path.append(abspath(dirname(dirname(__file__))))

from cUtil.apiUtil import calcSeqEntropy

def extractSeqEntropy(fastaPattern, window=100, namePattern='seq_entropy_%d', limit=0.2):

  name = namePattern % window
  outFileName = name + '.npz'
  fastaFiles = glob(fastaPattern)

  letters = 'CGATN'
  iDict = {}
  for i, x in enumerate(letters):
    iDict[x] = i

  counts = zeros(len(letters), int)
  kwArgs = {}

  for fastaFile in fastaFiles:
    print(fastaFile)
    fileObj = gzip.open(fastaFile, 'rt')

    seq = []
    append = seq.append
 
    for line in fileObj:
      if line[0] != '>':
        append(line[:-1])

    seq = ''.join(seq)
    chromo = fastaFile.split('.')[-3]
    n = len(seq)

    entropy = calcSeqEntropy(seq, window)
    
    m = int(n/window)
    entropy = entropy[:m*window]
    entropy = entropy.reshape(m, window)
    entropy = entropy.mean(axis=1)
     
    rKey = 'dataTrack/regions/%s/%s' % (name, chromo)
    vKey = 'dataTrack/values/%s/%s' % (name, chromo)
    
    regions = arange(0, m*window, window, dtype=uint32)
    regions = dstack([regions, regions+window])[0]
        
    # Clip outliers, e.g. polyN
    
    idx = (entropy < limit).nonzero()
    entropy = entropy[idx]
    regions = regions[idx]

    idx = (entropy > 0.0).nonzero()
    entropy = entropy[idx]
    regions = regions[idx]
        
    values= dstack([entropy, entropy])[0]
    
    kwArgs[rKey] = regions
    kwArgs[vKey] = array(values, float32)
 
  savez(outFileName, **kwArgs)
    

for window in [10000]:
  extractSeqEntropy('/data/genome/mm10_chromosomes/*.fa.gz', window)
