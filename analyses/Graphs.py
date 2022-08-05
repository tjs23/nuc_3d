from gui.qtgui.BasePopup import BasePopup
from gui.qtgui.Graph import GraphAxis, GraphDataSet, Graph

import numpy, math
from numpy.fft import fft
from numpy import abs, array, zeros, concatenate, hstack, log, log10, log2, histogram, arange
from math import ceil

from matplotlib import pyplot as plt


def graphChromoIntermingling(nuc, chromosomes, radii=arange(1.0, 3.2, 0.2)):

  xyVals = []
  dataSets = []
  models = nuc.getDisplayedModels()
  
  for r in radii:
    msg = 'Radius = %.2f' % r
    print(msg)
    frac = nuc.calcChromoIntermingling(chromosomes, radius=r, models=models,
                                       structure=None, label=None)
    xyVals.append((r, 100.0 * frac))
    
  popup = BasePopup(title='Graph popup')  
  xAxis = GraphAxis('Surface radius', labels=None, ticks=True, valRange=(max(radii), min(radii)))
  yAxis = GraphAxis('% intermingling', labels=None, ticks=True)
  
    
  ds = GraphDataSet(xyVals, '', '#A00000', plotType='line', symbol='circle', symbolSize=3)
  dataSets.append(ds)
    
  graph = Graph(popup, (xAxis, yAxis), dataSets, size=(800, 300),
                title='Chromosome intermingling') 


def graphPseudo4c(nuc, chromosome, position, groupNames, threshold=int(2e5), binSize=int(2e5)):

  popup = BasePopup(title='Graph popup')  
  xAxis = GraphAxis('Position', labels=None, ticks=True)
  yAxis = GraphAxis('log10(probability)', labels=None, ticks=True)
  yAxisB = GraphAxis('log10(observed/expected)', labels=None, ticks=True)
  
  
  dataSets = []
  dataSetsB = []
  
  for groupName in groupNames:
    contactDict = nuc.getCachedContacts(groupName)
    
    allCis = contactDict[chromosome][chromosome]
    deltas = abs(allCis[1]-allCis[0])
    numObs = allCis[2]
    
    nBins = deltas.max()/binSize
    allCounts, edges = histogram(deltas, bins=nBins, weights=numObs)
    pExp = allCounts / float(allCounts.sum())
    
    regionContacts = nuc.getRegionContacts(groupName, chromosome, position-threshold, position+threshold, cis=True, trans=False)
    
    positions = [x[1] for x in regionContacts]
    weights = [x[2] for x in regionContacts]
    minPos = min(positions)
    maxPos = max(positions)
    
    nBins = int((maxPos-minPos)/binSize)

    iRef = int((position-minPos)/binSize)
    seqSepBins = abs(arange(nBins) - iRef)  # [10,9,8,7,6,5,4,3,2,1,0,1,2,3 ] etc
    seqSepBins[iRef+1:] -= 1                # [10,9,8,7,6,5,4,3,2,1,0,0,1,2,] etc, accountys for rounding down of delta
    
    counts, edges = histogram(positions, bins=nBins, weights=weights)
    nz = counts.nonzero()
     
    expected = (pExp/2)[seqSepBins] # Half of the counts above/below ref position
    
    counts = counts[nz]
    
    observed = counts / float(counts.sum())
    expected = expected[nz]
    edges = edges[nz]
    
    xyValues = [(edges[i], log10(x)) for i, x in enumerate(observed) if x]   
    color = '#000080' #nuc.getContactGroupColor(groupName).rgbHex
    ds = GraphDataSet(xyValues, 'Observed '+groupName, color, plotType='line', symbolSize=1)
    dataSets.append(ds)

    xyValues = [(edges[i], log10(x)) for i, x in enumerate(expected) if x]   
    color = '#800000' #nuc.getContactGroupColor(groupName).rgbHex
    ds = GraphDataSet(xyValues, 'Expected', color, plotType='line', symbolSize=1)
    dataSets.append(ds)
  
    bias = log10(observed/expected)
    xyValues = [(edges[i], x) for i, x in enumerate(bias) if x]   
    color = '#008000' #nuc.getContactGroupColor(groupName).rgbHex
    ds = GraphDataSet(xyValues, groupName, color, plotType='bar', symbolSize=3)
    dataSetsB.append(ds)
  
  posA = position - threshold
  posB = position + threshold
  graph = Graph(popup, (xAxis, yAxis), dataSets, size=(1000, 300),
                title='Region {}:{:,}-{:,} bp'.format(chromosome, posA, posB)) 

  graph = Graph(popup, (xAxis, yAxisB), dataSetsB, size=(1000, 300),
                title='Region {}:{:,}-{:,} bp'.format(chromosome, posA, posB)) 


def graphContactFourierTransform(nuc, minContacts=130, binSize=1e5, nCols=3, nPoints=500, groupName='singleCell'):

  popup = BasePopup(title='Graph popup')

  contactDict = nuc.getContacts(groupName, cis=False, trans=True)
  
  nPoints = 2**int(1+ceil(log(200e6/binSize, 2)))
  
  print('Num points: %d' % nPoints)
  
  base = 10.0
  p1 = nPoints * binSize / 5e7
  p2 = nPoints * binSize / 1e6
  
  xStart = log((binSize*nPoints)/p2, base)
  xEnd = log((binSize*nPoints)/p1, base)
  
  
  g = 0
  for (chromoA, chromoB) in contactDict:
    key = (chromoA, chromoB)
    contacts = contactDict[key]
    nContacts = len(contacts[0]) 
    
    if nContacts < minContacts:
      continue
    
    if chromoA == chromoB:
      posA = concatenate([contacts[0], contacts[1]])
    
    else:
      posA = contacts[0]
      posB = contacts[1]
      
    dataSets = []
   
    pMax = posA.max()
    nBins = int(math.ceil(pMax/binSize))
    #print(nContacts)
    hist, edges = numpy.histogram(posA, nBins, (0.0, pMax))
   
    if nBins < nPoints:
      signal = zeros(nPoints, float)
      signal[:nBins] = hist
   
    else:
      signal = hist[:nPoints]
   
    fftA = (fft(signal)[p1:p2]**2).real
    xyValues = [(log(binSize*nPoints/(i+1.0+p1), base), val) for i, val in enumerate(fftA)]
    color = '#FF0000' # nuc.getChromoColor(chromoA, dType=str)
    ds = GraphDataSet(xyValues, 'Chr '+chromoA, color, plotType='line', symbolSize=0)
    dataSets.append(ds)
 
    if chromoA != chromoB:

      pMax = posB.max()
      nBins = int(math.ceil(pMax/binSize))
      hist, edges = numpy.histogram(posB, nBins, (0.0, pMax))
      
      if nBins < nPoints:
        signal = zeros(nPoints, float)
        signal[:nBins] = hist
      
      else:
        signal = hist[:nPoints]
      
      fftB = (fft(signal)[p1:p2]**2).real
      xyValues = [(log(binSize*nPoints/(i+1.0+p1),base), val) for i, val in enumerate(fftB)]
      color = '#0000FF' # nuc.getChromoColor(chromoB, dType=str)
      ds = GraphDataSet(xyValues, 'Chr '+chromoB, color, plotType='line', symbolSize=0)
      dataSets.append(ds)
      
    xAxis = GraphAxis('Wavelength Log10(bp)', labels=None, ticks=True,
                      valRange=(xStart, xEnd))
    yAxis = GraphAxis('Intensity', labels=None, ticks=True,
                      valRange=(-3e3, 3e3))
    graph = Graph(popup, (xAxis, yAxis), dataSets,
                  title='Contact Power Spectrum Chr. %s:%s' % (chromoA, chromoB),
                  size=(2000/nCols, 100), grid=(int(g//nCols), g % nCols))
   
    g += 1

def graphContactSeqSep(nuc, groupName='singleCell', regions=None):

  popup = BasePopup(title='Graph popup')  
  #xAxis = GraphAxis('Seq. sepratation Log10(bp)', labels=None, ticks=True, valRange=(4.0, 8.0))
  xAxis = GraphAxis('Seq. sepratation Log10(bp)', labels=None, ticks=True, valRange=(4.3, 8.1))
  yAxis = GraphAxis('Log10 proportion', labels=None, ticks=True)
  step = 0.1
  step = 0.05
  
  dataSets = []
  contactDict = nuc.getContacts(groupName, cis=True, trans=False)
  
  if not contactDict:
    contactDict = nuc.getContacts(groupName, cis=True, trans=False)
  
  chromosomes = nuc.getDisplayedChromosomes()
  allseps = []
  for chromo in chromosomes:
    key = (chromo, chromo)
  
    if key not in contactDict:
     continue
     
    if chromo not in chromosomes:
      continue
    
    contacts = array(contactDict[key], int)
    
    minVal = contacts[:1].min()
    maxVal = contacts[:1].max()
    dVal = float(maxVal-minVal+1)
    
    name = key[0]
    seps = abs(contacts[0]-contacts[1])
    counts = contacts[2]

    indices = seps.nonzero()
    seps    = seps[indices]
    counts  = counts[indices]
    n = float(counts.sum())
    
    fracs = (dVal-seps)/dVal # fraction of chromosome that could give rise to each separation
    #seps  = log10(seps)
    seps  = seps / 1e6
    
    histDict = {}
    for i, sep in enumerate(seps):
      bin = int(sep/step)
      c = counts[i] / fracs[i]
      
      if bin in histDict:
        histDict[bin] += c
      else:
        histDict[bin] = c
    
    xyValues = []
    for bin in sorted(histDict):
      
      x = (bin*step) + step/2.0
      x = log10(x * 1e6)
      y = log10(histDict[bin]/n)
      xyValues.append((x,y))
        
    if not xyValues:
      continue
    
    color = nuc.getChromoColor(chromo).rgbHex
    ds = GraphDataSet(xyValues, name, color, plotType='line', symbol='circle', symbolSize=1)
    dataSets.append(ds)
    
  graph = Graph(popup, (xAxis, yAxis), dataSets, title='Contact sequence separation')  
 

def graphContactDistances(nuc, groupName, box_color='#0000A0', min_num=1, structure=None, showTrans=False):
  print("graphContactDistances")
  
  chromos = nuc.getChromosomes()
  chromo_idx = {c:i for i, c in enumerate(chromos)}
  
  n_chromo = len(chromos)
  n_cols = n_chromo
  n_rows = n_chromo
 
  fig, axarr = plt.subplots(1, 2)  #plt.figure()
  fig.suptitle('Contact 3D distances (%d chromosomes)' % (len(chromos)))  
  
  #ax = fig.add_subplot(1,1,1)
  ax = axarr[0]
  
  print("Fetch contacts")
  
  contDict = nuc.getCachedContacts(groupName)
  models = range(nuc.getNumModels(structure))
  
  data = []
  labels = []
  d_max = 0.0
  particGroup = nuc._getParticleGroup(structure)
  bead_sep = None
  
  for chr_a in chromos:
    if chr_a not in contDict:
      continue
    
    if len(particGroup[chr_a]['positions']) < 3:
      continue
    
    if bead_sep is None:
      bead_pos = array(particGroup[chr_a]['positions'])
      bead_sep = min(bead_pos[1:]-bead_pos[:-1])
    
    for chr_b in chromos:
      if (not showTrans) and (chr_a != chr_b):
        continue
      
      if chr_b not in contDict[chr_a]:
        continue
      
      if len(particGroup[chr_b]['positions']) < 3:
        continue

      contacts = contDict[chr_a][chr_b]
      
      if chr_a == chr_b:
        bead_a = array(contacts[0]/bead_sep, int)
        bead_b = array(contacts[1]/bead_sep, int)
        diffs = abs(bead_a-bead_b)
        idx = (diffs > 2).nonzero()[0]
        
        contacts = contacts[:,idx]
        
      if contacts.shape[1] < min_num:
        continue
      
      chrPosA = [(chr_a, pos) for pos in contacts[0]]
      chrPosB = [(chr_b, pos) for pos in contacts[1]]

      dists = nuc.getPositionDistances(chrPosA, chrPosB, models, structure=structure)
      
      if dists is None:
        continue
      
      print('%2s - %2s %.3f %.4f' % (chr_a, chr_b, dists.mean(), dists.std()))
      
      d_max = max(d_max, dists.max()) 
      data.append(dists)
      
      if showTrans:
        labels.append('%s:%s' % tuple(sorted((chr_a, chr_b))))
      else:
        labels.append('%s' % chr_a)
  
  
  all_data = concatenate(data)

  boxprops   = {'linewidth':2, 'color':box_color}
  flierprops = {'marker':'.', 'color':'#606060', 'markersize':2,
                'linestyle':'none'}
  meanprops = dict(marker='x', linewidth=2, markeredgecolor='black',
                   markerfacecolor='black')
  whiskerprops = {'linestyle':'-', 'linewidth':2, 'color':box_color}
  medianprops = {'linewidth':1, 'color':'black'}
  capprops= {'linewidth':2, 'color':box_color}
  
  positions= 10 + arange(len(data)) * 10

  ax.vlines(positions, 0, 10, color='#C0C0FF', alpha=0.5)
  ax.boxplot(data, positions=positions, whis=[5,95], widths=5.0,
             showmeans=True, bootstrap=100, boxprops=boxprops,
             flierprops=flierprops, medianprops=medianprops, meanprops=meanprops,
             whiskerprops=whiskerprops, capprops=capprops)

  ax.set_xticklabels(labels, rotation=90)
  ax.xaxis.set_ticks(positions)
  ax.set_xlim((0, 10*len(data)+5))
  ax.set_ylim((0, 10))
 
  ax = axarr[1]
  ax.set_title('Contact distance distribution')
  ax.set_xlabel('Distance (bead radii)')
  ax.set_ylabel('Probability mass')
      
  hist, edges = histogram(all_data, bins=100, range=(0, 10), normed=True) 
  
  viol = all_data[(all_data > 4.0).nonzero()]
  pc = 100.0 * len(viol) / float(len(all_data))
  
  print('%% > 4 bead radii : %.2f' % pc)
  
  ax.plot(edges[:-1], hist)

  ax.set_xlim((0, 10))
  ax.set_ylim((0.0, 0.7))
  
  plt.show()


def graphViolations(nuc):

  popup = BasePopup(title='Graph popup')
  xAxis = GraphAxis('% Violation', labels=None, ticks=True)
  yAxis = GraphAxis('Fraction', labels=None, ticks=True,)
  
  dataSets = []
  
  allChromos = nuc.getChromosomes()    
  chromosomes = nuc.getDisplayedChromosomes()
  
  violDict = nuc.calcRestraintViolations(chromosomes, cis=True, trans=False, upperOnly=True)
  #nRest = nuc.getNumRestraints(chromosomes, cis=True, trans=False)
  
  keys = sorted(violDict.keys())
  
  for key in keys:
    chrA, chrB = key
    name = ' '.join(sorted(set(key)))
    vals = numpy.array([100.0* x[2]/x[3] for x in violDict[key] if x[3]]) # fractional 
    
    if not len(vals):
      continue
    
    nBins = 25 # int(vals.max()/10)
    hist, edges = numpy.histogram(vals, nBins, (0, 400))
    hist = numpy.array(hist, float)
    hist /= hist.sum()
    xyValues = [((edges[i]+edges[i+1])/2.0, hist[i]) for i in range(len(hist))]
    
    color = nuc.getChromoColor(chrA).rgbHex
    
    ds = GraphDataSet(xyValues, name, color, plotType='line', symbol='circle', symbolSize=3)
    dataSets.append(ds)
    
  graph = Graph(popup, (xAxis, yAxis), dataSets, size=(800, 300),
                          title='Contact violation') 

    
def graphTransSpacing(nuc, groupName='singleCell'): 

  popup = BasePopup(title='Graph popup')
  xAxis = GraphAxis('Log Separation (kb)', labels=None, ticks=True)
  yAxis = GraphAxis('Count', labels=None, ticks=True,)
  dataSets = []
  
  contactDict = nuc.getContacts(groupName, cis=False)
  posDict = {}
  
  for pairKey in contactDict:
    chrA, chrB = pairKey
    contacts = contactDict[pairKey]
  
    posA = contacts[:,0]
    posB = contacts[:,1]
    
    if chrA in posDict:
      posDict[chrA].update(posA)
    else:
      posDict[chrA] = set(posA)

    if chrB in posDict:
      posDict[chrB].update(posB)
    else:
      posDict[chrB] = set(posB)
  
  vals = numpy.array([])
  rVals = numpy.array([])
  for chromo in posDict:
    pos = numpy.array(sorted(posDict[chromo]), int)
    deltas = pos[1:] - pos[:-1]
    vals = numpy.append(vals, deltas)
    
    start, end = nuc.getChromosomeLimits(chromo)
    rPos = numpy.random.uniform(start, end, len(pos)).tolist()
    rPos = numpy.array(sorted(rPos))
    
    rDeltas = rPos[1:] - rPos[:-1]
    rVals = numpy.append(rVals, rDeltas)
       
  
  vals = numpy.array(vals, float) / 1000.0
  vals = numpy.log(vals)
  hist, edges = numpy.histogram(vals, 40)
  hist = numpy.array(hist, float) / hist.sum()
  xyValues = [((edges[i]+edges[i+1])/2.0, hist[i]) for i in range(len(hist))]
  
  ds = GraphDataSet(xyValues, 'Observed', '#008000', plotType='line', symbol='circle', symbolSize=2)
  dataSets.append(ds)
          

  rVals = numpy.array(rVals, float) / 1000.0 
  rVals = numpy.log(rVals)
  hist, edges = numpy.histogram(rVals, 40)
  hist = numpy.array(hist, float) / hist.sum()
  xyValues = [((edges[i]+edges[i+1])/2.0, hist[i]) for i in range(len(hist))]
  
  ds = GraphDataSet(xyValues, 'Random', '#800000', plotType='line', symbol='circle', symbolSize=2)
  dataSets.append(ds)

  graph = Graph(popup, (xAxis, yAxis), dataSets, title='Trans contact sequence spacing') 
  
  # Vs random - uniform?
  
  # TBD num trans partners
  # Chi-square on exp trans contacts

