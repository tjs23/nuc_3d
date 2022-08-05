
from NucApi import readContacts

petFile = 'data/tanayPairs/mouseTh1/TN_SC_1.fend_pairs.gz'

nuc = readContacts(petFile, 'Test.nuc')
print('ID:', nuc.getId())
print('Created:', nuc.getCreationTime())
print('Accessed:', nuc.getAccessTime())

chromos = nuc.getChromosomes()
print('Chromosomes:', chromos)
print('ChrLimits:', nuc.getChromosomeLimits(chromos[0]))
print('ChrPos:', nuc.getParticlePositions(chromosomes=chromos[:1]))
print('NModels:', nuc.getNumModels())
print('NConts All:', nuc.getNumContacts(chromos, cis=True, trans=True))
print('NConts Cis:', nuc.getNumContacts(chromos, cis=True, trans=False))
print('NConts Trans:', nuc.getNumContacts(chromos, cis=False, trans=True))

"""

#chromo = chromos[0]
restr = nuc.getRestraints(chromosomes=chromos, cis=True, trans=False)
for key in restr:
  print 'Cis Restraints', key, len(restr[key])
 
# Backbone getRestraints key not a tuple, check calcStruc


for chromo in chromos:
  mapa = nuc.getMapability(chromo)[:,1]
  print 'Mapability', chromo, mapa.min(), mapa.mean(), mapa.max()

conts = nuc.getContacts(chromos, cis=True, trans=True, singleCell=True)
for key in conts:
  print 'Contacts:', key, len(conts[key]), conts[key][0]


position = ('X', 73700000)
conts = nuc.getCloseContacts(position, threshold=int(1e6), cis=True, trans=True, singleCell=True)
for chromo, pos, nObs in conts:
  print 'Close', chromo, pos, nObs

# Slow ATM
#shuffDict = nuc.getShuffledContacts(chromos, cis=True, trans=True, exponent=-1.0)
#for key in shuffDict:
#  print 'Shuffled:', key, len(shuffDict[key])

mat = nuc.getContactMatrix('X', 'X',binSize=int(1e6))
print 'Contact matrix:', mat

# Slow ATM
#distMat = nuc.getDistanceMatrix('X', models=None, binSize=int(1e6))
#print 'Dist matrix:', distMat

position = ('X', 73700000)
dists = nuc.getCisDistances(position, models=None, step=int(1e6))
print 'Cis dists', dists[:10]

model = 0
dMax, dMin = nuc.getModelSize(model, chromosomes=None)
print 'Model min, max size:', dMax, dMin


coords = nuc.getModelCoords(model, chromosomes=None, backbone=None)
print 'Model coords', coords.shape

bbox = nuc.getModelBoundingBox(model, chromosomes=None)
print 'Model bbox', bbox

positionsA = [('4', i*10000000) for i in range(11)]
coords = nuc.getPositionCoords(model, positionsA)
print 'Position coords:', coords

positionsB = [('3', i*10000000) for i in range(11)]
posDists = nuc.getPositionDistances(positionsA, positionsB, models=None)
print 'Position dists:', posDists

ran = nuc.getRandomCoords(model, '1', 10)
print 'Random coords:', ran

bbSep = nuc.getBackboneSpacing()
print 'Backbone spacing:', bbSep

fileName = 'NucTest.pdb'
nuc.exportPdbFile(fileName, chromosomes=['1','4','7','10'], scale=0.1)

fileName = 'NucText.tsv'
nuc.exportContacts(fileName, 'PFE', 'singleCell', chromosomes=['1','4','7','10'], cis=True, trans=True)

fileName = 'NucTestContacts.png'
img = nuc.exportContactImage(fileName, format='PNG', chromosomes=['1','4','7','10'], binSize=int(5e6))
img.show()

#nuc.exportRmsdMatrixImage(self, fileName, format='PNG', models=None, backboneOnly=FalseexportDistanceMatrixImage(self, fileName, chromosome, format='PNG',
#                              models=None, binSize=int(1e6), minValue=True)
                              

#fileName = 'NucTestDists.png'
#img = nuc.exportDistanceMatrixImage(fileName, '1', format='PNG', models=None, binSize=int(2e6), useMinVal=False)
#img.show()

code = 'H3K4me3'
filePath = '/home/tjs23/chromoVista/data/paperCells/layers/GSM758726_CD4_ChIP_H3K4me3.dat'
nuc.importGenomeData(filePath, code, format=None)

dLayer = nuc.getGenomeData(code, chromosomes=None)
for chromo in dLayer:
  regionData, valueData, strands, annotations = dLayer[chromo]
  print 'Genome data "%s":' % chromo, len(valueData)

dRegions = nuc.getGenomeDataRegions(code, 'X')
print 'Denome data regions:', dRegions[:10]

fileName = 'NucTestGenomeData.tsv'
nuc.exportGenomeData(code, fileName, format='tsv', chromosomes=None)


# ! Make cRestraints not implemented !
#nuc.calcStructure(chromosomes=None, numModels=1, tMax=5000.0, tMin=10.0,
#                  tSteps=1000, mdSteps=10, randomStart=True, randomWalk=False,
#                  trans=True, callback=None, randomSeed=None, randomWidth=400.0,
#                  backboneSep=1.0/5e5)

#code = 'SurfA'
#nuc.calcSurface(model, 'X', code=code, sigma=1.4, cubeSize=1.0, level=0.001)

#voxCoords = nuc.getVoxelCoords(code)
#print 'Voxel coords:', voxCoords.shape

chromosomes=['1','4','7','10']
nuc.calcDepths([model,], chromosomes)
dLayer = nuc.getGenomeData('depth', chromosomes=None)
for chromo in dLayer:
  regionData, valueData, strands, annotations = dLayer[chromo]
  print 'Depths "%s":' % chromo, valueData.min(), valueData.mean(), valueData.max()

nuc.calcDensity(model, chromosomes)
dLayer = nuc.getGenomeData('density', chromosomes=None)
for chromo in dLayer:
  regionData, valueData, strands, annotations = dLayer[chromo]
  print 'Density "%s":' % chromo, valueData.min(), valueData.mean(), valueData.max()


violDict = nuc.calcRestraintViolations(chromosomes, models=None, cis=True, trans=True, upperOnly=False)
for key in violDict:
  print 'Viols:', key, len(violDict[key])

#nuc.calcModelRmsd(models=None, chromosomes=None, backbone=None, weightThreshold=1.0)

#nuc.calcModelRmsdMatrix(models=None, backbone=None, weightThreshold=1.0)

#nuc.calcModelClusters(k=4, backbone=None)

inCode = 'H3K4me3'
outCode = 'BinnedH3K4me3'
nuc.calcBinnedGenomeData(inCode, outCode, 2e6, start=None)
dLayer = nuc.getGenomeData(outCode, chromosomes=None)
for chromo in dLayer:
  regionData, valueData, strands, annotations = dLayer[chromo]
  print 'Binned data layer "%s":' % chromo, len(valueData)

#nuc.addGenomeData(code, regionDict, valueDict, annoDict=None, stranded=None)

nuc.removeGenomeData(inCode)
print 'Remove genome data', nuc.getGenomeData(code, chromosomes=None)

chromosomes=['1','4','7','10']
# # # # # # # # # # ! Below needs to be quicker !
#nuc.selectSupportedContacts(chromosomes, threshold=int(2e6))

nuc.selectCisContacts(chromosomes, replace=True, singleCell=True)

nuc.selectTransContacts(chromosomes, replace=True, singleCell=True)

nuc.selectContacts(chromosomes, singleCell=True)

nuc.selectAllContacts(singleCell=True)

centers = nuc.modelCentre(models=None, chromosomes=None)
print centers

models = [0]
chromosomes = ['X',]
nuc.modelAlign(models, chromosomes, backbone=None, weightThreshold=10.0)

nuc.modelAlignAxes(model, chromosomes)

factor = 2.3
nuc.modelScale(factor, models, chromosomes)

angle = 1.02
nuc.modelRotate(angle, axis=(0,0,1), models=models, chromosomes=chromosomes)

vector = (0.2, -9.7, 17.3)
nuc.modelTranslate(vector, models, chromosomes=chromosomes)

nuc.modelMirror(models=[0], chromosomes=None)

matrix = -1 * eye(3)
nuc.modelTransform(matrix, models=None, chromosomes=None)

nuc.setRandomCoords(models=[0], chromosomes=['1',], centre=(0.0,0.0,0.0),
                    randomWalk=True, randomSeed=None, maxStep=1.0)
#img = nuc.getContactImage(chromosomes=['1','X'], binSize=int(5e5))
#img.show()

nuc.setTestCoordsCoil(chromosomes=['1',], numRestraints=1024, ppTurn1=30, ppTurn2=7,
                      rad1=0.4, rad2=0.2, contact=0.06, scale=45.0)
#img = nuc.getContactImage(chromosomes=['1','X'], binSize=int(5e5))
#img.show()
                      
#nuc.setTestCoordsHilbert(chromosomes=['X',], centre=(0.0,0.0,0.0), scale=5.0)
#fileName = 'NucTestDists.png'
#img = nuc.exportDistanceMatrixImage(fileName, 'X', format='PNG', models=None, binSize=int(2e6), useMinVal=False)
#img.show()

model = 0
nuc.setRandomContacts(model, numRestraints=512, chromosomes=['X'])
img = nuc.getContactImage(chromosomes=['1','X'], binSize=int(5e5))
img.show()
# Also tests nuc.setContacts(contactDict, single=True)

#nuc.setModelChromosomeCoords(coords, chromosome, model=None)

#nuc.setAllCoords(coords, chromosomes=None)

infoDict = {'user':'tjs23', 'banana':'rama'}
nuc.setSampleInfo(infoDict)
print 'SampleInfo:', nuc.getSampleInfo()

nuc.removeChromosomes(['18','19'])

chrX = 'X'
regions = [(123,1200000), (2300000, 2761216)]
nuc.setDomains(chrX, regions)
print 'ChrDomains:', nuc.getDomains(chrX)

#nuc.setBackboneSpacing(bbSep=int(5e4)) # Words but breaks downstream

pos = nuc.getParticlePositions(chromosomes=[chrX,])
values = ones(len(pos[chrX]), float) * 0.5
nuc.setRefMapability('X', values)

nuc.setMapability(chromosomes=['X'])

model = 0
matrix = array([[2,0,0],[0,2,0],[0,0,2]])
nuc.setModelTransform(model, matrix)

tform = nuc.getModelTransform(model)
print 'Transform', tform

code = 'TestGrid'
data = array(range(125), float).reshape(5,5,5)
nuc.addVoxelGrid(code, data, origin=(0,0,0), gridSize=(1,1,1))
origin, gridSize, gridData = nuc.getVoxelGrid(code)
print 'Voxel grid A:', origin, gridSize, gridData.shape, gridData.mean()

data = array(range(64), float).reshape(4,4,4)
nuc.setVoxelGrid(code, data, origin=(0,0,0), gridSize=(2,2,1))
origin, gridSize, gridData = nuc.getVoxelGrid(code)
print 'Voxel grid B:', origin, gridSize, gridData.shape, gridData.mean()

nuc.removeVoxelGrid(code)
print 'Voxel grid C:', nuc.getVoxelGrid(code)

code = 'TestVoxCoords'
data = [(0,0,0,1.0), (2.0,7.3,-1.0,2.2), (2,2,2,2), (100.8,5.7,9.6,-1.9)]
nuc.addVoxelCoords(code, data)
coords = nuc.getVoxelCoords(code)
print 'Voxel Coords A:', coords.shape, coords.mean(axis=0)

data = [(2.0,7.3,-1.0,2.2), (100.8,5.7,9.6,-1.9)]
nuc.setVoxelCoords(code, data)
coords = nuc.getVoxelCoords(code)
print 'Voxel Coords B:', coords.shape, coords.mean(axis=0)

nuc.removeVoxelCoords(code)
coords = nuc.getVoxelCoords(code)
print 'Voxel Coords C:', coords

code = 'TestFoci'
positions = [(1.0,2.0,3.0),(6.6,7.7,8.8),(-1,-2,-3),(99,7,0.002)]
heights = (1.0, 2.0, -5.9, 3.0)
nuc.addFoci(code, positions, heights, volumes=None, widths=None, annotations=None)
data, annotations = nuc.getFoci(code)
print 'Volume foci A:', data, annotations

positions = [(1.0,2.0,3.0),(6.6,7.7,8.8),(-1,-2,-3)]
heights = (1.0, 2.0, -5.9)
nuc.setFoci(code, positions, heights, volumes=None, widths=None, annotations=None)
data, annotations = nuc.getFoci(code)
print 'Volume foci B:', data, annotations

#nuc.removeFoci(code)

"""
