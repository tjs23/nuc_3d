import gc
from numpy import cov, linalg, sqrt, zeros, random, dot, array, vstack, empty
from random import randint, sample

def kMeansSpread(data, k, thresh=1e-10, verbose=False):

  n = len(data)
  index = randint(0, n-1)
  indices = set([index])
  
  influence = zeros(n)
  while len(indices) < k:
    diff = data - data[index]
    sumSq = (diff * diff).sum(axis=1) + 1.0
    influence += 1.0 / sumSq
    index = influence.argmin()
    
    while index in indices:
      index = randint(0, n-1)
    
    indices.add(index)    
  
  centers = vstack([data[i] for i in indices])
    
  return kMeans(data, k, centers, thresh, verbose)


def kMeans(data, k, centers=None, thresh=1e-10, verbose=False):
  
  if centers is None:
    centers = array( sample(list(data), k) )  # list() not needed in Python 2

  labels = empty(len(data), float)
  change = 1.0
  prev = []

  j = 0
  while change > thresh:

    clusters = [[] for x in range(k)]
    for i, vector in enumerate(data):
      diffs = centers - vector
      dists = (diffs * diffs).sum(axis=1)
      closest = dists.argmin()
      labels[i] = closest
      clusters[closest].append(vector)
     
    change = 0
    for i, cluster in enumerate(clusters):
      cluster = array(cluster)
      center = cluster.sum(axis=0)/len(cluster)
      diff = center - centers[i]
      change += (diff * diff).sum()
      centers[i] = center
    
    j += 1
    
    if verbose:
      print j, change 
    
  return centers, clusters, labels
  
  
def dbScanCluster(data, threshold, minNeighbour, kdTree=None):     
  
  if kdTree is None:
    from scipy.spatial import cKDTree
    kdTree = cKDTree(data, 10)
    
  nList = kdTree.query_ball_point(data, threshold) # Euclidian L2 norm
  
  clusters = []
  clustersAppend = clusters.append
  noise = []
  noiseAppend = noise.append
  pool = set(range(len(data)))
  poolPop = pool.pop
  
  while pool:
    i = poolPop()
    neighbours = nList[i]
    
    if len(neighbours) < minNeighbour:
      noiseAppend(i)
    
    else:
      cluster = [i]

      pool2 = set(neighbours)
      while pool2:
        j = pool2.pop()
        
        if j in pool:
          pool.remove(j)
          neighbours2 = nList[i]
 
          if len(neighbours2) < minNeighbour:
            noiseAppend(j)
          
          else:  
            pool2.update(neighbours2)
            cluster.append(j)
 
      clustersAppend(cluster)
  
  nList = []
  del nList
  
  #noiseData = [data[i] for i in noise]
  #clusterData = []
  #for cluster in clusters:
  #  clusterData.append( [data[i] for i in cluster] )
         
  return clusters, noise
  
  
def extractPrincipleComponent(data, precision=1e-19):

  samples, features = data.shape
  meanVec = data.mean(axis=0)
  dataC = data - meanVec
  
  pc1 = random.random(features)
  pc0 = pc1 - 1.0
   
  while abs((pc0-pc1).sum()) > precision:

    t = zeros(features)    
    for datum in dataC:
      t += dot(datum, pc1) * datum

    pc0 = pc1
    pc1 = t / sqrt(dot(t,t))
    
  return pc1


def kMedioids(distMatrix, centers, k):
  
  n = len(distMatrix)

  change = 1e99

  while change > 1e-8:
    clusters = [[] for x in range(k)]
    for i in range(n):
      bestDist = None
      
      for j in range(k):
        center = centers[j]
        dist = distMatrix[i,center]
        
        if (bestDist is None) or (dist < bestDist):
          bestDist = dist
          closest = j
      
      clusters[closest].append(i)
 
    change = 0.0
    for j, cluster in enumerate(clusters):
      oldCenter = centers[j]
      
      clusterSim = []
      for x in cluster:
        sim = 0.0
        for y in cluster:
          sim += distMatrix[x,y] ** 4
        clusterSim.append(sim)
      
      idx = clusterSim.index(min(clusterSim))
      center = cluster[idx]
      
      change += distMatrix[oldCenter, center]
      centers[j] = center
     
  return centers, clusters


def principleComponentAnalysis(data, n=3):

  samples, features = data.shape

  meanVec = data.mean(axis=0)
  dataC = (data - meanVec).T

  covar = cov(dataC)
  evals, evecs = linalg.eig(covar)

  indices = evals.argsort()[::-1]

  evecs = evecs[:,indices]

  basis = evecs[:,:n]
  energy = evals[:n].sum()

  return basis, energy


def _getJoinPair(distMatrix):

  n = len(distMatrix)

  minQ = None
  joinPair = None

  for i in range(n-1):
    sumRow = sum(distMatrix[i])
  
    for j in range(i+1, n):
      sumCol = sum(distMatrix[j])
    
      dist = distMatrix[i][j]
      q = (n-2)*dist - sumRow - sumCol
      
      if (minQ is None) or (q < minQ):
        minQ = q
        joinPair = [i,j]
  
  joinPair.sort()

  return joinPair
  
def _getDistToJunction(distMatrix, i, j):
  
  n = len(distMatrix)
  row = distMatrix[i]
  column = distMatrix[j]
 
  dist = distMatrix[i][j] + (sum(row)-sum(column))/(n-2)
  dist *= 0.5

  return dist
  
def neighbourJoinTree(distMatrix):
 
  joinOrder = []
  n = len(distMatrix)
  tree = range(n)
  
  while n > 2:

    x, y = _getJoinPair(distMatrix)

    node = (tree[x], tree[y])
    joinOrder.append(node)
    tree.append(node)

    del tree[y]
    del tree[x]
     
    distX = _getDistToJunction(distMatrix, x, y)
    distY = _getDistToJunction(distMatrix, y, x)
  
    distMatrix.append([0] * (n+1))

    for i in range(n):
      if i not in (x,y):

        dist = (distMatrix[x][i]-distX) + (distMatrix[y][i]-distY)
        dist *= 0.5
  
        distMatrix[i].append(dist)
        distMatrix[n][i] = dist

    del distMatrix[y]
    del distMatrix[x]
  
    for row in distMatrix:
      del row[y]
      del row[x]

    n -= 1

  tree = tuple(tree)
  joinOrder.append(tree)
  
  return tree, joinOrder    

def _hierarchicalRowCluster(dataMatrix):
  
  n = len(dataMatrix)
  distanceMatrix = zeros((n, n), float)
  
  for i, row in enumerate(dataMatrix):
    diffs = dataMatrix - row
    sqDiffs = diffs * diffs
    dists = sqrt(sqDiffs.sum(axis=1))
    distanceMatrix[i,:] = dists
    
  tree, joinOrder = neighbourJoinTree(distanceMatrix.tolist())
   
  rowOrder = list(tree)
  
  i = 0
  while i < len(rowOrder):
    
    while not isinstance(rowOrder[i], int):
      rowOrder[i:i+1] = rowOrder[i]
    
    i += 1
  
  return rowOrder

def hierarchicalCluster(data):
   
  rows = _hierarchicalRowCluster(data)
  data = data[rows]
  
  dataT = data.T
  cols = _hierarchicalRowCluster(dataT)
  dataT = dataT[cols]
  data = dataT.T
  
  return data,  rows,  cols
