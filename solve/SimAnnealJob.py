# Tiny Python wrapper to run calculations for each model in separate processes
# using separate Python interpreters


def simAnnealjob(coords, posDicts, stageCounts, indices, dists, ambig, temps, repulsions,
                 timeSteps, dynSteps, minDist, repDist, printInterval, callback=None,
                 masses=None, radii=None):
                 
  import gc
  from numpy import ones, int32
  from solve.SimAnneal import runDynamics
  from cUtil.apiUtil import interpolateChromoModelCoords
  
  """
  Run the simulated annealing structure calculation for only one coordinate model
  """
  
  chromosomes = sorted(posDicts[0].keys())
  nStages = len(posDicts)
    
  n0 = 0
  p0 = 0
  for i in range(nStages):
    posDict = posDicts[i]
    n1 = n0 + stageCounts[i]
    
    # interpolate coordinates, as required
    if i > 0:
      nCoords = sum([len(posDict[c]) for c in chromosomes])
      if len(coords) != nCoords: # Change of scale
        coords = interpolateChromoModelCoords(posDict, posDicts[i-1], coords)
    
    p1 = p0 + len(coords)

    if masses is None:
      # Masses currently unity, but could be set for irregular particle sizes
      stageMasses = ones(len(coords), float)
    else:
      stageMasses = masses[p0:p1]
    
    if radii is None:
      stageRadii = ones(len(coords), float)
    else:
      stageRadii = radii[p0:p1]
    
    stageTemps = temps[i]
    stageRepuls = repulsions[i]
    stageIndices = indices[n0:n1]
    stageDists = dists[n0:n1]
    stageAmbig = ambig[n0:n1]    
    nRepMax = int32(0)
    
    # run the annealing schedule  
    for j, temp in enumerate(stageTemps):
      gc.collect()
      t, nRepMax = runDynamics(coords, stageMasses, stageRadii, stageIndices, stageDists, stageAmbig,
                               temp, timeSteps[i], dynSteps[i], stageRepuls[j],
                               minDist=minDist, repDist=repDist, printInterval=printInterval,
                               nRepMax=nRepMax)
      
      nRepMax = int32(1.05 * nRepMax)
      
      if callback:
        callback(coords, posDict, chromosomes, i, j)
    
    gc.collect()
    n0 = n1
    p0 = p1
  
  return coords


def simAnnealjobCluster(coords, posDicts, stageCounts, indices, dists, ambig, temps, repulsions,
                        timeSteps, dynSteps, minDist, repDist, printInterval, callback=None,
                        masses=None, radii=None):
                 
  import gc
  from numpy import ones, int32
  from SimAnneal import runDynamics, interpolateCoords
  
  """
  Run the simulated annealing structure calculation for only one coordinate model
  """
  
  chromosomes = sorted(posDicts[0].keys())
  nStages = len(posDicts)
    
  n0 = 0
  p0 = 0
  for i in range(nStages):
    posDict = posDicts[i]
    n1 = n0 + stageCounts[i]
    
    # interpolate coordinates, as required
    if i > 0:
      nCoords = sum([len(posDict[c]) for c in chromosomes])
      if len(coords) != nCoords: # Change of scale
        coords = interpolateCoords(posDict, posDicts[i-1], coords)
    
    p1 = p0 + len(coords)

    if masses is None:
      # Masses currently unity, but could be set for irregular particle sizes
      stageMasses = ones(len(coords), float)
    else:
      stageMasses = masses[p0:p1]
    
    if radii is None:
      stageRadii = ones(len(coords), float)
    else:
      stageRadii = radii[p0:p1]
    
    stageTemps = temps[i]
    stageRepuls = repulsions[i]
    stageIndices = indices[n0:n1]
    stageDists = dists[n0:n1]
    stageAmbig = ambig[n0:n1]    
    nRepMax = int32(0)
    
    # run the annealing schedule  
    for j, temp in enumerate(stageTemps):
      gc.collect()
      t, nRepMax = runDynamics(coords, stageMasses, stageRadii, stageIndices, stageDists, stageAmbig,
                               temp, timeSteps[i], dynSteps[i], stageRepuls[j],
                               minDist=minDist, repDist=repDist, printInterval=printInterval,
                               nRepMax=nRepMax)
      
      nRepMax = int32(1.05 * nRepMax)
      
      if callback:
        callback(coords, posDict, chromosomes, i, j)
    
    gc.collect()
    n0 = n1
    p0 = p1
  
  return coords

def chromo_anneal_file_job(in_file, out_file, model, minDist=1.5, repDist=1.5):

  import numpy as np
  
  data_dict = np.load(in_file)

  coords = data_dict['coords']
  stage_counts = data_dict['stageCounts']
  
  pos_dicts = [{} for i in range(len(stage_counts))]
  
  for key in data_dict:
    if key.startswith('pos_'):
      i, chromo = key.split('_')[1:]
      pos_dicts[int(i)][chromo] = data_dict[key]  
  
  indices = data_dict['restrIndices']
  dists = data_dict['restrDists']
  ambig = data_dict['restrAmbig']
  temps = data_dict['temps']
  repulsions = data_dict['repulsScales']
  timeSteps = data_dict['timeSteps']
  dynSteps = data_dict['dynSteps']
    
  coords = coords[model] 
  coords = simAnnealjobCluster(coords, pos_dicts, stage_counts, indices, dists, ambig, temps, repulsions,
                        timeSteps, dynSteps, minDist, repDist, printInterval=100, callback=None,
                        masses=None, radii=None)
  
  np.save(out_file, coords)
  
  # This could be run locally by Ruffus
  # This could be run by API using call()/Popen and regular multprocessing
  #  - API would collate results at end
  
  
if __name__ == '__main__':

  import sys, os
  
  from glob import glob
  
  #from os.path import abspath, dirname, join, splitext, split
  #sys.path.append(abspath(dirname(dirname(__file__))))

  args = sys.argv[1:]
  
  n_models = 20
  
  in_files = sorted(glob('P2-E8-PA_test/*_job.npz'))
  
  index = int(args[0]) - 1
  
  model = index % n_models
  
  in_file = in_files[int(index // n_models)]
  
  out_file = os.path.splitext(in_file)[0] + '_coords_%04d.npy' % model
  
  chromo_anneal_file_job(in_file, out_file, model)

"""
qsub -t 1-400 calc_chromo_struc.sh -M tstevens@lmb.internal
qsub -hold_jid 9601923 -t 401-800 calc_chromo_struc.sh -M tstevens@lmb.internal


qsub -t 1-800 calc_chromo_struc.sh -M tstevens@lmb.internal
qsub -hold_jid 9593106 -t 801-1600 calc_chromo_struc.sh -M tstevens@lmb.internal
qsub -hold_jid 9593107 -t 1601-2400 calc_chromo_struc.sh -M tstevens@lmb.internal

qsub -t 1-600 calc_chromo_struc.sh -M tstevens@lmb.internal
qsub -hold_jid  -t 601-1200 calc_chromo_struc.sh -M tstevens@lmb.internal
qsub -hold_jid  -t 1201-1800 calc_chromo_struc.sh -M tstevens@lmb.internal
qsub -hold_jid  -t 1801-2400 calc_chromo_struc.sh -M tstevens@lmb.internal

qsub -hold_jid  -t 2401-3000 calc_chromo_struc.sh -M tstevens@lmb.internal
qsub -hold_jid  -t 3001-3600 calc_chromo_struc.sh -M tstevens@lmb.internal
qsub -hold_jid  -t 3601-4200 calc_chromo_struc.sh -M tstevens@lmb.internal
qsub -hold_jid  -t 4201-4800 calc_chromo_struc.sh -M tstevens@lmb.internal

"""

# Convert a list of input .nuc files into a larger array of numbered job files, one per cell-model

# qsub handy if Python script can run off a job number

# make function work off of .npz file and generate a .npy file (coords) 
