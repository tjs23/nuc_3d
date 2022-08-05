from libc.math cimport exp, abs, sqrt
from numpy cimport ndarray, double_t, int_t

import numpy, cython, time, gc

#This is set such that the energy units from the previous version are preserved: 1 / (0.5 * 1007.0)
#k_B in kcal/mol:
BOLTZMANN_K = 0.0019872041

#ctypedef int_t   int
#ctypedef double_t double

class NucCythonError(Exception):

  def __init__(self, err):

    Exception.__init__(self, err)

def interpolateCoords(posDict, prevPosDict, ndarray[double, ndim=2] coords):
  """
  Interpolate x,y,z particle positions for an array of seq positions to a new  
  seq positions e.g. for a change in binned resolution.
  """
  
  cdef int i, j, i0, j0, p1, p2, n, m, d, dMin
  cdef double f, g
  cdef ndarray[int, ndim=1] positions
  cdef ndarray[int, ndim=1] prevPositions
  
  cdef int a, b
  
  # Get index offsets for each chromo

  #print(posDict, prevPosDict)
  
  i = 0
  offsets = {}
  for chrA in sorted(posDict):
    offsets[chrA] = i
    i += len(posDict[chrA])
  
  a = i
  cdef ndarray[double, ndim=2] newCoords = numpy.empty((i, 3), float)
  
  i = 0
  prevOffsets = {}
  for chrA in sorted(prevPosDict):
    prevOffsets[chrA] = i
    i += len(prevPosDict[chrA])  
  
  b = i
  
  for chrA in posDict:
    i0 = offsets[chrA]
    j0 = prevOffsets[chrA]
    positions = posDict[chrA]
    prevPositions = prevPosDict[chrA]
    n = len(positions)
    m = len(prevPositions)
    
    for i in range(n):

      #find closest old positions for coordinate interpolation
      p1 = 0
      dMin = positions[i]-prevPositions[0]
      
      for j in range(1,m):
        d = positions[i]-prevPositions[j]
        
        if abs(d) < abs(dMin):
          p1 = j
          dMin = d #closest pos
          
        elif abs(d) > abs(dMin): # Seq positions were in order
          break  
    
      if dMin == 0: #new position coincides with an old position
        p2 = p1
     
      elif dMin > 0: #new pos is above p1
        p2 = min(p1+1, m-1)
      
      else: #new pos is below p1
        p2 = p1
        p1 = max(0, p1-1)
        dMin = positions[i] - prevPositions[p1]
      
      #calculate coordinates
      if prevPositions[p2] == prevPositions[p1]:
        newCoords[i0+i, 0] = coords[j0+p1, 0]
        newCoords[i0+i, 1] = coords[j0+p1, 1]
        newCoords[i0+i, 2] = coords[j0+p1, 2]
 
      else: #interpolate
        f = <float>dMin/<float>(prevPositions[p2]-prevPositions[p1])
        g = 1.0 - f
        
        newCoords[i0+i, 0] = g * coords[j0+p1, 0] + f * coords[j0+p2, 0]
        newCoords[i0+i, 1] = g * coords[j0+p1, 1] + f * coords[j0+p2, 1]
        newCoords[i0+i, 2] = g * coords[j0+p1, 2] + f * coords[j0+p2, 2]
  
  return newCoords

#Collect the repulsion list, i.e. the list of particle pairs that are closer than a minimum distance.
#This is helpful for calculating the energy penalty term arising from an excluded volume.
#Return the number of clashing particle pairs.
# RepMax: array of the max number of things that could be repulsing.
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cdef int getRepulsionList(ndarray[int,   ndim=2] repList,
                          ndarray[double, ndim=2] coords,
                          int nCoords, int nRepMax,
                          double repDist,
                          ndarray[double, ndim=1] radii):
  
  cdef int i, j 
  cdef int n = 0
  cdef double dx, dy, dz, d2
  cdef double distLim
  cdef double distLim2 
  
  for i in range(nCoords-2):
    for j in range(i+2, nCoords):
      distLim = repDist + radii[i] + radii[j]
      distLim2 = distLim * distLim
      
      dx = coords[i,0] - coords[j,0]
      if abs(dx) > distLim:
        continue

      dy = coords[i,1] - coords[j,1]
      if abs(dy) > distLim:
        continue

      dz = coords[i,2] - coords[j,2]
      if abs(dz) > distLim:
        continue

      d2 = dx*dx + dy*dy + dz*dz
 
      if d2 > distLim2:
        continue

      repList[n,0] = i
      repList[n,1] = j
      
      n += 1
      
      #return if we reached max, something is seriously wrong anyway
      if n == nRepMax:
        return n
      
  return n


#Calculate the temperature from the kinetic energy.
#N_dof 1/2 k_B T = 1/2 Sum_i^N m_i v_i^2
#where N_dof is the number of degrees of freedom (3 * nCoords)
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cdef double getTemp(ndarray[double] masses,
                    ndarray[double, ndim=2] veloc,
                    int nCoords): 
  cdef int i
  cdef double kin = 0.0
  
  for i in range(nCoords):
    kin += masses[i] * (veloc[i,0]*veloc[i,0] + veloc[i,1]*veloc[i,1] + veloc[i,2]*veloc[i,2])

  #return 0.5 * 1007.0 * kin / (3 * nCoords)
  return kin / (3 * nCoords * BOLTZMANN_K)
  

#Count the number of violated restraints (lower or upper distances),
#and also return the root mean squared violation (in units of length):
#1/sqrt(N_r) sqrt( Sum_i^{N_r} viol^2).
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def getStats(ndarray[int,   ndim=2] restIndices,
             ndarray[double, ndim=2] restLimits,
             ndarray[double, ndim=2] coords,
             int nRest):      

  cdef int i, nViol = 0
  cdef int j, k
  cdef double viol, dmin, dmax, dx, dy, dz, r, s = 0
  
  for i in range(nRest):
    j = restIndices[i,0]
    k = restIndices[i,1]
    
    if j == k:
      continue

    dmin = restLimits[i,0]
    dmax = restLimits[i,1]
    
    dx = coords[j,0] - coords[k,0]
    dy = coords[j,1] - coords[k,1]
    dz = coords[j,2] - coords[k,2]
    r = sqrt(dx*dx + dy*dy + dz*dz)

    if r < dmin:
      viol = dmin - r
      nViol += 1
      
    elif r > dmax:
      viol = r - dmax
      nViol += 1
      
    else:
      viol = 0

    s += viol * viol

  return nViol, sqrt(s/nRest)


#integrator part 1
#given a(t)
#v(t) -> v(t+dt/2)
#x(t) -> x(t+dt)
#don't update accelerations here, because it would interfere with the thermostatting
#thermostatting: Langevin thermostat:
#dv/dt = F/m - gamma * p + R(T_ref)
#using the Liouville formulation for time reversibility
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cdef void velocity_verlet_1(ndarray[double] masses,
                       ndarray[double, ndim=2] forces,
                       ndarray[double, ndim=2] accel,
                       ndarray[double, ndim=2] veloc,
                       ndarray[double, ndim=2] coords,
                       int nCoords, double tRef,
                       double tStep, double Langevin_gamma):

  cdef int i
  cdef double rtStep, half_tStep, decay

  half_tStep = 0.5 * tStep
  rtStep = 0.5 * tStep * tStep

  #thermostatting (Langevin)
  #decay the velocities for dt/2
  decay = exp (- Langevin_gamma * half_tStep)
  for i in range(nCoords):
    for j in range(3):
      veloc[i,j] *= decay

  #v(t+dt/2) = v(t) + 1/2 dt a(t)
  for i in range(nCoords):
    veloc[i,0] += half_tStep * accel[i,0]
    veloc[i,1] += half_tStep * accel[i,1]
    veloc[i,2] += half_tStep * accel[i,2]

  #x(t+dt) = x(t) + dt v(t+dt/2)
  for i in range(nCoords):
    coords[i,0] += tStep * veloc[i,0]
    coords[i,1] += tStep * veloc[i,1]
    coords[i,2] += tStep * veloc[i,2]


#integrator part 2
#given F(t+dt)
#a(t+dt) = F(t+dt) / m
#v(t+dt/2) -> v(t+dt)
#thermostatting: Langevin thermostat:
#dv/dt = F/m - gamma * p + R(T_ref)
#using the Liouville formulation for time reversibility
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cdef void velocity_verlet_2(ndarray[double] masses,
                         ndarray[double, ndim=2] forces,
                         ndarray[double, ndim=2] accel,
                         ndarray[double, ndim=2] veloc,
                         int nCoords, double tRef,
                         double tStep, double Langevin_gamma):

  cdef int i
  cdef double half_tStep, R, decay
  
  # sample standard normal distribution in advance
  cdef ndarray[double, ndim=2] normSample = numpy.random.randn(nCoords,3)

  half_tStep = 0.5 * tStep

  #update acceleration from the new forces
  #a(t+dt) = F(t+dt) / m
  for i in range(nCoords):
    accel[i,0] = forces[i,0] / masses[i]
    accel[i,1] = forces[i,1] / masses[i]
    accel[i,2] = forces[i,2] / masses[i]

  #thermostatting (Langevin)
  #add the random acceleration
  R = Langevin_gamma * BOLTZMANN_K * tRef / half_tStep
  for i in range(nCoords):
    accel[i,0] += sqrt(R/masses[i]) * normSample[i,0]
    accel[i,1] += sqrt(R/masses[i]) * normSample[i,1]
    accel[i,2] += sqrt(R/masses[i]) * normSample[i,2]

  #v(t+dt) = v(t+dt/2) + 1/2 dt a(t+dt)
  for i in range(nCoords):
    veloc[i,0] += half_tStep * (accel[i,0])
    veloc[i,1] += half_tStep * (accel[i,1])
    veloc[i,2] += half_tStep * (accel[i,2])

  #thermostatting (Langevin)
  #decay the velocities for dt/2 again
  decay = exp (- Langevin_gamma * half_tStep)
  for i in range(nCoords):
    for j in range(3):
      veloc[i,j] *= decay


#integrator part 1
#F(t), a(t)
#v(t) -> v(t+dt/2)
#x(t) -> x(t+dt)
#old thermostatting: Berendsen thermostat (not canonical):
#dv/dt = F/m + 1/(2 tau) ( T_ref / T - 1 ) v
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cdef void updateMotion(ndarray[double] masses,
                       ndarray[double, ndim=2] forces,
                       ndarray[double, ndim=2] accel,
                       ndarray[double, ndim=2] veloc,
                       ndarray[double, ndim=2] coords,
                       int nCoords, double tRef,
                       double tStep, double beta):

  cdef int i
  cdef double r, rtStep, temp

  rtStep = 0.5 * tStep * tStep
  temp = getTemp(masses, veloc, nCoords)
  #avoid division by 0 temperature
  temp = max(temp, 0.001)
  r = beta * (tRef/temp-1.0)
  
  for i in range(nCoords):

    #TODO: what is the thermostat r*veloc 
    accel[i,0] = forces[i,0] / masses[i] + r * veloc[i,0]
    accel[i,1] = forces[i,1] / masses[i] + r * veloc[i,1]
    accel[i,2] = forces[i,2] / masses[i] + r * veloc[i,2]

    #x(t+dt) = x(t) + dt v(t) + 1/2 dt^2 a(t)   
    coords[i,0] += tStep * veloc[i,0] + rtStep * accel[i,0]
    coords[i,1] += tStep * veloc[i,1] + rtStep * accel[i,1]
    coords[i,2] += tStep * veloc[i,2] + rtStep * accel[i,2]

    #v(t+dt/2) = v(t) + 1/2 dt a(t)
    #actually, this is bigger than what we want by 1/2 dt F(t)/m    
    veloc[i,0] += tStep * accel[i,0]
    veloc[i,1] += tStep * accel[i,1]
    veloc[i,2] += tStep * accel[i,2]


#integrator part 2
#a(t) is not yet updated
#F(t+dt)
#v(t+dt/2) -> v(t+dt)
#old thermostatting: Berendsen thermostat (not canonical):
#dv/dt = F/m + 1/(2 tau) ( T_ref / T - 1 ) v
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cdef void updateVelocity(ndarray[double] masses,
                         ndarray[double, ndim=2] forces,
                         ndarray[double, ndim=2] accel,
                         ndarray[double, ndim=2] veloc,
                         int nCoords, double tRef,
                         double tStep, double beta):

  cdef int i
  cdef double r, temp

  temp = getTemp(masses, veloc, nCoords)
  #avoid division by 0 temperature
  temp = max(temp, 0.001)
  r = beta * (tRef/temp-1.0)

  #v(t+dt) = v(t+dt/2) + 1/2 dt F(t+dt)/m
  #actually, we subtract the 1/2 dt a(t), by which we overshot at the other part of the leapfog algorithm  
  for i in range(nCoords):
    veloc[i,0] += 0.5 * tStep * (forces[i,0] / masses[i] + r * veloc[i,0] - accel[i,0])
    veloc[i,1] += 0.5 * tStep * (forces[i,1] / masses[i] + r * veloc[i,1] - accel[i,1])
    veloc[i,2] += 0.5 * tStep * (forces[i,2] / masses[i] + r * veloc[i,2] - accel[i,2])


#Calculating repulsive forces for the clashing particle pair list:
#E(ij) = k ( d_lim^2 - d^2 )^2     with d = x_i - x_j
#F_i = 4 k ( d_lim^2 - d^2 ) d
#F_j = - 4 k ( d_lim^2 - d^2 ) d
#k is the force constant fConst
#d_lim is the restraint distance limit repDist
#i and j are k and j
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cdef double getRepulsiveForce(ndarray[int,   ndim=2] repList,
                              ndarray[double, ndim=2] forces,
                              ndarray[double, ndim=2] coords,
                              int nRep, double fConst,
                              ndarray[double, ndim=1] radii):

  cdef int i, j, k
  cdef double dx, dy, dz, d2, dr, rjk 
  cdef double force = 0
  cdef double repDist2 

  if fConst == 0:
    return force

  for i from 0 <= i < nRep:
    j = repList[i,0]
    k = repList[i,1]
    repDist = radii[j] + radii[k]
    repDist2 = repDist * repDist


    dx = coords[k,0] - coords[j,0]
    if abs(dx) > repDist:
      continue

    dy = coords[k,1] - coords[j,1]
    if abs(dy) > repDist:
      continue

    dz = coords[k,2] - coords[j,2]
    if abs(dz) > repDist:
      continue

    d2 = dx*dx + dy*dy + dz*dz
    if d2 > repDist2:
      continue

    dr = repDist2 - d2
    #energy contribution
    force += fConst * dr * dr
    rjk = 4 * fConst * dr

    dx *= rjk
    dy *= rjk
    dz *= rjk

    #force contributions
    forces[j,0] -= dx
    forces[k,0] += dx

    forces[j,1] -= dy
    forces[k,1] += dy

    forces[j,2] -= dz
    forces[k,2] += dz
      
  return force
    

#Calculating constraint forces for the restraint particle pair list:
#
#        k ( d_min - |r| )^2     |r| < d_min
#E(ij) = 0                       d_min <= |r| <= d_max
#        k ( |r| - d_max )^2     d_max < |r| <= d_a
#        k ( a + A|r| + b/|r| )  d_a < |r|
#
#where r = x_i - x_j  for unambiguous constraints
#      r = [ Sum (1 / | x_i - x_j |^4) ]^(-1/4)
#k is the force constant fConst
#d_min is the lower limit for the restraint
#d_max is the upper limit for the restraint
#d_a is the switching distance for the restraint function
#a, b and A are parameters of an asymptotic function
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cdef double getRestraintForce(ndarray[double, ndim=2] forces,
                              ndarray[double, ndim=2] coords,
                              ndarray[int,   ndim=2] restIndices,
                              ndarray[double, ndim=2] restLimits,
                              ndarray[int, ndim=1] restAmbig,
                              int nRest, double fConst, double exponent=2.0,
                              double distSwitch=0.5, double asymptote=1.0):
               
  cdef int i, j, k, n, nAmbig
  cdef double a, b, da, d, dmin, dmax, dx, dy, dz
  cdef double r, r2, s2, rjk, ujk, force = 0, t

  b = asymptote*distSwitch*distSwitch - exponent*distSwitch*distSwitch*distSwitch
  # b = 0
  a = distSwitch*distSwitch - asymptote*distSwitch - b/ distSwitch
  # a = -0.25
  
  i = 0
  
  while i < nRest:

    # 1. calculate constrained variable r^2
    nAmbig = restAmbig[i]
    r2 = 0.0

    #unambiguous constraint
    # r^2 = | x_i - x_j |^2
    if nAmbig == 1:
      j = restIndices[i,0]
      k = restIndices[i,1]
      
      if j != k:
        dx = coords[j,0] - coords[k,0]
        dy = coords[j,1] - coords[k,1]
        dz = coords[j,2] - coords[k,2]
        r2 = dx*dx + dy*dy + dz*dz
    
    #ambiguous constraint
    # r^2 = 1 / sqrt ( Sum (1 / | x_i - x_j |^4) )
    # if one i-j distance is much shorter, it will dominate and r^2 approx r_ij^2 = | x_i - x_j |^2
    else:
      for n in range(nAmbig):
        j = restIndices[i+n,0]
        k = restIndices[i+n,1]
        
        if j != k:
          dx = coords[j,0] - coords[k,0]
          dy = coords[j,1] - coords[k,1]
          dz = coords[j,2] - coords[k,2]
          r = dx*dx + dy*dy + dz*dz
          r2 += 1.0 / (r * r)

      if r2 > 0:
        r2 = 1.0 / sqrt(r2)

    #skip (ambiguous) constraint(s) for which (all) particle pairs coincide
    if r2 <= 0:
      i += nAmbig
      continue


    # 2. calculate the constraint force and energy contributions
    dmin = restLimits[i,0]
    dmax = restLimits[i,1]

    da = dmax + distSwitch

    if r2 < dmin*dmin:
      r2 = max(r2, 1e-8)
      r = sqrt(r2)
      d = dmin - r
      ujk = fConst * d * d
      rjk = fConst * 2 * d
    
    elif r2 > dmax*dmax:
      r = sqrt(r2)
      d = r - dmax

      if r <= da:
        ujk = fConst * d * d
        rjk = - fConst * 2 * d

      else:
        ujk = fConst * (a + asymptote*d + b/d)
        rjk = - fConst * (asymptote - b/(d*d))

    else:
      ujk = rjk = 0
      r = 1.0

    #the energy contribution
    force += ujk

    #the force contributions for unambiguous constraints
    if nAmbig == 1:
      j = restIndices[i,0]
      k = restIndices[i,1]

      if j == k:
        i += nAmbig
        continue
      
      t = rjk / r
      dx = coords[j,0] - coords[k,0]
      dy = coords[j,1] - coords[k,1]
      dz = coords[j,2] - coords[k,2]
      
      dx *= t
      dy *= t
      dz *= t
      
      forces[j,0] += dx
      forces[j,1] += dy
      forces[j,2] += dz
      
      forces[k,0] -= dx
      forces[k,1] -= dy
      forces[k,2] -= dz

    #the force contributions for ambiguous constraints
    else:

      for n in range(nAmbig):
        j = restIndices[i+n,0]
        k = restIndices[i+n,1]
        
        if j == k:
          continue
          
        dx = coords[j,0] - coords[k,0]
        dy = coords[j,1] - coords[k,1]
        dz = coords[j,2] - coords[k,2]
        
        s2 = dx*dx + dy*dy + dz*dz
        t = rjk * r2 * r2 * r / (s2 * s2 * s2)
        
        dx *= t
        dy *= t
        dz *= t
        
        forces[j,0] += dx
        forces[k,0] -= dx
        forces[j,1] += dy
        forces[k,1] -= dy
        forces[j,2] += dz
        forces[k,2] -= dz
      
    i += nAmbig
    
    
  return force

#Remove sanity checks for speed
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
#NVT MD simulation.
#This is the lower level function for annealStructure.
#Configuration parameters:
#  coord particle coordinates
#  masses particle masses
#Constraint parameters:
#  fConstR repulsive interaction force constant
#  fConstD distance constraint force constant
#  restIndices the index pairs of constrained particle pairs
#  restLimits the distance limits for all particle pair constraints separately
#  restAmbig ambiguity mark (1:unambiguous, >1 how many of this and the coming ones belong together)
#Repulsive interaction parameters:
#  repDist helper distance for repulsion list building
#  minDist distance lower limit for repulsive interactions
#Temperature parameters:
#  beta is 1/2 tau^{-1} of the Berendsen thermostat
#  tRef reference or target temperature
#  tot0 TODO ??
#Other parameters:
#  tStep time step
#  nSteps number of steps to be taken
#  tTaken run time of the function
#  printInterval ouput frequency
#  use_Langevin 1: use the Langevin thermostat with the new integrator (default, this should be used),
#              0: Berendsen thermostat with the old integrator
#  LangevinTau_per_tStep the time constant of the Langevin thermostat, as a multiple of tStep
def runDynamics(ndarray[double, ndim=2] coords,
                ndarray[double, ndim=1] masses,
                ndarray[double, ndim=1] radii,
                ndarray[int, ndim=2] restIndices,
                ndarray[double, ndim=2] restLimits,
                ndarray[int, ndim=1] restAmbig,
                double tRef=1000.0, double tStep=0.001, int nSteps=1000,
                double fConstR=1.0, double fConstD=25.0, double beta=10.0,
                double minDist=2.25, double repDist=2.0,
                double tTaken=0.0, int printInterval=10000,
                double tot0=20.458, int use_Langevin=1,
                double LangevinTau_per_tStep=100.0,
                int nRepMax=0):
  
  
  #Error checking
  
  cdef int nRest = len(restIndices)
  cdef int nCoords = len(coords)
    
  if nCoords < 2:
    raise NucCythonError('Too few coodinates')
     
  indices = set(restIndices.ravel())
  if min(indices) < 0:
    raise NucCythonError('Restraint index negative') 
  
  if max(indices) >= nCoords:
    data = (max(indices), nCoords)
    raise NucCythonError('Restraint index "%d" out of bounds (> %d)' % data) 
  
  if nCoords != len(masses):
    raise NucCythonError('Masses list size does not match coordinates') 
 
  if nRest != len(restLimits):
    raise NucCythonError('Number of restraint index pairs does not match number of restraint limits')   

  if len(restAmbig) != nRest:
    raise NucCythonError('Size of ambiguity list does not match number of restraints')   

  if use_Langevin != 1 and use_Langevin != 0:
    raise NucCythonError('Unknown value for use_Langevin, has to be 1 (Langevin) or 0 (Berendsen)')   
  
  cdef int i, j, n, step, nViol, nRep = 0
  
  if not nRepMax:
    nRepMax = nCoords*800
    
  cdef double d2, maxDelta, dx, dy, dz, ek, rmsd, tStep0, temp, fDist, fRep
  # cdef double distLim = repDist + minDist
  cdef double deltaLim = 0.25 * repDist * repDist
  cdef double Langevin_gamma

  #Thermostat initialisation
  tStep0 = tStep * tot0  # ~0.02
  if use_Langevin == 1:
    #the default time constant is ~ (100 * time step used), to be consistent with the default of the old version
    Langevin_gamma = 1.0 / (LangevinTau_per_tStep * tStep0)  # default gamma = 2.0
    #print "using Langevin thermostat with gamma = ", Langevin_gamma
  else:
    beta /= tot0  # 10.0 / 20.458 = ~0.5, beta^{-1} = 2.0
    #print "using Berendsen thermostat with beta = ", beta


  #Velocity initialization

  veloc = numpy.random.normal(0.0, 1.0, (nCoords, 3))
  veloc *= sqrt(tRef / getTemp(masses, veloc, nCoords))

  #if printInterval > 0:
  #  print('Temp: %7.3f Max rep: %d' % (tRef, nRepMax))

  #Other initialisation
  
  zeroArray = numpy.zeros((nCoords, 3))
  repList = numpy.zeros((nRepMax, 2), numpy.int32)
  
  cdef ndarray[double, ndim=2] coordsPrev = numpy.array(coords)  # typing this array important for speed 
  cdef ndarray[double, ndim=2] accel = numpy.array(zeroArray)    # this less so
  cdef ndarray[double, ndim=2] forces
  
  t0 = time.time()

  #Main loop
    
  for step in range(nSteps):

    #update repulsion list    
    if step == 0:
      #calculate repulsion list from scratch
      nRep = getRepulsionList(repList, coords, nCoords, nRepMax, repDist, radii) # Handle errors
      
      for i in range(nCoords):
        coordsPrev[i,0] = coords[i,0]
        coordsPrev[i,1] = coords[i,1]
        coordsPrev[i,2] = coords[i,2]

      #forces at step = 0
      forces = numpy.array(zeroArray)
      fRep = getRepulsiveForce(repList, forces, coords, nRep,  fConstR, radii)
      fDist = getRestraintForce(forces, coords, restIndices, restLimits, restAmbig, nRest, fConstD)
      #initial a(0) = F(0) / m
      for i in range(nCoords):
        accel[i,0] = forces[i,0] / masses[i]
        accel[i,1] = forces[i,1] / masses[i]
        accel[i,2] = forces[i,2] / masses[i]
   
    else:
      #calculate maximum displacement
      maxDelta = 0.0
      for i in range(nCoords): 
        dx = coords[i,0] - coordsPrev[i,0]
        dy = coords[i,1] - coordsPrev[i,1]
        dz = coords[i,2] - coordsPrev[i,2]
        d2 = dx*dx + dy*dy + dz*dz
        
        if d2 > maxDelta:
          maxDelta = d2
          
          if maxDelta > deltaLim:
            break            

      #recalculate repulsion list only if a particle has moved more than delta_lim
      if maxDelta > deltaLim:
        nRep = getRepulsionList(repList, coords, nCoords, nRepMax, repDist, radii) # Handle errors
        
        for i in range(nCoords):
          coordsPrev[i,0] = coords[i,0]
          coordsPrev[i,1] = coords[i,1]
          coordsPrev[i,2] = coords[i,2]
 
    #integrator 1
    #v(t) -> v(t+dt/2)
    #x(t) -> x(t+dt)
    if use_Langevin == 1:
      velocity_verlet_1(masses, forces, accel, veloc, coords, nCoords, tRef, tStep0, Langevin_gamma)
    else:
      updateMotion(masses, forces, accel, veloc, coords, nCoords, tRef, tStep0, beta)
    
    #update energy and forces
    for i in range(nCoords):      # Removed forces = numpy.array(zeroArray) - caused Python GC issues
      forces[i,0] = 0.0
      forces[i,1] = 0.0
      forces[i,2] = 0.0
      
    fRep  = getRepulsiveForce(repList, forces, coords, nRep, fConstR, radii)
    fDist = getRestraintForce(forces, coords, restIndices, restLimits, restAmbig, nRest, fConstD)

    #integrator 2
    #update a(t+dt)
    #v(t+dt/2) -> v(t+dt)
    if use_Langevin == 1:
      velocity_verlet_2(masses, forces, accel, veloc, nCoords, tRef, tStep0, Langevin_gamma)
    else:
      updateVelocity(masses, forces, accel, veloc, nCoords, tRef, tStep0,  beta)

    #output
    if (printInterval > 0) and step % printInterval == 0:
      temp = getTemp(masses, veloc, nCoords)
      nViol, rmsd = getStats(restIndices, restLimits, coords, nRest)

      #if step == 0:
      #  print('t[?] T[?] E_repulsive[?] E_constraint[?] RMSD[?] n_viol_constraint n_rep')
      #print time temperature E_repulsive E_constraint RMSD num_constraint_violations num_repulsive_violation     
      
      data = (temp, fRep, fDist, rmsd, nViol, nRep)
      print('temp:%7.2lf  fRep:%7.2lf  fDist:%7.2lf  rmsd:%7.2lf  nViol:%5d  nRep:%5d' % data)

    tTaken += tStep
  
  #if printInterval > 0:
  #  print('Time: %.3f' % (time.time() - t0))
  
  return tTaken, nRep

 


