import sys, os, traceback
from time import time, sleep
from collections import Sequence

#os.environ['OMPI_MCA_mpi_yield_when_idle'] = '1'

from numpy import array, fromstring, empty, ndarray
from Worker import Worker

dirName = os.path.dirname(__file__)
controllerPath = os.path.join(dirName, 'Controller.py')

# Send() - blocks until buffer is safe
# Ssend() - blocks until data recieved

PICKLE = 'pickle'

class Null():
  pass

# Null output value.
# Can't use None as that may legitimately be returned.
NULL = Null()

class Job():
  """
  This class represents a parallel processing job that has been submitted to
  an Engine. A Job object will be returned from Engine.run() to provide
  access to results and to interrogate background processes etc.
  """
  
  def __init__(self, engine, waitJob, sysPath, jobSpecs, diffNpyArgs,
               sameNpyArgs, combineFunc, callbackArg):
    
    self._tStart = time()
    self._tEnd = None
    self._output = NULL
    
    self.engine = engine
    self.waitJob = waitJob
    
    self.sysPath = sysPath
    self.jobSpecs = jobSpecs
    self.diffNpyArgs = diffNpyArgs
    self.sameNpyArgs = sameNpyArgs
    self.combineFunc = combineFunc
    self.callbackArg = callbackArg

    self.error = None
  
  
  def start(self):
    
    comm = self.engine.comm
     
    comm.send(2, dest=0, tag=0)             # Go
    comm.send(self.sysPath, dest=0, tag=1) # So modules can be imported as with Engine
    comm.send(self.jobSpecs, dest=0, tag=1) # Non-numpy args etc
    diffNpyKeys, sameNpyKeys = self.jobSpecs[-2:]
    
    for key in diffNpyKeys:
      data = self.diffNpyArgs[key]
      comm.send((list(data.shape), '%s' % data.dtype), dest=0, tag=1) # Size
      comm.Send(data, dest=0, tag=2)   # Data buffer

    for key in sameNpyKeys:
      data = self.sameNpyArgs[key]
      comm.send((list(data.shape), '%s' % data.dtype), dest=0, tag=1) # Size
      comm.Send(data, dest=0, tag=2)   # Data buffer
    
    if self.waitJob:
      localWorker = Worker(comm, self.callbackArg)
      localWorker.idle() # Will be released by Controller when done
      self.wait()
    
    
  def wait(self):
    
    if self._output is not NULL:
      return
   
    engine = self.engine
    while not engine.comm.Iprobe(source=0, tag=1):
      sleep(0.02)
    
    specs = engine.comm.recv(source=0, tag=1)
    
    if not specs: # Worker failed. Null data.
      self._output = tuple()
      print("[Parallel Engine Failed]")
      #sys.exit(0)
      return
    
    size, dtype = specs
 
    if dtype == PICKLE:
      result = engine.comm.recv(source=0, tag=2) # timeout?
 
    else:
      result = empty(size, dtype=dtype)
      engine.comm.Recv(result, source=0, tag=2) # timeout?
 
    self._output = self.combineFunc(result) 
     
  
  def waitResult(self):
  
    if self._output is NULL:
      self.wait()
      self._stopped()
    
    return self._output
  
  
  def getResult(self):
  
    if self.isComplete():
      if self._output is NULL:
        self.wait()
        self._stopped()
    
      return self._output
  
  
  def isComplete(self):

    if self._output is not NULL:
      return True
  
    self.engine.comm.send(1, dest=0, tag=0) 
    resp = self.engine.comm.recv(source=0, tag=0)
    
    if resp == 1:
      return True
    
    elif resp == 2:
      return False
      
  
  def stop(self):
    
    self.engine.comm.send(1, dest=0, tag=0)
    resp = self.engine.comm.recv(source=0, tag=0) 
    self._stopped()
    
    # Timeout
  
  
  def timeElapsed(self):
  
    if self._tEnd is None:
      return time() - self._tStart
    
    else:
      return self._tEnd - self._tStart
  
  
  def hasFailed(self):
  
    bool(self.error)
    

  def _stopped(self):
    
    self._tEnd = time()
    
    
class Engine():
  """
  This class starts a Controller, which in turn starts persistent workers
  The Engine.run() method is called for a given Job
  An Engine need not wait, but a Controller always does
  The Controller does the work of controlling the Worker pool
  and collating output data.
  """
  
  def __init__(self, numWorkers=None):
    from mpi4py import MPI
    
    if not numWorkers:
      import multiprocessing
      numWorkers = multiprocessing.cpu_count()

    self.numWorkers = numWorkers
    
    cmdArgs = [controllerPath, str(numWorkers)]
    self.comm = MPI.COMM_SELF.Spawn(sys.executable, args=cmdArgs, maxprocs=1)
    self.comm.send(1, dest=0, tag=0)         # Ping
    status = self.comm.recv(source=0, tag=0) # Wait until controller is idle
  
  
  def run(self, func, sameArgs=None, diffArgs=None, numProcs=None,
          combineFunc=None, wait=True, callbackArg=None):
    
    moduleName = func.__module__
    funcName = func.__name__
    
    if moduleName == "__main__":
      moduleObj = __import__(moduleName, globals(), locals(), [moduleName], -1)
      moduleName = os.path.splitext(moduleObj.__file__)[0]
    
    nDiffArgs = None   
    diffObjArgs = {}
    diffNpyArgs = {}
    sameObjArgs = {}
    sameNpyArgs = {}
    
    if diffArgs:
      if not isinstance(diffArgs, dict):
        msg = 'Process-specific arguments "diffArgs" must be a dictionary'
        raise Exception(msg)
      
      for name, val in diffArgs.iteritems():
        if not isinstance(val, (Sequence, ndarray)):
          msg = 'Process-specific argument "%s" must be a collection; a value for each process'
          raise Exception(msg % (name,))
        
        if nDiffArgs is None:
          nDiffArgs = len(val)
          
        elif len(val) != nDiffArgs:
          msg = 'Process-specific argument "%s" length (%d) doesn\'t match length (%d)'
          msg += ' of other argument(s)'
          raise Exception(msg % (name, len(val), nDiffArgs))
          
        if isinstance(val, ndarray):
          diffNpyArgs[name] = val
        else:
          diffObjArgs[name] = val
        
    if sameArgs:
      if not isinstance(sameArgs, dict):
        msg = 'Common arguments "sameArgs" must be a dictionary'
        raise Exception(msg)
      
      for name, val in sameArgs.iteritems():
        if name in diffArgs:
          msg = 'Argument "%s" appears twice; in process-specific and common arguments'
          raise Exception(msg % name)
      
        if isinstance(val, ndarray):
          sameNpyArgs[name] = val
        else:
          sameObjArgs[name] = val
    
    if combineFunc is None:
      combineFunc = tuple
    
    if numProcs is None:
      numProcs = nDiffArgs or self.numWorkers
      
    # TBD: Check Controller process still alive
    
    diffNpyKeys = sorted(diffNpyArgs.keys())
    sameNpyKeys = sorted(sameNpyArgs.keys())
    sysPath = list(sys.path)
    if not sysPath[0]:
      sysPath[0] = os.path.abspath('.')
    
    jobSpecs = [moduleName, funcName, numProcs, wait,
                diffObjArgs, sameObjArgs, diffNpyKeys, sameNpyKeys]
                
    job = Job(self, wait, sysPath, jobSpecs, diffNpyArgs,
              sameNpyArgs, combineFunc, callbackArg)
                
    try:
      job.start()     
       
    except Exception as err:
      job.stop()     
      tb = sys.exc_info()[2]
      traceback.print_tb(tb)
      del tb
      self.stop()
      raise(err)
    
    return job
 
 
  def stop(self):
    
    self.comm.isend(0, dest=0, tag=0) # Stop
    self.comm.Disconnect()
    
  
# Destructor?

def testFunc(n, m):
  
  s = 0.0 
  for i in xrange(n):
    s += (i/float(n))**2+m
    
  return n, s, s/2


def testFunc2(n, m):
  
  return n**m
 

if __name__ == '__main__':  
  
  vals = [5000000,6000000,7000000,8000000]
  
  engine = Engine(4)
  
  #"""
  job = engine.run(testFunc, sameArgs={'m':2}, diffArgs={'n':vals}, wait=False)

  print(job)
  print(job.timeElapsed())
  print(job.isComplete())
  print(job.getResult())

  job.wait()

  print(job.getResult())
  print(job.isComplete())

  print("Time: %.3f" % job.timeElapsed())
  
  t1 = time()
  for n in vals:
    print(testFunc(n, 2))
  t2 = time()
  
  print("Ref time: %.3f" % (t2-t1))
  #"""
  
  job = engine.run(testFunc, diffArgs={'n':[5,6,7,8,9,10,11], 'm':[2,3,4,5,6,7,8]})
  
  for i, r in enumerate(job.getResult()):
    print(i+1, r)
  
  engine.stop()
