import numpy, sys, os, traceback
#os.environ['OMPI_MCA_mpi_yield_when_idle'] = '1'

from numpy import array, empty, ndarray, issubdtype

from time import time, sleep

PICKLE = 'pickle'

dirName = os.path.dirname(__file__)
workerPath = os.path.join(dirName, 'Worker.py')
  
class Controller():
  """
  Handles comminication with Worker processes, issuing data from
  a pool and collecting the results.
  The Controller is not a daemon; it only lasts as long
  as the invoking process.
  """
  
  def __init__(self, nWorkers=None, queue=False):
    
    # Later accept details for remote workers
    
    if nWorkers is None:
      nWorkers = MPI.COMM_WORLD.Get_size()
    
    self.comm = None
    self.nWorkers = nWorkers
    self.workers = []
    self.queue = []
    self.pComm = MPI.Comm.Get_parent()
    self.startWorkers()
    self.checkWorkers()
    self._parentWait = False


  def startWorkers(self):
    
    cmd = sys.executable
    self.comm = MPI.COMM_SELF.Spawn(cmd, args=[workerPath],
                                    maxprocs=self.nWorkers)
    self.workers = list(range(self.nWorkers))

  
  def stopWorkers(self):
    
    for rank in self.workers:
      self.comm.isend(0, dest=rank, tag=0) # Stop, nonblocking
 
    self.comm = None
    self.workers = []


  def _getArraySpecs(self, obj):
 
    if isinstance(obj, ndarray):
      isArray = True
      shape = obj.shape
      dtype = str(obj.dtype)
    
    else:
      try:
        arry = array(obj)
        shape = obj.shape
         
        if issubdtype(dtype, int):
          isArray = True
          dtype = str(dtype)
          obj = arry
        
        elif issubdtype(dtype, float):
          isArray = True
          dtype = str(dtype)
          obj = arry
        
        else: # Do not use numpy arrays for strings, objects etc.
          isArray = False
          dtype = PICKLE
      
      except Exception:
        isArray = False
        shape = 1
        dtype = PICKLE
    
    return obj, isArray, shape, dtype


  def restartWorkers(self):
  
    # e.g. if failed    
    self.stopWorkers()
    self.startWorkers()
    
    
  def loadFunc(self, moduleName, funcName):
  
    try:
      module = __import__(moduleName, globals(), locals(), [moduleName], -1)
    
    except Exception as err:
      msg = 'Controller : Could not import module "%s"' % moduleName
      print(err)
      raise Exception(msg)
 
    try:
      func = getattr(module, funcName)
    
    except Exception as err:
      msg = 'Controller : Could not access function "%s" in module "%s"' % (funcName, moduleName)
      print(err)
      raise Exception(msg)

    return func
  
  
  def checkWorkers(self):
    
    
    # Make async
    for rank in self.workers:
      self.comm.send(1, dest=rank, tag=0)
      status = self.comm.recv(source=rank, tag=0) # Timeout...
      if not status:
        return False
      
    return True
  
  
  def _checkEngineSignal(self):
    
    comm = self.pComm
    
    if comm.Iprobe(source=0, tag=0):  # A control signal
      query = comm.recv(source=0, tag=0)
  
      if query == 0: # Quit
        self.quit()
 
      elif query == 1: # Ping
        self.pComm.send(2, dest=0, tag=0) # Busy
 
      elif query == 2: # Other
        pass # 
 
      elif query == 3: # Cancel job
        self.stopWorkers()
        self.pComm.send(1, dest=0, tag=0) # Done


  def _getWorkerData(self, comm, rank):
    
    wData = comm.recv(source=rank, tag=1)
 
    if wData == 0:
      self.quit() # Something is wrong. Null data. Workers are failing.
      
    elif wData == 1:
      self.quit() # Something is wrong. Null data. Workers are failing.
      
    else:
      objs, npyIdx, sizes, dtypes = wData
 
    # Fetch arrays
    for i, size in enumerate(sizes):
      data = empty(size,  dtypes[i])
      comm.Recv(data, source=rank, tag=2) # Array data
      objs.insert(i, data)
 
    if len(objs) == 1:
      objs = objs[0]
    else:
      objs = tuple(objs)
 
    return objs
    
 
  def run(self, mFile, fName, sysPath, nProc, waiting,
          diffObjArgs, sameObjArgs, diffNpyArgs, sameNpyArgs):
    
    # Check workers
    
    #if not self.checkWorkers():
    #  self.restartWorkers()
    
    results = [None] * nProc
    args = []
    for i in range(nProc):
      kwargs = {}
      kwargsNpy = {}
      
      for key in sameObjArgs:
        kwargs[key] = sameObjArgs[key]

      for key in sameNpyArgs:
        kwargsNpy[key] = sameNpyArgs[key]       
      
      for key in diffObjArgs:
        iMax = len(diffObjArgs[key])
        kwargs[key] = diffObjArgs[key][i % iMax]

      for key in diffNpyArgs:
        iMax = len(diffNpyArgs[key])
        kwargsNpy[key] = diffNpyArgs[key][i % iMax]
      
      keysNpy = sorted(kwargsNpy.keys())
      args.append((kwargs, kwargsNpy, keysNpy))
    
    alloc = {} # Which input a worker is allocated to
    pool = list(range(self.nWorkers))
    
    i = 0
    while args:
      kwargs, kwargsNpy, keysNpy = args.pop(0)
      jobSpecs = [mFile, fName, kwargs, keysNpy]
      rank = pool.pop(0)
      alloc[rank] = i
      i += 1
      
      if waiting and (rank == 0):
        comm = self.pComm
      else:
        comm = self.comm
      
      comm.send(2, dest=rank, tag=0) # Go. Timeout?
      comm.send(sysPath, dest=rank, tag=1)      
      comm.send(jobSpecs, dest=rank, tag=1)      
        
      # Send arrays
      for key in keysNpy:
        data = kwargsNpy[key]
        comm.send((data.shape, '%s' % data.dtype), dest=rank, tag=1) # Size
        comm.Isend(data, dest=rank, tag=2) # the data
   
      #print "Controller sent data", rank
      
      while not pool: # Nothing more to start, wait
        
        # Wait for any free worker
        for rank in alloc:
          if waiting and (rank == 0):
            comm = self.pComm
          else:
            comm = self.comm
          
          if comm.Iprobe(source=rank, tag=1):
            self._checkEngineSignal()
            objs = self._getWorkerData(comm, rank)
            
            # Collect result and make worker available again
            pool.append(rank)
            results[alloc[rank]] = objs            
        
          sleep(0.01)
        
        for rank in pool:
          del alloc[rank]
    
    # When all jobs sent args is empty
    # some workers may not be done
    while alloc:
      for rank in alloc:
        if waiting and (rank == 0):
          comm = self.pComm
        else:
          comm = self.comm
 
        if comm.Iprobe(source=rank, tag=1):
          self._checkEngineSignal()
          objs = self._getWorkerData(comm, rank)
          results[alloc[rank]] = objs
          del alloc[rank]
          break
        
        sleep(0.01)
      
    if waiting:
      # Release the root process worker
      self.pComm.send(0, dest=0, tag=0)
           
    return tuple(results)
    
    
  def idle(self):
    
    again = self.pComm.recv(source=0, tag=0)

    while again:
      if again == 1:
        self.pComm.send(1, dest=0, tag=0) # Ping, not busy
        
      else:
        sysPath = self.pComm.recv(source=0, tag=1)
        sys.path = sysPath # Next args may include objects for which the module must be known
        jobSpecs = self.pComm.recv(source=0, tag=1)
        mFile, fName, nProc, waiting = jobSpecs[:4]
        diffObjArgs, sameObjArgs = jobSpecs[4:6]
        diffNpyKeys, sameNpyKeys = jobSpecs[6:]
       
        diffNpyArgs = {}
        sameNpyArgs = {}
        
        for key in diffNpyKeys:
          size, dtype = self.pComm.recv(source=0, tag=1)
          data = empty(size, dtype=dtype)
          self.pComm.Recv(data, source=0, tag=2)
          diffNpyArgs[key] = data
 
        for key in sameNpyKeys:
          size, dtype = self.pComm.recv(source=0, tag=1)
          data = empty(size, dtype=dtype)
          self.pComm.Recv(data, source=0, tag=2)
          sameNpyArgs[key] = data
        
        self._parentWait = waiting
        
        try:
          func = self.loadFunc(mFile, fName)
        except Exception as err:
          #print("Controller exception", err)
          self.quit(err)
 
        try:
          result = self.run(mFile, fName, sysPath, nProc, waiting, diffObjArgs,
                            sameObjArgs, diffNpyArgs, sameNpyArgs)
          
        except Exception as err:
          #print("Controller exception", err)
          self.quit(err)
 
        result, isArray, shape, dtype = self._getArraySpecs(result)
        self.pComm.send([shape, dtype], dest=0, tag=1)
        
        if isArray:
          self.pComm.Send(result, dest=0, tag=2)  # Data buffer
        else:
          self.pComm.send(result, dest=0, tag=2)  # Pickle
       
        self._parentWait = False
     
      while not self.pComm.Iprobe(source=0, tag=0):
        sleep(0.01)
                  
      again = self.pComm.recv(source=0, tag=0) # Timeout? Check workers?
         
    #print("[Controller End]")
    self.quit()

  
  def quit(self, exception=None):
    
    self.stopWorkers()
    
    #print("[Controller Stop]")    

    if self._parentWait:
      self.pComm.send(0, dest=0, tag=0)  # Release any local worker, bocking before disconnect
    
    self.pComm.send(0, dest=0, tag=1)  # Tell Engine Controller is stopping

    if exception:
      tb = sys.exc_info()[2]
      traceback.print_tb(tb)
      del tb
      raise(exception)

    self.pComm.Disconnect()
    
    if self.comm:
      self.comm.Disconnect()
    
    sys.exit(0)  


if __name__ == '__main__':

  from mpi4py import MPI
  
  script, nWorkers = sys.argv

  controller = Controller(int(nWorkers))
  controller.idle()
  controller.quit()
  
