import sys, os, traceback
#os.environ['OMPI_MCA_mpi_yield_when_idle'] = '1'
from numpy import ndarray, empty
from time import sleep

QUIT = 0
PING = 1
RUN  = 2

class Worker():

  def __init__(self, comm=None, callbackArg=None):
    
    if comm:
      self.comm = comm
    else:
      self.comm = MPI.Comm.Get_parent()
    
    self.callbackArg = callbackArg
    self.rank = self.comm.Get_rank()
    

  def loadFunc(self, moduleName, funcName):
  
    try:
      module = __import__(moduleName, globals(), locals(), [moduleName], -1)
    
    except Exception as err:
      msg = 'Worker : Could not import module "%s"' % moduleName
      print(msg)
      self.quit(err)
 
    try:
      func = getattr(module, funcName)
    
    except Exception as err:
      msg = 'Worker : Could not access function "%s" in module "%s"' % (funcName, moduleName)
      print(msg)
      self.quit(err)

    return func
    
  
  def quit(self, exception=None):    
    
    print("[Worker %d Quit]" % self.rank)
    self.comm.send(0, dest=0, tag=1) # Null return data
    
    if exception is not None:
      tb = sys.exc_info()[2]
      traceback.print_tb(tb)
      del tb
      raise(exception)
    
    sys.exit(0)
  
  
  def runFunc(self):
  
    sysPath = self.comm.recv(source=0, tag=1)
    sys.path = sysPath # In case objs in the following args needs modeuls to be known
    mFile, fName, kwargs, keysNpy = self.comm.recv(source=0, tag=1)
    #print("Worker %d : run" % self.rank)
    
    try:
      func = self.loadFunc(mFile, fName)
    
    except Exception as err:
      self.quit(err)  
    
    for key in keysNpy:
      size, dtype = self.comm.recv(source=0, tag=1)
      data = empty(size, dtype=dtype)
      self.comm.Recv(data, source=0, tag=2)
      kwargs[key] = data
    
    if self.callbackArg:
      argName, callback = self.callbackArg
      kwargs[argName] = callback
      
    try:
      #print("Worker %d : wait" % self.rank)
      result = func(**kwargs)
      #print("Worker %d : done" % self.rank)
    
    except Exception as err:
      self.quit(err)
      return    
    
    if not isinstance(result, tuple):
      result = (result,)
    
    objs = []
    arrayIndices = []
    arrayShapes = []
    arrayDtypes = []
    arrays = []
    
    for i, val in enumerate(result):
      if isinstance(val, ndarray):
        arrayIndices.append(i)
        arrayShapes.append(val.shape)
        arrayDtypes.append('%s' % val.dtype)
        arrays.append(val)
      else:
        objs.append(val)
    
    self.comm.send((objs, arrayIndices, arrayShapes, arrayDtypes), dest=0, tag=1)
    
    for data in arrays:
      self.comm.Send(data, dest=0, tag=2)
    
    #print("Worker %d : sent" % self.rank)
  
  
  def idle(self):
    
    #print("Worker %d : started" % self.rank, self.comm)
    code = self.comm.recv(source=0, tag=0)  # Instructions
    
    while code != QUIT:
        
      if code == PING:
        #print("Worker %d : ping" % self.rank, self.comm)
        self.comm.send(1, dest=0, tag=0)
        #print("Worker %d : pong" % self.rank, self.comm)
        
      elif code == RUN:
        self.runFunc()
      
      #print("Worker %d : idle" % self.rank)
      while not self.comm.Iprobe(source=0, tag=0):
        sleep(0.01)
     
      code = self.comm.recv(source=0, tag=0) # Next instructions
      #print("Worker %d : code %d" % (self.rank, code))
    
    #print("[Worker %d End]" % self.rank)
    #self.comm.send(0, dest=0, tag=0)


if __name__ == '__main__':
  
  from mpi4py import MPI
  
  w = Worker()
  w.idle()

