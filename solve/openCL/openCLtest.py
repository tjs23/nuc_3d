import pyopencl as cl

for platform in cl.get_platforms():
  for paramName in dir(cl.platform_info):
    val = getattr(cl.platform_info, paramName)
    if type(val) == int:
      print('PLATFORM:', paramName, platform.get_info(val))
  
  for device in platform.get_devices():
    for paramName in dir(cl.device_info):
      val = getattr(cl.device_info, paramName)
      if type(val) == int:
        try:
          print('  DEVICE:', paramName, device.get_info(val))
        
        except cl.LogicError:
          pass
      
mf = cl.mem_flags

from numpy import random, float32, empty, array

kernelFile = 'nbody.cl'
nSteps = 100

pos = 100.0 * random.rand(1000, 4).astype(float32)
pos[:,3] = 1.0
result = empty(pos.shape, float32)

context = cl.create_some_context()
devices = context.get_info(cl.context_info.DEVICES)

print() 
print() 
print("Context:", context)
print("Using devices:", devices)

localBufSize =  devices[0].get_info(cl.device_info.LOCAL_MEM_SIZE)
queue =  cl.CommandQueue(context)

posBufA   =  cl.Buffer(context, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=pos)
posBufB   =  cl.Buffer(context, mf.READ_WRITE, pos.nbytes)
velBuf    =  cl.Buffer(context, mf.READ_WRITE, pos.nbytes)
tempBuf   =  cl.LocalMemory(localBufSize)

kernelSrc = open(kernelFile).read()

prog =  cl.Program(context, kernelSrc).build()

sizeGlob = pos.size()
sizeLocl = None # Host specific
timeStep = float32(0.01)
minDist = float32(0.01)

for i in range(nSteps):

  prog.nbody(queue, sizeGlob, sizeLocl,  timeStep, minDist,
             posBufA, posBufB, velBuf, tempBuf)

  prog.nbody(queue, sizeGlob, sizeLocl, timeStep, minDist,
             posBufB, posBufA, velBuf, tempBuf)
  
  cl.enqueue_barrier(queue) # See if this could be done less often

  #print i
  
cl.enqueue_copy(queue, result, posBufA).wait()

print(i, (pos-result).mean(axis=0))
