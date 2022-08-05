from PIL import Image
from numpy import array, uint8
import subprocess

def pixmapToImage(pixmap, mode='RGB'):
  
  if pixmap.max() > 255:
    pixmap *= 255.0 / pixmap.max()

  pixmap = array(pixmap, uint8)
  img = Image.fromarray(pixmap, mode)
  
  return img

def numpyFramesToMovie(frameGenerator, outFile, fps=30, ffopts=[], codec='rawvideo'):
  
  proc = None
  for frameArray in frameGenerator:
    if proc is None:
    
      cmd =['avconv', '-y', # No propt for overwrite
            '-s', '%dx%d' % (fr.shape[1], fr.shape[0]),
            '-r', str(fps),
            '-an', # Disable audio, -stats 
            '-codec:v', codec,
            '-f', 'rawvideo', # force i/o format (container not codec)
            '-pix_fmt', 'rgb24',
            '-i', '-'] + ffopts + [outFile]
            
      proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
      
    proc.stdin.write(frameArray.tostring())
    
  proc.stdin.close()
  proc.wait()
