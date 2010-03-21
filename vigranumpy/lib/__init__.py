import vigranumpycore
import arraytypes
import impex
import convolution
import tensor
import filters
import noise
import segmentation
import morphology
import edgedetection
import classification

try:
    import fourier
except:
    print "WARNING: Unable to load module 'vigra.fourier'"

from vigranumpycore import *
from arraytypes import *
from impex import *

# auto-generate code for  additional Kernel generators:
def genKernelFactories():
   for oldName in dir(convolution.Kernel1D):
      if not oldName.startswith('init'): 
        continue
      #remove init from beginning and start with lower case character
      newName = oldName[4].lower() + oldName[5:] + 'Kernel'
      code = '''def %(newName)s(*args):
      k = convolution.Kernel1D()
      k.%(oldName)s(*args)
      return k
%(newName)s.__doc__ = convolution.Kernel1D.%(oldName)s.__doc__
convolution.%(newName)s=%(newName)s
''' % {'oldName': oldName, 'newName': newName}
      exec code
      
   for oldName in dir(convolution.Kernel2D):
      if not oldName.startswith('init'): 
        continue
      newName = oldName[4].lower() + oldName[5:] + 'Kernel2D'
      code = '''def %(newName)s(*args):
      k = convolution.Kernel2D()
      k.%(oldName)s(*args)
      return k
%(newName)s.__doc__ = convolution.Kernel2D.%(oldName)s.__doc__
convolution.%(newName)s=%(newName)s
''' % {'oldName': oldName, 'newName': newName}
      exec code

genKernelFactories()
selfdict = globals()
def searchfor(searchstring):
   for attr in selfdict.keys():
      contents = dir(selfdict[attr])
      for cont in contents:
         if ( cont.upper().find(searchstring.upper()) ) >= 0:
            print attr+"."+cont

def imshow(image):
    import matplotlib
    
    if image.ndim == 3:
        if image.shape[2] != 3:
            raise RuntimeError("vigra.imshow(): Multi channel image must have 3 channels.")
        if image.dtype != uint8:
            mi, ma = image.min(), image.max()
            if mi >= ma:
                image = image.__class__(image.shape, dtype=uint8)
            else:
                image = (image - mi) / (ma - mi)
        return matplotlib.pyplot.imshow(image.swapaxes(0,1))
    elif image.ndim == 2:
        return matplotlib.pyplot.imshow(image.swapaxes(0,1), cmap=matplotlib.cm.gray, \
                                     norm=matplotlib.cm.colors.Normalize())
    else:
        raise RuntimeError("vigra.imshow(): ndim must be 2 or 3.")

