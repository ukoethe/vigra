from vigranumpycore import *
from arraytypes import *
import convolution
import tensor
import impex
import filters
import noise
import segmentation
import morphology
import edgedetection
import classification

try:
    from vigranumpyfourier import *
except:
    print "vigranumpyfourier not found"

# auto-generate code for  additional Kernel generators:
def genKernelFactories():
   for k in dir(convolution.Kernel1D):
      if not k.startswith('init'): continue
      kn=k[4].lower()+k[5:] #remove init from beginning and start with lower case character
      code = '''def %sKernel(*args):
      k = convolution.Kernel1D()
      k.%s(*args)
      return k
convolution.%sKernel=%sKernel
convolution.%sKernel.__doc__ = convolution.Kernel1D.%s.__doc__
''' % (k,k,kn,k,kn,k)
      exec code

genKernelFactories()
selfdict = globals()
def searchfor(searchstring):
   for attr in selfdict.keys():
      contents = dir(selfdict[attr])
      for cont in contents:
         if ( cont.upper().find(searchstring.upper()) ) >= 0:
            print attr+"."+cont
