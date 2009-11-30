from vigranumpycmodule import *
from arraytypes import *

try:
    from vigranumpyfourier import *
except:
    print "vigranumpyfourier not found"

# auto-generate code for  additional Kernel generators:
def genKernelFactories():
   for k in dir(Kernel1D):
      if not k.startswith('init'): continue;
      code = '''@staticmethod
def %sKernel(*args):
      k = Kernel1D()
      k.%s(*args)
      return k
Kernel1D.%sKernel=%sKernel
Kernel1D.%sKernel.__doc__ = Kernel1D.%s.__doc__
''' % (k,k,k,k,k,k)
      exec code

genKernelFactories()
