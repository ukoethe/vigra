import vigranumpycore
import arraytypes
import impex
import filters
import sampling
import analysis
import learning
import noise

try:
    import fourier
except:
    print "WARNING: Unable to load module 'vigra.fourier'"


# import most frequently used functions    
from vigranumpycore import registerPythonArrayType, listExportedArrayKeys
from arraytypes import *
from impex import readImage, readImageFromHDF5, readVolume, readVolumeFromHDF5
from filters import convolve, gaussianSmoothing
from sampling import resize

# import enums
CLOCKWISE = sampling.RotationDirection.CLOCKWISE
COUNTER_CLOCKWISE = sampling.RotationDirection.COUNTER_CLOCKWISE
UPSIDE_DOWN = sampling.RotationDirection.UPSIDE_DOWN
CompleteGrow = analysis.SRGType.CompleteGrow
KeepContours = analysis.SRGType.KeepContours
StopAtThreshold = analysis.SRGType.StopAtThreshold
 
_selfdict = globals()
def searchfor(searchstring):
   '''Scan all vigra modules to find classes and functions containing 
      'searchstring' in their name.
   '''
   for attr in _selfdict.keys():
      contents = dir(_selfdict[attr])
      for cont in contents:
         if ( cont.upper().find(searchstring.upper()) ) >= 0:
            print attr+"."+cont

def imshow(image):
    '''Display a scalar or RGB image by means of matplotlib.
       If the image does not have one or three channels, an exception is raised.
       The image will be automatically scaled to the range 0...255.
    '''
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

        
# auto-generate code for  additional Kernel generators:
def _genKernelFactories(name):
   for oldName in dir(eval('filters.'+name)):
      if not oldName.startswith('init'): 
        continue
      #remove init from beginning and start with lower case character
      newName = oldName[4].lower() + oldName[5:] + 'Kernel'
      if name == 'Kernel2D':
        newName += '2D'
      code = '''def %(newName)s(*args):
      k = filters.%(name)s()
      k.%(oldName)s(*args)
      return k
%(newName)s.__doc__ = filters.%(name)s.%(oldName)s.__doc__
filters.%(newName)s=%(newName)s
''' % {'oldName': oldName, 'newName': newName, 'name': name}
      exec code

_genKernelFactories('Kernel1D')
_genKernelFactories('Kernel2D')
del _genKernelFactories

# define watershedsUnionFind()
def _genWatershedsUnionFind():
    def watershedsUnionFind(image, neighborhood=None, out = None):
        '''Compute watersheds of an image using the union find algorithm.
           If 'neighborhood' is 'None', it defaults to 8-neighborhood for 2D inputs
           and 6-neighborhood for 3D inputs.
           
           Calls :func:`watersheds` with parameters::\n\n
                watersheds(image, neighborhood=neighborhood, method='UnionFind', out=out)
        '''
        if neighborhood is None:
            neighborhood = 8 if image.spatialDimensions == 2 else 6
                
        return analysis.watersheds(image, neighborhood=neighborhood, method='UnionFind', out=out)
    
    watershedsUnionFind.__module__ = 'vigra.analysis'
    analysis.watershedsUnionFind = watershedsUnionFind

_genWatershedsUnionFind()
del _genWatershedsUnionFind

# define tensor convenience functions
def _genTensorConvenienceFunctions():
    def hessianOfGaussianEigenvalues(image, scale, out = None):
        '''Compute the eigenvalues of the Hessian of Gaussian at the given scale
           for a scalar image or volume.
           
           Calls :func:`hessianOfGaussian` and :func:`tensorEigenvalues`.
        '''

        return filters.tensorEigenvalues(filters.hessianOfGaussian(image, scale), out=out)
    
    hessianOfGaussianEigenvalues.__module__ = 'vigra.filters'
    filters.hessianOfGaussianEigenvalues = hessianOfGaussianEigenvalues

    def structureTensorEigenvalues(image, innerScale, outerScale, out = None):
        '''Compute the eigenvalues of the structure tensor at the given scales
           for a scalar or multi-channel image or volume.
           
           Calls :func:`structureTensor` and :func:`tensorEigenvalues`.
        '''

        return filters.tensorEigenvalues(filters.structureTensor(image, innerScale, outerScale), out=out)
    
    structureTensorEigenvalues.__module__ = 'vigra.filters'
    filters.structureTensorEigenvalues = structureTensorEigenvalues

_genTensorConvenienceFunctions()
del _genTensorConvenienceFunctions
