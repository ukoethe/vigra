#######################################################################
#
#         Copyright 2009-2010 by Ullrich Koethe
#
#    This file is part of the VIGRA computer vision library.
#    The VIGRA Website is
#        http://hci.iwr.uni-heidelberg.de/vigra/
#    Please direct questions, bug reports, and contributions to
#        ullrich.koethe@iwr.uni-heidelberg.de    or
#        vigra@informatik.uni-hamburg.de
#
#    Permission is hereby granted, free of charge, to any person
#    obtaining a copy of this software and associated documentation
#    files (the "Software"), to deal in the Software without
#    restriction, including without limitation the rights to use,
#    copy, modify, merge, publish, distribute, sublicense, and/or
#    sell copies of the Software, and to permit persons to whom the
#    Software is furnished to do so, subject to the following
#    conditions:
#
#    The above copyright notice and this permission notice shall be
#    included in all copies or substantial portions of the
#    Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND
#    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
#    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
#    OTHER DEALINGS IN THE SOFTWARE.
#
#######################################################################

import sys, os

_vigra_path = os.path.abspath(os.path.dirname(__file__))
_vigra_doc_path = _vigra_path + '/doc/vigranumpy/index.html'

if sys.platform.startswith('win'):
    # On Windows, add subdirectory 'dlls' to the PATH in order to find
    # the DLLs vigranumpy depends upon. Since this directory appears
    # at the end of PATH, already installed DLLs are always preferred.
    _vigra_dll_path = _vigra_path + '/dlls'
    if os.path.exists(_vigra_dll_path):
        os.putenv('PATH', os.getenv('PATH') + os.pathsep + _vigra_dll_path)

def _fallbackModule(moduleName, message):
    '''This function installs a fallback module with the given 'moduleName'.
       All function calls into this module raise an ImportError with the
       given 'message' that hopefully tells the user why the real module
       was not available.
    '''
    import sys
    moduleClass = vigranumpycore.__class__
    class FallbackModule(moduleClass):
        def __init__(self, name):
            moduleClass.__init__(self, name)
            self.__name__ = name
        def __getattr__(self, name):
            if name.startswith('__'):
                return moduleClass.__getattribute__(self, name)
            try:
                return moduleClass.__getattribute__(self, name)
            except AttributeError:
                raise ImportError("""%s.%s: %s""" % (self.__name__, name, self.__doc__))

    module = FallbackModule(moduleName)
    sys.modules[moduleName] = module
    module.__doc__ = """Module '%s' is not available.\n%s""" % (moduleName, message)

if not os.path.exists(_vigra_doc_path):
    _vigra_doc_path = "http://hci.iwr.uni-heidelberg.de/vigra/doc/vigranumpy/index.html"

__doc__ = '''VIGRA Computer Vision Library

HTML documentation is available in

   %s

Help on individual functions can be obtained via their doc strings
as usual.

The following sub-modules group related functionality:

* impex
* colors
* filters
* sampling
* fourier
* analysis
* learning
* noise
''' % _vigra_doc_path
 
from __version__ import version
import vigranumpycore
import arraytypes
import impex
import filters
import sampling
import analysis
import learning
import colors
import noise

try:
    import fourier
except:
    print "WARNING: Unable to load module 'vigra.fourier'"
    _fallbackModule('vigra.fourier', "   Probably, the fftw3 libraries could not be found during compilation or import.")
    import fourier

# import most frequently used functions
from vigranumpycore import registerPythonArrayType, listExportedArrayKeys
from arraytypes import *
from filters import convolve, gaussianSmoothing
from sampling import resize
from impex import readImage, readVolume
try:
    from impex import readImageFromHDF5, readVolumeFromHDF5
except:
    pass

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
    import matplotlib.pylab
    
    if image.ndim == 3:
        if image.shape[2] != 3:
            raise RuntimeError("vigra.imshow(): Multi channel image must have 3 channels.")
        if image.dtype != uint8:
            image = colors.linearRangeMapping(image, newRange=(0.0, 255.0),\
                                              out=image.__class__(image.shape, dtype=uint8))
        return matplotlib.pyplot.imshow(image.swapaxes(0,1).view(numpy.ndarray))
    elif image.ndim == 2:
        return matplotlib.pyplot.imshow(image.swapaxes(0,1).view(numpy.ndarray), cmap=matplotlib.cm.gray, \
                                     norm=matplotlib.cm.colors.Normalize())
    else:
        raise RuntimeError("vigra.imshow(): ndim must be 2 or 3.")

        
# auto-generate code for additional Kernel generators:
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
