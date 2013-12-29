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
    module.__doc__ = """Import of module '%s' failed.\n%s""" % (moduleName, message)

if not os.path.exists(_vigra_doc_path):
    _vigra_doc_path = "http://hci.iwr.uni-heidelberg.de/vigra/doc/vigranumpy/index.html"

__doc__ = '''VIGRA Computer Vision Library

HTML documentation is available in

   %s

Help on individual functions can be obtained via their doc strings
as usual.

The following sub-modules group related functionality:

* arraytypes (VigraArray and axistags, automatically imported into 'vigra')
* ufunc      (improved array arithmetic, automatically used by VigraArray)
* impex      (image and array I/O)
* colors     (color space transformations)
* filters    (spatial filtering, e.g. smoothing)
* sampling   (image and array re-sampling and interpolation)
* fourier    (Fourier transform and Fourier domain filters)
* analysis   (image analysis and segmentation)
* learning   (machine learning and classification)
* noise      (noise estimation and normalization)
* geometry   (geometric algorithms, e.g. convex hull)
* clustering   (hierarchical clustering algorithms, e.g. convex UCM,Ward)
''' % _vigra_doc_path
 
from __version__ import version
import vigranumpycore
import arraytypes
import impex
import sampling
import filters
import analysis
import learning
import colors
import noise
import geometry
import optimization
import clustering
import cgp

sampling.ImagePyramid = arraytypes.ImagePyramid

try:
    import fourier
except Exception, e:
    _fallbackModule('vigra.fourier', 
    '''
    %s
    
    Make sure that the fftw3 libraries are found during compilation and import.
    They may be downloaded at http://www.fftw.org/.''' % str(e))
    import fourier

# import most frequently used functions
from arraytypes import *
standardArrayType = arraytypes.VigraArray 
defaultAxistags = arraytypes.VigraArray.defaultAxistags

from impex import readImage, readVolume

def readHDF5(filenameOrGroup, pathInFile, order=None):
    '''Read an array from an HDF5 file.
    
       'filenameOrGroup' can contain a filename or a group object
       referring to an already open HDF5 file. 'pathInFile' is the name 
       of the dataset to be read, including intermediate groups. If the 
       first argument is a group object, the path is relative to this 
       group, otherwise it is relative to the file's root group.
       
       If the dataset has an attribute 'axistags', the returned array
       will have type :class:`~vigra.VigraArray` and will be transposed 
       into the given 'order' ('vigra.VigraArray.defaultOrder'
       will be used if no order is given).  Otherwise, the returned 
       array is a plain 'numpy.ndarray'. In this case, order='F' will 
       return the array transposed into Fortran order.
       
       Requirements: the 'h5py' module must be installed.
    '''
    import h5py
    if isinstance(filenameOrGroup, h5py.highlevel.Group):
        file = None
        group = filenameOrGroup
    else:
        file = h5py.File(filenameOrGroup, 'r')
        group = file['/']
    try:
        dataset = group[pathInFile]
        if not isinstance(dataset, h5py.highlevel.Dataset):
            raise IOError("readHDF5(): '%s' is not a dataset" % pathInFile)
        data = dataset.value
        axistags = dataset.attrs.get('axistags', None)
        if axistags is not None:
            data = data.view(arraytypes.VigraArray)
            data.axistags = arraytypes.AxisTags.fromJSON(axistags)
            if order is None:
                order = arraytypes.VigraArray.defaultOrder
            data = data.transposeToOrder(order)
        else:
            if order == 'F':
                data = data.transpose()
            elif order not in [None, 'C', 'A']:
                raise IOError("readHDF5(): unsupported order '%s'" % order)
    finally:
        if file is not None:
            file.close()
    return data
        
def writeHDF5(data, filenameOrGroup, pathInFile, compression=None):
    '''Write an array to an HDF5 file.

       'filenameOrGroup' can contain a filename or a group object
       referring to an already open HDF5 file. 'pathInFile' is the name of the
       dataset to be written, including intermediate groups. If the first
       argument is a group object, the path is relative to this group,
       otherwise it is relative to the file's root group. If the dataset already
       exists, it will be replaced without warning.

       If 'data' has an attribute 'axistags', the array is transposed to
       numpy order before writing. Moreover, the axistags will be
       stored along with the data in an attribute 'axistags'.

       'compression' can be set to 'gzip', 'szip' or 'lzf'
       gzip (standard compression),
       szip (available if HDF5 is compiled with szip. Faster compression, limited types),
       lzf (very fast compression, all types).
       The 'lzf' compression filter is many times faster than 'gzip' 
       at the cost of a lower compresion ratio.

       Requirements: the 'h5py' module must be installed.
    '''
    import h5py
    if isinstance(filenameOrGroup, h5py.highlevel.Group):
        file = None
        group = filenameOrGroup
    else:
        file = h5py.File(filenameOrGroup)
        group = file['/']
    try:
        levels = pathInFile.split('/')
        for groupname in levels[:-1]:
            if groupname == '':
                continue
            g = group.get(groupname, default=None)
            if g is None:
                group = group.create_group(groupname)
            elif not isinstance(g, h5py.highlevel.Group):
                raise IOError("writeHDF5(): invalid path '%s'" % pathInFile)
            else:
                group = g
        dataset = group.get(levels[-1], default=None)
        if dataset is not None:
            if isinstance(dataset, h5py.highlevel.Dataset):
                del group[levels[-1]]
            else:
                raise IOError("writeHDF5(): cannot replace '%s' because it is not a dataset" % pathInFile)
        try:
            data = data.transposeToNumpyOrder()
        except:
            pass
        dataset = group.create_dataset(levels[-1], data=data, compression=compression)
        if hasattr(data, 'axistags'):
            dataset.attrs['axistags'] = data.axistags.toJSON()
    finally:
        if file is not None:
            file.close()
        
impex.readHDF5 = readHDF5
readHDF5.__module__ = 'vigra.impex'
impex.writeHDF5 = writeHDF5
writeHDF5.__module__ = 'vigra.impex'

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

# FIXME: use axistags here
def imshow(image):
    '''Display a scalar or RGB image by means of matplotlib.
       If the image does not have one or three channels, an exception is raised.
       The image will be automatically scaled to the range 0...255 when its dtype 
       is not already 'uint8'.
    '''
    import matplotlib.pylab
    
    if not hasattr(image, 'axistags'):
        return matplotlib.pyplot.imshow(image)
    
    image = image.transposeToNumpyOrder()
    if image.channels == 1:
        image = image.dropChannelAxis().view(numpy.ndarray)
        plot = matplotlib.pyplot.imshow(image, cmap=matplotlib.cm.gray, \
                                         norm=matplotlib.cm.colors.Normalize())
        matplotlib.pylab.show()
        return plot
    elif image.channels == 3:
        if image.dtype != numpy.uint8:
            out = image.__class__(image.shape, dtype=numpy.uint8, axistags=image.axistags)
            image = colors.linearRangeMapping(image, newRange=(0.0, 255.0), out=out)
        plot = matplotlib.pyplot.imshow(image.view(numpy.ndarray))
        matplotlib.pylab.show()
        return plot
    else:
        raise RuntimeError("vigra.imshow(): Image must have 1 or 3 channels.")

        
# auto-generate code for additional Kernel generators:
def _genKernelFactories(name):
    for oldName in dir(eval('filters.'+name)):
        if not oldName.startswith('init'):
            continue
        #remove init from beginning and start with lower case character
        newPrefix = oldName[4].lower() + oldName[5:]
        if newPrefix == "explicitly":
            newPrefix = "explict"
        newName = newPrefix + 'Kernel'
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
    def hessianOfGaussianEigenvalues(image, scale, out=None, 
                                     sigma_d=0.0, step_size=1.0, window_size=0.0, roi=None):
        '''Compute the eigenvalues of the Hessian of Gaussian at the given scale
           for a scalar image or volume.
           
           Calls :func:`hessianOfGaussian` and :func:`tensorEigenvalues`.
        '''
        
        hessian = filters.hessianOfGaussian(image, scale, 
                                            sigma_d=sigma_d, step_size=step_size, 
                                            window_size=window_size, roi=roi)
        return filters.tensorEigenvalues(hessian, out=out)
    
    hessianOfGaussianEigenvalues.__module__ = 'vigra.filters'
    filters.hessianOfGaussianEigenvalues = hessianOfGaussianEigenvalues

    def structureTensorEigenvalues(image, innerScale, outerScale, out=None, 
                                   sigma_d=0.0, step_size=1.0, window_size=0.0, roi=None):
        '''Compute the eigenvalues of the structure tensor at the given scales
           for a scalar or multi-channel image or volume.
           
           Calls :func:`structureTensor` and :func:`tensorEigenvalues`.
        '''

        st = filters.structureTensor(image, innerScale, outerScale, 
                                     sigma_d=sigma_d, step_size=step_size, 
                                     window_size=window_size, roi=roi)
        return filters.tensorEigenvalues(st, out=out)
    
    structureTensorEigenvalues.__module__ = 'vigra.filters'
    filters.structureTensorEigenvalues = structureTensorEigenvalues

_genTensorConvenienceFunctions()
del _genTensorConvenienceFunctions

# define feature convenience functions
def _genFeaturConvenienceFunctions():
    def supportedFeatures(array):
        '''Return a list of feature names that are available for the given array. These feature
           names are the valid inputs to a call of :func:`extractFeatures`. E.g., to compute 
           just the first two features in the list, use::
           
                f = vigra.analysis.supportedFeatures(array)
                print "Computing features:", f[:2]
                r = vigra.analysis.extractFeatures(array, features=f[:2])
        '''
        
        return analysis.extractFeatures(array, None).supportedFeatures()

    supportedFeatures.__module__ = 'vigra.analysis'
    analysis.supportedFeatures = supportedFeatures

    def supportedRegionFeatures(array, labels):
        '''Return a list of feature names that are available for the given array and label array. 
           These feature names are the valid inputs to a call of 
           :func:`extractRegionFeatures`. E.g., to compute just the first two features in the 
           list, use::
           
                f = vigra.analysis.supportedRegionFeatures(array, labels)
                print "Computing features:", f[:2]
                r = vigra.analysis.extractRegionFeatures(array, labels, features=f[:2])
        '''
        
        return analysis.extractRegionFeatures(array, labels, None).supportedFeatures()
    
    supportedRegionFeatures.__module__ = 'vigra.analysis'
    analysis.supportedRegionFeatures = supportedRegionFeatures
    
    # implement the read-only part of the 'dict' API in FeatureAccumulator and RegionFeatureAccumulator
    def __len__(self):
        return len(self.keys())
    def __iter__(self):
        return self.keys().__iter__()
    def has_key(self, key):
        try:
            return self.isActive(key)
        except:
            return False
    def values(self):
        return [self[k] for k in self.keys()]
    def items(self):
        return [(k, self[k]) for k in self.keys()]
    def iterkeys(self):
        return self.keys().__iter__()
    def itervalues(self):
        for k in self.keys():
            yield self[k]
    def iteritems(self):
        for k in self.keys():
            yield (k, self[k])
    
    for k in ['__len__', '__iter__', 'has_key', 'values', 'items', 'iterkeys', 'itervalues', 'iteritems']:
        setattr(analysis.FeatureAccumulator, k, eval(k))
        setattr(analysis.RegionFeatureAccumulator, k, eval(k))

_genFeaturConvenienceFunctions()
del _genFeaturConvenienceFunctions






######################################
# extend cgp
######################################
import vigra

import pylab
import numpy
from  matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import matplotlib.cm as cm













#def cgpFromLabels(labels):
#    tgrid=cgp.TopologicalGrid(labels)
#    cgp=cgp.Cgp(tgrid)
#    return cgp ,tgrid
  
metaclass_more_cgp   = cgp.Cgp.__class__ 
metaclass_more_cell0 = cgp.Cell0.__class__ 
metaclass_more_cell1 = cgp.Cell1.__class__ 
metaclass_more_cell2 = cgp.Cell2.__class__ 

class injector_more_cgp(object):
    class __metaclass__(metaclass_more_cgp):
        def __init__(self, name, bases, dict):
            for b in bases:
                if type(b) not in (self, type):
                    for k,v in dict.items():
                        setattr(b,k,v)
            return type.__init__(self, name, bases, dict)

class injector_more_cell0(object):
    class __metaclass__(metaclass_more_cell0):
        def __init__(self, name, bases, dict):
            for b in bases:
                if type(b) not in (self, type):
                    for k,v in dict.items():
                        setattr(b,k,v)
            return type.__init__(self, name, bases, dict)

class injector_more_cell1(object):
    class __metaclass__(metaclass_more_cell1):
        def __init__(self, name, bases, dict):
            for b in bases:
                if type(b) not in (self, type):
                    for k,v in dict.items():
                        setattr(b,k,v)
            return type.__init__(self, name, bases, dict)

class injector_more_cell2(object):
    class __metaclass__(metaclass_more_cell2):
        def __init__(self, name, bases, dict):
            for b in bases:
                if type(b) not in (self, type):
                    for k,v in dict.items():
                        setattr(b,k,v)
            return type.__init__(self, name, bases, dict)


class more_cgp(injector_more_cgp, cgp.Cgp):


    def labelGrid(self,cellType,useTopologicalShape=True,out=None):
      lgrid=self.tgrid.labelGrid(cellType=cellType,useTopologicalShape=useTopologicalShape,out=out)
      if (cellType !=2):
        assert lgrid.min()==0
      assert lgrid.max()==self.numCells(cellType)
      return lgrid


    def pixelHistToCell2Hist(self,pixelHist):
        """  calculate a histogram for each region based
        on a pixel wise histogram 
        """  

        shape = self.shapeLabeling
        nBins = pixelHist.shape[2]
        assert shape[0]==pixelHist.shape[0]
        assert shape[1]==pixelHist.shape[1]

        labels=self.labelGrid(cellType=2,useTopologicalShape=False)

        #print "shapeLabeling",shape,";labelsShape ",labels.shape

        assert shape[0]==labels.shape[0]
        assert shape[1]==labels.shape[1]
        assert labels.min()==1
        assert labels.max()==self.numCells(2)

        inputData = numpy.require(pixelHist,dtype=numpy.float32)
        featureDict = analysis.extractRegionFeatures(image=inputData,labels=labels,features='mean',ignoreLabel=0)
        regionHist  =featureDict['mean']
        assert regionHist.ndim == 2
        assert regionHist.shape[0]==self.numCells(2)+1
        assert regionHist.shape[1]==nBins
        returnArray = regionHist[1:self.numCells(2)+1,:].copy()
        assert returnArray.shape[0]==self.numCells(2)

        regionHist=None
        featureDict=None
        
        return returnArray

    
    def accumulateCellHistogram(self,cellType,image,histogramRange=None,binCount=64,useCellMinMax=False,sigma=None):
      assert image.ndim ==2 or image.ndim==3

      data=image.reshape([image.shape[0],image.shape[1],-1])
      nChannels = data.shape[2]
      nCells    = self.numCells(cellType)

      # result array 
      cellHisto = numpy.ones([nCells,nChannels,binCount],dtype=numpy.float32)


      if histogramRange is not None:
        histogramRange=numpy.require(histogramRange,dtype=numpy.float32)
        assert histogramRange.ndim==2
        assert histogramRange.shape[0]==nChannels
        assert histogramRange.shape[1]==2


      # iterate over every channel
      for c in range(nChannels):

        # get histogram range and channel of data
        if histogramRange is None:
          hr = None 
        else :
          hr = histogramRange[c,:]
        d = data[:,:,c]

        # accumulate histogram for one(!) channel
        featureDict,activeFeatures = self.accumulateCellFeatures( cellType=cellType,image=d,features=['Histogram','Count'],
                                                                  histogramRange=hr, binCount=binCount,
                                                                  useCellMinMax=useCellMinMax,histMode=True)
        # write channel result into result array
        assert 'Histogram' in activeFeatures

        #print activeFeatures

        channelHist = featureDict['Histogram']
        channelHistCount = featureDict['Count']

        channelHistCount2=numpy.repeat(channelHistCount, binCount)
        channelHistCount2=channelHistCount2.reshape(-1,binCount)

        #print "chshape",channelHist.shape
        #print "cpunt ",channelHistCount2.shape


        #channelHistCount = 
        #print channelHist.reshape(01,channelHistCount.shape

        #channelHist=channelHist.reshape(-1)
        channelHist/=channelHistCount2



        #channelHist=channelHist.reshape([-1,binCount])
        #print "chshape",channelHist.shape
        assert channelHist.ndim == 2
        assert channelHist.shape[0]==nCells
        assert channelHist.shape[1]==binCount

        cellHisto[:,c,:]=channelHist[:,:]

      if sigma is not None:
        cellHisto2d = cellHisto.reshape([-1,binCount])
        cellHisto2d = gaussianSmoothing1d(cellHisto2d,sigma)
        cellHisto   = cellHisto2d.reshape([nCells,nChannels,binCount])

      return cellHisto








    def accumulateCellFeatures(self,cellType,image,features='all',histogramRange=None,binCount=64,useCellMinMax=False,histMode=False):
      def spatialShapeMatch(shapeA,shapeB):
        return shapeA[0]==shapeB[0] and shapeA[1]==shapeB[1]

      # check for valid input
      dataShape = image.shape
      if spatialShapeMatch(dataShape,self.shape):
        useTopologicalShape=True
      elif spatialShapeMatch(dataShape,self.shapeLabeling):
        useTopologicalShape=False
      else :
        raise RuntimeError("image.shape does neither match cgp.shape nor cgp.shapeLabeling")

      image = numpy.require(image,dtype=numpy.float32)

      nCells=self.numCells(cellType)
      #labels=self.labelGrid(cellType)
      #labels=numpy.ones(self.shape,dtype=numpy.uint32)
      labels=self.labelGrid(cellType=cellType,useTopologicalShape=useTopologicalShape)
      

      if histMode :
        hFeatures = ['Histogram','Count']
        assert image.ndim == 2
        if histogramRange  is None :
          if(useCellMinMax==False):
            histogramRange=(float(image.min()),float(image.max()))
          else:
            histogramRange='globalminmax'
        else:
          histogramRange = (float(histogramRange[0]),float(histogramRange[1]))
        values=analysis.extractRegionFeatures(image=image ,labels=labels, features=hFeatures, histogramRange=histogramRange ,binCount=binCount) 

      else:
        values=analysis.extractRegionFeatures(image=image,labels=labels,features=features,ignoreLabel=0)
      activeFeatures=values.activeFeatures()
      #del values
      
      featureDict=dict()
      for fname in activeFeatures :
        featureVals=values[fname]
        if isinstance(featureVals, numpy.ndarray) or issubclass(featureVals.__class__,numpy.ndarray):
          shape=featureVals.shape
          dim=len(shape)
          if dim==1:
            featureDict[fname]=featureVals[1:nCells+1].copy()
          elif dim==2:
            featureDict[fname]=featureVals[1:nCells+1,:].copy()
          elif dim==3:
            featureDict[fname]=featureVals[1:nCells+1,:,:].copy()
        elif isinstance(featureVals,(int ,long,float)):
          featureDict[fname]=featureVals
        else :
          raise RuntimeError("internal error in accumulateCellFeatures")

      values=None
      return featureDict,activeFeatures
    
    def cells(self,cellType):
      if(cellType==0):
        return self.cells0
      elif(cellType==1):
        return self.cells1
      elif(cellType==2):
        return self.cells2
      else:
        raise NameError("cellType must be 0,1,or 2")

    def matchMergedCgpCells(self,coarse_cgp):
      # fine cells to coarse label(s)
      cell_to_coarse=[dict() ,dict(), dict() ]
      # coarse labels to fine cell labels
      cell_to_fine=[None]*3
      cell_to_fine[0]=[None]*coarse_cgp.numCells(0)
      cell_to_fine[1]=[ list() ]*coarse_cgp.numCells(1)
      cell_to_fine[2]=[ list() ]*coarse_cgp.numCells(2)

      coarseLabeling=numpy.ones(self.shape,dtype=numpy.uint32)
      for cellType in range(3):
        coarseLabeling=coarse_cgp.labelGrid(cellType,out=coarseLabeling)
        for cell in self.cells(cellType):
          label=cell.label
          aPoint=cell.points[0]
          coarseLabel=coarseLabeling(aPoint)
          if coarseLabel!=0:
            # cell is still active in coarse graph
            cell_to_coarse[cellType][ label - 1 ]=coarseLabel
            if cellType!=0 :
              cell_to_fine[cellType][coarseLabel-1].append(label)
            else:
              cell_to_fine[cellType][coarseLabel-1]=label

      return cell_to_coarse,cell_to_fine






class _cell_helper(object):
  
  @staticmethod
  def adjacencyGen(cell):

    cgp=cell.cgp
    cellType=cell.cellType
    # get own cell label
    cellLabel=cell.label
    #get cells below
    assert cellType!=0
    cellsBelow=cgp.cells(cellType-1)

    for boundedByCellLabel in cell.boundedBy:
      # index of junction
      boundedByCellIndex=boundedByCellLabel-1
      # get bounds of boundedByCell
      bounds = cellsBelow[boundedByCellIndex].bounds
      for otherCellLabel in bounds:
        if otherCellLabel != cellLabel:
          yield otherCellLabel , boundedByCellLabel

  @staticmethod
  def adjacentCellsGen(cell):
    cells=cell.cgp.cells(cell.cellType)
    for cellLabel in _cell_helper.adjacencyGen(cell):
      yield cells[cellLabel-1]

  @staticmethod
  def boundingCellsGen(cell):
    assert cell.cellType <=1
    higherCells=cell.cgp.cells(cell.cellType+1)
    for label in cell.bounds:
      yield higherCells[label-1]

  @staticmethod
  def boundingByCellsGen(cell):
    assert cell.cellType >=1
    lowerCells=cell.cgp.cells(cell.cellType-1)
    for label in cell.boundedBy:
      yield lowerCells[label-1]


class more_cell0(injector_more_cell0,cgp.Cell0):

  def boundingCellsGen(self):
    return _cell_helper.boundingCellsGen(self)

class more_cell1(injector_more_cell1,cgp.Cell1):

  def adjacencyGen(self):
    return _cell_helper.adjacencyGen(self)
  def adjacentCellsGen(self):
    return _cell_helper.adjacentCellsGen(self)
  def boundingCellsGen(self):
    return _cell_helper.boundingCellsGen(self)
  def boundedByCellsGen(self):
    return _cell_helper.boundedByCellsGen(self)

class more_cell2(injector_more_cell2,cgp.Cell2):

  def adjacencyGen(self):
    return _cell_helper.adjacencyGen(self)
  def adjacentCellsGen(self):
    return _cell_helper.adjacentCellsGen(self)
  def boundedByCellsGen(self):
    return _cell_helper.boundedByCellsGen(self)










def visualize(
    img_rgb,
    cgp,
    edge_data_in=None,
    show=True,
    cmap=cm.jet,
    title=None,
    black=False
):
    img_rgb_raw=img_rgb.copy()
    if edge_data_in is not None:
        edge_data=edge_data_in.copy()
    else:
        edge_data=None
    img_rgb=numpy.squeeze(img_rgb)
    if img_rgb.ndim == 2:
      img = numpy.ones([img_rgb.shape[0],img_rgb.shape[1],3 ])
      img[:,:,0] = img_rgb[:,:]
      img[:,:,1] = img_rgb[:,:]
      img[:,:,2] = img_rgb[:,:]
    else :
      img=img_rgb.copy()
    img-=img.min()
    img/=img.max()
    # get edge map
    edgeMarkers=cgp.labelGrid(1,True)
    whereEdges=numpy.where(edgeMarkers!=0)
    edgeMarkers[whereEdges]=1

    


    if edge_data is not None :

        #edge_data=numpy.sqrt(edge_data)

        resMin=numpy.min(edge_data)
        resMax=numpy.max(edge_data)

        #print "mi ma",resMin,resMax
        edge_data[:]=(edge_data[:]-resMin)/(resMax-resMin)


        resImg=cgp.featureToImage(cellType=1,features=edge_data,ignoreInactive=False,inactiveValue=0.0)


        edgeValues=resImg[whereEdges]

        # transform 
        mycm=cm.ScalarMappable(norm=None, cmap=cmap)
        mycm.set_array(edgeValues.reshape(-1))

        colorArray=mycm.to_rgba(edgeValues)

        #print " shape ",colorArray.shape
        #print colorArray
        #img*=255
        if black:
          img[:]=0.0
        img[whereEdges[0],whereEdges[1],0]=colorArray[:,0]
        img[whereEdges[0],whereEdges[1],1]=colorArray[:,1]
        img[whereEdges[0],whereEdges[1],2]=colorArray[:,2]
        """
        img[whereEdges[0],whereEdges[1],0]=(1.0-resImg[whereEdges[0],whereEdges[1]])
        img[whereEdges[0],whereEdges[1],1]=0.0#resImg[whereEdges[0],whereEdges[1]]*255.0
        img[whereEdges[0],whereEdges[1],2]=resImg[whereEdges[0],whereEdges[1]]
        """

    elif edge_data is None :
      labelImage=cgp.tgrid.labelGrid(2,False)

      cedge     = vigra.analysis.regionImageToCrackEdgeImage(numpy.require(labelImage,dtype=numpy.uint32))

      #cedge[cedge!=0]=0
      whereEdges=numpy.where(cedge==0)
      #img/=255
      img[whereEdges[0],whereEdges[1],0]=0.0
      img[whereEdges[0],whereEdges[1],1]=0.0
      img[whereEdges[0],whereEdges[1],2]=0.0

    else :
        #img#/=255
        #img[whereEdges[0],whereEdges[1],0]=0.0
        #img[whereEdges[0],whereEdges[1],1]=0.0
        #img[whereEdges[0],whereEdges[1],2]=0.0
        #edgeData=numpy.ones()
        #resImg=cgp.featureToImage(cellType=1,features=whereEdges.astype(numpy.float32),ignoreInactive=False,inactiveValue=0.0)
        resImg=vigra.filters.discDilation(edgeMarkers.astype(numpy.uint8),1)

        whereEdges=numpy.where(resImg!=0)

        img[whereEdges[0],whereEdges[1],0]=0.0
        img[whereEdges[0],whereEdges[1],1]=0.0
        img[whereEdges[0],whereEdges[1],2]=0.0



    f = pylab.figure()
    for n, iimg in enumerate([img,img_rgb_raw/255]):
        #f.add_subplot(2, 1, n)  # this line outputs images on top of each other
        f.add_subplot(1, 2, n)  # this line outputs images side-by-side
        plt.imshow(numpy.swapaxes(iimg,0,1))



    #plt.imshow(  numpy.flipud(numpy.rot90(img) )  ,interpolation=None)
    if title is not None:
      plt.title(title)
    if(show):
        plt.show()



