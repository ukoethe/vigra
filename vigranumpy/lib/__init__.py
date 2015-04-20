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

import sys, os, time

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
* histogram  (histograms and channel representation)
* graphs     (grid graphs / graphs / graph algorithms)
* utilities  (priority queues)
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
import histogram
import graphs
import utilities
import blockwise

sampling.ImagePyramid = arraytypes.ImagePyramid



class Timer:    
    def __init__(self, name, verbose=True):
        self.name = name 
        self.verbose = verbose

    def __enter__(self):
        if self.verbose:
            print self.name, "..."
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.interval = self.end - self.start
        if self.verbose  :
            print "... took ", self.interval, "sec"







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

from vigranumpycore import ChunkedArrayFull, ChunkedArrayLazy, ChunkedArrayCompressed, ChunkedArrayTmpFile, Compression
try:
    from vigranumpycore import ChunkedArrayHDF5, HDF5Mode
except:
    pass
    

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
def imshow(image,show=True):
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
        if show:
            matplotlib.pylab.show()
        return plot
    elif image.channels == 3:
        if image.dtype != numpy.uint8:
            out = image.__class__(image.shape, dtype=numpy.uint8, axistags=image.axistags)
            image = colors.linearRangeMapping(image, newRange=(0.0, 255.0), out=out)
        plot = matplotlib.pyplot.imshow(image.view(numpy.ndarray))
        if show:
            matplotlib.pylab.show()
        return plot
    else:
        raise RuntimeError("vigra.imshow(): Image must have 1 or 3 channels.")


def segShow(img,labels,edgeColor=(0,0,0),alpha=0.3,show=False,returnImg=False,r=0):

    img = numpy.squeeze(img)
    if img.ndim ==2:
        img = numpy.concatenate( [ img[:,:,None]]*3 ,axis=2).astype(numpy.float32)
        img = taggedView(img, 'xyc')

    labels = numpy.squeeze(labels)
    crackedEdges = analysis.regionImageToCrackEdgeImage(labels).squeeze()
    #print "cracked shape",crackedEdges.shape
    whereEdge    =  numpy.where(crackedEdges==0)
    whereNoEdge  =  numpy.where(crackedEdges!=0)
    crackedEdges[whereEdge] = 1
    crackedEdges[whereNoEdge] = 0

    if r>0 :
        res = filters.discDilation(crackedEdges.astype(numpy.uint8),int(r) )
        whereEdge  =  numpy.where(res==1)

    imgToDisplay = resize(img,numpy.squeeze(crackedEdges).shape)
    imgToDisplay-=imgToDisplay.min()
    imgToDisplay/=imgToDisplay.max()
    for c in range(3):
        ic = imgToDisplay[:,:,c]
        ic[whereEdge]=(1.0-alpha)*edgeColor[c] + alpha*ic[whereEdge]

    if returnImg:
        return imgToDisplay
    return imshow(imgToDisplay,show=show)

def nestedSegShow(img,labels,edgeColors=None,scale=1,show=False,returnImg=False):

    shape=(labels.shape[0]*scale,labels.shape[1]*scale)
    if scale!=1:
        img=vigra.resize(img,shape)
    



    assert numpy.squeeze(labels).ndim==3
    nSegs  = labels.shape[2]


    if edgeColors is None :
      edgeColors=numpy.ones([nSegs,4])

      a  =numpy.array([0,0,0.0,0.6],dtype=numpy.float32)
      b  =numpy.array([1,0,0,0.4],dtype=numpy.float32)

      for s in range(nSegs):
        f=float(s)/float(nSegs-1)
        edgeColors[s,:]=f*b + (1.0-f)*a

    tShape=(img.shape[0]*2-1,img.shape[1]*2-1)

    imgToDisplay = resize(img,tShape)
    imgToDisplay-=imgToDisplay.min()
    imgToDisplay/=imgToDisplay.max()

    imgIn = imgToDisplay.copy()

    for si in range(nSegs):
        l = labels[:,:,si].copy()
        if scale!=1:
            l=resize(l.astype(numpy.float32),shape,order=0).astype(numpy.uint32)

        crackedEdges = analysis.regionImageToCrackEdgeImage(l)
        whereEdge    = numpy.where(crackedEdges==0)

   
        if len(edgeColors[si])<4:
            alpha = 0.0
        else:
            alpha = edgeColors[si,3]
        for c in range(3):
            icI = imgIn[:,:,c]
            ic  = imgToDisplay[:,:,c]
            ic[whereEdge]=(1.0-alpha) * edgeColors[si,c] + alpha*icI[whereEdge]
    if returnImg:
        return imgToDisplay
    return imshow(imgToDisplay,show=show)


def show():
    import matplotlib.pylab
    matplotlib.pylab.show()

        
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



# define watershedsReoptimization)
def _genWatershedsReoptimization():
    def watershedsReoptimization(labels,edgeIndicator,shrinkN,out=None,visu=False):
        # do unseeding

        #if visu :
        #  import matplotlib,numpy
        #  import pylab
        #  # A random colormap for matplotlib
        #  cmap = matplotlib.colors.ListedColormap ( numpy.random.rand ( 256,3))
        #  pylab.imshow ( numpy.swapaxes(labels,0,1) , cmap = cmap)
        #  pylab.show()


        seeds=analysis.segToSeeds(labels,long(shrinkN))

        if visu :
          import matplotlib,numpy
          import pylab
          # A random colormap for matplotlib
          cmap = matplotlib.colors.ListedColormap ( numpy.random.rand ( 256,3))
          pylab.imshow ( numpy.swapaxes(seeds,0,1) , cmap = cmap)
          pylab.show()



        #if seeds.ndim==2:
        #    seeds=analysis.labelImageWithBackground(seeds)
        #elif seeds.ndim==3:
        #    seeds=analysis.labelVolumeWithBackground(seeds)
        #else :
        #    raise RuntimeError("only implemented for 2d and 3d")

        if visu :
          import matplotlib,numpy
          import pylab
          # A random colormap for matplotlib
          cmap = matplotlib.colors.ListedColormap ( numpy.random.rand ( 256,3))
          pylab.imshow ( numpy.swapaxes(seeds,0,1) , cmap = cmap)
          pylab.show()

        return analysis.watersheds(edgeIndicator,seeds=seeds,out=out)
        
    watershedsReoptimization.__module__ = 'vigra.analysis'
    analysis.watershedsReoptimization = watershedsReoptimization

_genWatershedsReoptimization()
del _genWatershedsReoptimization


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

    def supportedConvexHullFeatures(labels):
        '''Return a list of Convex Hull feature names that are available for the given 2D label array. 
           These Convex Hull feature names are the valid inputs to a call of 
           :func:`extractConvexHullFeatures`. E.g., to compute just the first two features in the 
           list, use::
           
                f = vigra.analysis.supportedConvexHullFeatures(labels)
                print "Computing Convex Hull features:", f[:2]
                r = vigra.analysis.extractConvexHullFeatures(labels, features=f[:2])
        '''
        try:
            return analysis.extractConvexHullFeatures(labels, list_features_only=True)
        except:
            return []

    supportedConvexHullFeatures.__module__ = 'vigra.analysis'
    analysis.supportedConvexHullFeatures = supportedConvexHullFeatures

    def supportedSkeletonFeatures(labels):
        '''Return a list of Skeleton feature names that are available for the given 2D label array. 
           These Skeleton feature names are the valid inputs to a call of 
           :func:`extractSkeletonFeatures`. E.g., to compute just the first two features in the 
           list, use::
           
                f = vigra.analysis.supportedSkeletonFeatures(labels)
                print "Computing Skeleton features:", f[:2]
                r = vigra.analysis.extractSkeletonFeatures(labels, features=f[:2])
        '''
        try:
            return analysis.extractSkeletonFeatures(labels, list_features_only=True)
        except:
            return []

    supportedSkeletonFeatures.__module__ = 'vigra.analysis'
    analysis.supportedSkeletonFeatures = supportedSkeletonFeatures

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


MetricType = graphs.MetricType
# define grid graph convenience functions
# and extend grid graph classes
def _genGridGraphConvenienceFunctions():

    def gridGraph(shape,directNeighborhood=True):
        '''Return a grid graph with certain shape. 

            Parameters:

                - shape -- shape of the image
                - directNeighborhood -- use  4 (True) or 8 (False) neighborhood (default: True)

            Returns:

                - grid graph

            use::

                >>> # 4-connected
                >>> g = vigra.graps.gridGraph(shape=[10,20])
                >>> g.nodeNum
                200
                >>> # 8-connected
                >>> g = vigra.graps.gridGraph(shape=[10,20],directNeighborhood=False)

        '''
        if(len(shape)==2):
            return graphs.GridGraphUndirected2d(shape,directNeighborhood)
        elif(len(shape)==3):
            return graphs.GridGraphUndirected3d(shape,directNeighborhood)
        else:
            raise RuntimeError("GridGraph is only implemented for 2d and 3d grids")

    gridGraph.__module__ = 'vigra.graphs'
    graphs.gridGraph = gridGraph

    # extend grid graph via meta classes
    for cls in [graphs.GridGraphUndirected2d, graphs.GridGraphUndirected3d] :

        metaCls = cls.__class__

        class gridGraphInjector(object):
            class __metaclass__(metaCls):
                def __init__(self, name, bases, dict):
                    for b in bases:
                        if type(b) not in (self, type):
                            for k,v in dict.items():
                                setattr(b,k,v)
                    return type.__init__(self, name, bases, dict)

        ##inject some methods in the point foo
        class moreGridGraph(gridGraphInjector, cls):
            
            @property
            def shape(self):
                """ shape of grid graph"""
                return self.intrinsicNodeMapShape()

            
            def nodeSize(self):
                """ node map filled with 1.0"""
                size = graphs.graphMap(self,item='node',dtype=numpy.float32)
                size[:]=1
                return size

            def edgeLengths(self):
                """ node map filled with 1.0"""
                size = graphs.graphMap(self,item='edge',dtype=numpy.float32)
                size[:]=1
                return size

            def mergeGraph(self):
                if  len(self.shape)==2:
                    mg = graphs.GridGraphUndirected2dMergeGraph(self)
                else:
                    mg =  graphs.GridGraphUndirected3dMergeGraph(self)
                return mg


    def isGridGraph(obj):
        """ check if obj is gridGraph"""
        return isinstance(obj,(graphs.GridGraphUndirected2d , graphs.GridGraphUndirected3d))
    def isGridGraph2d(obj):
        """ check if obj is gridGraph"""
        return isinstance(obj,graphs.GridGraphUndirected2d)

    isGridGraph.__module__ = 'vigra.graphs'
    graphs.isGridGraph = isGridGraph  

    isGridGraph2d.__module__ = 'vigra.graphs'
    graphs.isGridGraph2d = isGridGraph2d  


_genGridGraphConvenienceFunctions()
del _genGridGraphConvenienceFunctions



def _genGraphConvenienceFunctions():

    def listGraph(nodes=0,edges=0):
        ''' Return an empty directed graph

            Parameters :

                - nodes : number of nodes to reserveEdges
                - edges : number of edges to reserve

            Returns :

                - graph
        '''
        return graphs.AdjacencyListGraph(nodes,edges)
        
    listGraph.__module__ = 'vigra.graphs'
    graphs.listGraph = listGraph

    def intrinsicGraphMapShape(graph,item):
        """ Intrinsic shape of node/edge/arc-map for a given graph.

            Node edge and arc maps are stored in numpy arrays by default.
            The instric shape may not be confused with the number
            of nodes/edges/arcs. The instric shape is used to
            allocate a numpy are which can store data for nodes/arcs/edgeSizes
            of a given graph.

            Parameters:

                - graph : input graph to get the shape for

                - item  : item must be ``'node'`` , ``'edge'`` or ``'arc'``

            Returns:

                - shape as tuple
        """
        if   item=='edge':
            return graph.intrinsicEdgeMapShape()
        elif item=='node':
            return graph.intrinsicNodeMapShape()
        elif item=='arc': 
            return graph.intrinsicArcMapShape()
        else :
            raise RuntimeError("%s is not valid,must be 'edge','node' or 'arc' "%item)

    intrinsicGraphMapShape.__module__ = 'vigra.graphs'
    graphs.intrinsicGraphMapShape = intrinsicGraphMapShape


    def graphMap(graph,item,dtype=numpy.float32,channels=1,addChannelDim=False):
        """ Return a graph map for a given graph item (``'node'`` , ``'edge'`` or ``'arc'``).

            Parameters:

                - graph    : graph to get a graph map for
                - item     : ``'node'`` , ``'edge'`` or ``'arc'``
                - dtype    : desired dtype 
                - channels : number of channels (default: 1)
                - addChannelDim -- add an explicit channelDim :(default: False)
                    only useful if channels == 1 

            Returns:

                - graphmap as numpy.ndarray / VigraArray
        """
        s = intrinsicGraphMapShape(graph,item)
        intrDim = len(s)
        if(channels==1) and addChannelDim==False:
            a=numpy.zeros(shape=s,dtype=dtype)
            if intrDim == 1:
                return taggedView(a,'x')
            elif intrDim == 2:
                return taggedView(a,'xy')
            elif intrDim == 3:
                return taggedView(a,'xyz')
            elif intrDim == 4:
                return taggedView(a,'xyzt')
            else :
                raise RuntimeError("graphs with intrisic dimension >4 are not supported")
        else:
            s = s+(channels,)
            a=numpy.zeros(shape=s,dtype=dtype)
            if intrDim == 1:
                return taggedView(a,'xc')
            elif intrDim == 2:
                return taggedView(a,'xyc')
            elif intrDim == 3:
                return taggedView(a,'xyzc')
            elif intrDim == 4:
                return taggedView(a,'xyztc')
            else :
                raise RuntimeError("graphs with intrisic dimension >4 are not supported")



    def graphMap2(graph,item,dtype=numpy.float32,channels=1,addChannelDim=False):
        """ Return a graph map for a given graph item (``'node'`` , ``'edge'`` or ``'arc'``).

            Parameters:

                - graph    : graph to get a graph map for
                - item     : ``'node'`` , ``'edge'`` or ``'arc'``
                - dtype    : desired dtype 
                - channels : number of channels (default: 1)
                - addChannelDim -- add an explicit channelDim :(default: False)
                    only useful if channels == 1 

            Returns:

                - graphmap as numpy.ndarray / VigraArray
        """
        s = intrinsicGraphMapShape(graph,item)
        intrDim = len(s)
        if(channels==1) and addChannelDim==False:
            a=numpy.zeros(shape=s,dtype=dtype)
            if intrDim == 1:
                return taggedView(a,'x')
            elif intrDim == 2:
                return taggedView(a,'xy')
            elif intrDim == 3:
                return taggedView(a,'xyz')
            elif intrDim == 4:
                return taggedView(a,'xyzt')
            else :
                raise RuntimeError("graphs with intrisic dimension >4 are not supported")
        else:
            s = s+(channels,)
            a=numpy.zeros(shape=s,dtype=dtype)
            if intrDim == 1:
                return taggedView(a,'xc')
            elif intrDim == 2:
                return taggedView(a,'xyc')
            elif intrDim == 3:
                return taggedView(a,'xyzc')
            elif intrDim == 4:
                return taggedView(a,'xyztc')
            else :
                raise RuntimeError("graphs with intrisic dimension >4 are not supported")

    graphMap.__module__ = 'vigra.graphs'
    graphs.graphMap = graphMap


    def mergeGraph(graph):
        """ get a merge graph from input graph.

            A merge graph might be usefull for hierarchical clustering
        """
        #mg = graph.mergeGraph()
        mg = graphs.__mergeGraph(graph)
        #mg.__base_graph__=graph
        return mg

    mergeGraph.__module__ = 'vigra.graphs'
    graphs.mergeGraph = mergeGraph


    INVALID = graphs.Invalid()
    graphs.INVALID = INVALID



    class ShortestPathPathDijkstra(object):
        def __init__(self,graph):
            """ shortest path computer

                Keyword Arguments:

                    - graph : input graph

            """
            self.pathFinder =  graphs._shortestPathDijkstra(graph)
            self.graph=graph
            self.source = None
            self.target = None
        def run(self,weights,source,target=None):
            """ run shortest path search

                Keyword Arguments:

                   - weights : edge weights encoding distance from two adjacent nodes

                   - source : source node 

                   - target : target node (default: None)
                        If target node is None, the shortest path
                        to all nodes!=source is computed

            """
            self.source = source
            self.target = target
            if target is None:
                self.pathFinder.run(weights,source)
            else:
                self.pathFinder.run(weights,source,target)
            return self

        def runIgnoreLargeWeights(self,weights,source,val):
            """ run shortest path search, nodes with all edge weights larger than val will be ignored

                Keyword Arguments:

                   - weights : edge weights encoding distance from two adjacent nodes

                   - source : source node
                   
                   - val : upper bound

            """
            self.source = source
            self.target = None
            self.pathFinder.runIgnoreLargeWeights(weights,source,val)
            return self

        def path(self,target=None,pathType='coordinates'):
            """ get the shortest path from source to target

                Keyword Arguments:

                    - weights : edge weights encoding distance from two adjacent nodes

                    - source : source node 

                    - target : target node (default: None)
                        If target node is None, the target specified
                        by 'run' is used.

                    pathType : 'coordinates' or 'ids' path (default: 'coordinates')



            """
            if target is None:
                assert self.target is not None
                target=self.target

            if pathType=='coordinates':
                return self.pathFinder.nodeCoordinatePath(target)
            elif pathType == 'ids':
                return self.pathFinder.nodeIdPath(target)
        def distance(self,target=None):
            """ get distance from source to target

                Keyword Arguments:
                    - target : target node (default: None)
                        If target node is None, the target specified
                        by 'run' is used.
            """
            if target is None:
                assert self.target is not None
                target=self.target
            return self.pathFinder.distance(target)

        def distances(self,out=None):
            """ return the full distance map"""
            return self.pathFinder.distances(out)
        def predecessors(self,out=None):
            """ return the full predecessors map"""
            return self.pathFinder.predecessors(out)


    ShortestPathPathDijkstra.__module__ = 'vigra.graphs'
    graphs.ShortestPathPathDijkstra = ShortestPathPathDijkstra

_genGraphConvenienceFunctions()
del _genGraphConvenienceFunctions


def _genRegionAdjacencyGraphConvenienceFunctions():



    class RegionAdjacencyGraph(graphs.AdjacencyListGraph):
        def __init__(self,graph=None ,labels=None ,ignoreLabel=None,reserveEdges=0, maxLabel=None, isDense=None):
            """ Region adjacency graph

                Keyword Arguments :
                    - graph : the base graph, the region adjacency graph should be based on

                    - labels : label map for the graph

                    - ignoreLabel : ignore a label in the labels map (default: None)

                    - reserveEdges : reserve a certain number of Edges

                Attributes:

                    - labels : labels passed in constructor

                    - ignoreLabel  : ignoreLabel passed in constructor

                    - baseGraphLabels : labels passed in constructor 
                        (fixme,dublicated attribute (see labels) )

                    - baseGraph : baseGraph is the graph passed in constructor

                    - affiliatedEdges : for each edge in the region adjacency graph,
                        a vector of edges of the baseGraph is stored in affiliatedEdges


            """
            if(graph is not None and labels is not None):
                super(RegionAdjacencyGraph,self).__init__(long(labels.max()+1),long(reserveEdges))

                if ignoreLabel is None and isDense is not None and isDense == True:
                    if ignoreLabel is None:
                        ignoreLabel=-1

                    self.labels          = labels
                    self.ignoreLabel     = ignoreLabel
                    self.baseGraphLabels = labels
                    self.baseGraph       = graph
                    if maxLabel is None:
                        maxLabel = int(numpy.max(labels))
                    # set up rag
                    self.affiliatedEdges = graphs._regionAdjacencyGraphFast(graph,labels,self,maxLabel,int(reserveEdges))

                else:

                    if ignoreLabel is None:
                        ignoreLabel=-1

                    self.labels          = labels
                    self.ignoreLabel     = ignoreLabel
                    self.baseGraphLabels = labels
                    self.baseGraph       = graph
                    # set up rag
                    self.affiliatedEdges = graphs._regionAdjacencyGraph(graph,labels,self,self.ignoreLabel)   
            else :
                super(RegionAdjacencyGraph,self).__init__(0,0)

        def mergeGraph(self):
            return graphs.AdjacencyListGraphMergeGraph(self)

        def accumulateSeeds(self, seeds, out=None):
            graph = self.baseGraph
            labels = self.labels
            return graphs._pyAccNodeSeeds(self, graph, labels, seeds, out)

        def accumulateEdgeFeatures(self,edgeFeatures,acc='mean',out=None):
            """ accumulate edge features from base graphs edges features

                Keyword Argument:

                    - edgeFeatures : edge features of baseGraph
                    - acc : used accumulator (default: 'mean')
                        Currently only 'mean' and 'sum' are implemented
                    - out :  preallocated edge map

                Returns :
                    accumulated edge features
            """
            graph = self.baseGraph
            affiliatedEdges = self.affiliatedEdges

            if isinstance(edgeFeatures, (graphs.ImplicitMEanEdgeMap_2d_float_float, graphs.ImplicitMEanEdgeMap_3d_float_float)):


                if graphs.isGridGraph(graph)==False:
                    raise RuntimeError("implicit edge maps are only implemented for grid graphs")

                return graphs._ragEdgeFeatures(self, graph, affiliatedEdges, edgeFeatures,acc, out)

            else:
                if self.edgeNum == 0:
                    raise RuntimeError("self.edgeNum == 0  => cannot accumulate edge features")
                if acc == 'mean':
                    weights = self.baseGraph.edgeLengths()
                else:
                    weights = graphs.graphMap(self.baseGraph,'edge',dtype=numpy.float32)
                    weights[:] = 1
                if graphs.isGridGraph2d(graph) and edgeFeatures.ndim == 4 :
                    return graphs._ragEdgeFeaturesMb(self,graph,affiliatedEdges,edgeFeatures,weights,acc,out)
                else:
                    return graphs._ragEdgeFeatures(self,graph,affiliatedEdges,edgeFeatures,weights,acc,out)


        def accumulateNodeFeatures(self,nodeFeatures,acc='mean',out=None):
            """ accumulate edge features from base graphs edges features

                Keyword Argument:

                    - nodeFeatures : node features of baseGraph
                    - acc : used accumulator (default: 'mean')
                        Currently only 'mean' and 'sum' are implemented
                    - out :  preallocated node map (default: None)

                Returns : 
                    accumulated node features
            """
            if self.edgeNum == 0 :
              raise RuntimeError("self.edgeNum == 0  => cannot accumulate edge features")
            graph = self.baseGraph
            labels = self.baseGraphLabels
            ignoreLabel = self.ignoreLabel
            if acc == 'mean':
              #print "get node size..."
              weights = self.baseGraph.nodeSize()
              #print "weights == ", weights
            else :
              weights = graphs.graphMap(self.baseGraph,'node',dtype=numpy.float32)
              weights[:]=1

            return graphs._ragNodeFeatures(self,graph,labels,nodeFeatures,weights,acc,ignoreLabel,out)

        def projectNodeFeatureToBaseGraph(self,features,out=None):
            """ project node features from this graph, to the base graph of this graph.

                Keyword Arguments:

                    - features : node feautres for this graph
                    - out :  preallocated node map of baseGraph (default: None)

                Returns :
                    projected node features of base graph
            """
            out=graphs._ragProjectNodeFeaturesToBaseGraph(
                rag=self,
                baseGraph=self.baseGraph,
                baseGraphLabels=numpy.squeeze(self.baseGraphLabels),
                ragNodeFeatures=features,
                ignoreLabel=self.ignoreLabel,
                out=out
            )
            #print "out",out.shape,out.dtype
            return out

        def projectLabelsBack(self,steps,labels=None,_current=0):
            """  project labels from current graph to baseGraph and repeat this recursively

                Keyword  Arguments:

                    - steps : how often should the labels be projected back
                    - labels : labels for the current graph (default: None)
                        If labels is None, each node gets its own label
            """
            if labels is None :
                # identity segmentation on this level
                labels = self.nodeIdMap()

            if steps == current :
                return labels
            else :
                labels = self.projectLabelsToBaseGraph(labels)
                return self.baseGraph.projectLabelsBack(steps,labels,_current+1)


        def projectLabelsToBaseGraph(self,labels=None):
            """ project node labels from this graph, to the base graph of this graph.

                Keyword Arguments:

                    - labels : node labels for this graph (default: None)
                        If labels is None, each node gets its own label
                    - out :  preallocated node map of baseGraph (default: None)

                Returns :
            """
            if labels is None :
                # identity segmentation on this level
                labels = self.nodeIdMap()
            return self.projectNodeFeatureToBaseGraph(features=labels)

        def projectBaseGraphGt(self, baseGraphGt, gt=None, gtQuality=None):
            bggt = numpy.require(baseGraphGt,dtype=numpy.uint32)
            gt, gtQuality = graphs._ragProjectGroundTruth(rag=self, graph=self.baseGraph,
                                                          labels=self.baseGraphLabels, gt=bggt,
                                                          ragGt=gt, ragGtQuality=gtQuality)
            return gt, gtQuality


        def edgeUVCoordinates(self, edgeId):

            try :
                ei = int(edgeId)
            except:
                ei = edgeId.id

            affEdges = self.affiliatedEdges
            uvCoords = affEdges.getUVCoordinates(self.baseGraph, ei)
            dim = uvCoords.shape[1]/2
            uCoords = uvCoords[:,0:dim]
            vCoords = uvCoords[:,dim:2*dim]
            return (uCoords,vCoords)

        def edgeTopologicalCoordinates(self, edgeId):
            uc,vc = self.edgeUVCoordinates(edgeId)
            return uc+vc

        def edgeCoordinates(self, edgeId):
            uc,vc = self.edgeUVCoordinates(edgeId)
            return (uc+vc)/2.0

    RegionAdjacencyGraph.__module__ = 'vigra.graphs'
    graphs.RegionAdjacencyGraph = RegionAdjacencyGraph

    class GridRegionAdjacencyGraph(graphs.RegionAdjacencyGraph):
        def __init__(self,graph=None,labels=None,ignoreLabel=None,reserveEdges=0, maxLabel=None, isDense=None):
            """ Grid Region adjacency graph

                A region adjaceny graph,where the base graph should be
                a grid graph or a GridRegionAdjacencyGraph.


                Keyword Arguments :
                    - graph : the base graph, the region adjacency graph should be based on

                    - labels : label map for the graph

                    - ignoreLabel : ignore a label in the labels map (default: None)

                    - reserveEdges : reserve a certain number of Edges

                Attributes :

                    - labels : labels passed in constructor

                    - ignoreLabel  : ignoreLabel passed in constructor

                    - baseGraphLabels : labels passed in constructor 
                        (fixme,dublicated attribute (see labels) )

                    - baseGraph : baseGraph is the graph passed in constructor

                    - affiliatedEdges : for each edge in the region adjacency graph,
                        a vector of edges of the baseGraph is stored in affiliatedEdges

                    - shape : shape of the grid graph which is a base graph in the
                        complete graph chain.


            """
            if graph is not None and labels is not None:
                if not (graphs.isGridGraph(graph) or  isinstance(graph,GridRegionAdjacencyGraph)):
                    raise RuntimeError("graph must be a GridGraph or a GridRegionAdjacencyGraph")
                super(GridRegionAdjacencyGraph, self).__init__(graph, labels, ignoreLabel, reserveEdges, maxLabel, isDense)
            else:
                super(GridRegionAdjacencyGraph, self).__init__()

        @property
        def shape(self):
            """ shape of the underlying grid graph"""
            return self.baseGraph.shape

        def projectLabelsToGridGraph(self,labels=None):
            """project labels of this graph to the underlying grid graph.

                Keyword Arguments :

                    - labels : node labeling of this graph (default: None)
                        If labels is None, each node gets its own label

                Returns :
                    grid graph labeling

            """
            if labels is None :
                # identity segmentation on this level
                labels = self.nodeIdMap()

            if graphs.isGridGraph(self.baseGraph):
                return self.projectLabelsToBaseGraph(labels)
            else :
                labels = self.projectLabelsToBaseGraph(labels)
                return self.baseGraph.projectLabelsToGridGraph(labels)

        def projectNodeFeaturesToGridGraph(self,features):
            """ project features of this graph to the underlying grid graph.
                Therefore project the features to an image.

                Keyword Arguments :

                    - features : nodeFeatures of the current graph

                Returns :
                    grid graph labeling

            """
            if graphs.isGridGraph(self.baseGraph):
                return self.projectNodeFeatureToBaseGraph(features)
            else :
                features = self.projectNodeFeatureToBaseGraph(features)
                return self.baseGraph.projectNodeFeaturesToGridGraph(features)

        def showNested(self,img,labels=None,returnImg=False):
            """ show the complet graph chain  / hierarchy given an RGB image

                Keyword Arguments:
                    - img : RGB image 

                    - labels : node labeling of this graph (default: None)
                        If labels is None, each node gets its own label
            """
            ll=[]
            if labels is not None:
              ll.append( self.projectLabelsToGridGraph(labels) )
            ll.append( self.projectLabelsToGridGraph() )

            g=self.baseGraph
            while graphs.isGridGraph(g)==False:
              ll.append( g.projectLabelsToGridGraph() )
              g=g.baseGraph


            ll.reverse()
            gridLabels = [l[...,numpy.newaxis] for l in ll ]
            gridLabels = numpy.concatenate(gridLabels,axis=2)


            return nestedSegShow(img,gridLabels,returnImg=returnImg)


        def show(self,img,labels=None,edgeColor=(0,0,0),alpha=0.3,returnImg=False):
            """ show the graph given an RGB image

                Keyword Arguments:
                    - img : RGB image 

                    - labels : node labeling of this graph (default: None)
                        If labels is None, each node gets its own label

                    - edgeColor : RGB tuple of edge color (default: (0,0,0) ).
                        Do not use values bigger than 1 in edgeColor.

                    - alpha : make edges semi transparent (default: 0.3).
                        0 means no transparency,1 means full transparency.
            """
            pLabels = self.projectLabelsToGridGraph(labels)
            return segShow(img,numpy.squeeze(pLabels),edgeColor=edgeColor,alpha=alpha,returnImg=returnImg)

        def nodeSize(self):
            """ get the geometric size of the nodes """
            baseNodeSizes = self.baseGraph.nodeSize()
            return self.accumulateNodeFeatures(baseNodeSizes,acc='sum')
        def edgeLengths(self):
            """ get the geometric length of the edges"""
            baseNodeSizes = self.baseGraph.edgeLengths()
            return self.accumulateEdgeFeatures(baseNodeSizes,acc='sum')


        def writeHDF5(self, filename, dset):
            if(graphs.isGridGraph(self.baseGraph)):

                sGraph    = self.serialize()
                sAffEdges = graphs._serialzieGridGraphAffiliatedEdges(self.baseGraph, self, self.affiliatedEdges )
                sLabels   = self.labels


                writeHDF5(numpy.array([self.ignoreLabel]), filename, dset+'/ignore_label')
                writeHDF5(sLabels, filename, dset+'/labels')
                writeHDF5(sGraph, filename, dset+'/graph')
                writeHDF5(sAffEdges, filename, dset+'/affiliated_edges')
                

            else:
                raise RuntimeError("only RAGs of Grid graph can be serialized")


        #def readHdf5(self, filename, dset):
        #    labels = readHdf5(filename,  dset+'/labels')
        #    shape = labels.shape
        #    self.baseGraph  = graphs.gridGraph(shape)



    GridRegionAdjacencyGraph.__module__ = 'vigra.graphs'
    graphs.GridRegionAdjacencyGraph = GridRegionAdjacencyGraph



    def loadGridRagHDF5(filename , dset):

        #print "load labels and make grid graph"
        labels = readHDF5(filename,  dset+'/labels')
        shape = labels.shape
        gridGraph = graphs.gridGraph(shape)
        #print gridGraph


        #print "load graph serialization"
        graphSerialization = readHDF5(filename, dset+'/graph')

        #print "make empty grid rag"
        gridRag = GridRegionAdjacencyGraph()

        #print "deserialize"
        gridRag.deserialize(graphSerialization)


        #print "load affiliatedEdges"
        affEdgeSerialization = readHDF5(filename, dset+'/affiliated_edges')
        #print "deserialize"
        affiliatedEdges = graphs._deserialzieGridGraphAffiliatedEdges(gridGraph, gridRag, affEdgeSerialization)


        ignoreLabel =  readHDF5(filename, dset+'/ignore_label')

        gridRag.affiliatedEdges = affiliatedEdges
        gridRag.labels          = taggedView(labels,"xyz")
        gridRag.ignoreLabel     = int(ignoreLabel[0])
        gridRag.baseGraphLabels = taggedView(labels,"xyz")
        gridRag.baseGraph       = gridGraph

        return gridRag




    loadGridRagHDF5.__module__ = 'vigra.graphs'
    graphs.loadGridRagHDF5 = loadGridRagHDF5

    def regionAdjacencyGraph(graph,labels,ignoreLabel=None,reserveEdges=0, maxLabel=None, isDense=None):
        """ Return a region adjacency graph for a labeld graph.

            Parameters:

                - graph  -- input graph
                - lables -- node-map with labels for each nodeSumWeights
                - ignoreLabel -- label to ingnore (default: None)
                - reserveEdges -- reverse a certain number of edges (default: 0)

            Returns:
                - rag -- instance of RegionAdjacencyGraph or GridRegionAdjacencyGraph
                    If graph is a GridGraph or a GridRegionAdjacencyGraph, a GridRegionAdjacencyGraph 
                    will be returned.
                    Otherwise a RegionAdjacencyGraph will be returned
        """
        if isinstance(graph , graphs.GridRegionAdjacencyGraph) or graphs.isGridGraph(graph):
            return GridRegionAdjacencyGraph(graph=graph, labels=labels, ignoreLabel=ignoreLabel,
                                            reserveEdges=reserveEdges, maxLabel=maxLabel, isDense=isDense)
        else:
            return RegionAdjacencyGraph(graph=graph, labels=labels, ignoreLabel=ignoreLabel,
                                        reserveEdges=reserveEdges, maxLabel=maxLabel, isDense=isDense)


    regionAdjacencyGraph.__module__ = 'vigra.graphs'
    graphs.regionAdjacencyGraph = regionAdjacencyGraph


   



    def gridRegionAdjacencyGraph(labels,ignoreLabel=None,reserveEdges=0, maxLabel=None, isDense=None):
        """ get a region adjacency graph and a grid graph from a labeling.

            This function will call 'graphs.gridGraph' and 'graphs.regionAdjacencyGraph'

            Keyword Arguments:
                - labels : label image
                - ignoreLabel : label to ingnore (default: None)
                - reserveEdges : reserve a number of edges (default: 0)
        """
        _gridGraph=graphs.gridGraph(numpy.squeeze(labels).shape)
        rag=graphs.regionAdjacencyGraph(graph=_gridGraph, labels=labels, ignoreLabel=ignoreLabel,
                                        reserveEdges=reserveEdges, maxLabel=maxLabel, isDense=isDense)
        return _gridGraph, rag

    gridRegionAdjacencyGraph.__module__ = 'vigra.graphs'
    graphs.gridRegionAdjacencyGraph = gridRegionAdjacencyGraph

_genRegionAdjacencyGraphConvenienceFunctions()
del _genRegionAdjacencyGraphConvenienceFunctions


def _genGraphSegmentationFunctions():

    def getNodeSizes(graph):
        """ get size of nodes:

            This functions will try to call 'graph.nodeSize()' .
            If this fails, a node map filled with 1.0 will be 
            returned

            Keyword Arguments:

                - graph : input graph
        """
        try:
            return graph.nodeSize()
        except:
            size = graphs.graphMap(graph,'node',dtype=numpy.float32)
            size[:]=1
            return size
    getNodeSizes.__module__ = 'vigra.graphs'
    graphs.getNodeSizes = getNodeSizes

    def getEdgeLengths(graph):
        """ get lengths/sizes of edges:

            This functions will try to call 'graph.edgeLength()' .
            If this fails, an edge map filled with 1.0 will be 
            returned

            Keyword Arguments:

                - graph : input graph
        """
        try:
            return graph.edgeLengths()
        except:
            size = graphs.graphMap(graph,'edge',dtype=numpy.float32)
            size[:]=1
            return size
    getEdgeLengths.__module__ = 'vigra.graphs'
    graphs.getEdgeLengths = getEdgeLengths


    def felzenszwalbSegmentation(graph,edgeWeights,nodeSizes=None,k=1.0,nodeNumStop=None,out=None):
        """ felzenszwalbs segmentation method

        Keyword Arguments :

            - graph : input graph

            - edgeWeights : edge weights / indicators

            - nodeSizes : size of each node (default: None) 
                If nodeSizes is None, 'getNodeSizes' will be called

            - k : free parameter in felzenszwalbs algorithms (default : 1.0)
                (todo: write better docu)

            - nodeNumStop : stop the agglomeration at a given nodeNum (default :None)
                If nodeNumStop is None, the resulting number of nodes does depends on k.


            - backgroundBias : backgroundBias (default  : None)

        """
        if nodeNumStop is None :
            nodeNumStop=-1
        if nodeSizes is None :
            nodeSizes=graphs.getNodeSizes(graph)
        return graphs._felzenszwalbSegmentation(graph=graph,edgeWeights=edgeWeights,nodeSizes=nodeSizes,
                                                k=k,nodeNumStop=nodeNumStop,out=out)


    felzenszwalbSegmentation.__module__ = 'vigra.graphs'
    graphs.felzenszwalbSegmentation = felzenszwalbSegmentation


    def edgeWeightedWatersheds(graph,edgeWeights,seeds,backgroundLabel=None,backgroundBias=None,out=None):
        """ edge weighted seeded watersheds

        Keyword Arguments :

            - graph : input graph

            - edgeWeights : evaluation weights

            - seeds : node map with seeds .
                For at least one node, seeds must be nonzero

            - backgroundLabel : a specific backgroundLabel (default : None)

            - backgroundBias : backgroundBias (default  : None)

        """
        if backgroundLabel is None and backgroundBias is None:
            return graphs._edgeWeightedWatershedsSegmentation(graph=graph,edgeWeights=edgeWeights,seeds=seeds,
                                                                out=out)
        else :
            if backgroundLabel is None or backgroundBias is None:
                raise RuntimeError("if backgroundLabel or backgroundBias is not None, the other must also be not None")
            return graphs._carvingSegmentation(graph=graph,edgeWeights=edgeWeights,seeds=seeds,
                                                backgroundLabel=backgroundLabel,backgroundBias=backgroundBias,out=out)

    edgeWeightedWatersheds.__module__ = 'vigra.graphs'
    graphs.edgeWeightedWatersheds = edgeWeightedWatersheds

    def nodeWeightedWatershedsSeeds(graph,nodeWeights,out=None):
        """ generate watersheds seeds

        Keyword Arguments :

            - graph : input graph

            - nodeWeights : node height map

            - out : seed map

        """
        return graphs._nodeWeightedWatershedsSeeds(graph=graph,nodeWeights=nodeWeights,out=out)

    nodeWeightedWatershedsSeeds.__module__ = 'vigra.graphs'
    graphs.nodeWeightedWatershedsSeeds = nodeWeightedWatershedsSeeds


    def shortestPathSegmentation(graph, edgeWeights, nodeWeights, seeds=None, out=None):
        """ node weighted seeded watersheds

        Keyword Arguments :

            - graph : input graph

            - edgeWeights : edge weight map 

            - nodeWeights : node weight map

            - seeds : node map with seeds (default: None)
                If seeds are None, 'nodeWeightedWatershedsSeeds' will be called

        """

        if seeds  is None:
            seeds = graphs.nodeWeightedWatershedsSeeds(graph=graph,nodeWeights=nodeWeights)
        return graphs._shortestPathSegmentation(graph=graph, edgeWeights=edgeWeights, nodeWeights=nodeWeights,
                                                seeds=seeds, out=out)

    shortestPathSegmentation.__module__ = 'vigra.graphs'
    graphs.shortestPathSegmentation = shortestPathSegmentation

    def nodeWeightedWatersheds(graph,nodeWeights,seeds=None,method='regionGrowing',out=None):
        """ node weighted seeded watersheds

        Keyword Arguments :

            - graph : input graph

            - nodeWeights : node height map / evaluation weights

            - seeds : node map with seeds (default: None)
                If seeds are None, 'nodeWeightedWatershedsSeeds' will be called

        """

        if seeds  is None:
            seeds = graphs.nodeWeightedWatershedsSeeds(graph=graph,nodeWeights=nodeWeights)
        if method!='regionGrowing':
            raise RuntimeError("currently only 'regionGrowing' is supported")
        return graphs._nodeWeightedWatershedsSegmentation(graph=graph,nodeWeights=nodeWeights,seeds=seeds,method=method,out=out)

    nodeWeightedWatersheds.__module__ = 'vigra.graphs'
    graphs.nodeWeightedWatersheds = nodeWeightedWatersheds

    def agglomerativeClustering(graph,edgeWeights=None,edgeLengths=None,nodeFeatures=None,nodeSizes=None,
            nodeNumStop=None,beta=0.5,metric='l1',wardness=1.0,out=None):
        """ agglomerative hierarchicalClustering
        Keyword Arguments :

            - graph : input graph

            - edgeWeights : edge weights / indicators (default : None)

            - edgeLengths : length  / weight of each edge (default : None)
                Since we do weighted mean agglomeration, a length/weight
                is needed for each edge to merge 2 edges w.r.t. weighted mean.
                If no edgeLengths is given, 'getEdgeLengths' is called.

            - nodeFeatures : a feature vector for each node (default: None)
                A feature vector as RGB values,or a histogram for each node.
                Within the agglomeration, an additional edge weight will be
                computed from the "difference" between the features of two adjacent nodes.
                The metric specified in the keyword 'metric' is used to compute this
                difference

            - nodeSizes : size  / weight of each node (default : None)
                Since we do weighted mean agglomeration, a size / weight
                is needed for each node to merge 2 edges w.r.t. weighted mean.
                If no nodeSizes is given, 'getNodeSizes' is called.

            - nodeNumStop : stop the agglomeration at a given nodeNum (default : graph.nodeNum/2)

            - beta : weight between edgeWeights and nodeFeatures based edgeWeights (default:0.5) :
                    0.0 means only edgeWeights (from keyword edge weights) and 1.0 means only edgeWeights 
                    from nodeFeatures differences

            - metric : metric used to compute node feature difference (default : 'l1')

            - wardness : 0 means do not apply wards critrion, 1.0 means fully apply wards critrion (default : 1.0)

            - out : preallocated nodeMap for the resulting labeling (default : None)

        Returns:

            A node labele map encoding the segmentation

        """

        assert edgeWeights is not None or nodeFeatures is not None

        print "prepare "

        if nodeNumStop is None:
            nodeNumStop = max(graph.nodeNum/2,min(graph.nodeNum,2))

        
        if edgeLengths is None :
            print "get edge length"
            edgeLengths = graphs.getEdgeLengths(graph)

        
        if nodeSizes is None:
            print "get node size"
            nodeSizes = graphs.getNodeSizes(graph)

        
        if edgeWeights is None :
            print "get wegihts length"
            edgeWeights = graphs.graphMap(graph,'edge')
            edgeWeights[:]=0

        if nodeFeatures is None :
            print "get node feat"
            nodeFeatures = graphs.graphMap(graph,'node',addChannelDim=True)
            nodeFeatures[:]=0


        #import sys
        #print "graph refcout", sys.getrefcount(graph)
        mg = graphs.mergeGraph(graph)
        #print "graph refcout", sys.getrefcount(graph)
        #mg = []
        #del mg
        #import gc
        #gc.collect()

        #print "graph refcout", sys.getrefcount(graph)
        #sys.exit(0)



        clusterOp = graphs.minEdgeWeightNodeDist(mg,edgeWeights=edgeWeights,edgeLengths=edgeLengths,
                                                    nodeFeatures=nodeFeatures,nodeSizes=nodeSizes,
                                                    beta=float(beta),metric=metric,wardness=wardness)

        

        hc = graphs.hierarchicalClustering(clusterOp, nodeNumStopCond=nodeNumStop,
                                           buildMergeTreeEncoding=False)

        hc.cluster()
        labels = hc.resultLabels(out=out)
        #del hc
        #del clusterOp
        #del mg
        return labels


    agglomerativeClustering.__module__ = 'vigra.graphs'
    graphs.agglomerativeClustering = agglomerativeClustering



    def minEdgeWeightNodeDist(mergeGraph,edgeWeights=None,edgeLengths=None,nodeFeatures=None,nodeSizes=None,outWeight=None,
        beta=0.5,metric='squaredNorm',wardness=1.0):
            graph=mergeGraph.graph()
            assert edgeWeights is not None or nodeFeatures is not None

            if edgeLengths is None :
                edgeLengths = graphs.getEdgeLengths(graph,addChannelDim=True)
            if nodeSizes is None:
                nodeSizes = graphs.getNodeSizes(graph,addChannelDim=True)

            if edgeWeights is None :
                edgeWeights = graphs.graphMap(graph,'edge',addChannelDim=True)
                edgeWeights[:]=0

            if nodeFeatures is None :
                nodeFeatures = graphs.graphMap(graph,'node',addChannelDim=True)
                nodeFeatures[:]=0

            if outWeight is None:
                outWeight=graphs.graphMap(graph,item='edge',dtype=numpy.float32)


            if  metric=='squaredNorm':
                nd=graphs.MetricType.squaredNorm
            elif  metric=='norm':
                nd=graphs.MetricType.norm
            elif  metric=='chiSquared':
                nd=graphs.MetricType.chiSquared
            elif metric in ('l1','manhattan'):
                nd=graphs.MetricType.manhattan
            elif isinstance(metric,graphs.MetricType):
                nd=metric
            else :
                raise RuntimeError("'%s' is not a supported distance type"%str(metric))

            # call unsave c++ function and make it sav
            op = graphs.__minEdgeWeightNodeDistOperator(mergeGraph,edgeWeights,edgeLengths,nodeFeatures,nodeSizes,outWeight,
                float(beta),nd,float(wardness))


            op.__base_object__=mergeGraph
            op.__outWeightArray__=outWeight
            op.edgeLengths=edgeLengths
            op.nodeSizes=nodeSizes
            op.edgeWeights=edgeWeights
            op.nodeFeatures=nodeFeatures
            return op


    minEdgeWeightNodeDist.__module__ = 'vigra.graphs'
    graphs.minEdgeWeightNodeDist = minEdgeWeightNodeDist    



    def pythonClusterOperator(mergeGraph,operator,useMergeNodeCallback=True,useMergeEdgesCallback=True,useEraseEdgeCallback=True):
      #call unsave function and make it save
      op = graphs.__pythonClusterOperator(mergeGraph,operator,useMergeNodeCallback,useMergeEdgesCallback,useEraseEdgeCallback)
      #op.__dict__['__base_object__']=mergeGraph
      #op.__base_object__=mergeGraph
      return op

    pythonClusterOperator.__module__ = 'vigra.graphs'
    graphs.pythonClusterOperator = pythonClusterOperator    

    def hierarchicalClustering(clusterOperator,nodeNumStopCond,buildMergeTreeEncoding=True):
        # call unsave c++ function and make it save
        hc = graphs.__hierarchicalClustering(clusterOperator,long(nodeNumStopCond),bool(buildMergeTreeEncoding))
        #hc.__dict__['__base_object__']=clusterOperator
        hc.__base_object__ = clusterOperator
        return hc

    hierarchicalClustering.__module__ = 'vigra.graphs'
    graphs.hierarchicalClustering = hierarchicalClustering

_genGraphSegmentationFunctions()
del _genGraphSegmentationFunctions



def _genHistogram():
    def gaussianHistogram(image,minVals,maxVals,bins=30,
                     sigma=3.0,sigmaBin=2.0,out=None):
        """ 
        """
        spatialDim  = image.ndim - 1
        out = histogram.gaussianHistogram_(image=image, minVals=minVals, maxVals=maxVals,
                                           bins=bins, sigma=sigma, sigmaBin=sigmaBin, 
                                           out=out)

        out = out.reshape(image.shape[0:spatialDim]+(-1,))
        if spatialDim == 2:
            out /= numpy.sum(out,axis=spatialDim)[:,:, numpy.newaxis]
        elif spatialDim == 3:
            out /= numpy.sum(out,axis=spatialDim)[:,:,:, numpy.newaxis]
        elif spatialDim == 4:
            out /= numpy.sum(out,axis=spatialDim)[:,:,:, :,numpy.newaxis]
        return out

    gaussianHistogram.__module__ = 'vigra.histogram'
    histogram.gaussianHistogram = gaussianHistogram


_genHistogram()
del _genHistogram


def _genGraphSmoothingFunctions():
    def recursiveGraphSmoothing( graph,nodeFeatures,edgeIndicator,gamma,
                               edgeThreshold,scale=1.0,iterations=1,out=None):
        """ recursive graph smoothing to smooth node features.
            Each node feature is smoothed with the features of neighbor nodes.
            The strength of the smoothing is computed from:

                "edgeIndicator > edgeThreshold ? 0 : exp(-1.0*gamma*edgeIndicator)*scale"

            Therefore this filter is edge preserving.

            Keyword Arguments :

                - graph : input graph

                - nodeFeatures : node features which should be smoothed

                - edgeIndicator  : edge indicator 

                - gamma  : scale edgeIndicator by gamma bevore taking the negative exponent

                - scale  : how much should a node be mixed with its neighbours per iteration 

                - iteration : how often should recursiveGraphSmoothing be called recursively

            Returns : 
                smoothed nodeFeatures

        """
        return graphs._recursiveGraphSmoothing(graph=graph,nodeFeatures=nodeFeatures,edgeIndicator=edgeIndicator,
                              gamma=gamma,edgeThreshold=edgeThreshold,scale=scale,iterations=iterations,out=out)

    recursiveGraphSmoothing.__module__ = 'vigra.graphs'
    graphs.recursiveGraphSmoothing = recursiveGraphSmoothing

_genGraphSmoothingFunctions()
del _genGraphSmoothingFunctions




def _genGraphMiscFunctions():
    
    def nodeFeaturesToEdgeWeights(graph,nodeFeatures,metric='l1',out=None):
        """ compute an edge indicator from node features .

            Keyword Arguments :
                - graph : input graph
                - nodeFeatures : node map with feature vector for each node
                - metric : metric / distance used to convert 2 node features to
                    an edge weight

            Returns :
                edge indicator
        """
        return graphs._nodeFeatureDistToEdgeWeight(graph=graph,nodeFeatures=nodeFeatures,metric=metric,out=out)

    nodeFeaturesToEdgeWeights.__module__ = 'vigra.graphs'
    graphs.nodeFeaturesToEdgeWeights = nodeFeaturesToEdgeWeights

    def eccentricityTransform(labels, out=None):
        """ Compute eccentricity transform on labeled image.

            Keyword Arguments :
                - labels : input image (labeled)

            Returns :
                eccentricity transform
        """
        return graphs._eccentricityTransform(labels=labels, out=out)

    eccentricityTransform.__module__ = 'vigra.graphs'
    graphs.eccentricityTransform = eccentricityTransform

    def eccentricityCenters(labels, out=None):
        """ Compute eccentricity centers on labeled image.
            Output is of shape (m, N), where m is the largest label in labels and N is the dimension.
            Since labels start with 1, column i of the output contains the center of the region with label i+1.

            Keyword Arguments :
                - labels : input image (labeled)

            Returns :
                centers of the labeled regions
        """
        return graphs._eccentricityCenters(labels=labels, out=out)

    eccentricityCenters.__module__ = 'vigra.graphs'
    graphs.eccentricityCenters = eccentricityCenters

    def eccentricityTransformWithCenters(labels, ecc=None, centers=None):
        """ Compute eccentricity transform and centers on labeled image.
            Shape of centers is (m, N), where m is the largest label in labels and N is the dimension.
            Since labels start with 1, column i of the centers contains the center of the region with label i+1.

            Keyword Arguments :
                - labels : input image (labeled)

            Returns :
                2-tuple with eccentricity transform and centers
        """
        return graphs._eccentricityTransformWithCenters(labels=labels, ecc=ecc, centers=centers)

    eccentricityTransformWithCenters.__module__ = 'vigra.graphs'
    graphs.eccentricityTransformWithCenters = eccentricityTransformWithCenters

_genGraphMiscFunctions()
del _genGraphMiscFunctions




def loadBSDGt(filename):
    import scipy.io as sio
    matContents = sio.loadmat(filename)
    ngt = len(matContents['groundTruth'][0])
    gts = []
    for gti in range(ngt):
        gt =  matContents['groundTruth'][0][gti][0]['Segmentation'][0]
        gt = numpy.swapaxes(gt,0,1)
        gt = gt.astype(numpy.uint32)
        print gt.min(),gt.max()
        gts.append(gt[:,:,None])
    gtArray = numpy.concatenate(gts,axis=2)
    print gtArray.shape
    return gtArray
