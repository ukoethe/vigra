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


def segShow(img,labels,edgeColor=(1,0,0) ):
    labels = numpy.squeeze(labels)
    crackedEdges = analysis.regionImageToCrackEdgeImage(labels)
    whereEdge    =  numpy.where(crackedEdges==0)
    imgToDisplay = resize(img,numpy.squeeze(crackedEdges).shape)
    imgToDisplay-=imgToDisplay.min()
    imgToDisplay/=imgToDisplay.max()
    for c in range(3):
        ic = imgToDisplay[:,:,c]
        ic[whereEdge]=edgeColor[c]

    imshow(imgToDisplay)

def nestedSegShow(img,labels,edgeColors,scale=1):

    shape=(labels.shape[0]*scale,labels.shape[1]*scale)
    if scale!=1:
        img=vigra.resize(img,shape)
    



    assert numpy.squeeze(labels).ndim==3
    nSegs  = labels.shape[2]

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

        edgeColor = edgeColors[si][0:3]
        if len(edgeColors[si])<4:
            alpha = 0.0
        else:
            alpha = edgeColors[si][3]
        for c in range(3):
            icI = imgIn[:,:,c]
            ic  = imgToDisplay[:,:,c]
            ic[whereEdge]=(1.0-alpha) * edgeColors[si][c] + alpha*icI[whereEdge]

    imshow(imgToDisplay)


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







def _genGraphConvenienceFunctions():

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
                return self.intrinsicNodeMapShape()

            def nodeSize(self):
                size = graphs.graphMap(self,item='node',dtype=numpy.float32)
                size[:]=1
                return size

            def edgeLength(self):
                size = graphs.graphMap(self,item='edge',dtype=numpy.float32)
                size[:]=1
                return size



    def isGridGraph(obj):
        return isinstance(obj,(graphs.GridGraphUndirected2d , graphs.GridGraphUndirected3d))
    isGridGraph.__module__ = 'vigra.graphs'
    graphs.isGridGraph = isGridGraph  

    def listGraph(nodes=0,edges=0):
        ''' Return an empty directed graph

            Parameters :

                - nodes -- number of nodes to reserveEdges
                - edges -- number of edges to reserve

            Returns :

                - graph
        '''
        return graphs.AdjacencyListGraph(nodes,edges)
        
    listGraph.__module__ = 'vigra.graphs'
    graphs.listGraph = listGraph


    class RegionAdjacencyGraph(graphs.AdjacencyListGraph):
        def __init__(self,graph,labels,ignoreLabel=None,reserveEdges=0):
            super(RegionAdjacencyGraph,self).__init__(long(labels.max()+1),long(reserveEdges))

            if ignoreLabel is None:
                ignoreLabel=-1

            self.labels          = labels
            self.ignoreLabel     = ignoreLabel
            self.baseGraphLabels = labels
            self.baseGraph       = graph
            # set up rag
            self.affiliatedEdges = graphs._regionAdjacencyGraph(graph,labels,self,self.ignoreLabel)   

        def accumulateEdgeFeatures(self,edgeFeatures,acc='mean',out=None):
            graph = self.baseGraph
            affiliatedEdges = self.affiliatedEdges
            return graphs._ragEdgeFeatures(self,graph,affiliatedEdges,edgeFeatures,acc,out)

        def accumulateNodeFeatures(self,nodeFeatures,acc='mean',out=None):
            graph = self.baseGraph
            labels = self.baseGraphLabels
            ignoreLabel = self.ignoreLabel
            return graphs._ragNodeFeatures(self,graph,labels,nodeFeatures,acc,ignoreLabel,out)

        def projectNodeFeatureToBaseGraph(self,features,out=None):
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

        def projectLabelsBack(self,steps,labels=None,current=0):
            if labels is None :
                # identity segmentation on this level
                labels = self.nodeIdMap()

            if steps == current :
                return labels
            else :
                labels = self.projectLabelsToBaseGraph(labels)
                return self.baseGraph.projectLabelsBack(steps,labels,current+1)


        def projectLabelsToBaseGraph(self,labels=None):
            if labels is None :
                # identity segmentation on this level
                labels = self.nodeIdMap()
            return self.projectNodeFeatureToBaseGraph(features=labels)

 



    RegionAdjacencyGraph.__module__ = 'vigra.graphs'
    graphs.RegionAdjacencyGraph = RegionAdjacencyGraph

    class GridRegionAdjacencyGraph(graphs.RegionAdjacencyGraph):
        def __init__(self,graph,labels,ignoreLabel=None,reserveEdges=0):
            if not (isGridGraph(graph) or  isinstance(graph,GridRegionAdjacencyGraph)):
                raise RuntimeError("graph must be a GridGraph or a GridRegionAdjacencyGraph")
            super(GridRegionAdjacencyGraph,self).__init__(graph,labels,ignoreLabel,reserveEdges)

        @property
        def shape(self):
            return self.baseGraph.shape

        def projectLabelsToGridGraph(self,labels=None):
            if labels is None :
                # identity segmentation on this level
                labels = self.nodeIdMap()

            if isGridGraph(self.baseGraph):
                return self.projectLabelsToBaseGraph(labels)
            else :
                labels = self.projectLabelsToBaseGraph(labels)
                return self.baseGraph.projectLabelsToGridGraph(labels)


        def show(self,img,labels=None):
            pLabels = self.projectLabelsToGridGraph(labels)
            segShow(img,numpy.squeeze(pLabels))

        def nodeSize(self):
            baseNodeSizes = self.baseGraph.nodeSize()
            return self.accumulateNodeFeatures(baseNodeSizes,acc='sum')
        def edgeLength(self):
            baseNodeSizes = self.baseGraph.edgeLength()
            return self.accumulateEdgeFeatures(baseNodeSizes,acc='sum')

    GridRegionAdjacencyGraph.__module__ = 'vigra.graphs'
    graphs.GridRegionAdjacencyGraph = GridRegionAdjacencyGraph



    def regionAdjacencyGraph(graph,labels,ignoreLabel=None,reserveEdges=0):
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
        if isinstance(graph , graphs.GridRegionAdjacencyGraph) or isGridGraph(graph):
            return GridRegionAdjacencyGraph(graph=graph,labels=labels,ignoreLabel=ignoreLabel,reserveEdges=reserveEdges)
        else:
            return RegionAdjacencyGraph(graph=graph,labels=labels,ignoreLabel=ignoreLabel,reserveEdges=reserveEdges) 

       



    regionAdjacencyGraph.__module__ = 'vigra.graphs'
    graphs.regionAdjacencyGraph = regionAdjacencyGraph

    def gridRegionAdjacencyGraph(labels,ignoreLabel=None,reserveEdges=0):
      _gridGraph=graphs.gridGraph(numpy.squeeze(labels).shape)
      rag=graphs.regionAdjacencyGraph(_gridGraph,labels,ignoreLabel,reserveEdges)
      return _gridGraph,rag

    gridRegionAdjacencyGraph.__module__ = 'vigra.graphs'
    graphs.gridRegionAdjacencyGraph = gridRegionAdjacencyGraph


    def intrinsicGraphMapShape(graph,item):
        """ Intrinsic shape of node/edge/arc-map for a given graph.

            Node edge and arc maps are stored in numpy arrays by default.
            The instric shape may not be confused with the number
            of nodes/edges/arcs. The instric shape is used to
            allocate a numpy are which can store data for nodes/arcs/edgeSizes
            of a given graph.

            Parameters:

                - graph -- input graph to get the shape for

                - item  -- item must be ``'node'`` , ``'edge'`` or ``'arc'``

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


    def graphMap(graph,item,dtype,channels=1):
        """ Return a graph map for a given graph item (``'node'`` , ``'edge'`` or ``'arc'``).

            Parameters:

                - graph    -- graph to get a graph map for
                - item     -- ``'node'`` , ``'edge'`` or ``'arc'``
                - dtype    -- desired dtype 
                - channels -- number of channels (default: 1)

            Returns:

                - graphmap as numpy.ndarray
        """
        s = intrinsicGraphMapShape(graph,item)
        intrDim = len(s)
        if(channels==1):
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




    def minEdgeWeightNodeDist(mergeGraph,edgeWeights,edgeSizes,nodeFeatures,nodeSize,outWeight=None,
        beta=0.5,nodeDistType='squaredNorm',wardness=1.0):
            # call unsave c++ function and make it save

            if outWeight is None:
                outWeight=graphs.graphMap(mergeGraph,item='edge',dtype=numpy.float32)


            if    nodeDistType=='squaredNorm':
                nd=0
            elif  nodeDistType=='norm':
                nd=1
            elif  nodeDistType=='chiSquared':
                nd=2
            else :
                raise RuntimeError("'%s' is not a supported distance type"%str(nodeDistType))

            op = graphs.__minEdgeWeightNodeDistOperator(mergeGraph,edgeWeights,edgeSizes,nodeFeatures,nodeSize,outWeight,
                float(beta),long(nd),float(wardness))
            op.__dict__['__base_object__']=mergeGraph
            op.__dict__['__outWeightArray__']=outWeight
            return op

    minEdgeWeightNodeDist.__module__ = 'vigra.graphs'
    graphs.minEdgeWeightNodeDist = minEdgeWeightNodeDist    

    def pythonClusterOperator(mergeGraph,opertator,useMergeNodeCallback=True,useMergeEdgesCallback=True,useEraseEdgeCallback=True):
      #call unsave function and make it save
      op = graphs.__pythonClusterOperator(mergeGraph,opertator,useMergeNodeCallback,useMergeEdgesCallback,useEraseEdgeCallback)
      op.__dict__['__base_object__']=mergeGraph
      return op

    pythonClusterOperator.__module__ = 'vigra.graphs'
    graphs.pythonClusterOperator = pythonClusterOperator    

    def hierarchicalClustering(clusterOperator,nodeNumStopCond):
        # call unsave c++ function and make it save
        hc = graphs.__hierarchicalClustering(clusterOperator,long(nodeNumStopCond))
        hc.__dict__['__base_object__']=clusterOperator
        return hc

    hierarchicalClustering.__module__ = 'vigra.graphs'
    graphs.hierarchicalClustering = hierarchicalClustering

    def mergeGraph(graph):
        # call unsave c++ function and make it save
        mg = graphs.__mergeGraph(graph)
        mg.__dict__['__base_object__']=graph
        return mg

    mergeGraph.__module__ = 'vigra.graphs'
    graphs.mergeGraph = mergeGraph



    def hierarchicalSuperpixels(labels,edgeIndicatorImage,nodeFeaturesImage,nSuperpixels,
        beta=0.5,nodeDistType='squaredNorm',wardness=0.0,verbose=False):

        shape  = numpy.squeeze(labels).shape
        eShape = numpy.squeeze(edgeIndicatorImage).shape
        tShape = tuple( [ 2*s-1 for s in shape  ]  )
        if(eShape==shape):
            eImg = sampling.resize(numpy.squeeze(edgeIndicatorImage),tShape)
        elif(eShape==tShape):
            eImg = numpy.squeeze(edgeIndicatorImage)
        else :
            raise RuntimeError("edge indactor image has wrong shape")

        


        dimension = len(labels.shape)
        assert dimension == 2 or dimension ==3
        if nodeFeaturesImage.ndim == dimension:
            nodeFeatureChannels=1
        else:
            assert nodeFeaturesImage.ndim == dimension+1
            nodeFeatureChannels = nodeFeaturesImage.shape[dimension]
        

        #if verbose : print "gridGraph"
        gridGraph = graphs.gridGraph(labels.shape)
        
        eIndicator =  graphs.edgeFeaturesFromInterpolatedImage(gridGraph,eImg)

        rag = graphs.regionAdjacencyGraph(graph=gridGraph,labels=labels,ignoreLabel=0)
        if verbose :print "regionAdjacencyGraph",rag
        hyperEdgeSizes = rag.accumulateEdgeSize()
        hyperNodeSizes = rag.accumulateEdgeSize()


        edgeMinWeight  = graphs.graphMap(graph=rag,item='edge',dtype=numpy.float32,channels=1)

        print edgeIndicatorImage.shape,edgeIndicatorImage.dtype

        edgeIndicator  = rag.accumulateEdgeFeatures(numpy.squeeze(eIndicator),acc='mean')
        nodeFeatures   = rag.accumulateNodeFeatures(nodeFeaturesImage,acc='mean')


        mergeGraph = graphs.mergeGraph(rag)
        clusterOperator = graphs.minEdgeWeightNodeDist(mergeGraph,edgeIndicator,hyperEdgeSizes,
            nodeFeatures,hyperNodeSizes,edgeMinWeight,
            beta=beta,nodeDistType=nodeDistType,wardness=wardness)
        hc  = graphs.hierarchicalClustering(clusterOperator,nodeNumStopCond=nSuperpixels,)
        hc.cluster()

        newLabels = labels.copy()
        newLabels = newLabels.reshape(-1)
        # this is inplace
        hc.reprNodeIds(newLabels)
        newLabels = newLabels.reshape(labels.shape)

        if(labels.ndim==2):
            newLabels = analysis.labelImage(newLabels)
        if(labels.ndim==3):
            newLabels = analysis.labelVolume(newLabels)
        return newLabels

    hierarchicalSuperpixels.__module__ = 'vigra.graphs'
    graphs.hierarchicalSuperpixels = hierarchicalSuperpixels

    INVALID = graphs.Invalid()
    #hierarchicalSuperpixels.__module__ = 'vigra.graphs'
    graphs.INVALID = INVALID


    def nodeIdsLabels(graph,nodeIds,labels,out=None):
        nodeIdsShape = nodeIds.shape
        if out is not None:
            out=out.reshape(-1)
        out = graphs._nodeIdsLabels(graph=graph,nodeIds=nodeIds.reshape(-1) ,labels=labels,out=out)
        out = out.reshape(nodeIdsShape)
        return out

    nodeIdsLabels.__module__ = 'vigra.graphs'
    graphs.nodeIdsLabels = nodeIdsLabels


    def nodeIdsFeatures(graph,nodeIds,features,out=None):
        nodeMapShape = graph.intrinsicNodeMapShape()
        featureMapShape       = features.shape
        if len(nodeMapShape)+1==features.ndim:
            spatialInputShape = numpy.squeeze(nodeIds).shape
            numberOfChannels  = featureMapShape[-1]
            outShape = spatialInputShape + (numberOfChannels,)

            if out is not None:
                out = out.reshape([-1,numberOfChannels])
                out = taggedView(out,'xc')

            out = graphs._nodeIdsFeatures(graph=graph,nodeIds=nodeIds.reshape(-1) ,features=features,out=out)
            out = out.reshape(outShape)
            return out
        elif len(nodeMapShape)==features.ndim:
            raise RuntimeError("feature map has wrong dimension: feature map must have a channel dimension")
        else :
            raise RuntimeError("feature map has wrong dimension")
    
    nodeIdsFeatures.__module__ = 'vigra.graphs'
    graphs.nodeIdsFeatures = nodeIdsFeatures




    class ShortestPathPathDijkstra(object):
        def __init__(self,graph):
            self.pathFinder =  graphs._shortestPathDijkstra(graph)
            self.graph=graph
            self.source = None
            self.target = None
        def run(self,weights,source,target=None):
            self.source = source
            self.target = target
            if target is None:
                self.pathFinder.run(weights,source)
            else:
                self.pathFinder.run(weights,source,target)
            return self

        def path(self,target=None,pathType='coordinates'):
            if target is None:
                assert self.target is not None
                target=self.target

            if pathType=='coordinates':
                return self.pathFinder.nodeCoordinatePath(target)
            elif pathType == 'ids':
                return self.pathFinder.nodeIdPath(target)
        def distance(self,target):
          return self.pathFinder.distance(target)
        def distances(self,out=None):
            self.pathFinder.distances(out)
        def predecessors(self,out=None):
            self.pathFinder.predecessors(out)


    ShortestPathPathDijkstra.__module__ = 'vigra.graphs'
    graphs.ShortestPathPathDijkstra = ShortestPathPathDijkstra



_genGraphConvenienceFunctions()
del _genGraphConvenienceFunctions

