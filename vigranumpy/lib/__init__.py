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

import sys, os, time, math
from numbers import Number
from multiprocessing import cpu_count
try:
    import pylab
except Exception, e:
    pass


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
def imshow(image,show=True, **kwargs):
    '''Display a scalar or RGB image by means of matplotlib.
       If the image does not have one or three channels, an exception is raised.
       The image will be automatically scaled to the range 0...255 when its dtype
       is not already 'uint8' and neither 'cmap' nor 'norm' are specified in kwargs
    '''
    import matplotlib.pylab

    if not hasattr(image, 'axistags'):
        plot = matplotlib.pyplot.imshow(image, **kwargs)
        if show:
            matplotlib.pylab.show()
        return plot

    image = image.transposeToNumpyOrder()
    if image.channels == 1:
        image = image.dropChannelAxis().view(numpy.ndarray)
        if 'cmap' in kwargs.keys():
            cmap = kwargs['cmap']
        elif 'norm' in kwargs.keys():
            norm = kwargs['norm']
        else:
            cmap = matplotlib.cm.gray
            norm = matplotlib.cm.colors.Normalize()
        plot = matplotlib.pyplot.imshow(image, cmap=cmap, norm=norm, **kwargs)
        if show:
            matplotlib.pylab.show()
        return plot
    elif image.channels == 3:
        if image.dtype != numpy.uint8:
            out = image.__class__(image.shape, dtype=numpy.uint8, axistags=image.axistags)
            image = colors.linearRangeMapping(image, newRange=(0.0, 255.0), out=out)
        plot = matplotlib.pyplot.imshow(image.view(numpy.ndarray), **kwargs)
        if show:
            matplotlib.pylab.show()
        return plot
    else:
        raise RuntimeError("vigra.imshow(): Image must have 1 or 3 channels.")


def multiImshow(images,shape, show=True):
    nImg = len(images)
    f = pylab.figure()

    s = tuple(shape)
    for c,iname in enumerate(images.keys()):
        data,itype = images[iname]
        if itype == 'img':

            ax1 = f.add_subplot(s[0],s[1],c+1)
            imshow(data,show=False)
            ax1.set_title(iname)
            pylab.axis('off')
    if show :
        pylab.show()

def segShow(img,labels,edgeColor=(0,0,0),alpha=0.3,show=False,returnImg=False,r=0):

    img = numpy.squeeze(img)
    if img.ndim ==2:
        img = numpy.concatenate( [ img[:,:,None]]*3 ,axis=2).astype(numpy.float32)
        img = taggedView(img, 'xyc')

    labels = numpy.squeeze(labels)
    crackedEdges = analysis.regionImageToCrackEdgeImage(labels+1).squeeze()
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





# define tensor convenience functions
def _genDistanceTransformFunctions():

    def distanceTransform(array,background=True,norm=2,pixel_pitch=None, out=None):
        if array.squeeze().ndim == 2:
            return filters.distanceTransform2D(array,background=background,norm=norm,
                                               pixel_pitch=pixel_pitch, out=out)
        elif array.squeeze().ndim == 3:
            return filters.distanceTransform3D(array.astype('float32'),background=background,norm=2)
        else:
            raise RuntimeError("distanceTransform is only implemented for 2D and 3D arrays")

    distanceTransform.__module__ = 'vigra.filters'
    filters.distanceTransform = distanceTransform



_genDistanceTransformFunctions()
del _genDistanceTransformFunctions




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
                    #print "Weights",weights
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


        def showEdgeFeature(self, img, edgeFeature, cmap='jet', returnImg=False, labelMode=False):
            import matplotlib
            assert graphs.isGridGraph(self.baseGraph)
            imgOut = img.copy().squeeze()
            if imgOut.ndim == 2:
                imgOut = numpy.concatenate([imgOut[:,:,None]]*3,axis=2)
            imgOut = taggedView(imgOut,'xyc')
            imgOut-=imgOut.min()
            imgOut/=imgOut.max()

            if not labelMode:
                edgeFeatureShow = edgeFeature.copy()
                mi = edgeFeatureShow.min()
                ma = edgeFeatureShow.max()
                cm = matplotlib.cm.ScalarMappable(cmap=cmap)
                rgb = cm.to_rgba(edgeFeatureShow)[:,0:3]
                print rgb.shape

                if(ma > mi):
                    edgeFeatureShow -=mi
                    edgeFeatureShow /= edgeFeatureShow.max()
                else:
                    edgeFeatureShow[:] = 1

            for e in self.edgeIter():

                u,v = self.edgeUVCoordinates(e.id)

                if not labelMode:
                    showVal = rgb[e.id,:]
                else:
                    if edgeFeature[e.id] == 0:
                        showVal=[0,0,1]
                    elif edgeFeature[e.id] == 1:
                        showVal=[0,1,0]
                    elif edgeFeature[e.id] == -1:
                        showVal=[1,0,0]

                imgOut[u[:,0],u[:,1],:] = showVal
                imgOut[v[:,0],v[:,1],:] = showVal
                #print u.shape
            if returnImg:
                return imgOut
            imshow(imgOut)




        def nodeSize(self):
            """ get the geometric size of the nodes """
            if graphs.isGridGraph(self.baseGraph):
                return graphs._ragNodeSize(self, self.baseGraph, self.labels, self.ignoreLabel)
            else:
                baseNodeSizes = self.baseGraph.nodeSize()
                return self.accumulateNodeFeatures(baseNodeSizes,acc='sum')
        def edgeLengths(self):
            """ get the geometric length of the edges"""
            if graphs.isGridGraph(self.baseGraph):
                return graphs._ragEdgeSize(self,self.affiliatedEdges)
            else:
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


    class TinyEdgeLabelGui(object):
        def __init__(self, rag, img, edgeLabels = None, labelMode=True):

            if labelMode and isinstance(edgeLabels, numpy.ndarray):
                assert set(numpy.unique(edgeLabels)).issubset({-1, 0, 1}), 'if labelMode is true only label values of [-1, 0, 1] are permitted'

            self.press = None
            self.rag = rag
            self.img = img
            self.edgeLabels = edgeLabels
            self.dim = len(img.shape)
            self.zOffset = 0
            self.edgeRag2dToRag = None
            self.edgeRagToRag2d = None
            if self.dim == 3:
                self.zOffset = self.img.shape[2]/2

            self.visuImg = numpy.array(img, dtype=numpy.float32)
            self.visuImg -= self.visuImg.min()
            self.visuImg /= self.visuImg.max()

            self.rag2d = None
            self.visuImg2d = None

            self.labelMode = labelMode

            if self.edgeLabels is None :
                self.edgeLabels = numpy.zeros(self.rag.edgeNum, dtype=numpy.float32)
            self.edgeLabels2d = None

            self.slice2d()

            self.implot = None
            self.currentLabel  = 1

            self.brushSize = 1


        def startGui(self):
            from functools import partial
            import pylab as plt
            from matplotlib.widgets import Slider, Button, RadioButtons


            ax = plt.gca()
            fig = plt.gcf()

            imgWithEdges =self.rag2d.showEdgeFeature(self.visuImg2d, self.edgeLabels2d, returnImg=True, labelMode=self.labelMode)
            self.implot = ax.imshow(numpy.swapaxes(imgWithEdges,0,1))

            ff = partial(self.onclick, self)

            cid = fig.canvas.mpl_connect('button_press_event', self.onclick)

            fig.canvas.mpl_connect('key_press_event', self.press_event)

            fig.canvas.mpl_connect('scroll_event', self.scroll)

            fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
            fig.canvas.mpl_connect('button_release_event', self.on_release)

            if self.labelMode:
                axcolor = 'lightgoldenrodyellow'
                axamp  = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
                self.slideBrush = Slider(axamp, 'brush-size', 1, 20.0, valinit=2)

                self.slideBrush.on_changed(self.updateBrushSize)

            plt.show()


        def updateBrushSize(self, val):
            self.brushSize = int(val+0.5)



        def press_event(self, event):
            sys.stdout.flush()
            if event.key=='0' or event.key=='3':
                self.currentLabel = 0
            if event.key=='1':
                self.currentLabel = 1
            if event.key=='2':
                self.currentLabel = -1


        def slice2d(self):
            if self.dim==3:
                labels = self.rag.labels[:,:,self.zOffset].squeeze()
                gg = graphs.gridGraph(labels.shape)
                self.rag2d = graphs.regionAdjacencyGraph(gg, labels)
                # update edges 2d:
                self.edgeLabels2d = numpy.zeros(self.rag2d.edgeNum, dtype=numpy.float32)

                # update edge correlation
                self.edgeIdRag2dToRag = dict()
                self.edgeIdRagToRag2d = dict()
                for edge in self.rag2d.edgeIter():
                    edge3d = self.rag.findEdge(edge.u, edge.v)
                    self.edgeIdRag2dToRag[edge.id] = edge3d.id
                    self.edgeIdRagToRag2d[edge3d.id] = edge.id
                self.visuImg2d = self.visuImg[:,:,self.zOffset]

                # update edge 2d status:
                for i in numpy.arange(self.edgeLabels2d.shape[0]):
                    self.edgeLabels2d[i] = self.edgeLabels[self.edgeIdRag2dToRag[i]]

            elif self.dim==2:
                self.rag2d = self.rag
                self.visuImg2d = self.visuImg
                self.edgeIdRag2dToRag = dict()
                for edge in self.rag.edgeIter():
                    self.edgeIdRag2dToRag[edge.id] = edge.id

                self.edgeIdRagToRag2d = self.edgeIdRag2dToRag
                self.edgeLabels2d = self.edgeLabels

            else:
                print 'warning: bad dimension!'


        def scroll(self, event):
            import pylab as plt
            if self.dim==3:
                if event.button == 'up':
                    self.zOffset += 1
                else:
                    self.zOffset -= 1

                self.zOffset = self.zOffset % self.visuImg.shape[2]
                self.slice2d()
                imgWithEdges = self.rag2d.showEdgeFeature(self.visuImg2d, self.edgeLabels2d,returnImg=True, labelMode=self.labelMode)
                self.implot.set_data(numpy.swapaxes(imgWithEdges,0,1))
                plt.draw()
        def on_motion(self, event):

            if self.press is None:
                return

            print event.xdata, event.ydata
            self.handle_click(event)

        def on_release(self, event):
            self.press = None

        def onclick(self, event):
            self.press = event.xdata, event.ydata
            print event.xdata, event.ydata
            try:
                self.handle_click(event)
            except:
                pass
        def handle_click(self, event):

            import pylab as plt
            if event.button==1:
                self.currentLabel = 1
            if event.button==2:
                self.currentLabel = 0
            if event.button==3:
                self.currentLabel = -1


            img = self.img
            rag  = self.rag2d
            labels = rag.baseGraphLabels
            shape = img.shape
            if event.xdata != None and event.ydata != None:
                xRaw,yRaw = event.xdata,event.ydata
                if xRaw >=0.0 and yRaw>=0.0 and xRaw<img.shape[0] and yRaw<img.shape[1]:
                    x,y = long(math.floor(event.xdata)),long(math.floor(event.ydata))

                    #print "X,Y",x,y
                    l = labels[x,y]
                    others  = []

                    bs = self.brushSize
                    for xo in range(-1*bs, bs+1):
                        for yo in range(-1*bs, bs+1):
                            xx = x+xo
                            yy = y+yo
                            if xo is not 0 or yo is not 0:
                                if  xx >=0 and xx<shape[0] and \
                                    yy >=0 and yy<shape[0]:
                                    otherLabel = labels[xx, yy]
                                    if l != otherLabel:
                                        edge = rag.findEdge(long(l), long(otherLabel))
                                    #print edge
                                        others.append((xx,yy,edge))
                                        #break
                        #if other is not None:
                        #    pass

                    if self.labelMode:
                        for other in others:
                            eid = other[2].id
                            oldLabel  = self.edgeLabels[self.edgeIdRag2dToRag[eid]]

                            if self.currentLabel == oldLabel:
                                newLabel = oldLabel
                            else:
                                newLabel = self.currentLabel



                            self.edgeLabels[self.edgeIdRag2dToRag[eid]] = newLabel
                            self.edgeLabels2d[eid] = newLabel
                        imgWithEdges = rag.showEdgeFeature(self.visuImg2d, self.edgeLabels2d,returnImg=True, labelMode=self.labelMode)
                        self.implot.set_data(numpy.swapaxes(imgWithEdges,0,1))
                        plt.draw()


    TinyEdgeLabelGui.__module__ = 'vigra.graphs'
    graphs.TinyEdgeLabelGui = TinyEdgeLabelGui


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








    def seededSegmentation(graph, nodeMap=None, edgeMap=None, seeds=None, alg='ws',out=None,**kwargs):
        """
            alg:
                - 'ws' watershed
                - 'sp' shortest path
                - 'crf' crf/mrf method
                - 'hc' hierarchical-clustering method
        """

        if alg == 'ws':
            # "default" node weighted watershed
            if nodeMap is not None and edgeMap is None:
                seg = graphs.nodeWeightedWatersheds(graph=graph,
                                                         nodeWeights=nodeMap,
                                                         seeds=seeds,out=out)
            # edge weighted watershed
            elif nodeMap is None and edgeMap is not None:
                seg = graphs.edgeWeightedWatersheds(graph=graph,
                                                         edgeWeights=edgeMap,
                                                         seeds=seeds,out=out)
            # hybrid (not yet implemented)
            elif nodeMap is not None and edgeMap is not None:
                raise RuntimeError("Not Yet Implemented")
            else :
                # error
                raise RuntimeError("error")

        elif alg == 'sp':
            # "default" shortest path
            if nodeMap is None and edgeMap is  None:
                raise RuntimeError("Not Yet Implemented")
            elif nodeMap is not None or edgeMap is not None:
                if nodeMap is None:
                    nodeMap = graphs.graphMap(graph,'node',dtype='float32')
                    nodeMap[:] = 0
                if edgeMap is None:
                    edgeMap = graphs.graphMap(graph,'edge',dtype='float32')
                    edgeMap[:] = 0
                seg = graphs.shortestPathSegmentation(graph=graph,
                                                           edgeWeights=edgeMap,
                                                           nodeWeights=nodeMap,
                                                           seeds=seeds,out=out)

            else :
                # error
                raise RuntimeError("error")

        elif alg == 'crf':
            raise RuntimeError("Not Yet Implemented")


        return seg

    seededSegmentation.__module__ = 'vigra.graphs'
    graphs.seededSegmentation = seededSegmentation


    def wsDtSegmentation(pmap, pmin, minMembraneSize, minSegmentSize, sigmaMinima, sigmaWeights, cleanCloseSeeds=True):
        """A probability map 'pmap' is provided and thresholded using pmin.
        This results in a mask. Every connected component which has fewer pixel
        than 'minMembraneSize' is deleted from the mask. The mask is used to
        calculate the signed distance transformation.

        From this distance transformation the segmentation is computed using
        a seeded watershed algorithm. The seeds are placed on the local maxima
        of the distanceTrafo after smoothing with 'sigmaMinima'.

        The weights of the watershed are defined by the inverse of the signed
        distance transform smoothed with 'sigmaWeights'.

        'minSegmentSize' determines how small the smallest segment in the final
        segmentation is allowed to be. If there are smaller ones the corresponding
        seeds are deleted and the watershed is done again.

        If 'cleanCloseSeeds' is True, multiple seed points that are clearly in the
        same neuron will be merged with a heuristik that ensures that no seeds of
        two different neurons are merged.
        """

        def cdist(xy1, xy2):
            # influenced by: http://stackoverflow.com/a/1871630
            d = numpy.zeros((xy1.shape[1], xy1.shape[0], xy1.shape[0]))
            for i in numpy.arange(xy1.shape[1]):
                d[i,:,:] = numpy.square(numpy.subtract.outer(xy1[:,i], xy2[:,i]))
            d = numpy.sum(d, axis=0)
            return numpy.sqrt(d)

        def findBestSeedCloserThanMembrane(seeds, distances, distanceTrafo, membraneDistance):
            """ finds the best seed of the given seeds, that is the seed with the highest value distance transformation."""
            closeSeeds = distances <= membraneDistance
            numpy.zeros_like(closeSeeds)
            # iterate over all close seeds
            maximumDistance = -numpy.inf
            mostCentralSeed = None
            for seed in seeds[closeSeeds]:
                if distanceTrafo[seed[0], seed[1], seed[2]] > maximumDistance:
                    maximumDistance = distanceTrafo[seed[0], seed[1], seed[2]]
                    mostCentralSeed = seed
            return mostCentralSeed


        def nonMaximumSuppressionSeeds(seeds, distanceTrafo):
            """ removes all seeds that have a neigbour that is closer than the the next membrane

            seeds is a list of all seeds, distanceTrafo is array-like
            return is a list of all seeds that are relevant.

            works only for 3d
            """
            seedsCleaned = set()

            # calculate the distances from each seed to the next seeds.
            distances = cdist(seeds, seeds)
            for i in numpy.arange(len(seeds)):
                membraneDistance = distanceTrafo[seeds[i,0], seeds[i,1], seeds[i,2]]
                bestAlternative = findBestSeedCloserThanMembrane(seeds, distances[i,:], distanceTrafo, membraneDistance)
                seedsCleaned.add(tuple(bestAlternative))
            return numpy.array(list(seedsCleaned))


        def volumeToListOfPoints(seedsVolume, threshold=0.):
            return numpy.array(numpy.where(seedsVolume > threshold)).transpose()


        def placePointsInVolumen(points, shape):
            volumen = numpy.zeros(shape)
            points = numpy.maximum(points, numpy.array((0, 0, 0)))
            points = numpy.minimum(points, numpy.array(shape) - 1)
            for point in (numpy.floor(points)).astype(int):
                volumen[point[0], point[1], point[2]] = 1
            return volumen

        # get the thresholded pmap
        binary = numpy.zeros_like(pmap, dtype=numpy.uint32)
        binary[pmap >= pmin] = 1

        # delete small CCs
        labeled = analysis.labelVolumeWithBackground(binary)
        analysis.sizeFilterSegInplace(labeled, int(numpy.max(labeled)), int(minMembraneSize), checkAtBorder=True)

        # use cleaned binary image as mask
        mask = numpy.zeros_like(binary, dtype = numpy.float32)
        mask[labeled > 0] = 1.

        # perform signed dt on mask
        dt = filters.distanceTransform3D(mask)
        dtInv = filters.distanceTransform3D(mask, background=False)
        dtInv[dtInv>0] -= 1
        dtSigned = dt.max() - dt + dtInv

        dtSignedSmoothMinima = filters.gaussianSmoothing(dtSigned, sigmaMinima)
        dtSignedSmoothWeights = filters.gaussianSmoothing(dtSigned, sigmaWeights)

        seeds = analysis.localMinima3D(dtSignedSmoothMinima, neighborhood=26, allowAtBorder=True)

        if cleanCloseSeeds:
            seeds = nonMaximumSuppressionSeeds(volumeToListOfPoints(seeds), dt)
            seeds = placePointsInVolumen(seeds, mask.shape).astype(numpy.uint32)

        seedsLabeled = analysis.labelVolumeWithBackground(seeds)
        segmentation = analysis.watershedsNew(dtSignedSmoothWeights, seeds = seedsLabeled, neighborhood=26)[0]

        analysis.sizeFilterSegInplace(segmentation, int(numpy.max(segmentation)), int(minSegmentSize), checkAtBorder=True)

        segmentation = analysis.watershedsNew(dtSignedSmoothWeights, seeds = segmentation, neighborhood=26)[0]

        return segmentation


    wsDtSegmentation.__module__ = 'vigra.analysis'
    analysis.wsDtSegmentation = wsDtSegmentation



    def agglomerativeClustering(graph,edgeWeights=None,edgeLengths=None,nodeFeatures=None,nodeSizes=None,
            nodeLabels=None,nodeNumStop=None,beta=0.5,metric='l1',wardness=1.0,out=None):
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

        if nodeLabels is None:
            nodeLabels = graphs.graphMap(graph,'node',dtype='uint32')



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
                                                    nodeLabels=nodeLabels,
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



    def minEdgeWeightNodeDist(mergeGraph,edgeWeights=None,edgeLengths=None,nodeFeatures=None,nodeSizes=None,
        nodeLabels=None,outWeight=None,
        beta=0.5,metric='squaredNorm',wardness=1.0, gamma=10000000.0):
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

            if nodeLabels is None :
                nodeLabels = graphs.graphMap(graph,'node',dtype='uint32')
                nodeLabels[:]=0


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
            print "nodeLabels ",nodeLabels.shape, nodeLabels.dtype
            op = graphs.__minEdgeWeightNodeDistOperator(mergeGraph,edgeWeights,edgeLengths,nodeFeatures,nodeSizes,outWeight,nodeLabels,
                float(beta),nd,float(wardness),float(gamma))


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


    def gaussianRankOrder(image, minVal=None, maxVal=None,
                     bins=20, sigmas=None, ranks=[0.1,0.25,0.5,0.75,0.9],
                     out=None):
        image = numpy.require(image.squeeze(),dtype='float32')
        nDim = image.ndim
        if sigmas is None:
            sigmas =  (2.0,)*nDim + (float(bins)/10.0,)

        ranks = numpy.require(ranks,dtype='float32')
        sigmas = numpy.require(sigmas,dtype='float32')
        assert len(sigmas) == image.ndim + 1




        if minVal is None :
            minVal = image.min()
        if maxVal is None :
            maxVal = image.max()

        #print "image",image.shape,image.dtype
        #print "ranks",ranks.shape,ranks.dtype
        #print "sigmas",sigmas
        return histogram._gaussianRankOrder(image=image,
                                            minVal=float(minVal),
                                            maxVal=float(maxVal),
                                            bins=int(bins),
                                            sigmas=sigmas,ranks=ranks,
                                            out=out)

    gaussianRankOrder.__module__ = 'vigra.histogram'
    histogram.gaussianRankOrder = gaussianRankOrder


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






def _genBlockwiseFunctions():

    def makeTuple(val, ndim):
        tvals = None
        if isinstance(val, Number):
            tvals = (float(val),)*ndim
        else :
            tvals = tuple(val)
            if len(tvals) != ndim:
                raise RuntimeError("sigma/innerScale/outerScale must be as long as ndim, or must be a scalar")
        return tvals

    def getConvolutionOptionsClass(ndim):
        assert ndim >=2 and ndim <= 5
        if ndim == 2 :
            return blockwise.BlockwiseConvolutionOptions2D
        elif ndim == 3 :
            return blockwise.BlockwiseConvolutionOptions3D
        elif ndim == 4 :
            return blockwise.BlockwiseConvolutionOptions4D
        elif ndim == 5 :
            return blockwise.BlockwiseConvolutionOptions5D

    def convolutionOptions(blockShape, sigma=None,innerScale=None, outerScale=None, numThreads = cpu_count()):
        ndim = len(blockShape)
        options = getConvolutionOptionsClass(ndim)()
        options.blockShape = blockShape
        options.numThreads = numThreads

        if sigma is not None:
            sigma = makeTuple(sigma,ndim)
            options.stdDev = sigma

        if innerScale is not None:
            options.innerScale = makeTuple(innerScale,ndim)

        if outerScale is not None:
            options.outerScale = makeTuple(outerScale,ndim)

        return options

    convolutionOptions.__module__ = 'vigra.blockwise'
    blockwise.convolutionOptions = convolutionOptions
    blockwise.convOpts = convolutionOptions

    def gaussianSmooth(image,options,out=None):
        out = blockwise._gaussianSmooth(image,options,out)
        return out
    gaussianSmooth.__module__ = 'vigra.blockwise'
    blockwise.gaussianSmooth = gaussianSmooth

    def gaussianGradient(image,options,out=None):
        out = blockwise._gaussianGradient(image,options,out)
        return out
    gaussianGradient.__module__ = 'vigra.blockwise'
    blockwise.gaussianGradient = gaussianGradient

    def gaussianGradientMagnitude(image,options,out=None):
        out = blockwise._gaussianGradientMagnitude(image,options,out)
        return out
    gaussianGradientMagnitude.__module__ = 'vigra.blockwise'
    blockwise.gaussianGradientMagnitude = gaussianGradientMagnitude


    def hessianOfGaussianEigenvalues(image,options,out=None):
        out = blockwise._hessianOfGaussianEigenvalues(image,options,out)
        return out
    hessianOfGaussianEigenvalues.__module__ = 'vigra.blockwise'
    blockwise.hessianOfGaussianEigenvalues = hessianOfGaussianEigenvalues

    def hessianOfGaussianFirstEigenvalue(image,options,out=None):
        out = blockwise._hessianOfGaussianFirstEigenvalue(image,options,out)
        return out
    hessianOfGaussianFirstEigenvalue.__module__ = 'vigra.blockwise'
    blockwise.hessianOfGaussianFirstEigenvalue = hessianOfGaussianFirstEigenvalue

    def hessianOfGaussianLastEigenvalue(image,options,out=None):
        out = blockwise._hessianOfGaussianLastEigenvalue(image,options,out)
        return out
    hessianOfGaussianLastEigenvalue.__module__ = 'vigra.blockwise'
    blockwise.hessianOfGaussianLastEigenvalue = hessianOfGaussianLastEigenvalue


_genBlockwiseFunctions()
del _genBlockwiseFunctions


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





def pmapSeeds(pmap):
    pass
