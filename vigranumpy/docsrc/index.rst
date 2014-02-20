.. vigranumpy documentation master file, created by
   sphinx-quickstart on Thu Dec 03 15:51:41 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Vigranumpy Reference
====================

.. contents::

.. toctree::
   :maxdepth: 2

Introduction
------------

Vigranumpy exports the functionality of the C++ image processing library `VIGRA <../vigra/index.html>`_ to Python. It can be invoked by importing the vigra module::

    import vigra

Vigranumpy is based on the popular `numpy <http://numpy.scipy.org/>`_ module and uses its ndarray data structure to represent image and volume data. It introduces the C++ wrapper class NumpyArray_ to allow transparent execution of VIGRA C++ functions on numpy arrays. Thus, vigranumpy is fully interoperable with existing numpy functionality, including various tools for image display such as matplotlib. Since vigranumpy uses `boost_python <http://www.boost.org/doc/libs>`_, it is able to support function overloading (which plain Python does not provide), so that calling syntax is largely uniform, regardless of the type and dimension of the input arguments.

Basic calling syntax is similar to C++, with one important difference: Arguments for output images are optional. If no output image is provided, vigranumpy will allocate it as appropriate. In either case, the output image will be returned by the function, for example::

    # allocate new result image
    >>> smoothImage = vigra.gaussianSmoothing(inputImage, scale)
    
    # reuse and overwrite existing result image
    >>> smoothImage = vigra.gaussianSmoothing(inputImage, scale, out=smoothImage)
    
Unless otherwise noted, all functions expect and create arrays with dtype=numpy.float32.

Another important concern is the interpretation and ordering of the array axes. Numpy does not provide any means to attach semantics to axes, but relies purely on the convention that the most important axis is last, as in ``array[y, x]`` or ``array[z, y, x]`` ("C-order"). However, there is no way to enforce this convention in a program, since arrays can be transposed outside of the user's control (e.g. when saving data to a file). Moreover, many imaging libraries (e.g. `Image Magick <http://www.imagemagick.org/>`_, `OpenCV <http://opencv.willowgarage.com/>`_, `Qt <http://qt-project.org/>`_ and the C++ version of VIGRA) use the opposite convention where the x-axis comes first, as in ``array[x, y]`` or ``array[x, y, z]``. This makes it very difficult to convert solutions developed in Python into a fast C++ version, because one has to reverse all indices without making mistakes. Matters become even more complicated when multi-channel (e.g. color) images are considered -- should the color index now be first or last?

To solve these ambiguities in a clean way, vigranumpy introduces the concept of **axistags** which is realized in class :class:`vigra.AxisTags`. Every :class:`~vigra.VigraArray` (which is a subclass of numpy.ndarray) gets a new property ``array.axistags`` that describes axis semantics, and all vigranumpy functions account for and preserve axistag information. Unfortunately, this functionality cannot easily be retrofitted to numpy.ndarray itself. Therefore, we employ the following conversion rules between Python and C++ arrays:

* When the Python array has **no** ``array.axistags`` property, it is mapped to the C++ NumpyArray 
  **without** any change in axis ordering. Since most VIGRA functions can work on arbitrarily 
  transposed arrays, you will get correct results, but execution may be slower because the 
  processor cache is poorly utilized in certain axis orders.
  
  Moreover, this may lead to overload resolution ambiguities. For example, when the array has shape 
  ``(3, 60, 40)``, vigranumpy has no way to decide if this is a 2-dimensional RGB image or
  a 3-dimensional array that happens to have only 3 slices. Thus, vigranumpy may not always 
  execute the function you actually intended to call.
  
* When the Python array **has** the ``array.axistags`` property, it is transposed into a 
  **canonical** axis ordering before vigranumpy executes a function, and the results are 
  transposed back into the original ordering. Likewise, functions that change axis ordering
  (such as ``array.swapaxes(0,1)``) or reduce the number of axes (such as ``array.max(axis=1)``)
  as well as array arithmetic operations preserve axistags (see section :ref:`sec-dtype-coercion`). 
  Thus, you can work in any desired axis order without loosing control. Overload ambiguities 
  can no longer occur because a function cannot be called when the axistags are unsuitable.

Detailed information about the use of axistags is given in section :ref:`sec-vigraarray` below. Section :ref:`sec-own-modules` describes how you can take advantage of the axistags mechanism in your own C++ code.

.. _sec-vigraarray:
    
Axistags and the VigraArray Data Structure
------------------------------------------

While vigranumpy can directly work on numpy.ndarrays, this would not give us the advantages of axistags as described above. Therefore, vigranumpy introduces its own array class :class:`~vigra.VigraArray` which is a subclass of numpy.ndarray, but re-implements many of its methods so that axistags are respected. Arrays with a conforming ``axistags`` property are most easily constructed by one of the predefined :ref:`array factories <subsec-array-factories>`. A **view with axistags** can be created from an existing numpy.ndarray by means of the function :py:func:`~vigra.taggedView` (in contrast, factory functions create copies of the given arrays, not views). We illustrate the ideas by some examples::

    >>> width, height, depth = 300, 200, 3
    
    # create a 2-dimensional RGB image
    >>> rgb = vigra.RGBImage((width, height))
    >>> rgb.shape
    (300, 200, 3)
    >>> rgb.axistags             # short output: only axis keys
    x y c
    >>> print rgb.axistags       # long output
    AxisInfo: 'x' (type: Space)
    AxisInfo: 'y' (type: Space)
    AxisInfo: 'c' (type: Channels) RGB
    
    # create a 3-dimensional scalar volume
    >>> volume = vigra.ScalarVolume((width, height, depth))
    >>> volume.shape
    (300, 200, 3)        # same shape as before
    >>> volume.axistags
    x y z                # but different semantic interpretation
    >>> print volume.axistags
    AxisInfo: 'x' (type: Space)
    AxisInfo: 'y' (type: Space)
    AxisInfo: 'z' (type: Space)

It is also possible to attach additional information to the axistags, in particular the resolution of the axis, and a text comment. The resolution will be correctly adjusted when the image is resized::

    >>> rgb.axistags['x'].resolution = 1.2  # in some unit of length
    >>> rgb.axistags['y'].resolution = 1.4  # in some unit of length
    >>> rgb.axistags['c'].description = 'fluorescence microscopy, DAPI and GFP staining'
    >>> print rgb.axistags
    AxisInfo: 'x' (type: Space, resolution=1.2)
    AxisInfo: 'y' (type: Space, resolution=1.4)
    AxisInfo: 'c' (type: Channels) fluorescence microscopy, DAPI and GFP staining
    
    # interpolate the image to twice its original size
    >>> rgb2 = vigra.sampling.resize(rgb, shape=(2*width-1, 2*height-1))
    >>> print rgb2.axistags
    AxisInfo: 'x' (type: Space, resolution=0.6)
    AxisInfo: 'y' (type: Space, resolution=0.7)
    AxisInfo: 'c' (type: Channels) fluorescence microscopy, DAPI and GFP staining

When the array is transposed, the axistags are transposed accordingly. When axes are dropped from the array, the corresponding entries are dropped from the axistags property::

    # transpose the volume into reverse axis order
    >>> transposed_volume = volume.transpose()
    >>> transposed_volume.shape
    (3, 200, 300)
    >>> transposed_volume.axistags
    z y x
    
    # get a view to the first slice (z == 0)
    >>> first_slice = volume[..., 0]
    >>> first_slice.shape
    (300, 200)
    >>> first_slice.axistags
    x y
    
    # get the maximum of each slice
    >>> volume.max(axis=0).max(axis=0)
    VigraArray(shape=(3,), axistags=z, dtype=float32, data=
    [ 0.  0.  0.])
    
    # likewise, but specify axes by their keys
    >>> volume.max(axis='x').max(axis='y')
    VigraArray(shape=(3,), axistags=z, dtype=float32, data=
    [ 0.  0.  0.])

The initial ordering of the axes is controlled by the argument ``order`` that can optionally be passed to the VigraArray constuctor or the factory functions. If ``order`` is not explicitly provided, it is determined by the static property :attr:`VigraArray.defaultOrder` (which yields 'V' order). The following orders are currently supported:

.. _array-order-parameter:

    'C' order:
        Both strides and axes are arranged in descending order, as in a 
        plain numpy.ndarray. For example, axistags will be 'y x c' or 
        'z y x c'. array.flags['C_CONTIGUOUS'] will be true.

    'F' order:
        Both strides and axes are arranged in ascending order, i.e. 
        opposite to 'C' order. For example, axistags will be 'c x y' 
        or 'c x y z'. array.flags['F_CONTIGUOUS'] will be true.

    'V' order:
        VIGRA-order is an interleaved memory layout that simulates 
        vector-valued pixels or voxels: Channels will be the last axis 
        and have the smallest stride, whereas all other axes are arranged 
        in ascending order. For example, axistags will be 'x y c' or 
        'x y z c'.

    'A' order:
        Defaults to 'V' when a new array is created, and means
        'preserve order' when an existing array is copied.
        
The meaning of 'ascending' or 'descending' order is determined by two rules: the primary order is according to axis type (see :class:`vigra.AxisType`), where ``Channels < Space < Angle < Time < Frequency < Unknown``. The secondary order (between axes of the same type) is lexicographic, such that 'x' < 'y' < 'z'. Usage examples::

    >>> rgbv = vigra.RGBImage((width, height), order='V')
    >>> rgbv.shape
    (300, 200, 3)
    >>> rgbv.axistags
    x y c
    
    >>> rgbc = vigra.RGBImage((width, height), order='C')
    >>> rgbc.shape
    (200, 300, 3)
    >>> rgbc.axistags
    y x c
    
    >>> rgbf = vigra.RGBImage((width, height), order='F')
    >>> rgbf.shape
    (3, 300, 200)
    >>> rgbf.axistags
    c x y
    
Functions that reduce the array to a one-dimensional shape (``flatten()``, ``flat``, ``ravel()``, ``take()``) always transpose the array into 'C' order before flattening.

Axistags are stored in a list-like class :class:`vigra.AxisTags`, whose individual entries are of type :class:`vigra.AxisInfo`. The simplest way to attach axistags to a plain numpy.ndarray (by creating a view of type VigraArray) is via the convenience function :func:`vigra.taggedView`.

More On the Motivation and Use of Axistags
------------------------------------------

History of the problem
^^^^^^^^^^^^^^^^^^^^^^

A Fortran array declaration::

    real f(20,10)

is compiled such that the elements of the first index are consecutive in memory. In contrast, a C array declaration::

    float c[20][10];

places the elements of the last index consecutive in memory. The two possibilities are commonly referred to as "Fortran order" and "C order", but we will see that these terms are actually ambiguous because their meaning depends on what the array is used for.

Arrays as Matrices
^^^^^^^^^^^^^^^^^^

When the array represents a matrix, users require the syntax to conform to the mathematical syntax conventions (in order to avoid bugs when typing formulas). In mathematics, m\ :sub:`ij` means that the first index i refers to rows, and the second index j refers to columns. Thus, Fortran's ``f(i, j)`` and C's ``c[i][j]`` have exactly the same meaning, but the memory layouts of the two matrices differ: Since in Fortran the first index changes quickest, this is referred to as the "column major" format (the elements of a column are consecutive in memory). In contrast, C uses the "row major" format (the elements of a row are consecutive in memory).

In short: When the array is a matrix, the syntax is the same, but the memory layouts differ.

Arrays as Images
^^^^^^^^^^^^^^^^

When the array represents an image, users require the memory layout to be the same as in the image files which store the data (in order to ensure fast I/O). All image formats (JPEG, TIFF, PNG, ...) follow the convention of analog television to scan an image left to right, then top to bottom. Consequently, the horizontal pixels are consecutive in memory, and therefore Fortran (where the first axis changes fastest) must use the syntax ``f(x, y)``, whereas C (where the last axis changes fastest) must use ``c[y][x]``.

In short: When the array is an image, the memory layout is the same, but the syntax differs.

Thus, we have four basic conventions: FM and CM for matrices, FI and CI for images. The meaning of the terms "Fortran order" and "C order" depends on whether we are talking about matrices (where they refer to FM and CM, i.e. a difference in memory layout) or images (where they refer to FI and CI, i.e. a difference in index order).

Multi-Dimensional Arrays
^^^^^^^^^^^^^^^^^^^^^^^^

When we add more dimensions, the confusion increases, because there are no universally accepted memory and indexing conventions. For example, an RGB image can be stored in "interleaved format"::

    RGB RGB RGB ...
    RGB RGB RGB ...
    :
    :

where the color values of each pixel are consecutive in memory, or in "banded format"::

    R R R ...
    R R R ...
    :
    G G G ...
    G G G ...
    :
    B B B ...
    B B B ...
    :
    
where we have a separate scalar image for each color. In Fortran, interleaved and banded images must be indexed as ``f(color, x, y)`` and ``f(x, y, color)`` respectively, whereas in C we must use ``c[y][x][color]`` or ``c[color][y][x]``.

VIGRA and numpy
^^^^^^^^^^^^^^^

From the beginning, VIGRA adopted Fortran conventions, i.e. its default behavior is according to FM and FI (this is possible because VIGRA uses array classes, where the mapping from indices to memory is encapsulated in the appropriate way).

In contrast, numpy adopted C conventions, i.e. its default behavior is CM and CI. 

In addition, both packages provide array views which keep the memory layout intact, but change the index order. Thus, VIGRA also supports the CI convention, and numpy also supports FI. Note that changing the index order is only allowed for images. Matrices always use the fixed index order dictated by mathematics where transpose(m) is a well-defined mathematical operation (which just happens to revert the index order). Therefore, the existence of array views does not imply that VIGRA supports CM or numpy supports FM.

However, numpy's array constructor provides an argument 'order' which can take the values 'C' (default) and 'F' to construct arrays with C or Fortran memory order. By this mechanism, numpy also supports the FM convention (and thus all four basic possibilities).

But here is the catch: When you get a numpy array view, you have no way to tell which convention it adheres to. It simply doesn't contain this information.

The usual solution to this problem is to enforce a fixed axis order in the entire application, but this workaround breaks down when the application must simultaneously handle arrays with different meaning (e.g. sequences 'xyt' vs. volumes 'xyz' vs. multi-spectral images 'xyc') or when the application uses modules with conflicting requirements (e.g. numpy's 'yx' vs. PyQt's 'xy').

Vigranumpy Axistags
^^^^^^^^^^^^^^^^^^^

This is precisely where axistags enter: They attach information to array views that allows the user to figure out which convention applies to a given view. Thus, the primary purpose of axistags is entirely passive - they just keep track of how users manipulate the axis order when creating new array views. The main use case of axistags is therefore to re-adjust the axis order whenever an algorithm has specific order requirements. Using axistags, one can easily ensure that arrays conform to the desired order at the beginning of every algorithm. Consider the following example: Let 'image' be a scalar 2D image with axistags. Then we can ask for the image's width and height independently of the current axis order::

    width  = image.width    # this works for any axis order!
    height = image.height
    
Now suppose we want to execute a numpy algorithm which expects the [y, x] ordering. We simply transpose the array before calling the algorithm like this::

    # adjust the axis order
    numpy_image = image.transposeToNumpyOrder()
    
    # execute the algorithm
    for y in xrange(height):
        for x in xrange(width):
            numpy_image[y, x] = ...   # note the index order
            
When we want to execute a VIGRA algorithm which expects the [x, y] ordering, we do::

    # adjust the axis order
    vigra_image = image.transposeToVigraOrder()
    
    # execute the algorithm
    for y in xrange(height):
        for x in xrange(width):
            vigra_image[x, y] = ...   # note the index order

Notice that the order of the loops (the inner loop runs over x) is the same in both cases: To maximize cache locality and therefore speed, the inner loop should operate on consecutive memory cells. Since image memory order derives from file memory order (where the x-axis is consecutive), and array views can never change the memory order, this loop ordering is always preferable, regardless of the index order.

Vigranumpy Conventions
^^^^^^^^^^^^^^^^^^^^^^

To handle axis meaning in a well-defined way, vigranumpy adopts the following conventions, which are designed such that The Right Thing (TM) should happen automatically, at least most of the time:

1. When the array represents a matrix, no axistags are allowed because the index order has a fixed semantic meaning and must not be messed around with. In vigranumpy, this requirement is enforced by an assertion::

    vigra_precondition( !matrix.axistags(), "matrix must not have axistags");
    
   in the C++ gluecode functions. This applies, for example, to the feature matrices passed to a random forest and to unsupervised learning algorithms. If desired, we can introduce additional axistags for features and samples in the future because this is a common use case.

2. When arrays represent image data with up to five dimensions, axistags should be used. To sort indices according to the requirements of the next algorithm to be executed, the appropriate convenience function should be called (many more convenience functions are documented in :py:class:`vigra.VigraArray`)::

    numpy_array   = array.transposeToNumpyOrder()    # gives 'yx', 'zyx' etc. 
    vigra_array   = array.transposeToVigraOrder()    # gives 'xy', 'xyz' etc.
    ilastik_array = array.view5D()                   # gives 'txyzc' (inserts singleton axes if necessary)
    user_defined  = array.withAxes('y', 'x', 't')    # specify order explicitly (inserts singleton axes if necessary)
    
   Algorithms with special order requirements can then check for the correct order in an assertion. 

3. The function ``vigra.taggedView()`` allows to attach axistags to an array very conveniently. For example, when you know from the context that the axes of a given array are to be interpreted as 'xyzt' in that order, you can make this knowledge explicit by calling::

    tagged_array = vigra.taggedView(array, 'xyzt')

4. When you call a vigranumpy function that executes a C++ VIGRA function, ``image.transposeToVigraOrder()`` will be invoked automatically on the input arrays, and the original axis order will be restored in the output arrays.

5. When you call a vigranumpy function that is forwarded to a numpy function (in particular, a ufunc like ``+``, ``-``, ``sqrt`` etc.), ``image.transposeToNumpyOrder()`` will be invoked automatically on the input arrays, and the original axis order will be restored in the output arrays.

6. When vigranumpy writes arrays to a file, it will always order the axes such that the memory order conforms to the established file conventions (e.g. values along the x-axis are consecutive). In particular, when you use ``vigra.impex.writeHDF5()`` to create HDF5 datasets, ``array.transposeToNumpyOrder()`` will be called before writing the data (this is a consequence of item 5, because ``writeHDF5()`` eventually forwards the actual work to h5py). In addition, the axistags (in numpy order) will be stored in a dataset attribute ``axistags``.

7. When vigranumpy reads data from a file, it preserves the file's memory order and attaches the appropriate axistags. In case of images, the axis order follows from the usual file conventions. If you call ``vigra.impex.readHDF5()`` to read HDF5, the axistags will be read from the attribute ``axistags`` (if present). Upon return, the read functions automatically call ``array.transposeToVigraOrder()``, but this only changes the index order, not the memory layout. This latter convention was adopted to ensure that the default index order is the same in the Python and C++ versions of VIGRA. 

8. When you display an image via ``image.imshow()`` or ``image.show()``, the axes are re-ordered automatically such that the image is displayed upright (i.e. x goes to the right, y goes down). If you want to override this (i.e. want to enforce transposed display), you can remove the axistags by calling ``image.view(numpy.ndarray)``.

Item 4. has two important consequences one should be aware of:

* When a function needs parameters that depend on the axis order (such as the new shape in ``vigra.sampling.resize()`` or axis-dependent sigmas in ``vigra.filters.gaussianSmoothing()``), these parameters must be re-ordered on the C++ side in the same way as the axes. This is achieved by a call to ``NumpyArray::permuteLikewise()`` in the C++ gluecode functions.

* When a vigranumpy function computes a gradient image, the gradient vector elements will always be stored in the order (dx, dy, dz), regardless of the array's original axis ordering. The same applies to any function producing vector-valued elements whose interpretation depends on the axis meaning (e.g. the Hessian matrix and the structure tensor). For example, the output of the Hessian is ordered as (dxx, dxy, dxz, dyy, dyz, dzz).

Axistag Reference
-----------------

.. autoclass:: vigra.AxisType

----------------

.. autoclass:: vigra.AxisInfo
    :members: key, typeFlags, resolution, description, isSpatial, isTemporal, isChannel, isFrequency, isAngular, isType, compatible
    
----------------

.. autoclass:: vigra.AxisTags
    :members: index, channelIndex, innerNonchannelIndex, axisTypeCount, setChannelDescription, toJSON, fromJSON

VigraArray Reference
--------------------

.. autoclass:: vigra.VigraArray
    :show-inheritance:
    :members: defaultAxistags, channelIndex, innerNonchannelIndex, channels, spatialDimensions, width, height, depth, duration, dropChannelAxis, insertChannelAxis, withAxes, view5D, asRGB, __getitem__, subarray, bindAxis, channelIter, sliceIter, spaceIter, timeIter, copyValues, swapaxes, transpose, T, transposeToOrder, transposeToDefaultOrder, transposeToNormalOrder, transposeToNumpyOrder, transposeToVigraOrder, permutationToOrder, permutationToNormalOrder, permutationFromNormalOrder, permutationToNumpyOrder, permutationFromNumpyOrder, permutationToVigraOrder, permutationFromVigraOrder, writeHDF5, writeImage, writeSlices, receiveSocket, sendSocket, qimage, show
    
    .. attribute:: VigraArray.axistags
    
      The :class:`~vigra.AxisTags` object of this array. 

    .. attribute:: VigraArray.defaultOrder
    
      Get the default axis ordering, currently 'V' (:ref:`VIGRA order <array-order-parameter>`).

-------------

.. autofunction:: vigra.newaxis(axisinfo=vigra.AxisInfo())
.. autofunction:: vigra.taggedView
.. autofunction:: vigra.dropChannelAxis

-------------

.. _subsec-array-factories:

.. autofunction:: vigra.Image
.. autofunction:: vigra.ScalarImage
.. autofunction:: vigra.Vector2Image
.. autofunction:: vigra.Vector3Image
.. autofunction:: vigra.Vector4Image
.. autofunction:: vigra.RGBImage

-------------

.. autofunction:: vigra.Volume
.. autofunction:: vigra.ScalarVolume
.. autofunction:: vigra.Vector2Volume
.. autofunction:: vigra.Vector3Volume
.. autofunction:: vigra.Vector4Volume
.. autofunction:: vigra.Vector6Volume
.. autofunction:: vigra.RGBVolume
   

Import and Export Functions
---------------------------

The module vigra.impex defines read and write functions for image and volume data. Note
that the contents of this module are automatically imported into the vigra module, so
you may call 'vigra.readImage(...)' instead of 'vigra.impex.readImage(...)' etc.

.. automodule:: vigra.impex
   :members:

   
.. _sec-dtype-coercion:

Mathematical Functions and Type Coercion
----------------------------------------

vigranumpy supports all arithmetic and algebraic functions defined in  
`numpy.ufunc <http://docs.scipy.org/doc/numpy/reference/ufuncs.html#available-ufuncs>`_, but re-implements them in module `vigra.ufunc` to take full advantage of axistags. 

.. automodule:: vigra.ufunc


Color and Intensity Manipulation
--------------------------------

The module vigra.colors provides functions to adjust image brightness and contrast, 
and to transform between different color spaces. 
See `Color Conversions <../vigra/group__ColorConversions.html>`_ in the C++ documentation
for more information.

.. automodule:: vigra.colors
   :members:


Filters
-------

The module vigra.filters provides operators that consider a window around each pixel, compute
one or several numbers from the values in the window, and store the results in the
corresponding pixel of the output image. This includes convolution, non-linear diffusion, 
morphological operators, feature detectors (such as the structure tensor) etc.

.. automodule:: vigra.filters
   :members:


Sampling: Image Resizing and Image Pyramids
-------------------------------------------

The module vigra.sampling contains methods to change the number and/or location of
the image sampling points, such as resizing, rotation, and interpolation.

.. automodule:: vigra.sampling
   :members:
   
---------------------------------------------

Spline image views implement an interpolated view for an image which can be accessed 
at real-valued coordinates (in contrast to the plain image, which can only be
accessed at integer coordinates). Module vigra.sampling defines::

    SplineImageView0
    SplineImageView1
    SplineImageView2
    SplineImageView3
    SplineImageView4
    SplineImageView5
    
The number denotes the spline interpolation order of the respective classes. 
Below, we describe SplineImageView3 in detail, but the other classes work 
analogously. See SplineImageView_ in the C++ documentation for more detailed information.

.. autoclass:: vigra.sampling.SplineImageView3
   :members:

-------------

.. autoclass:: vigra.sampling.ImagePyramid
   :show-inheritance:
   :members:


Fourier Transforms
------------------

The module vigra.fourier contains functions for Fourier transforms, Cosine/Sine 
transforms, and Fourier-domain filters.

.. automodule:: vigra.fourier
   :members:

   
Image Analysis
--------------

The module vigra.analysis contains segmentation algorithms (e.g. watershed), edge and 
corner detection, localization of maxima and minima etc.

.. automodule:: vigra.analysis
   :members:

   
Geometry
--------

The module vigra.geometry contains geometric primitives (such as polygons) and related algorithms.

.. automodule:: vigra.geometry
   :members:


Optimization
------------

The module vigra.optimization provides functions for constrained and unconstrained linear regression.

.. automodule:: vigra.optimization
   :members:


Machine Learning
----------------

The module vigra.learning will eventually provide a wide range of machine learning 
tools. Right now, it only contains an implementation of the random forest classifier
and probabilistic latent semantic analysis (pLSA) as an example for unsupervised learning.

.. automodule:: vigra.learning
   :members:

.. autoclass:: vigra.learning.RandomForest
   :members:
   
For more information, refer to RandomForest_ in the C++ documentation.

.. autoclass:: vigra.learning.RandomForestOld
   :members:


Noise Estimation and Normalization
----------------------------------

The module vigra.noise provides noise estimation and normalization according to a
method proposed by Foerstner.

.. automodule:: vigra.noise
   :members:


Histogram and Channel Representation
----------------------------------

The module vigra.histogram provides histograms and channel representation

.. automodule:: vigra.histogram
   :members:



Graphs and Algorithms on Graphs
----------------------------------

The module vigra.graphs provides graphs and graph algorithms

.. automodule:: vigra.graphs
   :members:

.. autoclass:: vigra.graphs.GridGraphUndirected2d
   :members:

.. autoclass:: vigra.graphs.GridGraphUndirected3d
   :members:

.. autoclass:: vigra.graphs.AdjacencyListGraph
    :members:



Utilities 
----------------------------------

The module vigra.utilities provides  utilities and tools
like priority queues with changeable priorities

.. automodule:: vigra.utilities
   :members:
    





.. _sec-own-modules:

Writing Your Own C++ Modules
----------------------------

When you want to write your own vigranumpy extension modules, first make sure that you compile and link with the same versions of numpy and boost_python that your current vigranumpy installation uses. Otherwise, communication between new and existing modules will not work (and even crash). Then follow these steps:

1. Create the main module source file. This file contains the module's 'init' function. Let's assume that the module will be called 'my_module', and the file is 'my_module.cxx'. A stub for 'my_module.cxx' typically looks like this::

        // define PY_ARRAY_UNIQUE_SYMBOL (required by the numpy C-API)
        #define PY_ARRAY_UNIQUE_SYMBOL my_module_PyArray_API

        // include the vigranumpy C++ API
        #include <Python.h>
        #include <boost/python.hpp>
        #include <vigra/numpy_array.hxx>
        #include <vigra/numpy_array_converters.hxx>
        
        ... // your includes
        
        ... // implementation of your wrapper functions and classes
        
        using namespace boost::python;
        
        // the argument of the init macro must be the module name
        BOOST_PYTHON_MODULE_INIT(my_module)
        {
            // initialize numpy and vigranumpy
            vigra::import_vigranumpy();
            
            // export a function
            def("my_function", &my_function, 
                (arg("arg1"), arg("arg2"), ...),
                "Documentation");

            // export a class and its member functions
            class_<MyClass>("MyClass",
                "Documentation")
                .def("foo", &MyClass::foo,
                     (arg("arg1"), arg("arg2"), ...),
                     "Documentation")
            ;
                     
            ... // more module functionality (refer to boost_python documentation)
        }
    
2. When your module uses additional C++ source files, they should start with the following defines::

        // this must define the same symbol as the main module file (numpy requirement)
        #define PY_ARRAY_UNIQUE_SYMBOL my_module_PyArray_API
        #define NO_IMPORT_ARRAY

3. Implement your wrapper functions. Numpy ndarrays are passed to C++ via the wrapper classes NumpyArray_ and NumpyAnyArray_. You can influence the conversion from Python to C++ by using different instantiations of NumpyArray, as long as the Python array supports the axistags attribute (refer to :ref:`axis order definitions <array-order-parameter>` for the meaning of the term 'ascending order')::

            // We add a 'using' declaration for brevity of our examples.
            // In actual code, you should probably prefer explicit namespace qualification.
        using namespace vigra;

            // Accept any array type and return an arbitrary array type.
            // Returning NumpyAnyArray is always safe, because at that point
            // C++ no longer cares about the particular type of the array.
        NumpyAnyArray foo(NumpyAnyArray array);
        
            // Accept a 3-dimensional float32 array and transpose it 
            // into ascending axis order ('F' order).
        void foo(NumpyArray<3, float> array);
        
            // Accept a 2-dimensional float32 array with an arbitrary number of channels and
            // transpose the axes into VIGRA ('V') order (channels are last, other axes ascending).
            // Note that the NumpyArray dimension is 3 to account for the channel dimension.
            // If the original numpy array has no channel axis, vigranumpy will automatically
            // insert a singleton axis.
        void foo(NumpyArray<3, Multiband<float> > array);
        
            // Accept a 2-dimensional float32 array that has only a single channel
            // (that is, 'array.channels == 1' must hold on the Python side).
            // Non-channel axes are transposed into ascending order.
            // Note that the NumpyArray dimension is now 2.
        void foo(NumpyArray<2, Singleband<float> > array);
        
            // Accept a float32 array that has 2 non-channel dimensions and 
            // exactly 3 channels (i.e. 'array.channels == 3' on the Python side). 
            // Non-channel axes are transposed into ascending order.
            // Note that the NumpyArray dimension is again 2, but the pixel type is 
            // now a vector.
            // The conversion will only succeed if the channel axis is unstrided on
            // the Python side (that is, the following expression is True:
            //      array.strides[array.channelIndex] == array.dtype.itemsize).
        void foo(NumpyArray<2, TinyVector<float, 3> > array);
        void foo(NumpyArray<2, RGBValue<float> > array);
    
   Or course, these functions can also be templated. 
   
   When your functions return newly allocated arrays, it is usually desirable to transfer the input's axistags to the output (otherwise, vigranumpy will use :meth:`~vigra.VigraArray.defaultAxistags` as a fallback). There is a standard vigranumpy idiom for this task which assumes that the wrapped function has an optional parameter 'output' for a possibly pre-allocated output array. The axistags are then transferred by reshaping the output array with a ``taggedShape()`` (which is a combination of a shape and axistags)::
   
        NumpyAnyArray
        foo(NumpyArray<3, Multiband<float32> > input, 
            NumpyArray<3, Multiband<float32> > output = boost::python::object())
        {
            // Reshape only if the output array was not explicitly passed in.
            // Otherwise, use the output array as is.
            output.reshapeIfEmpty(input.taggedShape(), 
                      "error message when shape is unsuitable.");
                      
            ... // your algorithm
        }
        
   It is also possible to modify the tagged shape before it is applied to the output array::
   
        input.taggedShape()
             .resize(Shape2(new_width, new_height))
             .setChannelCount(new_channel_count)
             .setChannelDescription("a description")
             
   The C++ code can be multi-threaded when you unlock Python's global interpreter lock. After unlocking, your wrapper code must not call any Python functions, so the unlock statement should go after ``output.reshapeIfEmpty()``::
   
        NumpyAnyArray
        foo(NumpyArray<3, Multiband<float32> > input, 
            NumpyArray<3, Multiband<float32> > output = boost::python::object())
        {
            output.reshapeIfEmpty(input.taggedShape(), "Message.");
            
                // Allow parallelization from here on. The destructor of
                // _pythread will automatically regain the global interpreter lock
                // just before this function returns to Python.
            PyAllowThreads _pythread;
            
            ... // your algorithm
        }

4. Export your wrapped functions. ``boost::python::def`` is called in its usual way, with one simple extension: Since vigranumpy does not know which NumpyArray variants you are going to use, appropriate converter functions between Python and C++ must be registered on demand. You do this by enclosing your function pointer into a call to the 'registerConverters()' function::

        // in the module's init function
        def("my_function", vigra::registerConverters(&my_function),
           (arg("arg1"), ...),
           "Documentation");
           
If you need more information, it is always a good idea to look at the source code of the existing vigranumpy modules.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
