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

Vigranumpy is based on the popular `numpy <http://numpy.scipy.org/>`_ module and uses its ndarray data structure to represent image and volume data. Thus, it is fully interoperable with existing numpy functionality, including various tools for image display such as matplotlib. Since vigranumpy uses `boost_python <http://www.boost.org/doc/libs>`_, it is able to support function overloading (which plain Python does not provide), so that calling syntax is largely uniform, regardless of the type and dimension of the input arguments.

Basic calling syntax is similar to C++, with one important difference: Arguments for output images are optional. If no output image is provided, vigranumpy will allocate it as appropriate. In either case, the output image will be returned by the function, for example::

    # allocate new result image
    smoothImage = vigra.gaussianSmoothing(inputImage, scale)
    
    # reuse and overwrite existing result image
    smoothImage = vigra.gaussianSmoothing(inputImage, scale, out=smoothImage)
    
Unless otherwise noted, all functions expect and create arrays with dtype=numpy.float32.

Another important property is vigranumpy's indexing convention. In order to be compatible with the index order of the VIGRA C++ version and many other libraries (e.g. `Qt <http://qt.nokia.com/>`_ and `Image Magick <http://www.imagemagick.org/>`_), and with standard mathematical notation, images are indexed in the following order::

    value = scalarImage[x, y]
    value = multibandImage[x, y, channel]
    
    value = scalarVolume[x, y, z]
    value = multibandVolume[x, y, z, channel]

where x is the horizontal axis (increasing left to right), and y is the vertical axis (increasing top to bottom). This convention differs from the `Python Imaging Library <http://www.pythonware.com/products/pil/>`_ and Matlab, where the spatial indices must be given in reverse order (e.g. scalarImage[y, x]). Either convention has advantages and disadvantages. In the end, we considered compatibility between the Python and C++ versions of VIGRA to be a very critical feature in order to prevent subtle errors when porting from one language to the other, so we went with the [x, y] order described. Note that this convention does *not* change the internal representation of the data in memory. It only changes the indexing order, so that one can switch between the different conventions by simply swapping axes, for example::

    vigraImage  = array2D.swapaxes(0, 1).view(vigra.ScalarImage)
    array2D     = vigraImage.swapaxes(0, 1).view(numpy.ndarray)
   
    vigraVolume = array3D.swapaxes(0, 2).view(vigra.ScalarVolume)
    array3D     = vigraVolume.swapaxes(0, 2).view(numpy.ndarray)

In order to turn your own C++ VIGRA functions into Python modules, look at the VIGRA wrapper class
NumpyArray_ and the source code of the existing vigranumpy modules.

    
Image and Volume Data Structures
--------------------------------

Vigranumpy can work directly on numpy.ndarrays. However, plain ndarrays do not carry
any information about the semantics of the different coordinate axes. For example,
one cannot distinguish a 2-dimensional RGB image from a scalar volume data set that
happens to contain only three slices. In order to distinguish between arrays that 
have the same structure but different interpretation, vigra.arraytypes provides the 
following array classes::

    numpy.ndarray
        Image
            ScalarImage
            Vector2Image
            Vector3Image
                RGBImage
            Vector4Image
        Volume
            ScalarVolume
            Vector2Volume
            Vector3Volume
                RGBVolume
            Vector4Volume
            Vector6Volume
    
    list 
        ImagePyramid

where indentation encodes inheritance. Below, we describe :class:`~vigra.Image`, 
:class:`~vigra.ScalarImage`, :class:`~vigra.RGBImage`, and  :class:`~vigra.ImagePyramid` 
in detail, the other classes work analogously. The new array classes serve several purposes:

* Semantic interpretation improves code readability.

* vigra.arraytypes maximize compatibility with corresponding VIGRA C++ types. In
  particular, vigra.arraytype constructors ensure that arrays are created with 
  the most appropriate memory layout on the Python side. For example, RGBImage 
  (with dtype=numpy.float32) can be mapped to MultiArrayView<2, RGBValue<float> >. 

* The distinction of different array types allows for more fine-grained overload
  resolution and thus improves mapping from Python to C++. For example, 
  gaussianSmoothing() works differently for 2D RGBImages and 3D ScalarVolumes, 
  although both kinds of arrays have the same 3-dimensional memory layout.

* The array types help simplify use of the vigranumpy indexing convention::

    image[x, y, channel]
    volume[x, y, z, channel]
  
  In particular, they overload '__str__' and '__repr__' (used in print), 'flatten', 
  and 'imshow' (for matplotlib-based image display) so that these functions work in the 
  expected way (i.e. images are printed in horizontal scan order and are displayed 
  upright). Note that other Python imaging modules (such as PIL) use a different
  indexing convention (namely image[y, x, channel]).

* vigra.arraytypes and vigra.ufunc overload numpy.ufunc (i.e. basic mathematical 
  functions for arrays). See :ref:`sec-dtype-coercion` for details.
  
Mapping between C++ types and Python types is controlled by the following two functions:

.. autofunction:: vigra.registerPythonArrayType

.. autofunction:: vigra.listExportedArrayKeys

----------------

.. autoclass:: vigra.Image
   :show-inheritance:
   :members:
  
   .. attribute:: channels
   
      the number of channels of this image (e.g. RGB has three channels)

   .. attribute:: spatialDimensions
   
      number of spatial dimensions (always 2 for images). This is useful for 
      distinguishing RGBImage from ScalarVolume in overload resolution.

.. autofunction:: vigra.imshow

-------------

.. autoclass:: vigra.ScalarImage
   :show-inheritance:
   :members:

-------------

.. autoclass:: vigra.RGBImage
   :show-inheritance:
   :members:

-------------

.. autoclass:: vigra.ImagePyramid
   :show-inheritance:
   :members:
   

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

Vigra images and volumes support all arithmetic and algebraic functions defined in  
`numpy.ufunc <http://docs.scipy.org/doc/numpy/reference/ufuncs.html#available-ufuncs>`_. 

.. automodule:: vigra.ufunc

As usual, these functions are applied independently at each pixel.
Vigranumpy overloads the numpy-versions of these functions in order to make their
behavior more suitable for image analysis. In particular, we changed two aspects:

* The memory layout of the input arrays is preserved in the result arrays. 
  In contrast, plain numpy.ufuncs always create C-order arrays, even if 
  the inputs have a different order (e.g. as Fortran-order). 

* The value types of result arrays (i.e. their 'dtype') are determined in a way 
  that is more suitable for image processing than the original numpy conversion rules.  

Array dtype conversion (aka coercion) is implemented by the function 
vigra.ufunc.Function.common_type according to the following coercion rules:

.. automethod::  vigra.ufunc.Function.common_type(in_type, out_type)


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


Sampling
--------

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


Machine Learning
----------------

The module vigra.learning will eventually provide a wide range of machine learning 
tools. Right now, it only contains an implementation of the random forest classifier
and probabilistic latent sementic analysis (pLSA) as an example for unsupervised learning.

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


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
