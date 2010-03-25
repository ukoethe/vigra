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

Vigranumpy exports the functionality of the C++ image processing library `VIGRA <../vigra/index.html>`_ to Python. It is based on the popular `numpy <http://numpy.scipy.org/>`_ module and uses its ndarray data structure to represent image and volume data. Thus, it is fully interoperable with existing numpy functionality, including various tools for image display such as matplotlib. Since vigranumpy uses `boost_python <http://www.boost.org/doc/libs>`_, it is able to use function overloading (which plain Python does not support), so that calling syntax is largely uniform, regardless of the type and dimension of the input arguments.

Basic calling syntax is similar to C++, with one important difference: Arguments for output images are optional. If no output image is provided, vigranumpy will allocate it as appropriate. In either case, the output image will be returned by the function, for example::

    # allocate new result image
    smoothImage = gaussianSmoothing(inputImage, scale)
    
    # reuse and overwrite existing result image
    smoothImage = gaussianSmoothing(inputImage, scale, out=smoothImage)

Another important property is vigranumpy's indexing convention. In order to be compatible with the index order of the VIGRA C++ version and many other libraries (e.g. `Qt <http://qt.nokia.com/>`_ and `Image Magick <http://www.imagemagick.org/>`_), and with standard mathematical notation, images are indexed in the following order::

    value = scalarImage[x, y]
    value = multibandImage[x, y, channel]
    
    value = scalarVolume[x, y, z]
    value = multibandVolume[x, y, z, channel]

where x is the horizontal axis (increasing left to right), and y is the vertical axis (increasing top to bottom). This convention differs from the `Python Imaging Library <http://www.pythonware.com/products/pil/>`_ and Matlab, where the spatial indices must be given in reverse order (e.g. scalarImage[y, x]). Either convention has advantages and disadvantages. In the end, we considered compatibility between the Python and C++ versions of VIGRA to be critical in order to prevent subtle errors when porting from one language to the other, so we went with the convention described.
   
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
   

.. _sec-dtype-coercion:

Point Operators and Type Coercion
---------------------------------

Vigra images and volumes support all arithmetic and algebraic functions defined in  
`numpy.ufunc <http://docs.scipy.org/doc/numpy/reference/ufuncs.html#available-ufuncs>`_. 
As usual, these operations are applied independently at each pixel.
However, vigranumpy overloads numpy.ufunc in a way that makes operator behavior 
more suitable for image analysis. In particular, we changed two aspects:

* The memory layout of the input arrays is preserved in the result arrays. 
  In contrast, plain numpy.ufuncs always create C-order arrays, even if 
  the inputs have a different order (e.g. as Fortran-order). 

* The value types of result arrays (i.e. their 'dtype') are determined in a way 
  that is more suitable for image processing than the original numpy conversion rules.  

Array dtype conversion (aka coercion) is implemented by the function 
vigra.ufunc.Function.common_type according to the following coercion rules:

.. automethod::  vigra.ufunc.Function.common_type


Import and Export Functions
---------------------------

The module vigra.impex defines read and write functions for image and volume data. Note
that the contents of this module are automatically imported into the vigra module, so
you may call 'vigra.readImage(...)' instead of 'vigra.impex.readImage(...)' etc.

.. automodule:: vigra.impex
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
accessed at integer coordinates). Below, we describe SplineImageView3, but the module 
also defines SplineImageView0 ... SplineImageView5 (representing spline interpolation
of the respective order) which work analogously. See SplineImageView_ 
in the C++ documentation for mode detailed information.

.. autoclass:: vigra.sampling.SplineImageView3
   :members:
.. automethod::  vigra.sampling.SplineImageView3.__init__


Image Analysis
--------------

The module vigra.analysis contains segmentation algorithms (e.g. watershed), edge and 
corner detection, localization of maxima and minima etc.

.. automodule:: vigra.analysis
   :members:


Machine Learning
----------------

The module vigra.learning will eventually provide a wide range of machine learning 
tools. Right now, it only contains an implementation of the random forest classifier.

.. automodule:: vigra.learning
   :members:

.. autoclass:: vigra.learning.RandomForest
   :members:
   
For more information, refer to RandomForest_ in the C++ decomentation.

.. autoclass:: vigra.learning.RandomForest_new
   :members:

Noise Normalization
-------------------

The module vigra.noise provides noise estimation and normalization according to a
method proposed by Förstner.

.. automodule:: vigra.noise
   :members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
