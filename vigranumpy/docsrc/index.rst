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

Vigranumpy exports the functionality of the C++ image processing library `VIGRA <../vigra/index.html>`_ to Python. It is based on the popular `numpy <http://numpy.scipy.org/>`_ module and uses its ndarray data structure to represent image and volume data. Thus, it is fully interoperable with existing numpy functionality, including various tools for image display such as matplotlib. 

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

This convention differs from the `Python Imaging Library <http://www.pythonware.com/products/pil/>`_ and Matlab, where the spatial indices must be given in reverse order (e.g. scalarImage[y, x]). Either convention has advantages and disadvantages. In the end, we considered compatibility between the Python and C++ versions of VIGRA to be critical in order to prevent subtle errors when porting from one language to the other, so we went with the convention described.
   
Image and Volume Data Structures
--------------------------------

Vigranumpy can work directly on numpy.ndarrays. However, plain ndarrays do not carry
any information about the semantics of the different coordinate axes. For example,
one cannot distinguish a 2-dimensional RGB image from a scalar volume data set that
happens to contain only three slices. In order to distinguish between arrays that 
have the same structure but different interpretation, vigra.arraytypes provides the 
following array classes::

    Image
        ScalarImage
        Vector2Image
        Vector3Image
        Vector4Image
        RGBImage
    Volume
        ScalarVolume
        Vector2Volume
        Vector3Volume
        Vector4Volume
        Vector6Volume
        RGBVolume

with the obvious inheritance relationships. Below, we describe Image, ScalarImage, 
and RGBImage in detail, the other classes work analogously. The new array classes 
serve several purposes:

* The semantic interpretation improves code readability.

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
  and 'imshow' (for matplotlib-based image display) so that they work in the 
  expected way (i.e. images are printed in horizontal scan order and are displayed 
  upright). Note that other Python imaging modules (such as PIL) use a different
  indexing convention (namely image[y, x, channel]).

* vigra.arraytypes and vigra.ufunc overload numpy.ufunc (i.e. basic mathematical 
  functions for arrays) so that the memory layout of the input arrays is preserved 
  in the result (whereas plain numpy.ufuncs always create C-order arrays, even if 
  the inputs have Fortran-order). vigra.ufunc also implements array dtype coercion 
  in a way that is more suitable for image processing than the original coercion. 
  See :ref:`sec-dtype-coercion` for details.

----------------

.. autoclass:: vigra.Image
   :show-inheritance:
   :members:
  
   .. attribute:: channels
   
      the number of channels of this image (e.g. RGB has three channels)

   .. attribute:: spatialDimensions
   
      number of spatial dimensions (always 2 for images). This is useful for 
      distinguishing RGBImage from ScalarVolume in overload resolution.

-------------

.. autoclass:: vigra.ScalarImage
   :show-inheritance:
   :members:

-------------

.. autoclass:: vigra.RGBImage
   :show-inheritance:
   :members:
   

.. _sec-dtype-coercion:

Type Coercion in Point Operators
--------------------------------

When an arithmetic or algebraic function is called for an image (or
set of images), it is applied to each pixel separately. The value
type of the results (i.e. the result array's 'dtype') is automatically
determined by the function vigra.ufunc.Function.common_type according 
to the following coercion rules:

.. automethod::  vigra.ufunc.Function.common_type

Core Image Processing and Analysis Functions
--------------------------------------------

.. automodule:: vigra
   :members:

Common Filters
--------------

.. automodule:: vigra.filters
   :members:

Convolution Functions
---------------------

.. automodule:: vigra.convolution
   :members:

Basic Morphological Operations
------------------------------

.. automodule:: vigra.morphology
   :members:

Tensor Image Processing
-----------------------

.. automodule:: vigra.tensor
   :members:

Region Segmentation Algorithms
------------------------------

.. automodule:: vigra.segmentation
   :members:

Noise Normalization
-------------------

.. automodule:: vigra.noise
   :members:

Edge and Corner Detection
-------------------------

.. automodule:: vigra.edgedetection
   :members:

Classification Functions
------------------------

.. automodule:: vigra.classification
   :members:

Import and export Functions
---------------------------

.. automodule:: vigra.impex
   :members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
