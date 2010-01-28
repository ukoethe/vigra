.. vigranumpy documentation master file, created by
   sphinx-quickstart on Thu Dec 03 15:51:41 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to vigranumpy's documentation!
======================================

Contents:

.. toctree::
   :maxdepth: 2
   
Image and Volume Data Structures
--------------------------------

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

Core Image Processing and Analysis Functions
---------------------------------------

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

