VIGRA Computer Vision Library
=============================

[![Build Status](https://travis-ci.org/ukoethe/vigra.png?branch=master)](https://travis-ci.org/ukoethe/vigra)

                Copyright 1998-2013 by Ullrich Koethe

    This file is part of the VIGRA computer vision library.
    You may use, modify, and distribute this software according
    to the terms stated in the LICENSE.txt file included in
    the VIGRA distribution.

    The VIGRA Website is
        http://hci.iwr.uni-heidelberg.de/vigra/                       
    Please direct questions, bug reports, and contributions to        
        ullrich.koethe@iwr.uni-heidelberg.de    or                    
        vigra@informatik.uni-hamburg.de                               


    THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR
    IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
    WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.


Installation
------------

Installation instructions can be found in the file 
  $VIGRA_PATH/doc/vigra/Installation.html
If the documentation has not yet been generated (e.g. when you build from a development 
snapshot), you find these instructions in
  $VIGRA_PATH/docsrc/installation.dxx
or online at
  http://hci.iwr.uni-heidelberg.de/vigra/doc/vigra/Installation.html

Documentation
-------------

If you downloaded an official release, the documentation can be found in $VIGRA_PATH/doc/vigra/, the start file 
is $VIGRA_PATH/doc/vigra/index.html. Online documentation for the latest release is at 
http://hci.iwr.uni-heidelberg.de/vigra/doc/vigra/index.html

When you use the development version from github, you can generate documentation by `make doc`. Up-to-date 
online documentation for the 'master' branch is at http://ukoethe.github.io/vigra/doc/vigra/index.html

Download
--------

VIGRA can be downloaded at http://hci.iwr.uni-heidelberg.de/vigra/#download The official development 
repository is at https://github.com/ukoethe/vigra

What is VIGRA
-------------

VIGRA is a computer vision library that puts its main emphasis on flexible algorithms, because algorithms represent the principle know-how of this field. The library was consequently built using generic programming as introduced by Stepanov and Musser and exemplified in the C++ Standard Template Library. By writing a few adapters (image iterators and accessors) you can use VIGRA's algorithms on top of your data structures, within your environment. Alternatively, you can also use the data structures provided within VIGRA, which can be easily adapted to a wide range of applications. VIGRA's flexibility comes almost for free: Since the design uses compile-time polymorphism (templates), performance of the compiled program approaches that of a traditional, hand tuned, inflexible, solution.



