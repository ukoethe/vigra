#######################################################################
#                                                                      
#         Copyright 2009-2013 by Ullrich Koethe and Thorben Kroeger
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

import sys
print >> sys.stderr, "\nexecuting test file", __file__
execfile('set_paths.py')

from nose.tools import assert_equal, raises, assert_raises, assert_true
import vigra
import numpy

_have_qt = False
try:
    from PyQt4.QtGui import QImage
    from qimage2ndarray import byte_view
    _have_qt = True
except:
    pass

def test_gray2qimage():
    if not _have_qt:
        return
    
    a = numpy.random.random((100,200)).astype(numpy.float32)-0.5
    a[0,0] =-0.5 #make sure we get the correct bounds
    a[0,1] = 0.5
    vigra.impex.writeImage(a.swapaxes(0,1), "tmp1.png")

    img = QImage(a.shape[1], a.shape[0], QImage.Format_ARGB32_Premultiplied)
    n = numpy.asarray([-0.5, 0.5], dtype=numpy.float32)
    vigra.colors.gray2qimage_ARGB32Premultiplied(a, byte_view(img), n)
    img.save("tmp2.png")
    
    tmp1 = vigra.impex.readImage("tmp1.png")
    tmp2 = vigra.impex.readImage("tmp2.png")
    
    for i in range(3):
        assert_true( (tmp1 == tmp2[:,:,i]).all() )
    assert_true( (tmp2[:,:,3] == 255).all() )

def test_alphamodulated2qimage():
    if not _have_qt:
        return
    
    a = numpy.random.random((100,200)).astype(numpy.float32)-0.5
    a[0,0] =-0.5
    a[0,1] = 0.5
    vigra.impex.writeImage(a.swapaxes(0,1), "tmp1.png")
    
    tintColor = numpy.asarray([0.0, 1.0, 0.0], dtype=numpy.float32)
    
    img = QImage(a.shape[1], a.shape[0], QImage.Format_ARGB32_Premultiplied)
    n = numpy.asarray([-0.5, 0.5], dtype=numpy.float32)
    
    vigra.colors.alphamodulated2qimage_ARGB32Premultiplied(a, byte_view(img), tintColor, n)
    
    tmp1 = vigra.impex.readImage("tmp1.png").view(numpy.ndarray).squeeze().astype(numpy.uint8)
    tmp2 = byte_view(img)
    tmp2 = tmp2.swapaxes(0,1)
    
    assert_true( tmp1.shape[0:2] == tmp2.shape[0:2] )
    assert_true( tmp2.ndim == 3 )
    assert_true( tmp2.shape[2] == 4 )
    
    assert_true( (tmp2[:,:,0] == 0).all() )
    assert_true( (tmp2[:,:,2] == 0).all() )
    assert_true( (tmp2[:,:,3] == tmp1).all() )
    assert_true( (tmp2[:,:,1] == tmp1).all() )

if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__file__)
