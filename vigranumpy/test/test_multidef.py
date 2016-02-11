#######################################################################
#                                                                      
#         Copyright 2015-2016 by Ullrich Koethe and Philip Schill      
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

from __future__ import division, print_function
import sys
print("\nexecuting test file", __file__, file=sys.stderr)
exec(compile(open('set_paths.py', "rb").read(), 'set_paths.py', 'exec'))

from nose.tools import assert_equal, raises, assert_raises
import vigra
import numpy as np

def checkAboutSame(i1,i2):
    assert(i1.shape==i2.shape)
    difference=np.sum(np.abs(i1-i2))/float(np.size(i1))
    assert(difference<5)

def checkEqual(i1, i2):
    assert(i1.shape==i2.shape)
    assert((i1==i2).all())

def test_convexHull():
    points = np.array([[0, 0], [2, 0], [2, 1], [0, 1], [1, 0.5]], dtype=np.float32)
    res = np.array([[0, 0], [0, 1], [2, 1], [2, 0], [0, 0]], dtype=np.float32)
    res = vigra.taggedView(res)
    hull = vigra.geometry.convexHull(points)
    checkAboutSame(hull, res)
    hull = vigra.geometry.convexHull(vigra.taggedView(points))
    checkAboutSame(res, hull)
    assert_raises(ValueError, vigra.geometry.convexHull, points.transpose())
    assert_raises(ValueError, vigra.geometry.convexHull, points.astype(np.uint8))
    assert_raises(ValueError, vigra.geometry.convexHull, 0, "a")

def test_convolveOneDimension():
    im = np.array([[0, 0, 0, 0, 0],
                   [0, 0, 1, 0, 0],
                   [0, 0, 0, 0, 0]], dtype=np.float64)
    imc = np.array([[0, 0, 0, 0, 0],
                    [0, 1, 2, 3, 0],
                    [0, 0, 0, 0, 0]], dtype=np.float64)
    k = vigra.filters.explictKernel(-1, 1, np.array([1, 2, 3], dtype=np.float64))
    res = vigra.filters.convolveOneDimension(im, 0, k)
    checkAboutSame(res, imc)
    assert_raises(ValueError, vigra.filters.convolveOneDimension, im)
    assert_raises(ValueError, vigra.filters.convolveOneDimension, im, 0)
    assert_raises(ValueError, vigra.filters.convolveOneDimension, im.astype(np.uint8), 0, k)
    assert_raises(ValueError, vigra.filters.convolveOneDimension, [0, 1], 0, k)

def test_convolve():
    # Test convolve with a 2D kernel.
    im = np.array([[0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0],
                   [0, 0, 1, 0, 0],
                   [0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0]], dtype=np.float64)
    imc0 = np.array([[0, 0, 0, 0, 0],
                     [0, 1, 2, 3, 0],
                     [0, 4, 5, 6, 0],
                     [0, 7, 8, 9, 0],
                     [0, 0, 0, 0, 0]], dtype=np.float64)
    k = vigra.filters.Kernel2D()
    k.initExplicitly((-1, -1), (1, 1), np.array([[1, 2, 3],
                                                 [4, 5, 6],
                                                 [7, 8, 9]], dtype=np.float64))
    res = vigra.filters.convolve(im, k)
    checkAboutSame(res, imc0)
    assert_raises(ValueError, vigra.filters.convolve, im)
    assert_raises(ValueError, vigra.filters.convolve, im.astype(np.uint8), k)
    assert_raises(ValueError, vigra.filters.convolve, [0, 1], k)

    # Test convolve with a 1D kernel.
    imc1 = np.array([[0, 0, 0, 0, 0],
                     [0, 1, 2, 3, 0],
                     [0, 2, 4, 6, 0],
                     [0, 3, 6, 9, 0],
                     [0, 0, 0, 0, 0]], dtype=np.float64)
    k = vigra.filters.Kernel1D()
    k.initExplicitly(-1, 1, np.array([1, 2, 3], dtype=np.float64))
    res = vigra.filters.convolve(im, k)
    checkAboutSame(res, imc1)

    # Test convolve with two 1D kernels.
    imc2 = np.array([[0, 0, 0, 0, 0],
                     [0, 0.5, 1, 2, 0],
                     [0, 1, 2, 4, 0],
                     [0, 1.5, 3, 6, 0],
                     [0, 0, 0, 0, 0]], dtype=np.float64)
    k0 = vigra.filters.explictKernel(-1, 1, np.array([1, 2, 3], dtype=np.float64))
    k1 = vigra.filters.explictKernel(-1, 1, np.array([0.5, 1, 2], dtype=np.float64))
    res = vigra.filters.convolve(im, (k0, k1))
    checkAboutSame(res, imc2)

def test_gaussianSmoothing():
    im = np.array([[0, 0, 0],
                   [0, 1, 0],
                   [0, 0, 0]], dtype=np.float64)
    imc = np.array([[ 0.04532707,  0.16757448,  0.04532707],
                    [ 0.16757448,  0.61952398,  0.16757448],
                    [ 0.04532707,  0.16757448,  0.04532707]])
    res = vigra.filters.gaussianSmoothing(im, 0.5)
    checkAboutSame(res, imc)
    assert_raises(ValueError, vigra.filters.gaussianSmoothing, im)
    assert_raises(ValueError, vigra.filters.gaussianSmoothing, im.astype(np.int8), 0.5)

def test_laplacianOfGaussian():
    im = np.zeros((6, 5), dtype=np.float64)
    im[2, 2] = 1
    imc = np.array([[ 0.07002893,  0.08373281,  0.08667994,  0.08373281,  0.07002893],
                    [ 0.08373281,  0.0174964 , -0.08323153,  0.0174964 ,  0.08373281],
                    [ 0.08646211, -0.0837286 , -0.31618431, -0.0837286 ,  0.08646211],
                    [ 0.07846283,  0.00859282, -0.09564048,  0.00859282,  0.07846283],
                    [ 0.0352323 ,  0.04236348,  0.04414477,  0.04236348,  0.0352323 ],
                    [ 0.01053995,  0.01780715,  0.0248179 ,  0.01780715,  0.01053995]])
    res = vigra.filters.laplacianOfGaussian(im)
    checkAboutSame(res, imc)
    assert_raises(ValueError, vigra.filters.laplacianOfGaussian, im.astype(np.uint8))

def test_gaussianDivergence():
    im = np.zeros((5, 5, 2), dtype=np.float64)
    im[2, 2] = [0, 1]
    imc = np.array([[  0.,   2.47013247e-02,   5.96311195e-19, -2.47013247e-02,   0.],
                    [  0.,   5.63656325e-02,   1.35525272e-18, -5.63656325e-02,   0.],
                    [  0.,   9.12597369e-02,   2.16840434e-18, -9.12597369e-02,   0.],
                    [  0.,   5.63656325e-02,   1.35525272e-18, -5.63656325e-02,   0.],
                    [  0.,   2.47013247e-02,   5.96311195e-19, -2.47013247e-02,   0.]])
    res = vigra.filters.gaussianDivergence(im)
    checkAboutSame(res, imc)
    assert_raises(ValueError, vigra.filters.gaussianDivergence, np.zeros((5, 6, 3), dtype=np.float64))
    assert_raises(ValueError, vigra.filters.gaussianDivergence, np.zeros((10, 11, 12, 2), dtype=np.float64))

def test_multiBinaryErosionDilation():
    im = np.array([[0, 0, 0, 0, 0],
                   [0, 1, 1, 0, 0],
                   [0, 0, 1, 0, 0],
                   [0, 0, 0, 0, 0]], dtype=np.uint8)
    imc = np.array([[0, 1, 1, 0, 0],
                    [1, 1, 1, 1, 0],
                    [0, 1, 1, 1, 0],
                    [0, 0, 1, 0, 0]], dtype=np.uint8)
    imc2 = np.array([[1, 1, 1, 0, 0],
                    [1, 1, 1, 1, 0],
                    [0, 1, 1, 1, 0],
                    [0, 0, 1, 0, 0]], dtype=np.uint8)
    res = vigra.filters.multiBinaryErosion(imc, 1)
    checkEqual(res, im)
    assert_raises(ValueError, vigra.filters.multiBinaryErosion, imc)
    assert_raises(ValueError, vigra.filters.multiBinaryErosion, imc.astype(np.float64))
    res = vigra.filters.multiBinaryDilation(im, 1)
    checkEqual(res, imc)
    assert_raises(ValueError, vigra.filters.multiBinaryDilation, im)
    assert_raises(ValueError, vigra.filters.multiBinaryDilation, im.astype(np.float64))

    res = vigra.filters.multiBinaryOpening(im, 1)
    checkEqual(res, np.zeros((4, 5), dtype=np.uint8))
    res = vigra.filters.multiBinaryOpening(imc, 1)
    checkEqual(res, imc)
    assert_raises(ValueError, vigra.filters.multiBinaryOpening, im)
    assert_raises(ValueError, vigra.filters.multiBinaryOpening, im.astype(np.float64))
    res = vigra.filters.multiBinaryClosing(im, 1)
    checkEqual(res, im)
    res = vigra.filters.multiBinaryClosing(imc, 1)
    checkEqual(res, imc2)
    assert_raises(ValueError, vigra.filters.multiBinaryClosing, im)
    assert_raises(ValueError, vigra.filters.multiBinaryClosing, im.astype(np.float64))

def test_multiGrayscaleErosionDilation():
    im = np.array([[0, 0, 0, 0, 0],
                   [0, 128, 100, 0, 0],
                   [0, 0, 128, 0, 0],
                   [0, 0, 0, 0, 0]], dtype=np.uint8)
    imc0 = np.array([[0, 0, 0, 0, 0],
                     [0, 81, 81, 0, 0],
                     [0, 0, 81, 0, 0],
                     [0, 0, 0, 0, 0]], dtype=np.uint8)
    imc1 = np.array([[ 56,  92,  64,  28,   0],
                     [ 92, 128, 100,  64,   0],
                     [ 56,  92, 128,  92,   0],
                     [  0,  56,  92,  56,   0]], dtype=np.uint8)
    imc2 = np.array([[0, 47, 19, 0, 0],
                     [47, 128, 100, 19, 0],
                     [0, 47, 128, 47, 0],
                     [0, 0, 47, 0, 0]], dtype=np.uint8)
    res = vigra.filters.multiGrayscaleErosion(im, 9)
    checkEqual(res, imc0)
    assert_raises(ValueError, vigra.filters.multiGrayscaleErosion, im)
    assert_raises(ValueError, vigra.filters.multiGrayscaleErosion, im.astype(np.int16), 9)
    res = vigra.filters.multiGrayscaleDilation(im, 6)
    checkEqual(res, imc1)
    assert_raises(ValueError, vigra.filters.multiGrayscaleDilation, im)
    assert_raises(ValueError, vigra.filters.multiGrayscaleDilation, im.astype(np.int16), 6)

    res = vigra.filters.multiGrayscaleOpening(im, 9)
    checkEqual(res, imc0)
    assert_raises(ValueError, vigra.filters.multiGrayscaleOpening, im)
    assert_raises(ValueError, vigra.filters.multiGrayscaleOpening, im.astype(np.int16), 9)
    res = vigra.filters.multiGrayscaleClosing(im, 9)
    checkEqual(res, imc2)
    assert_raises(ValueError, vigra.filters.multiGrayscaleClosing, im)
    assert_raises(ValueError, vigra.filters.multiGrayscaleClosing, im.astype(np.int16), 9)

def test_distanceTransform():
    # Test distanceTransform and distanceTransform2D.
    im = np.array([[0, 0, 0, 0, 0],
                   [0, 1, 1, 0, 0],
                   [0, 0, 1, 0, 0],
                   [0, 0, 0, 0, 0]], dtype=np.uint32)
    s2 = 1.41421354
    s5 = 2.23606801
    imd = np.array([[s2, 1, 1, s2, s5],
                    [1, 0, 0, 1, 2],
                    [s2, 1, 0, 1, 2],
                    [s5, s2, 1, s2, s5]], dtype=np.float32)
    res = vigra.filters.distanceTransform(im)
    checkAboutSame(res, imd)
    assert_raises(ValueError, vigra.filters.distanceTransform, np.zeros((5, 6, 7, 8), dtype=np.float32))

    # Test vectorDistanceTransform.
    im = im.astype(np.float32)
    imd = np.array([[[1, 1], [1, 0], [1, 0], [1, -1], [1, -2]],
                    [[0, 1], [0, 0], [0, 0], [0, -1], [0, -2]],
                    [[-1, 1], [0, 1], [0, 0], [0, -1], [0, -2]],
                    [[-1, 2], [-1, 1], [-1, 0], [-1, -1], [-1, -2]]], dtype=np.float32)
    res = vigra.filters.vectorDistanceTransform(im)
    checkAboutSame(res, imd)
    assert_raises(ValueError, vigra.filters.vectorDistanceTransform, np.zeros((5, 6, 7, 8), dtype=np.float32))
    
def test_boundaryDistanceTransform():
    # Test boundaryDistanceTransform.
    im = np.array([[2, 2, 2, 2, 2, 2],
                   [2, 2, 2, 2, 2, 2],
                   [1, 1, 1, 1, 2, 2],
                   [1, 1, 1, 1, 1, 2],
                   [1, 1, 1, 1, 1, 1],
                   [1, 1, 1, 1, 1, 1]], dtype=np.float32)
    s2 = 0.91421354
    s3 = 1.73606801
    s5 = 2.32842708
    s6 = 3.10555124
    imd = np.array([[1.5, 1.5, 1.5, 1.5, s3, s5],
                    [0.5, 0.5, 0.5, 0.5, s2, s3],
                    [0.5, 0.5, 0.5, 0.5, 0.5, s2],
                    [1.5, 1.5, 1.5, s2, 0.5, 0.5],
                    [2.5, 2.5, s5, s3, s2, 0.5],
                    [3.5, 3.5, s6, s5, s3, 1.5]], dtype=np.float32)
    res = vigra.filters.boundaryDistanceTransform(im)
    checkAboutSame(res, imd)
    assert_raises(ValueError, vigra.filters.boundaryDistanceTransform, im.astype(np.uint8))
    assert_raises(ValueError, vigra.filters.boundaryDistanceTransform, np.zeros((5, 6, 7, 8), dtype=np.float32))

    # Test boundaryVectorDistanceTransform.
    imd2 = np.array([[[ 1.5,  0. ], [ 1.5,  0. ], [ 1.5,  0. ], [ 1.5,  0. ], [ 1.5, -1. ], [ 1.5, -2. ]],
                     [[ 0.5,  0. ], [ 0.5,  0. ], [ 0.5,  0. ], [ 0.5,  0. ], [ 0.5, -1. ], [ 1. , -1.5]],
                     [[-0.5,  0. ], [-0.5,  0. ], [-0.5,  0. ], [-0.5,  0. ], [ 0. , -0.5], [ 0.5, -1. ]],
                     [[-1.5,  0. ], [-1.5,  0. ], [-1.5,  0. ], [-1. ,  0.5], [-0.5,  0. ], [ 0. , -0.5]],
                     [[-2.5,  0. ], [-2.5,  0. ], [-2. ,  1.5], [-1.5,  1. ], [-1. ,  0.5], [-0.5,  0. ]],
                     [[-3.5,  0. ], [-3.5,  0. ], [-2.5,  2. ], [-2. ,  1.5], [-1.5,  1. ], [-1.5,  0. ]]],
                    dtype=np.float32)
    res = vigra.filters.boundaryVectorDistanceTransform(im)
    checkAboutSame(res, imd2)
    assert_raises(ValueError, vigra.filters.boundaryVectorDistanceTransform, im.astype(np.uint8))
    assert_raises(ValueError, vigra.filters.boundaryVectorDistanceTransform, np.zeros((5, 6, 7, 8), dtype=np.float32))

def test_eccentricityTransform():
    # Test eccentricityTransform.
    im = np.array([[2, 2, 2, 3, 3, 3, 3],
                   [2, 1, 1, 1, 1, 1, 3],
                   [2, 1, 1, 1, 1, 1, 3],
                   [2, 1, 1, 1, 1, 1, 3],
                   [2, 2, 2, 3, 3, 3, 3]], dtype=np.float32)
    s2 = 1.41421354
    imd = np.array([[2, s2+1, s2+2, s2+3, s2+2, s2+1, 2],
                    [1, s2+1, s2, 1, s2, s2+1, 1],
                    [0, 2, 1, 0, 1, 2, 0],
                    [1, s2+1, s2, 1, s2, s2+1, 1],
                    [2, s2+1, s2+2, s2+3, s2+2, s2+1, 2]])
    res = vigra.filters.eccentricityTransform(im)
    checkAboutSame(res, imd)
    assert_raises(ValueError, vigra.filters.eccentricityTransform, im.astype(np.int32))
    assert_raises(ValueError, vigra.filters.eccentricityTransform, np.zeros((5, 6, 7, 8), dtype=np.float32))

    # Test eccentricityCenters.
    c = vigra.filters.eccentricityCenters(im)
    assert(c[1] == (2, 3) and c[2] == (2, 0) and c[3] == (2, 6))
    assert_raises(ValueError, vigra.filters.eccentricityCenters, im.astype(np.int32))
    assert_raises(ValueError, vigra.filters.eccentricityCenters, np.zeros((5, 6, 7, 8), dtype=np.float32))
    
    # Test eccentricityTransformWithCenters.
    res, c = vigra.filters.eccentricityTransformWithCenters(im)
    checkAboutSame(res, imd)
    assert(c[1] == (2, 3) and c[2] == (2, 0) and c[3] == (2, 6))
    assert_raises(ValueError, vigra.filters.eccentricityTransformWithCenters, im.astype(np.int32))
    assert_raises(ValueError, vigra.filters.eccentricityTransformWithCenters, np.zeros((5, 6, 7, 8), dtype=np.float32))

def test_gaussianGradient():
    im = np.array([[0, 0, 0, 0, 0],
                   [0, 1, 1, 0, 0],
                   [0, 1, 1, 1, 0],
                   [0, 0, 1, 1, 0],
                   [0, 0, 0, 0, 0]], dtype=np.float64)
    imd0 = np.array([[[ -8.01406117e-18 ,  8.02309608e-18],
                     [ -9.18662066e-18 ,  8.61150404e-02],
                     [ -7.83858551e-18 , -1.21313252e-01],
                     [ -2.95644367e-18 , -1.88750139e-01],
                     [ -5.94896474e-19 ,  0.00000000e+00]],
                    [[  8.61150404e-02 , -1.30104261e-17],
                     [  1.29233242e-01 ,  1.29233242e-01],
                     [  2.03991002e-01 , -1.01871715e-01],
                     [  2.15420146e-01 , -2.15420146e-01],
                     [  1.88750139e-01 ,  1.73472348e-18]],
                    [[ -1.21313252e-01 , -1.73472348e-18],
                     [ -1.01871715e-01 ,  2.03991002e-01],
                     [  2.77555756e-17 , -2.58040117e-17],
                     [  1.01871715e-01 , -2.03991002e-01],
                     [  1.21313252e-01 , -2.94902991e-17]],
                    [[ -1.88750139e-01 , -1.38777878e-17],
                     [ -2.15420146e-01 ,  2.15420146e-01],
                     [ -2.03991002e-01 ,  1.01871715e-01],
                     [ -1.29233242e-01 , -1.29233242e-01],
                     [ -8.61150404e-02 ,  8.67361738e-19]],
                    [[ -5.94896474e-19 ,  0.00000000e+00],
                     [ -2.95644367e-18 ,  1.88750139e-01],
                     [ -7.83858551e-18 ,  1.21313252e-01],
                     [ -9.18662066e-18 , -8.61150404e-02],
                     [ -8.01406117e-18 ,  8.02309608e-18]]], dtype=np.float64)
    imd1 = np.array([[  1.13399844e-17 , 8.61150404e-02 , 1.21313252e-01 , 1.88750139e-01 , 5.94896474e-19],
                     [  8.61150404e-02 , 1.82763404e-01 , 2.28013542e-01 , 3.04650092e-01 , 1.88750139e-01],
                     [  1.21313252e-01 , 2.28013542e-01 , 3.78974801e-17 , 2.28013542e-01 , 1.21313252e-01],
                     [  1.88750139e-01 , 3.04650092e-01 , 2.28013542e-01 , 1.82763404e-01 , 8.61150404e-02],
                     [  5.94896474e-19 , 1.88750139e-01 , 1.21313252e-01 , 8.61150404e-02 , 1.13399844e-17]], dtype=np.float64)
    res = vigra.filters.gaussianGradient(im, 1)
    checkAboutSame(res, imd0)
    assert_raises(ValueError, vigra.filters.gaussianGradient, im)
    assert_raises(ValueError, vigra.filters.gaussianGradient, im.astype(np.uint8), 1)
    assert_raises(ValueError, vigra.filters.gaussianGradient, np.zeros((5, 6, 7, 8, 9), dtype=np.float64), 1)

    res = vigra.filters.gaussianGradientMagnitude(im, 1)
    checkAboutSame(res, imd1)
    assert_raises(ValueError, vigra.filters.gaussianGradientMagnitude, im)
    assert_raises(ValueError, vigra.filters.gaussianGradientMagnitude, im.astype(np.uint8), 1)
    assert_raises(ValueError, vigra.filters.gaussianGradientMagnitude, np.zeros((5, 6, 7, 8, 9, 10), dtype=np.float64), 1)

def test_hessianOfGaussian():
    im = np.array([[0, 0, 0, 0, 0],
                   [0, 1, 1, 0, 0],
                   [0, 1, 1, 1, 0],
                   [0, 0, 1, 1, 0],
                   [0, 0, 0, 0, 0]], dtype=np.float64)
    imd = np.array([[[  2.03117070e-01,  4.84491752e-34,  2.03117070e-01],
                     [  2.65611971e-01, -1.09103007e-18, -9.96424231e-02],
                     [  3.34984252e-01,  3.98564636e-18, -2.12308959e-01],
                     [  2.93927033e-01,  4.46302120e-18,  9.19749106e-02],
                     [  2.36835873e-01, -2.40741243e-35,  2.36835873e-01]],
                    [[ -9.96424231e-02,  3.46944695e-18,  2.65611971e-01],
                     [ -8.92026671e-02,  7.56238735e-02, -8.92026671e-02],
                     [ -1.41639200e-02,  5.51365585e-02, -2.62273105e-01],
                     [  7.17062698e-02, -2.89764871e-02,  7.17062698e-02],
                     [  9.19749106e-02,  7.80625564e-18,  2.93927033e-01]],
                    [[ -2.12308959e-01,  0.00000000e+00,  3.34984252e-01],
                     [ -2.62273105e-01,  5.51365585e-02, -1.41639200e-02],
                     [ -3.06656412e-01,  1.30341282e-01, -3.06656412e-01],
                     [ -2.62273105e-01,  5.51365585e-02, -1.41639200e-02],
                     [ -2.12308959e-01, -2.60208521e-18,  3.34984252e-01]],
                    [[  9.19749106e-02, -7.58941521e-18,  2.93927033e-01],
                     [  7.17062698e-02, -2.89764871e-02,  7.17062698e-02],
                     [ -1.41639200e-02,  5.51365585e-02, -2.62273105e-01],
                     [ -8.92026671e-02,  7.56238735e-02, -8.92026671e-02],
                     [ -9.96424231e-02,  0.00000000e+00,  2.65611971e-01]],
                    [[  2.36835873e-01, -2.40741243e-35,  2.36835873e-01],
                     [  2.93927033e-01, -4.46302120e-18,  9.19749106e-02],
                     [  3.34984252e-01, -3.98564636e-18, -2.12308959e-01],
                     [  2.65611971e-01,  1.09103007e-18, -9.96424231e-02],
                     [  2.03117070e-01,  4.84491752e-34,  2.03117070e-01]]], dtype=np.float64)
    res = vigra.filters.hessianOfGaussian(im, 1)
    checkAboutSame(res, imd)
    assert_raises(ValueError, vigra.filters.hessianOfGaussian, im)
    assert_raises(ValueError, vigra.filters.hessianOfGaussian, im.astype(np.uint8), 1)
    assert_raises(ValueError, vigra.filters.hessianOfGaussian, np.zeros((5, 6, 7, 8, 9), dtype=np.float32), 1)

def test_structureTensor():
    im = np.array([[0, 0, 0],
                   [0, 1, 1]], dtype=np.float64)
    imd = np.array([[[0, 0, 0.00163116],
                     [0, 0, 0.01856249],
                     [0, 0, 0.00163116]],
                    [[0, 0, 0.01856249],
                     [0, 0, 0.21124014],
                     [0, 0, 0.01856249]]], dtype=np.float64)
    res = vigra.filters.structureTensor(im, 0.2, 0.4)
    checkAboutSame(imd, res)
    assert_raises(ValueError, vigra.filters.structureTensor, im)
    assert_raises(ValueError, vigra.filters.structureTensor, im, 0.2)
    assert_raises(ValueError, vigra.filters.structureTensor, im.astype(np.uint8), 0.2, 0.4)
    assert_raises(ValueError, vigra.filters.structureTensor, np.zeros((5, 6, 7, 8, 9, 10), dtype=np.float64), 0.5, 0.75)

def test_vectorToTensor():
    im = np.zeros((2, 3, 2), dtype=np.float32)
    im[0, 1] = [4, 5]
    im[1, 1] = [3, 7]
    imd = np.array([[[0, 0, 0],
                     [16, 20, 25],
                     [0, 0, 0]],
                    [[0, 0, 0],
                     [9, 21, 49],
                     [0, 0, 0]]], dtype=np.float32)
    res = vigra.filters.vectorToTensor(im)
    checkAboutSame(res, imd)
    vigra.filters.vectorToTensor(np.zeros((3, 4, 5, 3), dtype=np.float32))
    assert_raises(ValueError, vigra.filters.vectorToTensor, im.astype(np.uint32))
    assert_raises(ValueError, vigra.filters.vectorToTensor, np.zeros((2, 3, 4), dtype=np.float32))
    assert_raises(ValueError, vigra.filters.vectorToTensor, np.zeros((3, 4, 5, 4), dtype=np.float32))

def test_tensorTraceDetEigen():
    # Test tensorTrace.
    im = np.array([[0, 0, 0],
                   [0, 1, 1]], dtype=np.float64)
    im = vigra.filters.structureTensor(im, 0.2, 0.4)
    imd = np.array([[ 0.00163116,  0.01856249,  0.00163116],
                    [ 0.01856249,  0.21124014,  0.01856249]], dtype=np.float64)
    res = vigra.filters.tensorTrace(im)
    checkAboutSame(res, imd)
    assert_raises(ValueError, vigra.filters.tensorTrace, im.astype(np.uint32))
    assert_raises(ValueError, vigra.filters.tensorTrace, np.zeros((3, 4, 5), dtype=np.float64))

    # Test tensorDeterminant.
    vigra.filters.tensorDeterminant(im)
    assert_raises(ValueError, vigra.filters.tensorDeterminant, im.astype(np.uint32))
    assert_raises(ValueError, vigra.filters.tensorDeterminant, np.zeros((3, 4, 5), dtype=np.float64))

    # Test tensorEigenvalues.
    imd = np.array([[[ 0.00163116, 0 ],
                     [ 0.01856249, 0 ],
                     [ 0.00163116, 0 ]],
                    [[ 0.01856249, 0 ],
                     [ 0.21124014, 0 ],
                     [ 0.01856249, 0 ]]], dtype=np.float64)
    res = vigra.filters.tensorEigenvalues(im)
    checkAboutSame(res, imd)
    assert_raises(ValueError, vigra.filters.tensorEigenvalues, im.astype(np.uint32))
    assert_raises(ValueError, vigra.filters.tensorEigenvalues, np.zeros((3, 4, 5), dtype=np.float64))

def test_applyColortable():
    vigra.colors.applyColortable(np.zeros((10, 11), dtype=np.uint8), np.zeros((10, 4), dtype=np.uint8))
    assert_raises(ValueError, vigra.colors.applyColortable, np.zeros((10, 11), dtype=np.uint8))
    assert_raises(ValueError, vigra.colors.applyColortable, np.zeros((10, 11), dtype=np.uint8), np.zeros((10, 4), dtype=np.uint32))
    assert_raises(ValueError, vigra.colors.applyColortable, np.zeros((10, 11, 12), dtype=np.uint8), np.zeros((10, 4), dtype=np.uint8))

def test_linearRangeMapping():
    im = np.zeros((10, 11), dtype=np.uint8)
    im[5, 5] = 128
    vigra.colors.linearRangeMapping(im)
    im = np.zeros((10, 11, 12, 13), dtype=np.uint8)
    im[5, 5, 5, 6] = 128
    assert_raises(ValueError, vigra.colors.linearRangeMapping, im)

def test_labelImage():
    im = np.array([[1, 1, 2, 2],
                   [1, 0, 2, 2],
                   [0, 0, 1, 1]], dtype=np.uint8)
    imd0 = np.array([[1, 1, 3, 3],
                    [1, 2, 3, 3],
                    [2, 2, 4, 4]], dtype=np.uint8)
    imd1 = np.array([[1, 1, 2, 2],
                     [1, 0, 2, 2],
                     [0, 0, 3, 3]], dtype=np.uint8)
    res = vigra.analysis.labelImage(im)
    checkEqual(res, imd0)
    vigra.analysis.labelImage(im, "direct")
    assert_raises(ValueError, vigra.analysis.labelImage, im.astype(np.float64))
    assert_raises(ValueError, vigra.analysis.labelImage, np.zeros((3, 4, 5), dtype=np.uint8))
    res = vigra.analysis.labelImageWithBackground(im)
    checkEqual(res, imd1)
    vigra.analysis.labelImageWithBackground(im, "direct")
    assert_raises(ValueError, vigra.analysis.labelImageWithBackground, im.astype(np.float64))
    assert_raises(ValueError, vigra.analysis.labelImageWithBackground, np.zeros((3, 4, 5), dtype=np.uint8))

def test_labelVolume():
    im = np.zeros((4, 5, 6), dtype=np.uint8)
    vigra.analysis.labelVolume(im)
    vigra.analysis.labelVolume(im, "direct")
    assert_raises(ValueError, vigra.analysis.labelVolume, im.astype(np.float64))
    assert_raises(ValueError, vigra.analysis.labelVolume, np.zeros((3, 4, 5, 6), dtype=np.uint8))
    vigra.analysis.labelVolumeWithBackground(im)
    vigra.analysis.labelVolumeWithBackground(im, "direct")
    assert_raises(ValueError, vigra.analysis.labelVolumeWithBackground, im.astype(np.float64))
    assert_raises(ValueError, vigra.analysis.labelVolumeWithBackground, np.zeros((3, 4, 5, 6), dtype=np.uint8))

def test_labelMultiArray():
    im = np.zeros((4, 5, 6, 7), dtype=np.uint8)
    vigra.analysis.labelMultiArray(im)
    vigra.analysis.labelMultiArray(im, "direct")
    assert_raises(ValueError, vigra.analysis.labelMultiArray, im, "direct", 0)
    assert_raises(ValueError, vigra.analysis.labelMultiArray, im.astype(np.float64))
    vigra.analysis.labelMultiArrayWithBackground(im)
    vigra.analysis.labelMultiArrayWithBackground(im, "direct")
    vigra.analysis.labelMultiArrayWithBackground(im, "direct", 0)
    assert_raises(ValueError, vigra.analysis.labelMultiArrayWithBackground, im, "direct", 0, 1)
    assert_raises(ValueError, vigra.analysis.labelMultiArrayWithBackground, im.astype(np.float64))

def test_extendedLocalMinima():
    im = np.array([[5, 4, 5, 6],
                   [6, 3, 4, 5],
                   [7, 4, 5, 6]], dtype=np.uint8)
    imd = np.array([[0, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 0]], dtype=np.uint8)
    res = vigra.analysis.extendedLocalMinima(im, 1)
    checkEqual(res, imd)
    vigra.analysis.extendedLocalMinima(im, 1, 8)
    vigra.analysis.extendedLocalMinima(im.astype(np.float32))
    vigra.analysis.extendedLocalMinima(im.astype(np.float32), 1.5)
    assert_raises(ValueError, vigra.analysis.extendedLocalMinima, im)
    assert_raises(ValueError, vigra.analysis.extendedLocalMinima, im, 1, 8, 1)
    assert_raises(ValueError, vigra.analysis.extendedLocalMinima, im, 1.0)

def test_extendedLocalMinima3D():
    im = np.zeros((3, 4, 5), dtype=np.uint8)
    vigra.analysis.extendedLocalMinima3D(im)
    vigra.analysis.extendedLocalMinima3D(im, 5)
    vigra.analysis.extendedLocalMinima3D(im, 5, 6)
    assert_raises(ValueError, vigra.analysis.extendedLocalMinima3D, im, 1.5)
    assert_raises(ValueError, vigra.analysis.extendedLocalMinima3D, im, 5, 6, 1)

def test_watersheds():
    im_f = vigra.arraytypes.ScalarImage(np.random.rand(100, 200)*255, dtype=np.float32)
    im_i = vigra.arraytypes.ScalarImage(np.random.rand(100, 200)*255, dtype=np.uint32)

    vigra.analysis.watersheds(im_f)
    vigra.analysis.watersheds(im_f, seeds=im_i)
    vigra.analysis.watersheds(im_f, method="RegionGrowing")
    vigra.analysis.watersheds(im_f, seeds=im_i, method="RegionGrowing")
    assert_raises(ValueError, vigra.analysis.watersheds, im_f, seeds=im_i.astype(np.float32), method="RegionGrowing")
    assert_raises(ValueError, vigra.analysis.watersheds, im_f, seeds=np.zeros((100, 200, 3), dtype=np.uint32), method="RegionGrowing")
    assert_raises(ValueError, vigra.analysis.watersheds, np.zeros((100, 200, 3), dtype=np.float32), seeds=np.zeros((100, 200), dtype=np.uint32), method="RegionGrowing")

    vigra.analysis.watershedsNew(im_f)
    vigra.analysis.watershedsNew(im_f, seeds=im_i)
    vigra.analysis.watershedsNew(im_f, method="RegionGrowing")
    vigra.analysis.watershedsNew(im_f, seeds=im_i, method="RegionGrowing")
    assert_raises(ValueError, vigra.analysis.watershedsNew, im_f, seeds=im_i.astype(np.float32), method="RegionGrowing")
    assert_raises(ValueError, vigra.analysis.watershedsNew, im_f, seeds=np.zeros((100, 200, 3), dtype=np.uint32), method="RegionGrowing")
    assert_raises(ValueError, vigra.analysis.watershedsNew, np.zeros((100, 200, 3), dtype=np.float32), seeds=np.zeros((100, 200), dtype=np.uint32), method="RegionGrowing")

def test_superpixels():
    im_f = vigra.arraytypes.ScalarImage(np.random.rand(100, 200)*255, dtype=np.float32)
    im_i = vigra.arraytypes.ScalarImage(np.random.rand(100, 200)*255, dtype=np.uint32)

    vigra.analysis.slicSuperpixels(im_f, 0.5, 1)
    assert_raises(ValueError, vigra.analysis.slicSuperpixels, im_f)
    assert_raises(ValueError, vigra.analysis.slicSuperpixels, im_f, 0.5)
    assert_raises(ValueError, vigra.analysis.slicSuperpixels, np.zeros((100, 200, 10, 20), dtype=np.float32), 0.5, 1)
