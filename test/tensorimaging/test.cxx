/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */                
/*                                                                      */
/************************************************************************/

/* tensorimaging tests */

#include <iostream>
#include <functional>
#include <cmath>
#include "unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/impex.hxx"
#include "vigra/tensorutilities.hxx"
#include "vigra/orientedtensorfilters.hxx"
#include "vigra/boundarytensor.hxx"
#include "vigra/gradient_energy_tensor.hxx"

using namespace vigra;

struct TensorUtilityTest
{
    typedef vigra::DImage Image;
    typedef vigra::DVector2Image V2Image;
    typedef vigra::DVector3Image V3Image;

    TensorUtilityTest()
    : w(5), h(5),
      img2(w, h), img3(w, h)
    {
        for(int y = 0; y < h; ++y)
        {
            for(int x = 0; x < w; ++x)
            {
                img2(x,y)[0] = x;
                img2(x,y)[1] = y;
            }
        }
        vectorToTensor(srcImageRange(img2), destImage(img3));
    }

    void vector2TensorTest()
    {
        for(int y = 0; y < h; ++y)
        {
            for(int x = 0; x < w; ++x)
            {
                shouldEqual(img3(x,y)[0], x*x);
                shouldEqual(img3(x,y)[1], -x*y);
                shouldEqual(img3(x,y)[2], y*y);
            }
        }
    }
    
    void tensorEigenRepresentationTest()
    {
        V3Image res(img3.size());
        tensorEigenRepresentation(srcImageRange(img3), destImage(res));
        
        for(int y = 0; y < h; ++y)
        {
            for(int x = 0; x < w; ++x)
            {
                shouldEqual(res(x,y)[0]+res(x,y)[1], img3(x,y)[0] + img3(x,y)[2]);
                shouldEqual(res(x,y)[0]*res(x,y)[1], 
                            img3(x,y)[0]*img3(x,y)[2]-sq(img3(x,y)[1]));
                if(x == 0 && y == 0)
                {
                    shouldEqual(res(x,y)[2], 0.0);
                }
                else
                {
                    shouldEqualTolerance(res(x,y)[2], 
                                    VIGRA_CSTD::atan2((double)-y, (double)x), 1e-12);
                }
            }
        }
    }
    
    void tensorTraceTest()
    {
        Image res(img3.size());
        tensorTrace(srcImageRange(img3), destImage(res));
        
        for(int y = 0; y < h; ++y)
        {
            for(int x = 0; x < w; ++x)
            {
                shouldEqual(res(x,y), img3(x,y)[0] + img3(x,y)[2]);
            }
        }
    }
    
    
    void tensorToEdgeCornerTest()
    {
        V2Image edge(img3.size());
        Image corner(img3.size());
        tensorToEdgeCorner(srcImageRange(img3), destImage(edge), destImage(corner));
        V3Image eigen(img3.size());
        tensorEigenRepresentation(srcImageRange(img3), destImage(eigen));
        
        for(int y = 0; y < h; ++y)
        {
            for(int x = 0; x < w; ++x)
            {
                shouldEqual(edge(x,y)[0], eigen(x,y)[0]-eigen(x,y)[1]);
                shouldEqual(edge(x,y)[1], eigen(x,y)[2]);
                shouldEqual(corner(x,y), 2.0*eigen(x,y)[1]);
            }
        }
    }
    
    int w, h;
    V2Image img2;
    V3Image img3;    
};

#define defImpulseResponseTest(xorder, yorder) \
    void rieszTransform##xorder##yorder##Test() \
    { \
        Image res(img1.size()), ref(img1.size()); \
        ImageImportInfo iref("riesz" #xorder #yorder ".xv"); \
        importImage(iref, destImage(ref)); \
        \
        rieszTransformOfLOG(srcImageRange(img1), destImage(res), 2.0, xorder, yorder); \
        \
        shouldEqualSequenceTolerance(res.begin(), res.end(), ref.begin(), 1e-12); \
    } \
    


struct EdgeJunctionTensorTest
{
    typedef vigra::DImage Image;
    typedef vigra::DVector2Image V2Image;
    typedef vigra::DVector3Image V3Image;

    EdgeJunctionTensorTest()
    : img1(17, 17)
    {
        img1.init(0.0);
        img1(8,8) = 1.0;

        ImageImportInfo i2("l2.xv");
        img2.resize(i2.size());
        importImage(i2, destImage(img2));
    }

    defImpulseResponseTest(0, 0)
    defImpulseResponseTest(1, 0)
    defImpulseResponseTest(0, 1)
    defImpulseResponseTest(2, 0)
    defImpulseResponseTest(1, 1)
    defImpulseResponseTest(0, 2)

    void boundaryTensorTest0()
    {
        V3Image res(img1.size()), ref(img1.size());
        ImageImportInfo iref("boundaryTensor.xv");
        importImage(iref, destImage(ref));
        
        boundaryTensor(srcImageRange(img1), destImage(res), 2.0);
        
        for(V3Image::iterator i = res.begin(), j = ref.begin(); i < res.end(); ++i, ++j)
        {
            shouldEqualTolerance(i->magnitude(), j->magnitude(), 1e-12);
        }
    }

    void boundaryTensorTest1()
    {
        V3Image bt(img2.size());
        Image res(img2.size()), ref(img2.size());
        ImageImportInfo iref("l2_boundary1.xv");
        importImage(iref, destImage(ref));

        boundaryTensor1(srcImageRange(img2), destImage(bt), 2.0);
        tensorTrace(srcImageRange(bt), destImage(res));
        
        shouldEqualSequenceTolerance(res.begin(), res.end(), ref.begin(), 1e-12);
    }

    void boundaryTensorTest2()
    {
        V3Image bt(img2.size());
        Image res(img2.size()), ref(img2.size());
        ImageImportInfo iref("l2_boundary.xv");
        importImage(iref, destImage(ref));

        boundaryTensor(srcImageRange(img2), destImage(bt), 2.0);
        tensorTrace(srcImageRange(bt), destImage(res));
        
        shouldEqualSequenceTolerance(res.begin(), res.end(), ref.begin(), 1e-12);
    }

    void boundaryTensorTest3()
    {
        // does not produce the correct result
        Image res(img2.size());
        {
            V3Image bt(img2.size()), bt2(img2.size());

            boundaryTensor3(srcImageRange(img2), destImage(bt), destImage(bt2), 1.0);
            combineTwoImages(srcImageRange(bt), srcImage(bt2), destImage(bt2), 
                         std::plus<V3Image::value_type>());
            tensorTrace(srcImageRange(bt2), destImage(res));
        }
        Image ref(img2.size());
        ImageImportInfo iref("l2_boundary3.xv");
        importImage(iref, destImage(ref));

        shouldEqualSequenceTolerance(res.begin(), res.end(), ref.begin(), 1e-12);
    }
    
    void hourglassTest()
    {
        V2Image gradient(img2.size());
        V3Image tensor(img2.size()), smoothedTensor(img2.size());
        Image res(img2.size()), ref(img2.size());
        ImageImportInfo iref("l2_hourglass.xv");
        importImage(iref, destImage(ref));

        gaussianGradient(srcImageRange(img2), destImage(gradient), 0.7);
        vectorToTensor(srcImageRange(gradient), destImage(tensor));
        hourGlassFilter(srcImageRange(tensor), destImage(smoothedTensor), 2.8, 0.4);
        tensorTrace(srcImageRange(smoothedTensor), destImage(res));
  
        shouldEqualSequenceTolerance(res.begin(), res.end(), ref.begin(), 1e-12);
    }

    void energyTensorTest()
    {
        using namespace functor;
        
        V3Image get(img2.size());
        Image res(img2.size()), ref(img2.size());
        ImageImportInfo iref("l2_get.xv");
        importImage(iref, destImage(ref));

        Kernel1D<double> smooth, grad;
        smooth.initGaussian(1.0);
        grad.initGaussianDerivative(1.0, 1);
        gradientEnergyTensor(srcImageRange(img2), destImage(get), grad, smooth);
        tensorTrace(srcImageRange(get), destImage(res));

        combineTwoImages(srcImageRange(res), srcImage(ref), destImage(res),
                         Arg1() - Arg2());

        Image::iterator i = res.begin(), end = res.end();        
        for(; i != end; ++i)
            shouldEqualTolerance(*i, 0.0, 1e-12);
    }

    Image img1, img2;
};

struct TensorTestSuite
: public vigra::test_suite
{
    TensorTestSuite()
    : vigra::test_suite("TensorTestSuite")
    {
        add( testCase( &TensorUtilityTest::vector2TensorTest));
        add( testCase( &TensorUtilityTest::tensorEigenRepresentationTest));
        add( testCase( &TensorUtilityTest::tensorTraceTest));
        add( testCase( &TensorUtilityTest::tensorToEdgeCornerTest));

        add( testCase( &EdgeJunctionTensorTest::rieszTransform00Test));
        add( testCase( &EdgeJunctionTensorTest::rieszTransform10Test));
        add( testCase( &EdgeJunctionTensorTest::rieszTransform01Test));
        add( testCase( &EdgeJunctionTensorTest::rieszTransform20Test));
        add( testCase( &EdgeJunctionTensorTest::rieszTransform11Test));
        add( testCase( &EdgeJunctionTensorTest::rieszTransform02Test));
        add( testCase( &EdgeJunctionTensorTest::boundaryTensorTest0));
        add( testCase( &EdgeJunctionTensorTest::boundaryTensorTest1));
        add( testCase( &EdgeJunctionTensorTest::boundaryTensorTest2));
        add( testCase( &EdgeJunctionTensorTest::hourglassTest));
        add( testCase( &EdgeJunctionTensorTest::energyTensorTest));
    }
};

int main()
{
    TensorTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

