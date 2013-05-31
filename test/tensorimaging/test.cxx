/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe                     */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
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
#include "vigra/multi_array.hxx"
#include "vigra/multi_math.hxx"

using namespace vigra;

struct TensorUtilityTest
{
    typedef vigra::MultiArray<2, double> Image;
    typedef vigra::MultiArray<2, TinyVector<double, 2> > V2Image;
    typedef vigra::MultiArray<2, TinyVector<double, 3> > V3Image;

    TensorUtilityTest()
    : w(5), h(5),
      img2(Shape2(w, h)), img3(Shape2(w, h))
    {
        for(int y = 0; y < h; ++y)
        {
            for(int x = 0; x < w; ++x)
            {
                img2(x,y)[0] = x;
                img2(x,y)[1] = y;
            }
        }
        vectorToTensor(img2, img3);
    }

    void vector2TensorTest()
    {
        for(int y = 0; y < h; ++y)
        {
            for(int x = 0; x < w; ++x)
            {
                shouldEqual(img3(x,y)[0], x*x);
                shouldEqual(img3(x,y)[1], x*y);
                shouldEqual(img3(x,y)[2], y*y);
            }
        }
    }
    
    void tensorEigenRepresentationTest()
    {
        V3Image res(img3.shape());
        tensorEigenRepresentation(img3, res);
        
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
                                    VIGRA_CSTD::atan2((double)y, (double)x), 1e-12);
                }
            }
        }
    }
    
    void tensorTraceTest()
    {
        Image res(img3.shape());
        tensorTrace(img3, res);
        
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
        V2Image edge(img3.shape());
        Image corner(img3.shape());
        tensorToEdgeCorner(img3, edge, corner);
        V3Image eigen(img3.shape());
        tensorEigenRepresentation(img3, eigen);
        
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
        \
        res.init(0.0); \
        rieszTransformOfLOG(View(img1), View(res), 2.0, xorder, yorder); \
        \
        shouldEqualSequenceTolerance(res.begin(), res.end(), ref.begin(), 1e-12); \
    } \
    


struct EdgeJunctionTensorTest
{
    typedef vigra::DImage Image;
    typedef vigra::MultiArrayView<2, double> View;
    typedef vigra::DVector2Image V2Image;
    typedef vigra::DVector3Image V3Image;
    typedef vigra::MultiArrayView<2, TinyVector<double, 3> > View3;

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
        V3Image res(img1.size()), res1(img1.size()), ref(img1.size());
        ImageImportInfo iref("boundaryTensor.xv");
        importImage(iref, destImage(ref));
        
        boundaryTensor(srcImageRange(img1), destImage(res), 2.0);
        boundaryTensor(View(img1), View3(res1), 2.0);
        
        for(V3Image::iterator i = res.begin(), j = res1.begin(), k = ref.begin(); i < res.end(); ++i, ++j, ++k)
        {
            shouldEqualTolerance(i->magnitude(), k->magnitude(), 1e-12);
            shouldEqualTolerance(j->magnitude(), k->magnitude(), 1e-12);
        }
    }

    void boundaryTensorTest1()
    {
        V3Image bt(img2.size()), bt1(img2.size());
        Image res(img2.size()), res1(img2.size()), ref(img2.size());
        ImageImportInfo iref("l2_boundary1.xv");
        importImage(iref, destImage(ref));

        boundaryTensor1(srcImageRange(img2), destImage(bt), 2.0);
        tensorTrace(srcImageRange(bt), destImage(res));
        
        shouldEqualSequenceTolerance(res.begin(), res.end(), ref.begin(), 1e-12);

        boundaryTensor1(View(img2), View3(bt1), 2.0);
        tensorTrace(View3(bt1), View(res1));
        
        shouldEqualSequenceTolerance(res1.begin(), res1.end(), ref.begin(), 1e-12);
    }

    void boundaryTensorTest2()
    {
        V3Image bt(img2.size());
        Image res(img2.size()), ref(img2.size());
        ImageImportInfo iref("l2_boundary.xv");
        importImage(iref, destImage(ref));

        boundaryTensor(View(img2), View3(bt), 2.0);
        tensorTrace(View3(bt), View(res));
        
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
        hourGlassFilter(View3(tensor), View3(smoothedTensor), 2.8, 0.4);
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
        gradientEnergyTensor(View(img2), View3(get), grad, smooth);
        tensorTrace(View3(get), View(res));

        using namespace multi_math;
        should(all(View(res) - View(ref) < 1e-12));
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

int main(int argc, char ** argv)
{
    TensorTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

