/* tensorimaging tests */

#include <iostream>
#include <functional>
#include <cmath>
#include "unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/impex.hxx"
#include "vigra/tensorutilities.hxx"
#include "vigra/boundarytensor.hxx"

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
                shouldEqual(res(x,y)[2], VIGRA_CSTD::atan2((double)-y, (double)x));
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
        shouldEqualSequence(res.begin(), res.end(), ref.begin()); \
    } \
    


struct BoundaryTensorTest
{
    typedef vigra::DImage Image;
    typedef vigra::DVector3Image VImage;

    BoundaryTensorTest()
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
        VImage res(img1.size()), ref(img1.size());
        ImageImportInfo iref("boundaryTensor.xv");
        importImage(iref, destImage(ref));
        
        boundaryTensor(srcImageRange(img1), destImage(res), 2.0);
        
        shouldEqualSequence(res.begin(), res.end(), ref.begin());
    }

    void boundaryTensorTest1()
    {
        VImage bt(img2.size());
        Image res(img2.size()), ref(img2.size());
        ImageImportInfo iref("l2_boundary1.xv");
        importImage(iref, destImage(ref));

        boundaryTensor1(srcImageRange(img2), destImage(bt), 2.0);
        tensorTrace(srcImageRange(bt), destImage(res));
        
        shouldEqualSequence(res.begin(), res.end(), ref.begin());
    }

    void boundaryTensorTest2()
    {
        VImage bt(img2.size());
        Image res(img2.size()), ref(img2.size());
        ImageImportInfo iref("l2_boundary.xv");
        importImage(iref, destImage(ref));

        boundaryTensor(srcImageRange(img2), destImage(bt), 2.0);
        tensorTrace(srcImageRange(bt), destImage(res));
        
        shouldEqualSequence(res.begin(), res.end(), ref.begin());
    }

    void boundaryTensorTest3()
    {
        // does not produce the correct result
        Image res(img2.size());
        {
            VImage bt(img2.size()), bt2(img2.size());

            boundaryTensor3(srcImageRange(img2), destImage(bt), destImage(bt2), 1.0);
            combineTwoImages(srcImageRange(bt), srcImage(bt2), destImage(bt2), 
                         std::plus<VImage::value_type>());
            tensorTrace(srcImageRange(bt2), destImage(res));
        }
        Image ref(img2.size());
        ImageImportInfo iref("l2_boundary3.xv");
        importImage(iref, destImage(ref));

        shouldEqualSequence(res.begin(), res.end(), ref.begin());
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

        add( testCase( &BoundaryTensorTest::rieszTransform00Test));
        add( testCase( &BoundaryTensorTest::rieszTransform10Test));
        add( testCase( &BoundaryTensorTest::rieszTransform01Test));
        add( testCase( &BoundaryTensorTest::rieszTransform20Test));
        add( testCase( &BoundaryTensorTest::rieszTransform11Test));
        add( testCase( &BoundaryTensorTest::rieszTransform02Test));
        add( testCase( &BoundaryTensorTest::boundaryTensorTest0));
        add( testCase( &BoundaryTensorTest::boundaryTensorTest1));
        add( testCase( &BoundaryTensorTest::boundaryTensorTest2));
    }
};

int main()
{
    TensorTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

