#include <iostream>
#include <functional>
#include "unittest.h"
#include "vigra/stdimage.hxx"
#include "vigra/labelimage.hxx"
#include "vigra/edgedetection.hxx"
#include "vigra/distancetransform.hxx"
#include "vigra/localminmax.hxx"
#include "vigra/seededregiongrowing.hxx"
#include "vigra/cornerdetection.hxx"

struct LabelingTest
{
    typedef DImage Image;

    LabelingTest()
    : img1(5,5), img2(5,5), img3(9,5), img4(11,11)
    {
        static const double in1[] = { 0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.0, 1.0, 1.0, 1.0, 0.0,
                                      0.0, 1.0, 1.0, 1.0, 0.0,
                                      0.0, 1.0, 1.0, 1.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0};
        
        Image::ScanOrderIterator i = img1.begin();
        Image::ScanOrderIterator end = img1.end();
        Image::Accessor acc = img1.accessor();
        const double * p = in1;
        
        for(; i != end; ++i, ++p)
        {
            acc.set(*p, i);
        }

        static const double in2[] = { 0.0, 1.0, 0.0, 1.0, 0.0,
                                      1.0, 0.0, 1.0, 0.0, 1.0,
                                      0.0, 1.0, 0.0, 1.0, 0.0,
                                      1.0, 0.0, 1.0, 0.0, 1.0,
                                      0.0, 1.0, 0.0, 1.0, 0.0};
        
        i = img2.begin();
        end = img2.end();
        p = in2;
        
        for(; i != end; ++i, ++p)
        {
            acc.set(*p, i);
        }
        
        static const int in3[] = { 
            0, 1, 0, 1, 0, 1, 0, 1, 0,
            0, 1, 0, 1, 0, 1, 0, 0, 0,
            0, 1, 0, 1, 0, 0, 0, 1, 0,
            0, 1, 0, 0, 0, 1, 1, 1, 0,
            0, 0, 0, 1, 1, 1, 0, 0, 0
        };
        
        i = img3.begin();
        end = img3.end();
        const int * p1 = in3;
        
        for(; i != end; ++i, ++p1)
        {
            acc.set(*p1, i);
        }

        static const int spiral[] = { 
            1,1,1,1,1,1,1,1,1,1,1,
            1,2,2,2,2,2,2,2,2,2,1,
            1,2,1,1,1,1,1,1,1,2,1,
            1,2,1,2,2,2,2,2,1,2,1,
            1,2,1,2,1,1,1,2,1,2,1,
            1,2,1,2,1,2,1,2,1,2,1,
            1,2,1,2,1,2,2,2,1,2,1,
            1,2,1,2,1,1,1,1,1,2,1,
            1,2,1,2,2,2,2,2,2,2,1,
            1,2,1,1,1,1,1,1,1,1,1,
            1,2,2,2,2,2,2,2,2,2,2
        };
        
        i = img4.begin();
        end = img4.end();
        p1 = spiral;
        
        for(; i != end; ++i, ++p1)
        {
            acc.set(*p1, i);
        }
    }
    
    void labelingFourTest1()
    {
        Image res(img1);
        
        should(2 == labelImage(srcImageRange(img1), destImage(res), false));
        
        Image::ScanOrderIterator i1 = img1.begin();
        Image::ScanOrderIterator i1end = img1.end();
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = img1.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(acc(i1) == acc(i2) - 1.0);
        }
    }
    
    void labelingFourTest2()
    {
        Image res(img3);
        
        should(6 == labelImage(srcImageRange(img3), destImage(res), false));
        
        static const int target[] = { 
            1, 2, 1, 3, 1, 4, 1, 5, 1,
            1, 2, 1, 3, 1, 4, 1, 1, 1,
            1, 2, 1, 3, 1, 1, 1, 6, 1,
            1, 2, 1, 1, 1, 6, 6, 6, 1,
            1, 1, 1, 6, 6, 6, 1, 1, 1
        };
        
        Image::ScanOrderIterator i = res.begin();
        Image::ScanOrderIterator iend = res.end();
        const int * p = target;
        Image::Accessor acc = res.accessor();
        
        for(; i != iend; ++i, ++p)
        {
            should(acc(i) == *p);
        }
    }
    
    void labelingFourTest3()
    {
        Image res(img4.size());
        
        should(2 == labelImage(srcImageRange(img4), destImage(res), false));
        
        Image::ScanOrderIterator i = res.begin();
        Image::ScanOrderIterator iend = res.end();
        Image::ScanOrderIterator id = img4.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i != iend; ++i, ++id)
        {
            should(acc(i) == acc(id));
        }
    }
    
    void labelingFourTest4()
    {
        static const int data[] = {
            1,1,1,1,1,1,1,1,1,2,
            2,1,1,1,1,1,1,2,1,2,
            2,1,1,1,1,2,1,2,1,2,
            2,1,1,2,1,2,1,2,1,2,
            2,1,2,2,2,2,2,2,2,2,
            2,2,2,3,3,3,3,3,3,3
            
        };
        
        int w=10;
        int h=6;
        Image img(w,h), res(w,h);
        
        std::copy(data, data+w*h, img.begin());
        
        should(3 == labelImage(srcImageRange(img), destImage(res), false));
        
        Image::ScanOrderIterator i = res.begin();
        Image::ScanOrderIterator iend = res.end();
        Image::ScanOrderIterator id = img.begin();
        Image::Accessor acc = res.accessor();
        
        for(int c=0; i != iend; ++i, ++id, ++c)
        {
            should(acc(i) == acc(id));
        }
    }
    
    void labelingToCellGridTest()
    {
        Image tmp(img1);
        Image res(9, 9);
        
        should(2 == labelImage(srcImageRange(img1), destImage(tmp), false));
        
        regionImageToCellGridImage(srcImageRange(tmp), destImage(res), 0.0);
        
        static const double desired[] = { 
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  
               1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,  
               1.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 1.0,  
               1.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 1.0,  
               1.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 1.0,  
               1.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 1.0,  
               1.0, 0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 1.0,  
               1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,  
               1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        
        const double * i1 = desired;
        const double * i1end = i1 + 81;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = img1.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }
    }
    
    void labelingEightTest1()
    {
        Image res(img2);
        
        should(2 == labelImage(srcImageRange(img2), destImage(res), true));
        
        Image::ScanOrderIterator i1 = img2.begin();
        Image::ScanOrderIterator i1end = img2.end();
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = img2.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(acc(i1) == acc(i2) - 1.0);
        }
    }
    
    void labelingEightTest2()
    {
        static const int data[] = {
            1,1,1,1,1,1,1,1,1,1,
            1,1,1,2,3,3,2,1,1,1,
            1,1,2,3,3,2,4,2,1,1,
            1,2,3,3,2,4,4,4,2,1,
            1,2,3,2,4,4,2,3,2,1,
            1,2,3,3,2,2,3,3,2,1,
            1,1,2,3,3,3,3,2,1,1,
            1,1,1,2,2,2,2,1,1,1,
            1,1,1,1,1,1,1,1,1,1
        };
        
        Image img(10,9), res(10,9);
        
        std::copy(data, data+90, img.begin());
        
        should(4 == labelImage(srcImageRange(img), destImage(res), true));
        
        Image::ScanOrderIterator i = res.begin();
        Image::ScanOrderIterator iend = res.end();
        Image::Accessor acc = res.accessor();
        const int * p = data;
        
        for(; i != iend; ++i, ++p)
        {
            should(acc(i) == *p);
        }
    }
    
    void labelingFourWithBackgroundTest1()
    {
        Image res(img1);
        res = 0.0;
        
        should(1 == labelImageWithBackground(srcImageRange(img1), 
                                             destImage(res), 
                                             false, 0.0));
        
        Image::ScanOrderIterator i1 = img1.begin();
        Image::ScanOrderIterator i1end = img1.end();
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = img1.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(acc(i1) == acc(i2));
        }
    }
    
    void labelingFourWithBackgroundTest2()
    {
        Image res(img4);
        
        should(1 == labelImageWithBackground(srcImageRange(img4), 
                                            destImage(res), 
                                            false, 2.0));
        
        Image::ScanOrderIterator i = res.begin();
        Image::ScanOrderIterator iend = res.end();
        Image::ScanOrderIterator id = img4.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i != iend; ++i, ++id)
        {
            should(acc(i) == acc(id));
        }
    }
    
    void labelingEightWithBackgroundTest()
    {
        Image res(img2);
        
        should(1 == labelImageWithBackground(srcImageRange(img2), 
                                             destImage(res), 
                                             true, 0.0));
        
        Image::ScanOrderIterator i1 = img2.begin();
        Image::ScanOrderIterator i1end = img2.end();
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = img2.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(acc(i1) == acc(i2));
        }
    }
    
    Image img1, img2, img3, img4;
};
    
struct EdgeDetectionTest
{
    typedef DImage Image;

    EdgeDetectionTest()
    : img1(5,5), img2(9,11)
    {
        static const double in1[] = { 0.0, 0.0, 0.0, 0.0, 0.0,
                                      0.0, 1.0, 1.0, 1.0, 0.0,
                                      0.0, 1.0, 1.0, 1.0, 0.0,
                                      0.0, 1.0, 1.0, 1.0, 0.0,
                                      0.0, 0.0, 0.0, 0.0, 0.0};
        
        Image::ScanOrderIterator i = img1.begin();
        Image::ScanOrderIterator end = img1.end();
        Image::Accessor acc = img1.accessor();
        const double * p = in1;
        
        for(; i != end; ++i, ++p)
        {
            acc.set(*p, i);
        }

        static const double in2[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        
        i = img2.begin();
        end = img2.end();
        p = in2;
        
        for(; i != end; ++i, ++p)
        {
            acc.set(*p, i);
        }
    }
    
    void edgeDetectionTest()
    {
        Image res(img1);
        res = 0.0;
        
        differenceOfExponentialEdgeImage(srcImageRange(img1), destImage(res), 
            0.7, 0.1, 1.0);
        
        static const double desired[] = {0.0, 1.0, 1.0, 1.0, 0.0,
                                         1.0, 0.0, 0.0, 0.0, 1.0,
                                         1.0, 0.0, 0.0, 0.0, 1.0,
                                         1.0, 0.0, 0.0, 0.0, 1.0,
                                         0.0, 1.0, 1.0, 1.0, 0.0};
                                         
        const double * i1 = desired;
        const double * i1end = i1 + 25;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }
    }
    
    void edgeToCellGridTest()
    {
        img1.upperLeft()[Diff2D(2,2)] = 0.0;
        
        Image res(9,9);
        res = 0.0;
        
        differenceOfExponentialCellGridImage(srcImageRange(img1), 
                                                destImage(res), 
                                                0.7, 0.1, 1.0);
        
        static const double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        
        const double * i1 = desired;
        const double * i1end = i1 + 81;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }
    }
    
    void removeShortEdgesTest()
    {
        img1.upperLeft()[Diff2D(2,2)] = 0.0;
        
        Image res(9,9);
        res = 0.0;
        
        differenceOfExponentialCellGridImage(srcImageRange(img1), 
                                                destImage(res), 
                                                0.7, 0.1, 1.0);
        removeShortEdges(srcImageRange(res), 9, 0.0);
        
        static const double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        
        const double * i1 = desired;
        const double * i1end = i1 + 81;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }
    }
    
    void beautifyCellGridTest()
    {
        img1.upperLeft()[Diff2D(2,2)] = 0.0;
        
        Image res(9,9);
        res = 0.0;
        
        differenceOfExponentialCellGridImage(srcImageRange(img1), 
                                                destImage(res), 
                                                0.7, 0.1, 1.0);
        beautifyCellGridImage(srcImageRange(res), 1.0, 0.0);
        
        static const double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        
        const double * i1 = desired;
        const double * i1end = i1 + 81;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }
    }
    
    void closeGapsInCellGridTest()
    {
        closeGapsInCellGridImage(srcImageRange(img2), 1.0);
        
        static const double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        
        const double * i1 = desired;
        const double * i1end = i1 + 99;
        Image::ScanOrderIterator i2 = img2.begin();
        Image::Accessor acc = img2.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }
    }
    
    Image img1, img2;
};
    
struct DistanceTransformTest
{
    typedef DImage Image;

    DistanceTransformTest()
    : img(7,7)
    {
        static const double in[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        
        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator end = img.end();
        Image::Accessor acc = img.accessor();
        const double * p = in;
        
        for(; i != end; ++i, ++p)
        {
            acc.set(*p, i);
        }
    }
    
    void distanceTransformL1Test()
    {
        double epsilon = 0.000001;
        
        Image res(img);
        
        distanceTransform(srcImageRange(img), destImage(res), 0.0, 1);
        
        Image::Iterator i = res.upperLeft();
        Image::Accessor acc = res.accessor();
        int x,y;
        
        for(y=0; y<7; ++y)
        {
            for(x=0; x<7; ++x)
            {
                double dist = acc(i, Diff2D(x,y));
                double dist1 = abs(2.0 - x) + abs(2.0 - y);
                double dist2 = abs(5.0 - x) + abs(5.0 - y);
                double desired = (dist1 < dist2) ? dist1 : dist2;
                
                should(abs(dist - desired) < epsilon);
            }
        }
    }
    
    void distanceTransformL2Test()
    {
        double epsilon = 0.000001;
        
        Image res(img);
        
        distanceTransform(srcImageRange(img), destImage(res), 0.0, 2);
        
        Image::Iterator i = res.upperLeft();
        Image::Accessor acc = res.accessor();
        int x,y;
        
        for(y=0; y<7; ++y)
        {
            for(x=0; x<7; ++x)
            {
                double dist = acc(i, Diff2D(x,y));
                double dist1 = sqrt((2.0 - x)*(2.0 - x) +
                                         (2.0 - y)*(2.0 - y));
                double dist2 = sqrt((5.0 - x)*(5.0 - x) +
                                         (5.0 - y)*(5.0 - y));
                double desired = (dist1 < dist2) ? dist1 : dist2;
                
                should(abs(dist - desired) < epsilon);
            }
        }
    }
    
    void distanceTransformLInfTest()
    {
        double epsilon = 0.000001;
        
        Image res(img);
        
        distanceTransform(srcImageRange(img), destImage(res), 0.0, 0);
        
        Image::Iterator i = res.upperLeft();
        Image::Accessor acc = res.accessor();
        int x,y;
        
        for(y=0; y<7; ++y)
        {
            for(x=0; x<7; ++x)
            {
                double dist = acc(i, Diff2D(x,y));
                double dist1 = abs(2.0 - x) < abs(2.0 - y) ?
                                      abs(2.0 - y) : abs(2.0 - x);
                double dist2 = abs(5.0 - x) < abs(5.0 - y) ?
                                      abs(5.0 - y) : abs(5.0 - x);
                double desired = (dist1 < dist2) ? dist1 : dist2;
                
                should(abs(dist - desired) < epsilon);
            }
        }
    }
    
    
    Image img;
};

struct LocalMinMaxTest
{
    typedef DImage Image;

    LocalMinMaxTest()
    : img(9,9)
    {
        static const double in[] = {
            0.0,  0.1,  0.1,  0.3,  0.5,  0.0,  0.0,  0.0, 0.0,
            0.0, -0.1,  0.1,  0.0,  1.0,  0.0,  0.0,  0.0, 0.0,
            0.0,  0.5,  2.0,  0.0,  2.0,  2.0,  2.0,  0.0, 0.0,
            0.0,  0.0,  1.0,  1.5,  1.0,  1.0,  0.0,  0.0, 0.0,
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,
            0.0,  0.0,  0.0,  0.0, -1.0, -1.5, -1.0,  0.0, 0.0,
            0.0,  0.0, -2.0, -2.0, -2.0,  0.0, -2.0, -0.5, 0.0,
            0.0,  0.0,  0.0,  0.0, -1.0,  0.0, -0.1,  0.1, 0.0,
            0.0,  0.0,  0.0,  0.0, -0.5, -0.3, -0.1, -0.1, 0.0};
        
        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator end = img.end();
        Image::Accessor acc = img.accessor();
        const double * p = in;
        
        for(; i != end; ++i, ++p)
        {
            acc.set(*p, i);
        }
    }
    
    void localMinimumTest()
    {
        Image res(img);
        res = 0;
        
        localMinima(srcImageRange(img), destImage(res), 1.0);
        
        static const double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        
        const double * i1 = desired;
        const double * i1end = i1 + 81;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }
    }

    void localMaximumTest()
    {
        Image res(img);
        res = 0;
        
        localMaxima(srcImageRange(img), destImage(res), 1.0);
        
        static const double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        
        const double * i1 = desired;
        const double * i1end = i1 + 81;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }
    }

    void extendedLocalMinimumTest()
    {
        Image res(img);
        res = 0;
        
        extendedLocalMinima(srcImageRange(img), destImage(res), 1.0);
        
        static const double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        
        const double * i1 = desired;
        const double * i1end = i1 + 81;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }
    }

    void extendedLocalMaximumTest()
    {
        Image res(img);
        res = 0;
        
        extendedLocalMaxima(srcImageRange(img), destImage(res), 1.0);
        
        static const double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        
        const double * i1 = desired;
        const double * i1end = i1 + 81;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }
   }

    Image img;
};

struct RegionGrowingTest
{
    typedef DImage Image;

    RegionGrowingTest()
    : img(7,7), seeds(7,7)
    {
        static const double in[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        
        Image tmp(7,7);
        
        Image::ScanOrderIterator i = tmp.begin();
        Image::ScanOrderIterator end = tmp.end();
        Image::Accessor acc = tmp.accessor();
        const double * p = in;
        
        for(; i != end; ++i, ++p)
        {
            acc.set(*p, i);
        }
        
        distanceTransform(srcImageRange(tmp), destImage(img), 0.0, 2);
        
        seeds = 0;
        labelImageWithBackground(srcImageRange(tmp), destImage(seeds), 
                                 false, 0.0);
    }
    
    struct DirectCostFunctor
    {
        typedef double value_type;
        typedef double cost_type;

        void operator()(double const &) {}

        double const & cost(double const & v) const
        {
            return v;
        }
    };
    
    void voronoiTest()
    {
        
        double epsilon = 0.000001;
        
        Image res(img);
        
        ArrayOfRegionStatistics<DirectCostFunctor> cost(2);
        seededRegionGrowing(srcImageRange(img), srcImage(seeds), 
                            destImage(res), cost);
        
        Image::Iterator i = res.upperLeft();
        Image::Accessor acc = res.accessor();
        int x,y;
        
        for(y=0; y<7; ++y)
        {
            for(x=0; x<7; ++x)
            {
                double dist = acc(i, Diff2D(x,y));
                double dist1 = sqrt((2.0 - x)*(2.0 - x) +
                                         (2.0 - y)*(2.0 - y));
                double dist2 = sqrt((5.0 - x)*(5.0 - x) +
                                         (5.0 - y)*(5.0 - y));
                double desired = (dist1 <= dist2) ? 1 : 2;
                
                should(abs(dist - desired) < epsilon);
            }
        }
    }
    
    Image img, seeds;
};

struct CornerDetectionTest
{
    typedef DImage Image;

    CornerDetectionTest()
    : img(7,7)
    {
        static const double in[] = {
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0,
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
        
        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator end = img.end();
        Image::Accessor acc = img.accessor();
        const double * p = in;
        
        for(; i != end; ++i, ++p)
        {
            acc.set(*p, i);
        }
    }
    
    void cornerResponseFunctionTest()
    {
        Image tmp(img);
        Image res(img);
        res = 0.0;
        
        cornerResponseFunction(srcImageRange(img), destImage(tmp), 1.0);
        localMaxima(srcImageRange(tmp), destImage(res), 1.0);
        
        static const double desired[] = {
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                                         
        const double * i1 = desired;
        const double * i1end = i1 + 49;
        Image::ScanOrderIterator i2 = res.begin();
        Image::Accessor acc = res.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(*i1 == acc(i2));
        }
    }
    
    Image img;
};

struct SimpleAnalysisTestSuite
: public TestSuite
{
    SimpleAnalysisTestSuite()
    : TestSuite("SimpleAnalysisTestSuite")
    {
        add( testCase( &LabelingTest::labelingFourTest1));
        add( testCase( &LabelingTest::labelingFourTest2));
        add( testCase( &LabelingTest::labelingFourTest3));
        add( testCase( &LabelingTest::labelingFourTest4));
        add( testCase( &LabelingTest::labelingToCellGridTest));
        add( testCase( &LabelingTest::labelingEightTest1));
        add( testCase( &LabelingTest::labelingEightTest2));
        add( testCase( &LabelingTest::labelingFourWithBackgroundTest1));
        add( testCase( &LabelingTest::labelingFourWithBackgroundTest2));
        add( testCase( &LabelingTest::labelingEightWithBackgroundTest));
        add( testCase( &EdgeDetectionTest::edgeDetectionTest));
        add( testCase( &EdgeDetectionTest::edgeToCellGridTest));
        add( testCase( &EdgeDetectionTest::removeShortEdgesTest));
        add( testCase( &EdgeDetectionTest::beautifyCellGridTest));
        add( testCase( &EdgeDetectionTest::closeGapsInCellGridTest));
        add( testCase( &DistanceTransformTest::distanceTransformL1Test));
        add( testCase( &DistanceTransformTest::distanceTransformL2Test));
        add( testCase( &DistanceTransformTest::distanceTransformLInfTest));
        add( testCase( &LocalMinMaxTest::localMinimumTest));
        add( testCase( &LocalMinMaxTest::localMaximumTest));
        add( testCase( &LocalMinMaxTest::extendedLocalMinimumTest));
        add( testCase( &LocalMinMaxTest::extendedLocalMaximumTest));
        add( testCase( &RegionGrowingTest::voronoiTest));
        add( testCase( &CornerDetectionTest::cornerResponseFunctionTest));
    }
};

int main()
{
    SimpleAnalysisTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

