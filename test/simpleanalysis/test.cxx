#include <iostream>
#include <functional>
#include <cmath>
#include "unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/labelimage.hxx"
#include "vigra/edgedetection.hxx"
#include "vigra/distancetransform.hxx"
#include "vigra/localminmax.hxx"
#include "vigra/seededregiongrowing.hxx"
#include "vigra/cornerdetection.hxx"
#include "vigra/symmetry.hxx"

using namespace vigra;

struct LabelingTest
{
    typedef vigra::DImage Image;

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
    
    void labelingToCrackEdgeTest()
    {
        Image tmp(img1);
        Image res(9, 9);
        
        should(2 == labelImage(srcImageRange(img1), destImage(tmp), false));
        
        regionImageToCrackEdgeImage(srcImageRange(tmp), destImage(res), 0.0);
        
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
    typedef vigra::DImage Image;

    EdgeDetectionTest()
    : img1(5,5), img2(9,11), imgCanny(40, 40)
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
        
        for(int y=0; y<40; ++y)
        {
            for(int x=0; x<40; ++x)
            {
                if(x==y)
                    imgCanny(x,y) = 0.0;
                else if(x > y) 
                    imgCanny(x,y) = 2.0;
                else
                    imgCanny(x,y) = -2.0;
            }
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
    
    void edgeToCrackEdgeTest()
    {
        img1.upperLeft()[vigra::Diff2D(2,2)] = 0.0;
        
        Image res(9,9);
        res = 0.0;
        
        differenceOfExponentialCrackEdgeImage(srcImageRange(img1), 
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
        img1.upperLeft()[vigra::Diff2D(2,2)] = 0.0;
        
        Image res(9,9);
        res = 0.0;
        
        differenceOfExponentialCrackEdgeImage(srcImageRange(img1), 
                                                destImage(res), 
                                                0.7, 0.1, 1.0);
        removeShortEdges(destImageRange(res), 9, 0.0);
        
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
    
    void beautifyCrackEdgeTest()
    {
        img1.upperLeft()[vigra::Diff2D(2,2)] = 0.0;
        
        Image res(9,9);
        res = 0.0;
        
        differenceOfExponentialCrackEdgeImage(srcImageRange(img1), 
                                                destImage(res), 
                                                0.7, 0.1, 1.0);
        beautifyCrackEdgeImage(destImageRange(res), 1.0, 0.0);
        
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
    
    void closeGapsInCrackEdgeTest()
    {
        closeGapsInCrackEdgeImage(destImageRange(img2), 1.0);
        
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
    
    void cannyEdgelListTest()
    {
        std::vector<vigra::Edgel> edgels;
        cannyEdgelList(srcImageRange(imgCanny), edgels, 1.0);
        int count = 0;
        for(unsigned int i=0; i<edgels.size(); ++i)
        {
            if (edgels[i].strength < 1.0e-10)
                continue;  // ignore edgels that result from round off error during convolution
            ++count;
            should(edgels[i].x == edgels[i].y);
            should(fabs(edgels[i].orientation-M_PI*0.75) < 0.1);
        }
        should(count == 75);
    }
    
    void cannyEdgeImageTest()
    {
        vigra::BImage result(40, 40);
        result = 0;
        
        cannyEdgeImage(srcImageRange(imgCanny), destImage(result), 1.0, 0.1, 1);
        
        for(int y=1; y<39; ++y)
        {
            for(int x=1; x<39; ++x)
            {
                if(x == y)
                    should(result(x,y) == 1);
                else
                    should(result(x,y) == 0);
            }
        }   
    }
    
    Image img1, img2, imgCanny;
};
    
struct DistanceTransformTest
{
    typedef vigra::DImage Image;

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
        Image res(img);
        
        distanceTransform(srcImageRange(img), destImage(res), 0.0, 1);
        
        Image::Iterator i = res.upperLeft();
        Image::Accessor acc = res.accessor();
        int x,y;
        
        for(y=0; y<7; ++y)
        {
            for(x=0; x<7; ++x)
            {
                double dist = acc(i, vigra::Diff2D(x,y));
                double dist1 = std::abs(2.0 - x) + std::abs(2.0 - y);
                double dist2 = std::abs(5.0 - x) + std::abs(5.0 - y);
                double desired = (dist1 < dist2) ? dist1 : dist2;
                
                shouldEqualTolerance(dist, desired, 1e-14);
            }
        }
    }
    
    void distanceTransformL2Test()
    {
        
        Image res(img);
        
        distanceTransform(srcImageRange(img), destImage(res), 0.0, 2);
        
        Image::Iterator i = res.upperLeft();
        Image::Accessor acc = res.accessor();
        int x,y;
        
        for(y=0; y<7; ++y)
        {
            for(x=0; x<7; ++x)
            {
                double dist = acc(i, vigra::Diff2D(x,y));
                double dist1 = VIGRA_CSTD::sqrt((2.0 - x)*(2.0 - x) +
                                         (2.0 - y)*(2.0 - y));
                double dist2 = VIGRA_CSTD::sqrt((5.0 - x)*(5.0 - x) +
                                         (5.0 - y)*(5.0 - y));
                double desired = (dist1 < dist2) ? dist1 : dist2;
                
                shouldEqualTolerance(dist, desired, 1e-7);
            }
        }
    }
    
    void distanceTransformLInfTest()
    {
        
        Image res(img);
        
        distanceTransform(srcImageRange(img), destImage(res), 0.0, 0);
        
        Image::Iterator i = res.upperLeft();
        Image::Accessor acc = res.accessor();
        int x,y;
        
        for(y=0; y<7; ++y)
        {
            for(x=0; x<7; ++x)
            {
                double dist = acc(i, vigra::Diff2D(x,y));
                double dist1 = std::abs(2.0 - x) < std::abs(2.0 - y) ?
                                      std::abs(2.0 - y) : std::abs(2.0 - x);
                double dist2 = std::abs(5.0 - x) < std::abs(5.0 - y) ?
                                      std::abs(5.0 - y) : std::abs(5.0 - x);
                double desired = (dist1 < dist2) ? dist1 : dist2;
                
                shouldEqualTolerance(dist, desired, 1e-14);
            }
        }
    }
    
    
    Image img;
};

struct LocalMinMaxTest
{
    typedef vigra::DImage Image;

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
    typedef vigra::DImage Image;

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
        typedef double argument_type;
        typedef double result_type;
        typedef double cost_type;

        void operator()(double const &) {}

        double const & cost(double const & v) const
        {
            return v;
        }
    };
    
    void voronoiTest()
    {
        Image res(img);
        
        vigra::ArrayOfRegionStatistics<DirectCostFunctor> cost(2);
        seededRegionGrowing(srcImageRange(img), srcImage(seeds), 
                            destImage(res), cost);
        
        Image::Iterator i = res.upperLeft();
        Image::Accessor acc = res.accessor();
        int x,y;
        
        for(y=0; y<7; ++y)
        {
            for(x=0; x<7; ++x)
            {
                double dist = acc(i, vigra::Diff2D(x,y));
                double dist1 = VIGRA_CSTD::sqrt((2.0 - x)*(2.0 - x) +
                                         (2.0 - y)*(2.0 - y));
                double dist2 = VIGRA_CSTD::sqrt((5.0 - x)*(5.0 - x) +
                                         (5.0 - y)*(5.0 - y));
                double desired = (dist1 <= dist2) ? 1 : 2;
                
                shouldEqualTolerance(dist, desired, 1e-14);
            }
        }
    }
    
    Image img, seeds;
};

struct InterestOperatorTest
{
    typedef vigra::DImage Image;

    InterestOperatorTest()
    : img(9,9)
    {
        static const double in[] = {
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
        
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
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
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
    
    void foerstnerCornerTest()
    {
        Image tmp(img);
        Image res(img);
        res = 0.0;
        
        foerstnerCornerDetector(srcImageRange(img), destImage(tmp), 1.0);
        localMaxima(srcImageRange(tmp), destImage(res), 1.0);
        
        static const double desired[] = {
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
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
    
    void rohrCornerTest()
    {
        Image tmp(img);
        Image res(img);
        res = 0.0;
        
        rohrCornerDetector(srcImageRange(img), destImage(tmp), 1.0);
        localMaxima(srcImageRange(tmp), destImage(res), 1.0);
        
        static const double desired[] = {
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
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
    
     
    void beaudetCornerTest()
    {
        Image tmp(img);
        Image res(img);
        res = 0.0;
        
        beaudetCornerDetector(srcImageRange(img), destImage(tmp), 1.0);
        localMaxima(srcImageRange(tmp), destImage(res), 1.0);
        
        static const double desired[] = {
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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
   
    void radialSymmetryTest()
    {
        Image tmp(img);
        Image res(img);
        res = 0.0;
        
        radialSymmetryTransform(srcImageRange(img), destImage(tmp), 1.0);
        localMaxima(srcImageRange(tmp), destImage(res), 1.0);
        localMinima(srcImageRange(tmp), destImage(res), -1.0);
        
        static const double desired[] = {
                   0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0,  0.0, 1.0,  0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 1.0, -1.0, 0.0, -1.0, 1.0, 0.0, 0.0,
                   0.0, 0.0, 0.0,  0.0, 1.0,  0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0, 0.0, 0.0};
                                         
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

struct SimpleAnalysisTestSuite
: public vigra::test_suite
{
    SimpleAnalysisTestSuite()
    : vigra::test_suite("SimpleAnalysisTestSuite")
    {
        add( testCase( &LabelingTest::labelingFourTest1));
        add( testCase( &LabelingTest::labelingFourTest2));
        add( testCase( &LabelingTest::labelingFourTest3));
        add( testCase( &LabelingTest::labelingFourTest4));
        add( testCase( &LabelingTest::labelingToCrackEdgeTest));
        add( testCase( &LabelingTest::labelingEightTest1));
        add( testCase( &LabelingTest::labelingEightTest2));
        add( testCase( &LabelingTest::labelingFourWithBackgroundTest1));
        add( testCase( &LabelingTest::labelingFourWithBackgroundTest2));
        add( testCase( &LabelingTest::labelingEightWithBackgroundTest));
        add( testCase( &EdgeDetectionTest::edgeDetectionTest));
        add( testCase( &EdgeDetectionTest::edgeToCrackEdgeTest));
        add( testCase( &EdgeDetectionTest::removeShortEdgesTest));
        add( testCase( &EdgeDetectionTest::beautifyCrackEdgeTest));
        add( testCase( &EdgeDetectionTest::closeGapsInCrackEdgeTest));
        add( testCase( &EdgeDetectionTest::cannyEdgelListTest));
        add( testCase( &EdgeDetectionTest::cannyEdgeImageTest));
        add( testCase( &DistanceTransformTest::distanceTransformL1Test));
        add( testCase( &DistanceTransformTest::distanceTransformL2Test));
        add( testCase( &DistanceTransformTest::distanceTransformLInfTest));
        add( testCase( &LocalMinMaxTest::localMinimumTest));
        add( testCase( &LocalMinMaxTest::localMaximumTest));
        add( testCase( &LocalMinMaxTest::extendedLocalMinimumTest));
        add( testCase( &LocalMinMaxTest::extendedLocalMaximumTest));
        add( testCase( &RegionGrowingTest::voronoiTest));
        add( testCase( &InterestOperatorTest::cornerResponseFunctionTest));
        add( testCase( &InterestOperatorTest::foerstnerCornerTest));
        add( testCase( &InterestOperatorTest::rohrCornerTest));
        add( testCase( &InterestOperatorTest::beaudetCornerTest));
        add( testCase( &InterestOperatorTest::radialSymmetryTest));
    }
};

int main()
{
    SimpleAnalysisTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

