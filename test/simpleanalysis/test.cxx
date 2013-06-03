/************************************************************************/
/*                                                                      */
/*              Copyright 2004-2010 by Ullrich Koethe                   */
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

#include <iostream>
#include <fstream>
#include <functional>
#include <cmath>
#include "unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/labelimage.hxx"
#include "vigra/edgedetection.hxx"
#include "vigra/distancetransform.hxx"
#include "vigra/localminmax.hxx"
#include "vigra/multi_localminmax.hxx"
#include "vigra/seededregiongrowing.hxx"
#include "vigra/cornerdetection.hxx"
#include "vigra/symmetry.hxx"
#include "vigra/watersheds.hxx"
#include "vigra/multi_watersheds.hxx"
#include "vigra/noise_normalization.hxx"
#include "vigra/affinegeometry.hxx"
#include "vigra/affine_registration.hxx"
#include "vigra/impex.hxx"

#ifdef HasFFTW3
# include "vigra/slanted_edge_mtf.hxx"
#endif

using namespace vigra;

struct LabelingTest
{
    typedef vigra::DImage Image;
    typedef vigra::MultiArrayView<2, double> View;

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
        int count = labelImage(srcImageRange(img1), destImage(res), false);
        
        should(2 == count);

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
        int count = labelImage(srcImageRange(img3), destImage(res), false);

        should(6 == count);

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

    void labelingToEdgeTest()
    {
        Image tmp(img1);
        Image res(9, 9), res2(img1.size());

        should(2 == labelImage(View(img1), View(tmp), false));

        regionImageToCrackEdgeImage(View(tmp), View(res), 0.0);

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

        regionImageToEdgeImage(View(tmp), View(res2), 1.0);

        static const double edges[] = {
            0, 1, 1, 1, 0, 
            1, 0, 0, 1, 0, 
            1, 0, 0, 1, 0, 
            1, 1, 1, 1, 0, 
            0, 0, 0, 0, 0};

        shouldEqualSequence(res2.begin(), res2.end(), edges);
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
        Image res(img2), res2(img2);

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

        should(1 == labelImageWithBackground(View(img2), View(res2),
                                             true, 0.0));
        should(View(res) == View(res2));
    }

    Image img1, img2, img3, img4;
};

struct EdgeDetectionTest
{
    typedef vigra::DImage Image;
    typedef vigra::MultiArrayView<2, double> View;

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

        differenceOfExponentialEdgeImage(View(img1), View(res),
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

        differenceOfExponentialCrackEdgeImage(View(img1), View(res),
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

        differenceOfExponentialCrackEdgeImage(View(img1), View(res),
                                              0.7, 0.1, 1.0);
        removeShortEdges(View(res), 9, 0.0);

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

        differenceOfExponentialCrackEdgeImage(View(img1), View(res),
                                              0.7, 0.1, 1.0);
        beautifyCrackEdgeImage(View(res), 1.0, 0.0);

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
        closeGapsInCrackEdgeImage(View(img2), 1.0);

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
        {
            std::vector<vigra::Edgel> edgels;
            cannyEdgelList(View(imgCanny), edgels, 1.0);
            int count = 0;
            for(unsigned int i=0; i<edgels.size(); ++i)
            {
                if (edgels[i].strength < 1.0e-10)
                    continue;  // ignore edgels that result from round off error during convolution
                ++count;
                should(edgels[i].x == edgels[i].y);
                should(VIGRA_CSTD::fabs(edgels[i].orientation-M_PI*0.25) < 0.1);
            }
            should(count == 75);

            std::vector<vigra::Edgel> edgelsThresh;
            double threshold = 1.25;
            cannyEdgelListThreshold(View(imgCanny), edgelsThresh, 1.0, threshold);
            count = 0;
            for(unsigned int i=0; i<edgels.size(); ++i)
            {
                if (edgels[i].strength <= threshold)
                    continue;  // ignore edgels below threshold
                should(edgels[i].x == edgelsThresh[count].x);
                should(edgels[i].y == edgelsThresh[count].y);
                ++count;
            }
            should(count == 38);
        }
        {
            std::vector<vigra::Edgel> edgels;
            MultiArray<2, TinyVector<double, 2> > grad(imgCanny.width(), imgCanny.height());
            gaussianGradient(View(imgCanny), grad, 1.0);

            cannyEdgelList(grad, edgels);
            int count = 0;
            for(unsigned int i=0; i<edgels.size(); ++i)
            {
                if (edgels[i].strength < 1.0e-10)
                    continue;  // ignore edgels that result from round off error during convolution
                ++count;
                should(edgels[i].x == edgels[i].y);
                should(VIGRA_CSTD::fabs(edgels[i].orientation-M_PI*0.25) < 0.1);
            }
            should(count == 75);
        }
    }

    void cannyEdgelList3x3Test()
    {
        std::vector<vigra::Edgel> edgels;
        cannyEdgelList3x3(View(imgCanny), edgels, 1.0);
        int count = 0;
        for(unsigned int i=0; i<edgels.size(); ++i)
        {
            if (edgels[i].strength < 1.0e-10)
                continue;  // ignore edgels that result from round off error during convolution
            ++count;
            should(edgels[i].x == edgels[i].y);
            should(VIGRA_CSTD::fabs(edgels[i].orientation-M_PI*0.25) < 0.1);
        }
        should(count == 38);

        std::vector<vigra::Edgel> edgelsThresh;
        double threshold = 1.3;
        cannyEdgelList3x3Threshold(View(imgCanny), edgelsThresh, 1.0, threshold);
        count = 0;
        for(unsigned int i=0; i<edgels.size(); ++i)
        {
            if (edgels[i].strength <= threshold)
                continue;  // ignore edgels below threshold
            should(edgels[i].x == edgelsThresh[count].x);
            should(edgels[i].y == edgelsThresh[count].y);
            ++count;
        }
        should(count == 36);
    }

    void cannyEdgeImageTest()
    {
        vigra::BImage result(40, 40);
        result = 0;

        cannyEdgeImage(View(imgCanny), MultiArrayView<2, unsigned char>(result), 1.0, 0.1, 1);

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

    void cannyEdgeImageWithThinningTest()
    {
        vigra::BImage result(40, 40);
        result = 0;

        cannyEdgeImageWithThinning(View(imgCanny), MultiArrayView<2, unsigned char>(result), 1.0, 0.1, 1, false);

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
    typedef vigra::MultiArrayView<2, double> View;

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
        Image res1(img);

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

        distanceTransform(View(img), View(res1), 0.0, 1);
        should(View(res) == View(res1));
    }

    void distanceTransformL2Test()
    {

        Image res(img);

        distanceTransform(View(img), View(res), 0.0, 2);

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

        distanceTransform(View(img), View(res), 0.0, 0);

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

template <class T>
struct EqualWithToleranceFunctor
{
    EqualWithToleranceFunctor(T tolerance = 2.0 * NumericTraits<T>::epsilon())
    : t(tolerance)
    {}
    
    bool operator()(T l, T r) const
    {
        return abs(l-r) <= t;
    }
    
    T t;
};

struct LocalMinMaxTest
{
    typedef vigra::MultiArray<2, double> Image;
    typedef vigra::MultiArray<3, double> Volume;
    typedef MultiArrayShape<3>::type Shp3D;

    LocalMinMaxTest()
    : img(Shape2(9,9)), vol()
    {
        static const double in[] = {
            0.2,  0.1,  0.1,  0.3,  0.5,  0.3,  0.0,  0.0, -0.1,
            0.0, -0.1,  0.1,  0.0,  1.0,  0.0,  0.3,  0.0,  0.0,
            0.0,  0.5,  2.0,  0.0,  2.0,  2.0,  2.0,  0.0, -1.1,
            0.1,  0.0,  1.0,  1.5,  1.0,  1.0,  0.0,  0.0,  0.0,
            0.0,  0.1,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
            0.0,  0.0,  0.0,  0.0, -1.0, -1.5, -1.0,  0.0, -0.3,
            0.0,  0.0, -2.0, -2.0, -2.0,  0.0, -2.0, -0.5,  0.0,
            1.0,  0.0,  0.0,  0.0, -1.0,  0.0, -0.1,  0.1,  0.0,
            0.0,  0.0,  0.0,  0.0, -0.5, -0.3, -0.1, -0.1,  0.0};

        Image::iterator i = img.begin();
        Image::iterator end = img.end();
        const double * p = in;

        for(; i != end; ++i, ++p)
        {
            *i = *p;
        }

        //prepare the multiarray
        vol.reshape(Shp3D(10,20,50),0);

        vol(1,1,1)=10;
        vol(5,5,5)=350;

        vol(8,3,5)=9; //plateaux
        vol(8,4,5)=9;
        vol(8,5,5)=9;

        vol(6,7,7)=-0.5;
        vol(7,7,7)=-1;
        vol(7,1,15)=-100;
        vol(7,1,19)=-20;

        vol(3,15,26)=-1; //plateaux
        vol(3,15,27)=-1;
        vol(3,15,28)=-1;
        vol(3,16,26)=-1;

        vol(9,18,35)=-100; //on the border is skipped
        vol(0,1,49)=100; //on the border is skipped

    }


    void localMinimum3DTest()
    {
        Volume res(vol.shape()), res2(vol.shape());

        localMinima3D(srcMultiArrayRange(vol), destMultiArray(res), 1, NeighborCode3DSix());

        Volume desired(vol);
        desired.init(0);

        desired(7,7,7)=1;
        desired(7,1,15)=1;
        desired(7,1,19)=1;


        for(int z=0; z<vol.shape(2); ++z)
            for(int y=0; y<vol.shape(1); ++y)
                for(int x=0; x<vol.shape(0); ++x)
                    shouldEqual(res(x,y,z), desired(x,y,z));

        should(3 == localMinima(vol, res2, LocalMinmaxOptions().neighborhood(0).markWith(1)));
        should(res == res2);

        res2 = 0;
        localMinima3D(vol, res2, 1, NeighborCode3DSix());
        should(res == res2);
    }

    void extendedLocalMinimum3DTest()
    {
        Volume res(vol.shape()), res2(vol.shape());

        extendedLocalMinima3D(srcMultiArrayRange(vol), destMultiArray(res),1,NeighborCode3DSix());

        Volume desired(vol);
        desired.init(0);

        desired(7,7,7)=1;
        desired(7,1,15)=1;
        desired(7,1,19)=1;

        desired(3,15,26)=1; //plateau
        desired(3,15,27)=1;
        desired(3,15,28)=1;
        desired(3,16,26)=1;

        for(int z=0; z<vol.shape(2); ++z)
            for(int y=0; y<vol.shape(1); ++y)
                for(int x=0; x<vol.shape(0); ++x)
                    shouldEqual(res(x,y,z), desired(x,y,z));

        should(4 == localMinima(vol, res2, LocalMinmaxOptions().neighborhood(0).markWith(1).allowPlateaus()));
        should(res == res2);
    }
    
    void extendedLocalMinimum3DTest2()
    {
        Volume res(vol.shape()), res2(vol.shape());

        extendedLocalMinima3D(srcMultiArrayRange(vol), destMultiArray(res),1,NeighborCode3DTwentySix());

        Volume desired(vol);
        desired.init(0);

        desired(7,7,7)=1;
        desired(7,1,15)=1;
        desired(7,1,19)=1;

        desired(3,15,26)=1; //plateaux
        desired(3,15,27)=1;
        desired(3,15,28)=1;
        desired(3,16,26)=1;

        for(int z=0; z<vol.shape(2); ++z)
            for(int y=0; y<vol.shape(1); ++y)
                for(int x=0; x<vol.shape(0); ++x)
                    shouldEqual(res(x,y,z), desired(x,y,z));

        should(4 == localMinima(vol, res2, LocalMinmaxOptions().neighborhood(1).markWith(1).allowPlateaus()));
        should(res == res2);
    }

    void localMaximum3DTest()
    {
        Volume res(vol.shape()), res2(vol.shape());

        localMaxima3D(srcMultiArrayRange(vol), destMultiArray(res), 1, NeighborCode3DSix());

        Volume desired(vol);
        desired.init(0);

        desired(1,1,1)=1;
        desired(5,5,5)=1;

        for(int z=0; z<vol.shape(2); ++z)
            for(int y=0; y<vol.shape(1); ++y)
                for(int x=0; x<vol.shape(0); ++x)
                    shouldEqual(res(x,y,z), desired(x,y,z));

        should(2 == localMaxima(vol, res2, LocalMinmaxOptions().neighborhood(0).markWith(1)));
        should(res == res2);

        res2 = 0;
        localMaxima3D(vol, res2, 1, NeighborCode3DSix());
        should(res == res2);
    }

    void extendedLocalMaximum3DTest()
    {
        Volume res(vol.shape()), res2(vol.shape());

        extendedLocalMaxima3D(srcMultiArrayRange(vol), destMultiArray(res),1,NeighborCode3DSix());

        Volume desired(vol);
        desired.init(0);

        desired(1,1,1)=1;
        desired(5,5,5)=1;

        desired(8,3,5)=1;
        desired(8,4,5)=1;
        desired(8,5,5)=1;

        for(int z=0; z<vol.shape(2); ++z)
            for(int y=0; y<vol.shape(1); ++y)
                for(int x=0; x<vol.shape(0); ++x)
                    shouldEqual(res(x,y,z), desired(x,y,z));

        should(3 == localMaxima(vol, res2, LocalMinmaxOptions().neighborhood(0).markWith(1).allowPlateaus()));
        should(res == res2);
    }
    
    void extendedLocalMaximum3DTest2()
    {
        Volume res(vol.shape()), res2(vol.shape());

        extendedLocalMaxima3D(srcMultiArrayRange(vol), destMultiArray(res),1,NeighborCode3DTwentySix());

        Volume desired(vol);
        desired.init(0);

        desired(1,1,1)=1;
        desired(5,5,5)=1;

        desired(8,3,5)=1;
        desired(8,4,5)=1;
        desired(8,5,5)=1;

        for(int z=0; z<vol.shape(2); ++z)
            for(int y=0; y<vol.shape(1); ++y)
                for(int x=0; x<vol.shape(0); ++x)
                    shouldEqual(res(x,y,z), desired(x,y,z));

        should(3 == localMaxima(vol, res2, LocalMinmaxOptions().neighborhood(1).markWith(1).allowPlateaus()));
        should(res == res2);
    }

    void localMinimumTest()
    {
        Image res(img.shape()), res2(img.shape());

        localMinima(srcImageRange(img), destImage(res));

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        shouldEqualSequence(res.begin(), res.end(), desired);

        should(2 == localMinima(img, res2));
        should(res == res2);

        res.init(0);
        localMinima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().allowAtBorder());
        desired[8] = 1.0;
        desired[26] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);

        res2.init(0);
        should(4 == localMinima(img, res2, LocalMinmaxOptions().allowAtBorder()));
        should(res == res2);
    }

    void localMinimum4Test()
    {
        Image res(img.shape()), res2(img.shape());

        localMinima(srcImageRange(img), destImage(res), LocalMinmaxOptions().neighborhood(4));

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        shouldEqualSequence(res.begin(), res.end(), desired);

        should(5 == localMinima(img, res2, LocalMinmaxOptions().neighborhood(4)));
        should(res == res2);

        res.init(0);
        localMinima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().neighborhood(0).allowAtBorder());
        desired[8] = 1.0;
        desired[26] = 1.0;
        desired[53] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);

        res2.init(0);
        should(8 == localMinima(img, res2, LocalMinmaxOptions().allowAtBorder().neighborhood(4)));
        should(res == res2);
    }

    void localMinimumTestThr()
    {
        Image res(img.shape()), res2(img.shape());

        localMinima(srcImageRange(img), destImage(res),
                    LocalMinmaxOptions().neighborhood(8).markWith(1.0).threshold(-1.0));

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        shouldEqualSequence(res.begin(), res.end(), desired);

        should(1 == localMinima(img, res2,
                                LocalMinmaxOptions().neighborhood(8).markWith(1.0).threshold(-1.0)));
        should(res == res2);

        res.init(0);
        localMinima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().neighborhood(8).threshold(-1.0).allowAtBorder());
        desired[26] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);

        res2.init(0);
        should(2 == localMinima(img, res2, 
                                LocalMinmaxOptions().neighborhood(8).threshold(-1.0).allowAtBorder()));
        should(res == res2);
    }

    void localMaximumTest()
    {
        Image res(img.shape()), res2(img.shape());

        localMaxima(srcImageRange(img), destImage(res));

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        shouldEqualSequence(res.begin(), res.end(), desired);

        should(2 == localMaxima(img, res2));
        should(res == res2);

        res.init(0);
        localMaxima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().allowAtBorder());
        desired[0] = 1.0;
        desired[63] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);

        res2.init(0);
        should(4 == localMaxima(img, res2, 
                                LocalMinmaxOptions().allowAtBorder()));
        should(res == res2);
    }

    void localMaximum4Test()
    {
        Image res(img.shape()), res2(img.shape());

        localMaxima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().neighborhood(4));

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        shouldEqualSequence(res.begin(), res.end(), desired);

        should(4 == localMaxima(img, res2,
                                LocalMinmaxOptions().neighborhood(4)));
        should(res == res2);

        res.init(0);
        localMaxima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().neighborhood(DirectNeighborhood).allowAtBorder());
        desired[0] = 1.0;
        desired[27] = 1.0;
        desired[63] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);

        res2.init(0);
        should(7 == localMaxima(img, res2, 
                                LocalMinmaxOptions().neighborhood(DirectNeighborhood).allowAtBorder()));
        should(res == res2);
    }

    void localMaximumTestThr()
    {
        Image res(img.shape()), res2(img.shape());

        localMaxima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().markWith(1.0).neighborhood(8).threshold(0.2));

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        shouldEqualSequence(res.begin(), res.end(), desired);

        should(1 == localMaxima(img, res2,
                                LocalMinmaxOptions().neighborhood(8).markWith(1.0).threshold(0.2)));
        should(res == res2);

        res.init(0);
        localMaxima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().allowAtBorder().threshold(0.2));
        desired[63] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);

        res2.init(0);
        should(2 == localMaxima(img, res2, 
                                LocalMinmaxOptions().neighborhood(IndirectNeighborhood).threshold(0.2).allowAtBorder()));
        should(res == res2);
    }

    void extendedLocalMinimumTest()
    {
        Image res(img.shape()), res2(img.shape());

        extendedLocalMinima(srcImageRange(img), destImage(res), 1.0);

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        shouldEqualSequence(res.begin(), res.end(), desired);

        res.init(0);
        localMinima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().allowPlateaus());
        shouldEqualSequence(res.begin(), res.end(), desired);

        should(4 == localMinima(img, res2, 
                                LocalMinmaxOptions().allowPlateaus()));
        should(res == res2);

        res.init(0);
        localMinima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().allowAtBorder().allowPlateaus());
        desired[8] = 1.0;
        desired[26] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);

        res2.init(0);
        should(6 == localMinima(img, res2, 
                                LocalMinmaxOptions().allowAtBorder().allowPlateaus()));
        should(res == res2);
    }

    void extendedLocalMinimum4Test()
    {
        Image res(img.shape()), res2(img.shape());

        extendedLocalMinima(srcImageRange(img), destImage(res), 1.0, FourNeighborCode());

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        shouldEqualSequence(res.begin(), res.end(), desired);
 
        res.init(0);
        localMinima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().allowPlateaus().neighborhood(4));
        shouldEqualSequence(res.begin(), res.end(), desired);

        should(7 == localMinima(img, res2, 
                                LocalMinmaxOptions().allowPlateaus().neighborhood(0)));
        should(res == res2);

        res.init(0);
        localMinima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().neighborhood(4).allowAtBorder().allowPlateaus());
        desired[8] = 1.0;
        desired[26] = 1.0;
        desired[53] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);

        res2.init(0);
        should(10 == localMinima(img, res2, 
                                LocalMinmaxOptions().allowAtBorder().allowPlateaus().neighborhood(4)));
        should(res == res2);
   }

    void extendedLocalMaximumTest()
    {
        Image res(img.shape()), res2(img.shape());

        extendedLocalMaxima(srcImageRange(img), destImage(res), 1.0);

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        shouldEqualSequence(res.begin(), res.end(), desired);
        res.init(0);
        localMaxima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().allowPlateaus());
        shouldEqualSequence(res.begin(), res.end(), desired);

        should(4 == localMaxima(img, res2, 
                                LocalMinmaxOptions().allowPlateaus()));
        should(res == res2);

        res.init(0);
        localMaxima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().allowAtBorder().allowPlateaus());
        desired[0] = 1.0;
        desired[63] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);

        res2.init(0);
        should(6 == localMaxima(img, res2, 
                                LocalMinmaxOptions().allowAtBorder().allowPlateaus()));
        should(res == res2);
   }

    void extendedLocalMaximum4Test()
    {
        Image res(img.shape()), res2(img.shape());

        extendedLocalMaxima(srcImageRange(img), destImage(res), 1.0, FourNeighborCode());

        double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        shouldEqualSequence(res.begin(), res.end(), desired);

        res.init(0);
        localMaxima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().allowPlateaus().neighborhood(4));
        shouldEqualSequence(res.begin(), res.end(), desired);

        should(6 == localMaxima(img, res2, 
                                LocalMinmaxOptions().allowPlateaus().neighborhood(DirectNeighborhood)));
        should(res == res2);

        res.init(0);
        localMaxima(srcImageRange(img), destImage(res), 
                    LocalMinmaxOptions().neighborhood(4).allowAtBorder().allowPlateaus());
        desired[0] = 1.0;
        desired[27] = 1.0;
        desired[63] = 1.0;
        shouldEqualSequence(res.begin(), res.end(), desired);

        res2.init(0);
        should(9 == localMaxima(img, res2, 
                                LocalMinmaxOptions().allowAtBorder().allowPlateaus().neighborhood(4)));
        should(res == res2);
   }

    void plateauWithHolesTest()
    {
        static const double in[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 2.0, 3.0, 2.0, 2.0, 1.0, 0.0, 0.0,
            0.0, 1.0, 3.0, 4.0, 4.1, 3.0, 1.0, 0.0, 0.0,
            0.0, 1.0, 2.0, 2.0, 3.0, 2.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        std::copy(in, in+81, img.begin());
        Image res(img.shape()), res2(img.shape());

        extendedLocalMaxima(srcImageRange(img), destImage(res), 1.0,
                            EightNeighborCode(),
                            EqualWithToleranceFunctor<Image::value_type>(0.2));
        
        static const double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        shouldEqualSequence(res.begin(), res.end(), desired);
    }

    Image img;
    Volume vol;
};

struct WatershedsTest
{
    typedef vigra::MultiArray<2, double> Image;
    typedef vigra::MultiArray<2, int> IntImage;

    WatershedsTest()
    : img(Shape2(9,9))
    {
        static const double in[] = {
            0.0,  0.1,  0.1,  0.3,  0.5,  0.3,  0.0,  0.0, 0.0,
            0.0, -0.1,  0.1,  0.0,  1.0,  0.0,  0.3,  0.0, 0.0,
            0.0,  0.5,  2.0,  0.0,  2.0,  2.0,  2.0,  0.0, 0.0,
            0.0,  0.0,  1.0,  1.5,  1.0,  1.0,  3.0,  4.0, 0.0,
            0.3,  0.1,  1.5,  0.0,  0.0,  0.0,  0.0,  3.0, 3.0,
            0.0,  0.0,  0.0,  0.0, -1.0, -1.5, -1.0,  0.0, 0.0,
            0.0,  0.0, -2.0, -2.0, -2.0,  0.0, -2.1, -0.5, 0.0,
            0.0,  0.0,  0.0,  0.0, -1.0,  0.0, -0.1,  0.1, 0.0,
            0.0,  0.0,  0.0,  0.0, -0.5, -0.3, -0.1, -0.1, 0.0};

        Image::iterator i = img.begin();
        Image::iterator end = img.end();
        const double * p = in;

        for(; i != end; ++i, ++p)
        {
            // transform data to a range suitable for BucketQueue (in the turbo algorithm)
            *i = *p*10.0 + 30.0;
        }
    }

    void watershedsTest()
    {
        IntImage res(img.shape()), res2(img.shape());

        /*******************************************************************/
        
        static const double desired[] = {
            1.0,  1.0,  1.0,  2.0,  3.0,  3.0,  3.0,  3.0,  4.0,
            1.0,  1.0,  1.0,  2.0,  2.0,  3.0,  3.0,  3.0,  4.0,
            1.0,  1.0,  1.0,  2.0,  2.0,  3.0,  3.0,  3.0,  4.0,
            1.0,  1.0,  1.0,  5.0,  6.0,  6.0,  6.0,  3.0,  4.0,
            7.0,  5.0,  5.0,  5.0,  6.0,  6.0,  6.0,  6.0,  6.0,
            7.0,  5.0,  5.0,  5.0,  5.0,  6.0,  6.0,  6.0,  6.0,
            7.0,  5.0,  5.0,  5.0,  5.0,  6.0,  6.0,  6.0,  6.0,
            7.0,  5.0,  5.0,  5.0,  5.0,  6.0,  6.0,  6.0,  6.0,
            7.0,  7.0,  7.0,  5.0,  5.0,  5.0,  5.0,  5.0,  5.0};

        int count = watershedsUnionFind(srcImageRange(img), destImage(res));
        
        shouldEqual(7, count);
        shouldEqualSequence(res.begin(), res.end(), desired);

        /*******************************************************************/
        
        // break ties explicitly to make the test independent of tie breaking rules
        img(3,1) -= 0.01;
        img(3,2) -= 0.01;
        img(5,4) += 0.01;
        img(6,4) += 0.01;
        img(7,5) += 0.01;
        img(8,5) += 0.01;

        static const double desiredSeeds[] = {
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  
            0.0,  2.0,  0.0,  3.0,  0.0,  1.0,  0.0,  1.0,  1.0,  
            0.0,  0.0,  0.0,  3.0,  0.0,  0.0,  0.0,  1.0,  1.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  
            0.0,  0.0,  4.0,  4.0,  4.0,  0.0,  5.0,  0.0,  0.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};


        res.init(0);
        count = generateWatershedSeeds(srcImageRange(img), destImage(res),
                                       SeedOptions().extendedMinima());
        shouldEqual(5, count);
        shouldEqualSequence(res.begin(), res.end(), desiredSeeds);

        res2.init(0);
        should(5 == generateWatershedSeeds(img, res2, IndirectNeighborhood, SeedOptions().extendedMinima()));
        shouldEqualSequence(res2.begin(), res2.end(), desiredSeeds);

        /*******************************************************************/
        
        static const double desiredRG[] = {
            2.0,  2.0,  2.0,  3.0,  3.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  3.0,  3.0,  1.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  3.0,  3.0,  3.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  3.0,  3.0,  3.0,  3.0,  1.0,  1.0,  1.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  1.0,  1.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0};

        count = watershedsRegionGrowing(srcImageRange(img), destImage(res));

        shouldEqual(5, count);
        shouldEqualSequence(res.begin(), res.end(), desiredRG);

        /*******************************************************************/
        
        static const double desiredRGC[] = {
            2.0,  2.0,  0.0,  3.0,  0.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  0.0,  3.0,  0.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  0.0,  3.0,  0.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  
            0.0,  0.0,  0.0,  4.0,  4.0,  0.0,  5.0,  0.0,  0.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  0.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  0.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  0.0,  0.0,  0.0,  0.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  4.0};

        res.init(0);
        count = watershedsRegionGrowing(srcImageRange(img), destImage(res),
                                        WatershedOptions().keepContours()
                                          .seedOptions(SeedOptions().extendedMinima()));
        shouldEqual(5, count);
        shouldEqualSequence(res.begin(), res.end(), desiredRGC);

        res2.init(0);
        generateWatershedSeeds(img, res2, IndirectNeighborhood, SeedOptions().extendedMinima());
        should(5 == watershedsMultiArray(img, res2, IndirectNeighborhood, WatershedOptions().regionGrowing().keepContours()));
        should(res == res2);

        /*******************************************************************/
        
        static const double desiredTRG[] = {
            2.0,  2.0,  2.0,  3.0,  3.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  2.0,  3.0,  3.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  2.0,  3.0,  3.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  3.0,  3.0,  3.0,  5.0,  1.0,  1.0,  1.0,  
            4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0,  5.0,  1.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0,  5.0,  5.0,  
            4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  4.0,  5.0,  5.0};

        res.init(0);
        count = watershedsRegionGrowing(srcImageRange(img), destImage(res),
                                        WatershedOptions().turboAlgorithm()
                                          .seedOptions(SeedOptions().extendedMinima()));
        shouldEqual(5, count);
        shouldEqualSequence(res.begin(), res.end(), desiredTRG);

        res2.init(0);
        generateWatershedSeeds(img, res2, IndirectNeighborhood, SeedOptions().extendedMinima());
        should(5 == watershedsMultiArray(img, res2, IndirectNeighborhood, WatershedOptions().regionGrowing()));
        should(res == res2);

        res2.init(1);  // check that this is overridden by explicit seed computation
        should(5 == watershedsMultiArray(img, res2, IndirectNeighborhood, WatershedOptions().regionGrowing().seedOptions(SeedOptions().extendedMinima())));
        should(res == res2);

#if 0
        std::cerr << count << "\n";
        for(int y=0;y<9;++y)
        {
            std::cerr << "            ";
            for(int x=0;x<9;++x)
                std::cerr << res(x,y) << ".0,  ";
            std::cerr << "\n\n";
        }
#endif /* #if 0 */
    }

    void watersheds4Test()
    {
        IntImage res(img.shape()), res2(img.shape());

        /*******************************************************************/
        
        static const double desired[] = {
            1.0,  1.0,  1.0,  2.0,  2.0,  3.0,  4.0,  4.0,  5.0,
            1.0,  1.0,  1.0,  2.0,  2.0,  3.0,  3.0,  4.0,  5.0,
            6.0,  1.0,  2.0,  2.0,  2.0,  3.0,  4.0,  4.0,  5.0,
            6.0,  6.0,  6.0,  7.0,  7.0,  8.0,  9.0,  4.0,  5.0,
           10.0,  7.0,  7.0,  7.0,  7.0,  8.0,  9.0,  9.0,  9.0,
           10.0,  7.0,  7.0,  7.0,  7.0,  8.0,  9.0,  9.0,  9.0,
           10.0,  7.0,  7.0,  7.0,  7.0,  9.0,  9.0,  9.0,  9.0,
           10.0, 10.0,  7.0,  7.0,  7.0,  7.0,  9.0,  9.0,  7.0,
           10.0, 10.0, 10.0,  7.0,  7.0,  7.0,  7.0,  7.0,  7.0};

        int count = watershedsUnionFind(srcImageRange(img), destImage(res), FourNeighborCode());

        should(10 == count);
        shouldEqualSequence(res.begin(), res.end(), desired);

        /*******************************************************************/
        
        // break ties explicitly to make the test independent of tie breaking rules
        img(3,0) += 0.01;
        img(3,4) += 0.01;
        img(5,1) += 0.01;
        img(0,3) += 0.01;
        img(1,3) += 0.01;
        img(8,5) += 0.01;
        img(6,8) += 0.01;
        img(7,8) += 0.01;
        img(8,8) += 0.01;

        static const double desiredSeeds[] = {
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  
            0.0,  2.0,  0.0,  3.0,  0.0,  4.0,  0.0,  1.0,  1.0,  
            0.0,  0.0,  0.0,  3.0,  0.0,  0.0,  0.0,  1.0,  1.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  5.0,  0.0,  0.0,  0.0,  
            0.0,  0.0,  6.0,  6.0,  6.0,  0.0,  7.0,  0.0,  0.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  
            0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};

        res.init(0);
        count = generateWatershedSeeds(srcImageRange(img), destImage(res), FourNeighborCode(),
                                       SeedOptions().extendedMinima());
        shouldEqual(7, count);
        shouldEqualSequence(res.begin(), res.end(), desiredSeeds);

        res2.init(0);
        should(7 == generateWatershedSeeds(img, res2, DirectNeighborhood, SeedOptions().extendedMinima()));
        shouldEqualSequence(res2.begin(), res2.end(), desiredSeeds);

        /*******************************************************************/
        
        static const double desiredRG[] = {
            2.0,  2.0,  2.0,  3.0,  3.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  3.0,  3.0,  4.0,  4.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  3.0,  3.0,  3.0,  4.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  2.0,  3.0,  5.0,  5.0,  1.0,  1.0,  1.0,  
            6.0,  6.0,  6.0,  6.0,  5.0,  5.0,  5.0,  1.0,  1.0,  
            6.0,  6.0,  6.0,  6.0,  5.0,  5.0,  5.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  7.0,  7.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  7.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  7.0,  7.0,  7.0};

        count = watershedsRegionGrowing(srcImageRange(img), destImage(res), FourNeighborCode());

        shouldEqual(7, count);
        shouldEqualSequence(res.begin(), res.end(), desiredRG);

        /*******************************************************************/
        
        static const double desiredRGC[] = {
            2.0,  2.0,  2.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  0.0,  3.0,  0.0,  4.0,  0.0,  1.0,  1.0,  
            2.0,  2.0,  0.0,  3.0,  0.0,  0.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  2.0,  0.0,  5.0,  5.0,  0.0,  1.0,  1.0,  
            0.0,  0.0,  0.0,  0.0,  5.0,  5.0,  5.0,  0.0,  0.0,  
            6.0,  6.0,  6.0,  6.0,  0.0,  5.0,  0.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  0.0,  7.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  0.0,  7.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  0.0,  7.0,  7.0};

        res.init(0);
        count = watershedsRegionGrowing(srcImageRange(img), destImage(res), FourNeighborCode(),
                                        WatershedOptions().keepContours()
                                          .seedOptions(SeedOptions().extendedMinima()));
        shouldEqual(7, count);
        shouldEqualSequence(res.begin(), res.end(), desiredRGC);

        res2.init(0);
        generateWatershedSeeds(img, res2, DirectNeighborhood, SeedOptions().extendedMinima());
        should(7 == watershedsMultiArray(img, res2, DirectNeighborhood, WatershedOptions().regionGrowing().keepContours()));
        should(res == res2);

        /*******************************************************************/
        
        static const double desiredTRG[] = {
            2.0,  2.0,  2.0,  3.0,  1.0,  1.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  2.0,  3.0,  3.0,  4.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  3.0,  3.0,  3.0,  4.0,  1.0,  1.0,  1.0,  
            2.0,  2.0,  2.0,  3.0,  6.0,  5.0,  7.0,  1.0,  1.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  5.0,  7.0,  7.0,  1.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  5.0,  7.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  7.0,  7.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  7.0,  7.0,  7.0,  
            6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  6.0,  6.0};

        res.init(0);
        count = watershedsRegionGrowing(srcImageRange(img), destImage(res), FourNeighborCode(),
                                        WatershedOptions().turboAlgorithm()
                                          .seedOptions(SeedOptions().extendedMinima()));
        shouldEqual(7, count);
        shouldEqualSequence(res.begin(), res.end(), desiredTRG);

        res2.init(0);
        generateWatershedSeeds(img, res2, DirectNeighborhood, SeedOptions().extendedMinima());
        should(7 == watershedsMultiArray(img, res2, DirectNeighborhood, WatershedOptions().regionGrowing()));
        should(res == res2);

        res2.init(1);  // check that this is overridden by explicit seed computation
        should(7 == watershedsMultiArray(img, res2, DirectNeighborhood, WatershedOptions().regionGrowing().seedOptions(SeedOptions().extendedMinima())));
        should(res == res2);

#if 0
        std::cerr << count << "\n";
        for(int k=0; k<res.size(); ++k)
        {
            std::cerr << res[k] << (res[k] == res2[k] ? " " : "*");
            if(k%res.shape(0) == res.shape(0)-1)
                std::cerr << "\n";
        }
#endif /* #if 0 */
    }

    Image img;
};

struct RegionGrowingTest
{
    typedef vigra::DImage Image;
    typedef MultiArrayView<2, double> View;

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
        seededRegionGrowing(View(img), View(seeds), View(res), cost);

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

                if(VIGRA_CSTD::fabs(dist1 - dist2) > 1e-10)
                    shouldEqual(dist, desired);
            }
        }

        vigra::IImage wres(img.size());

        watershedsRegionGrowing(srcImageRange(img), destImage(wres),
                                WatershedOptions().completeGrow()
                                 .seedOptions(SeedOptions().minima().threshold(1.0)));

        shouldEqualSequence(res.begin(), res.end(), wres.begin());
    }

    void voronoiWithBorderTest()
    {
        Image res(img);
        Image::value_type reference[] = {
            1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 0,
            1, 1, 1, 1, 1, 0, 2,
            1, 1, 1, 1, 0, 2, 2,
            1, 1, 1, 0, 2, 2, 2,
            1, 1, 0, 2, 2, 2, 2,
            1, 0, 2, 2, 2, 2, 2
        };

        vigra::ArrayOfRegionStatistics<DirectCostFunctor> cost(2);
        seededRegionGrowing(srcImageRange(img), srcImage(seeds),
                            destImage(res), cost, KeepContours);

        shouldEqualSequence(res.begin(), res.end(), reference);
    }

    Image img, seeds;
};

struct InterestOperatorTest
{
    typedef vigra::DImage Image;
    typedef vigra::MultiArrayView<2, double> View;

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
        Image tmp1(img);
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

        cornerResponseFunction(View(img), View(tmp1), 1.0);
        should(View(tmp) == View(tmp1));
    }

    void foerstnerCornerTest()
    {
        Image tmp(img);
        Image tmp1(img);
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

        foerstnerCornerDetector(View(img), View(tmp1), 1.0);
        should(View(tmp) == View(tmp1));
    }

    void rohrCornerTest()
    {
        Image tmp(img);
        Image tmp1(img);
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

        rohrCornerDetector(View(img), View(tmp1), 1.0);
        should(View(tmp) == View(tmp1));
    }


    void beaudetCornerTest()
    {
        Image tmp(img);
        Image tmp1(img);
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

        beaudetCornerDetector(View(img), View(tmp1), 1.0);
        should(View(tmp) == View(tmp1));
    }

    void radialSymmetryTest()
    {
        Image tmp(img);
        Image res(img);
        res = 0.0;

        radialSymmetryTransform(View(img), View(tmp), 1.0);
        localMaxima(srcImageRange(tmp), destImage(res), 1.0);
        localMinima(srcImageRange(tmp), destImage(res), -1.0);

        static const double desired[] = {
                   0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0,  0.0, 0.0,  0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0,  0.0, 1.0,  0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0,
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

struct NoiseNormalizationTest
{
    typedef vigra::BImage U8Image;
    typedef vigra::DImage GrayImage;
    typedef vigra::DRGBImage RGBImage;
    typedef MultiArrayView<2, double> View;
    U8Image u8image;
    GrayImage image;
    RGBImage rgb;

    NoiseNormalizationTest()
    {
        vigra::ImageImportInfo info("noiseNormalizationTest.xv");
        vigra_precondition(info.width() == 400 && info.height() == 20,
           "NoiseNormalizationTest: input image has wrong size.");
           
        u8image.resize(info.size());
        importImage(info, destImage(u8image));
        image.resize(info.size());
        importImage(info, destImage(image));
        rgb.resize(info.size());
        for(unsigned int band = 0; band < 3; ++band)
        {
            vigra::VectorElementAccessor<RGBImage::Accessor> dband(band, rgb.accessor());
            importImage(info, destImage(rgb, dband));
        }
    }
    
    template <class Iterator, class Accessor>
    void checkVariance(Iterator ul, Accessor const & a, double tolerance)
    {
        for(unsigned int k = 0; k < 20; ++k)
        {
            double sum = 0.0, sum2 = 0.0;
            for(unsigned int y = 0; y < 20; ++y)
            {
                for(unsigned int x = 20*k; x < 20*(k+1); ++x)
                {
                    sum += a(ul, Diff2D(x, y));
                    sum2 += sq(a(ul, Diff2D(x, y)));
                }
            }
            
            sum /= 400.0;
            sum2 /= 400.0;
            
            shouldEqualTolerance(VIGRA_CSTD::sqrt(sum2 - sq(sum))-1.0, 0.0, tolerance);
       }
    }
    
    void testParametricNoiseNormalizationU8()
    {
        GrayImage res(image.size());
        linearNoiseNormalization(srcImageRange(u8image), destImage(res), 1.0, 0.02);
        checkVariance(res.upperLeft(), res.accessor(), 0.1);
        linearNoiseNormalization(srcImageRange(u8image), destImage(res));
        checkVariance(res.upperLeft(), res.accessor(), 0.1);
        quadraticNoiseNormalization(srcImageRange(u8image), destImage(res), 1.0, 0.02, 0.0);
        checkVariance(res.upperLeft(), res.accessor(), 0.1);
        quadraticNoiseNormalization(srcImageRange(u8image), destImage(res));
        checkVariance(res.upperLeft(), res.accessor(), 0.1);
    }
  
    void testParametricNoiseNormalization()
    {
        GrayImage res(image.size());
        linearNoiseNormalization(View(image), View(res), 1.0, 0.02);
        checkVariance(res.upperLeft(), res.accessor(), 0.1);
        linearNoiseNormalization(View(image), View(res));
        checkVariance(res.upperLeft(), res.accessor(), 0.1);
        quadraticNoiseNormalization(View(image), View(res), 1.0, 0.02, 0.0);
        checkVariance(res.upperLeft(), res.accessor(), 0.1);
        quadraticNoiseNormalization(View(image), View(res));
        checkVariance(res.upperLeft(), res.accessor(), 0.1);
   }
  
    void testParametricNoiseNormalizationRGB()
    {
        RGBImage res(rgb.size());
        linearNoiseNormalization(srcImageRange(rgb), destImage(res), 1.0, 0.02);
        for(unsigned int band = 0; band < 3; ++band)
        {
            vigra::VectorElementAccessor<RGBImage::Accessor> dband(band, res.accessor());
            checkVariance(res.upperLeft(), dband, 0.1);
        }
        linearNoiseNormalization(srcImageRange(rgb), destImage(res));
        for(unsigned int band = 0; band < 3; ++band)
        {
            vigra::VectorElementAccessor<RGBImage::Accessor> dband(band, res.accessor());
            checkVariance(res.upperLeft(), dband, 0.1);
        }
        quadraticNoiseNormalization(srcImageRange(rgb), destImage(res), 1.0, 0.02, 0.0);
        for(unsigned int band = 0; band < 3; ++band)
        {
            vigra::VectorElementAccessor<RGBImage::Accessor> dband(band, res.accessor());
            checkVariance(res.upperLeft(), dband, 0.1);
        }
        quadraticNoiseNormalization(srcImageRange(rgb), destImage(res));
        for(unsigned int band = 0; band < 3; ++band)
        {
            vigra::VectorElementAccessor<RGBImage::Accessor> dband(band, res.accessor());
            checkVariance(res.upperLeft(), dband, 0.1);
        }
    }
  
    void testNonparametricNoiseNormalizationU8()
    {
        GrayImage res(image.size());        
        nonparametricNoiseNormalization(srcImageRange(u8image), destImage(res));
        checkVariance(res.upperLeft(), res.accessor(), 0.1);
    }
  
    void testNonparametricNoiseNormalization()
    {
        GrayImage res(image.size());        
        nonparametricNoiseNormalization(srcImageRange(image), destImage(res));
        checkVariance(res.upperLeft(), res.accessor(), 0.1);
   }
  
    void testNonparametricNoiseNormalizationRGB()
    {
        RGBImage res(rgb.size());        
        nonparametricNoiseNormalization(srcImageRange(rgb), destImage(res));
        for(unsigned int band = 0; band < 3; ++band)
        {
            vigra::VectorElementAccessor<RGBImage::Accessor> dband(band, res.accessor());
            checkVariance(res.upperLeft(), dband, 0.1);
        }
    }
};

#ifdef HasFFTW3
struct SlantedEdgeMTFTest
{
    typedef vigra::MultiArray<2, double> Image;
    typedef vigra::ArrayVector<vigra::TinyVector<double, 2> > Result;
    typedef Result::value_type Pair;
    
    Image image;
    Result reference;

    SlantedEdgeMTFTest()
    {
        vigra::ImageImportInfo info("slantedEdgeMTF.xv");
           
        image.reshape(info.shape());
        importImage(info, destImage(image));
        
        reference.push_back(Pair(0, 1));
        reference.push_back(Pair(0.0564351, 0.981739));
        reference.push_back(Pair(0.11287, 0.929577));
        reference.push_back(Pair(0.169305, 0.850509));
        reference.push_back(Pair(0.22574, 0.754266));
        reference.push_back(Pair(0.282175, 0.651081));
        reference.push_back(Pair(0.33861, 0.549492));
        reference.push_back(Pair(0.395045, 0.454718));
        reference.push_back(Pair(0.45148, 0.368628));
        reference.push_back(Pair(0.507915, 0.291512));
        reference.push_back(Pair(0.564351, 0.223585));
        reference.push_back(Pair(0.620786, 0.165506));
        reference.push_back(Pair(0.677221, 0.117569));
        reference.push_back(Pair(0.733656, 0.0788625));
        reference.push_back(Pair(0.790091, 0.0475224));
        reference.push_back(Pair(0.846526, 0.021419));
        reference.push_back(Pair(0.902961, 0));
    }
    
    void testSlantedEdgeMTF()
    {
        Result res;
        slantedEdgeMTF(image, res);
        
        shouldEqual(res.size(), reference.size());
        
        for(unsigned int k = 0; k < res.size(); ++k)
        {
            shouldEqualTolerance(res[k][0], reference[k][0], 1e-5);
            shouldEqualTolerance(res[k][1], reference[k][1], 1e-5);
        }
        
        shouldEqualTolerance(mtfFitGaussian(res), 0.5, 1e-2);
    }
};
#endif // HasFFTW3

struct AffineRegistrationTest
{
    typedef vigra::DImage Image;
    typedef vigra::MultiArrayView<2, double> View;
    typedef vigra::TinyVector<double, 2> Vector2;
    typedef vigra::ArrayVector<Vector2> PointList;
    
    Image image;

    AffineRegistrationTest()
    {
        ImageImportInfo info("lenna128.xv");
        image.resize(info.size());
        importImage(info, destImage(image));
    }
    
    void testCorrespondingPoints()
    {
        Matrix<double> point(3,1), res(3,1);
        point(2,0) = 1.0;
        
        PointList src(3), dest(3);
        src[0] = Vector2(1.6, 2.9);     
        src[1] = Vector2(-1.3, -4.0);     
        src[2] = Vector2(7.1, -0.4);     
        dest[0] = Vector2(-2.9, 1.6);     
        dest[1] = Vector2(12.6, -3.4);     
        dest[2] = Vector2(-3.3, 4.2);
        
        for(int k=1; k<=3; ++k)
        {     
            Matrix<double> a = affineMatrix2DFromCorrespondingPoints(src.begin(), src.begin()+k, dest.begin());
            for(int i=0; i<k; ++i)
            {
                point(0,0) = src[i][0];
                point(1,0) = src[i][1];
                res = a * point;
                shouldEqualTolerance(res(0,0), dest[i][0], 1e-14);
                shouldEqualTolerance(res(1,0), dest[i][1], 1e-14);
                shouldEqual(res(2,0), 1.0);
            }
        }
    }

    void testTranslationRegistration()
    {
        Matrix<double> m = translationMatrix2D(Vector2(5.0, 10.0));
        Image timg(image.size());
        affineWarpImage(SplineImageView<2, double>(srcImageRange(image)), destImageRange(timg), m);
        
        Matrix<double> estimated = identityMatrix<double>(3);
        estimateTranslation(srcImageRange(image), srcImageRange(timg), estimated,
                            AffineMotionEstimationOptions<1>().highestPyramidLevel(3));
        
        for(int i=0; i<9; ++i)
            shouldEqualTolerance(m.data()[i] - estimated.data()[i], 0.0, 1e-6);
        
        estimated = identityMatrix<double>(3);
        estimateTranslation(View(image), View(timg), estimated,
                            AffineMotionEstimationOptions<1>().highestPyramidLevel(3));
        
        for(int i=0; i<9; ++i)
            shouldEqualTolerance(m.data()[i] - estimated.data()[i], 0.0, 1e-6);
    }

    void testSimilarityRegistration()
    {
        Matrix<double> m = translationMatrix2D(Vector2(5.0, 10.0)) * 
                             rotationMatrix2DDegrees(5.0)* scalingMatrix2D(0.9);
        Image timg(image.size());
        affineWarpImage(SplineImageView<2, double>(srcImageRange(image)), destImageRange(timg), m);

        Matrix<double> estimated = identityMatrix<double>(3);
        estimateSimilarityTransform(srcImageRange(image), srcImageRange(timg), estimated,
                            AffineMotionEstimationOptions<>().useLaplacianPyramid(false));
        
        for(int i=0; i<9; ++i)
            shouldEqualTolerance(m.data()[i] - estimated.data()[i], 0.0, 1e-6);

        estimated = identityMatrix<double>(3);
        estimateSimilarityTransform(srcImageRange(image), 
                            srcIterRange(timg.upperLeft(), timg.lowerRight()-Diff2D(20,20)), estimated,
                            AffineMotionEstimationOptions<>().useLaplacianPyramid(true));
        
        for(int i=0; i<9; ++i)
            shouldEqualTolerance(m.data()[i] , estimated.data()[i], 1e-2);

        estimated = identityMatrix<double>(3);
        estimateSimilarityTransform(View(image),View(timg).subarray(Shape2(), Shape2(-20)), estimated,
                            AffineMotionEstimationOptions<>().useLaplacianPyramid(true));
        
        for(int i=0; i<9; ++i)
            shouldEqualTolerance(m.data()[i] , estimated.data()[i], 1e-2);
    }

    void testAffineRegistration()
    {
        Matrix<double> m = translationMatrix2D(Vector2(5.0, 10.0)) * 
                             rotationMatrix2DDegrees(5.0)* scalingMatrix2D(1.0, 0.9);
        Image timg(image.size());
        affineWarpImage(SplineImageView<2, double>(srcImageRange(image)), destImageRange(timg), m);

        Matrix<double> estimated = identityMatrix<double>(3);
        estimateAffineTransform(srcImageRange(image), srcImageRange(timg), estimated);
        
        for(int i=0; i<9; ++i)
            shouldEqualTolerance(m.data()[i] - estimated.data()[i], 0.0, 1e-6);

        estimated = identityMatrix<double>(3);
        estimateAffineTransform(View(image), View(timg), estimated);
        
        for(int i=0; i<9; ++i)
            shouldEqualTolerance(m.data()[i] - estimated.data()[i], 0.0, 1e-6);
    }
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
        add( testCase( &LabelingTest::labelingToEdgeTest));
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
        add( testCase( &EdgeDetectionTest::cannyEdgelList3x3Test));
        add( testCase( &EdgeDetectionTest::cannyEdgeImageTest));
        add( testCase( &EdgeDetectionTest::cannyEdgeImageWithThinningTest));
        add( testCase( &DistanceTransformTest::distanceTransformL1Test));
        add( testCase( &DistanceTransformTest::distanceTransformL2Test));
        add( testCase( &DistanceTransformTest::distanceTransformLInfTest));

        add( testCase( &LocalMinMaxTest::localMinimumTest));
        add( testCase( &LocalMinMaxTest::localMinimum4Test));
        add( testCase( &LocalMinMaxTest::localMinimumTestThr));
        add( testCase( &LocalMinMaxTest::localMaximumTest));
        add( testCase( &LocalMinMaxTest::localMaximum4Test));
        add( testCase( &LocalMinMaxTest::localMaximumTestThr));
        add( testCase( &LocalMinMaxTest::extendedLocalMinimumTest));
        add( testCase( &LocalMinMaxTest::extendedLocalMinimum4Test));
        add( testCase( &LocalMinMaxTest::extendedLocalMaximumTest));
        add( testCase( &LocalMinMaxTest::extendedLocalMaximum4Test));

        add( testCase( &LocalMinMaxTest::extendedLocalMaximum3DTest));
        add( testCase( &LocalMinMaxTest::extendedLocalMinimum3DTest));
        add( testCase( &LocalMinMaxTest::extendedLocalMaximum3DTest2));
        add( testCase( &LocalMinMaxTest::extendedLocalMinimum3DTest2));
        add( testCase( &LocalMinMaxTest::localMaximum3DTest));
        add( testCase( &LocalMinMaxTest::localMinimum3DTest));

        add( testCase( &LocalMinMaxTest::plateauWithHolesTest));
        add( testCase( &WatershedsTest::watershedsTest));
        add( testCase( &WatershedsTest::watersheds4Test));
        add( testCase( &RegionGrowingTest::voronoiTest));
        add( testCase( &RegionGrowingTest::voronoiWithBorderTest));
        add( testCase( &InterestOperatorTest::cornerResponseFunctionTest));
        add( testCase( &InterestOperatorTest::foerstnerCornerTest));
        add( testCase( &InterestOperatorTest::rohrCornerTest));
        add( testCase( &InterestOperatorTest::beaudetCornerTest));
        add( testCase( &InterestOperatorTest::radialSymmetryTest));
        add( testCase( &NoiseNormalizationTest::testParametricNoiseNormalization));
        add( testCase( &NoiseNormalizationTest::testNonparametricNoiseNormalization));
        add( testCase( &NoiseNormalizationTest::testParametricNoiseNormalizationU8));
        add( testCase( &NoiseNormalizationTest::testNonparametricNoiseNormalizationU8));
        add( testCase( &NoiseNormalizationTest::testParametricNoiseNormalizationRGB));
        add( testCase( &NoiseNormalizationTest::testNonparametricNoiseNormalizationRGB));
        add( testCase( &AffineRegistrationTest::testCorrespondingPoints));
        add( testCase( &AffineRegistrationTest::testTranslationRegistration));
        add( testCase( &AffineRegistrationTest::testSimilarityRegistration));
        add( testCase( &AffineRegistrationTest::testAffineRegistration));
#ifdef HasFFTW3
        add( testCase( &SlantedEdgeMTFTest::testSlantedEdgeMTF));
#endif // HasFFTW3
    }
};

int main(int argc, char ** argv)
{
    SimpleAnalysisTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

