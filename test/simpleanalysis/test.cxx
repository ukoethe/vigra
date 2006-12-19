/************************************************************************/
/*                                                                      */
/*                 Copyright 2004 by Ullrich Koethe                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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
#include "vigra/seededregiongrowing.hxx"
#include "vigra/cornerdetection.hxx"
#include "vigra/symmetry.hxx"
#include "vigra/watersheds.hxx"
#include "vigra/noise_normalization.hxx"
#include "vigra/slanted_edge_mtf.hxx"
#include "vigra/affinegeometry.hxx"
#include "vigra/affine_registration.hxx"
#include "vigra/impex.hxx"

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
        std::vector<vigra::Edgel> edgels3x3;
        cannyEdgelList3x3(srcImageRange(imgCanny), edgels3x3, 1.0);
        int count = 0;
        for(unsigned int i=0; i<edgels.size(); ++i)
        {
            if (edgels[i].strength < 1.0e-10)
                continue;  // ignore edgels that result from round off error during convolution
            ++count;
            should(edgels[i].x == edgels[i].y);
            should(VIGRA_CSTD::fabs(edgels[i].orientation-M_PI*0.75) < 0.1);
        }
        should(count == 75);
    }

    void cannyEdgelList3x3Test()
    {
        std::vector<vigra::Edgel> edgels;
        cannyEdgelList3x3(srcImageRange(imgCanny), edgels, 1.0);
        int count = 0;
        for(unsigned int i=0; i<edgels.size(); ++i)
        {
            if (edgels[i].strength < 1.0e-10)
                continue;  // ignore edgels that result from round off error during convolution
            ++count;
            should(edgels[i].x == edgels[i].y);
            should(VIGRA_CSTD::fabs(edgels[i].orientation-M_PI*0.75) < 0.1);
        }
        should(count == 38);
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

    void cannyEdgeImageWithThinningTest()
    {
        vigra::BImage result(40, 40);
        result = 0;

        cannyEdgeImageWithThinning(srcImageRange(imgCanny), destImage(result), 1.0, 0.1, 1, false);

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
    typedef vigra::DImage Image;

    LocalMinMaxTest()
    : img(9,9)
    {
        static const double in[] = {
            0.0,  0.1,  0.1,  0.3,  0.5,  0.3,  0.0,  0.0, 0.0,
            0.0, -0.1,  0.1,  0.0,  1.0,  0.0,  0.3,  0.0, 0.0,
            0.0,  0.5,  2.0,  0.0,  2.0,  2.0,  2.0,  0.0, 0.0,
            0.0,  0.0,  1.0,  1.5,  1.0,  1.0,  0.0,  0.0, 0.0,
            0.0,  0.1,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 0.0,
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
        res.init(0);

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

        shouldEqualSequence(res.begin(), res.end(), desired);
    }

    void localMinimum4Test()
    {
        Image res(img);
        res.init(0);

        localMinima(srcImageRange(img), destImage(res), 1.0, FourNeighborCode());

        static const double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        shouldEqualSequence(res.begin(), res.end(), desired);
    }

    void localMaximumTest()
    {
        Image res(img);
        res.init(0);

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

        shouldEqualSequence(res.begin(), res.end(), desired);
    }

    void localMaximum4Test()
    {
        Image res(img);
        res.init(0);

        localMaxima(srcImageRange(img), destImage(res), 1.0, FourNeighborCode());

        static const double desired[] = {
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
    }

    void extendedLocalMinimumTest()
    {
        Image res(img);
        res.init(0);

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

        shouldEqualSequence(res.begin(), res.end(), desired);
    }

    void extendedLocalMinimum4Test()
    {
        Image res(img);
        res.init(0);

        extendedLocalMinima(srcImageRange(img), destImage(res), 1.0, FourNeighborCode());

        static const double desired[] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        shouldEqualSequence(res.begin(), res.end(), desired);
    }

    void extendedLocalMaximumTest()
    {
        Image res(img);
        res.init(0);

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

        shouldEqualSequence(res.begin(), res.end(), desired);
   }

    void extendedLocalMaximum4Test()
    {
        Image res(img);
        res.init(0);

        extendedLocalMaxima(srcImageRange(img), destImage(res), 1.0, FourNeighborCode());

        static const double desired[] = {
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
        Image res(img.size(), 0.0);

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
};

struct WatershedsTest
{
    typedef vigra::DImage Image;

    WatershedsTest()
    : img(9,9)
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

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator end = img.end();
        Image::Accessor acc = img.accessor();
        const double * p = in;

        for(; i != end; ++i, ++p)
        {
            acc.set(*p, i);
        }
    }

    void watershedsTest()
    {
        Image res(img);
        res.init(0);

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

        int count = watersheds(srcImageRange(img), destImage(res));

#if 0
        std::cerr << count << "\n";
        for(int y=0;y<9;++y)
        {
            for(int x=0;x<9;++x)
                std::cerr << res(x,y) << "     ";
            std::cerr << "\n";
        }
#endif /* #if 0 */

        
        should(7 == count);

        shouldEqualSequence(res.begin(), res.end(), desired);
    }

    void watersheds4Test()
    {
        Image res(img);
        res.init(0);

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

        int count = watersheds(srcImageRange(img), destImage(res), FourNeighborCode());

#if 0
        std::cerr << count << "\n";
        for(int y=0;y<9;++y)
        {
            for(int x=0;x<9;++x)
                std::cerr << res(x,y) << "     ";
            std::cerr << "\n";
        }
#endif /* #if 0 */

        should(10 == count);

        shouldEqualSequence(res.begin(), res.end(), desired);
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

                if(VIGRA_CSTD::fabs(dist1 - dist2) > 1e-10)
                    shouldEqual(dist, desired);
            }
        }
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

struct NoiseNormalizationTest
{
    typedef vigra::BImage U8Image;
    typedef vigra::DImage GrayImage;
    typedef vigra::DRGBImage RGBImage;
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
        linearNoiseNormalization(srcImageRange(image), destImage(res), 1.0, 0.02);
        checkVariance(res.upperLeft(), res.accessor(), 0.1);
        linearNoiseNormalization(srcImageRange(image), destImage(res));
        checkVariance(res.upperLeft(), res.accessor(), 0.1);
        quadraticNoiseNormalization(srcImageRange(image), destImage(res), 1.0, 0.02, 0.0);
        checkVariance(res.upperLeft(), res.accessor(), 0.1);
        quadraticNoiseNormalization(srcImageRange(image), destImage(res));
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

struct SlantedEdgeMTFTest
{
    typedef vigra::DImage Image;
    typedef vigra::ArrayVector<vigra::TinyVector<double, 2> > Result;
    typedef Result::value_type Pair;
    
    Image image;
    Result reference;

    SlantedEdgeMTFTest()
    {
        vigra::ImageImportInfo info("slantedEdgeMTF.xv");
           
        image.resize(info.size());
        importImage(info, destImage(image));
        
        reference.push_back(Pair(0, 1));
        reference.push_back(Pair(0.0564351, 0.961336));
        reference.push_back(Pair(0.11287, 0.905812));
        reference.push_back(Pair(0.169305, 0.831278));
        reference.push_back(Pair(0.225741, 0.740566));
        reference.push_back(Pair(0.282176, 0.641225));
        reference.push_back(Pair(0.338611, 0.541988));
        reference.push_back(Pair(0.395046, 0.448508));
        reference.push_back(Pair(0.451481, 0.363593));
        reference.push_back(Pair(0.507916, 0.287531));
        reference.push_back(Pair(0.564351, 0.220531));
        reference.push_back(Pair(0.620787, 0.163246));
        reference.push_back(Pair(0.677222, 0.115963));
        reference.push_back(Pair(0.733657, 0.0777855));
        reference.push_back(Pair(0.790092, 0.0468734));
        reference.push_back(Pair(0.846527, 0.0211264));
        reference.push_back(Pair(0.902962, 0));
    }
    
    void testSlantedEdgeMTF()
    {
        Result res;
        slantedEdgeMTF(srcImageRange(image), res);
        
        shouldEqual(res.size(), reference.size());
        
        for(unsigned int k = 0; k < res.size(); ++k)
        {
            shouldEqualTolerance(res[k][0], reference[k][0], 1e-5);
            shouldEqualTolerance(res[k][1], reference[k][1], 1e-5);
        }
        
        shouldEqualTolerance(mtfFitGaussian(res), 0.5, 1e-2);
    }
};

struct AffineRegistrationTest
{
    typedef vigra::DImage Image;
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
        add( testCase( &EdgeDetectionTest::cannyEdgelList3x3Test));
        add( testCase( &EdgeDetectionTest::cannyEdgeImageTest));
        add( testCase( &EdgeDetectionTest::cannyEdgeImageWithThinningTest));
        add( testCase( &DistanceTransformTest::distanceTransformL1Test));
        add( testCase( &DistanceTransformTest::distanceTransformL2Test));
        add( testCase( &DistanceTransformTest::distanceTransformLInfTest));
        add( testCase( &LocalMinMaxTest::localMinimumTest));
        add( testCase( &LocalMinMaxTest::localMinimum4Test));
        add( testCase( &LocalMinMaxTest::localMaximumTest));
        add( testCase( &LocalMinMaxTest::localMaximum4Test));
        add( testCase( &LocalMinMaxTest::extendedLocalMinimumTest));
        add( testCase( &LocalMinMaxTest::extendedLocalMinimum4Test));
        add( testCase( &LocalMinMaxTest::extendedLocalMaximumTest));
        add( testCase( &LocalMinMaxTest::extendedLocalMaximum4Test));
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
        add( testCase( &SlantedEdgeMTFTest::testSlantedEdgeMTF));
        add( testCase( &AffineRegistrationTest::testCorrespondingPoints));
        add( testCase( &AffineRegistrationTest::testTranslationRegistration));
        add( testCase( &AffineRegistrationTest::testSimilarityRegistration));
        add( testCase( &AffineRegistrationTest::testAffineRegistration));
    }
};

int main()
{
    SimpleAnalysisTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

