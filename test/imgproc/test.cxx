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

#include <iostream>
#include <functional>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include "unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/stdimagefunctions.hxx"
#include "vigra/splineimageview.hxx"
#include "vigra/basicgeometry.hxx"
#include "vigra/affinegeometry.hxx"
#include "vigra/impex.hxx"
#include "vigra/meshgrid.hxx"
#include "vigra/multi_array.hxx"
#include "vigra/multi_math.hxx"
#include "vigra/functorexpression.hxx"

using namespace vigra;

struct Counter   // for initImageWithFunctor
{
    Counter() : count(0) {}
        
    int operator()() const { return count++; }
    
    mutable int count;
};


struct ImageFunctionsTest
{
    typedef vigra::DImage Image;
    typedef vigra::MultiArrayView<2, double> View;
    typedef vigra::DRGBImage RGBImage;
    typedef vigra::MultiArrayView<2, RGBValue<double>> RGBView;
    typedef Image::value_type GrayValue;
    typedef RGBImage::value_type RGBValue;

    ImageFunctionsTest()
    : img(3,3), mask(Shape2(3,3)), rgb(3,3), col(1.0, 2.0, 3.0)
    {
        Image::Accessor acc = img.accessor();
        Image::ScanOrderIterator i = img.begin();

        acc.set(1.1, i);
        ++i;
        acc.set(2.2, i);
        ++i;
        acc.set(3.3, i);
        ++i;
        acc.set(4.4, i);
        ++i;
        acc.set(5.5, i);
        ++i;
        acc.set(6.6, i);
        ++i;
        acc.set(7.7, i);
        ++i;
        acc.set(8.8, i);
        ++i;
        acc.set(9.9, i);
        ++i;
        should(i == img.end());

        mask.init(1);
        mask.begin()[0] = 0;
        mask.begin()[8] = 0;

        rgb.init(col);
    }

    void copyImageTest()
    {
        Image img1(3,3);

        copyImage(srcImageRange(img), destImage(img1));

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i, ++i1)
        {
            should(acc(i) == acc(i1));
        }

        img1 = 0.0;
        copyImage(View(img), View(img1));
        should(View(img) == View(img1));
    }

    void copyImageIfTest()
    {
        Image img1(3,3);
        img1 = 42.0;
        Image img2(3,3);
        img2 = 42.0;

        copyImageIf(srcImageRange(img), maskImage(mask), destImage(img1));
        copyImageIf(View(img), mask, View(img2));

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::ScanOrderIterator i1end = img1.end();
        Image::Accessor acc = img.accessor();

        should(acc(i1) == 42.0);
        ++i;
        ++i1;
        --i1end;
        should(acc(i1end) == 42.0);

        for(; i1 != i1end; ++i, ++i1)
        {
            should(acc(i) == acc(i1));
        }

        should(View(img1) == View(img2));
    }

    void copyRedBandTest()
    {
        vigra::RedAccessor<RGBValue> red;

        Image img1(3,3);
        copyImage(srcImageRange(rgb, red), destImage(img));
        copyImage(RGBView(rgb).expandElements(0).bindInner(0), View(img1));

        Image::ScanOrderIterator i = img.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i)
        {
            should(acc(i) == 1.0);
        }

        should(View(img) == View(img1));

        img = 42.0;
        col.setRed(42.0);

        copyImage(srcImageRange(img), destImage(rgb, red));

        RGBImage::ScanOrderIterator ri = rgb.begin();
        RGBImage::Accessor racc = rgb.accessor();

        for(; ri != rgb.end(); ++ri)
        {
            should(racc(ri) == col);
        }

        img1 = 12.0;
        col.setBlue(12.0);

        copyImage(View(img1), RGBView(rgb).expandElements(0).bindInner(2));
        using namespace multi_math;
        should(all(RGBView(rgb) == col));
    }

    void copyGreenBandTest()
    {
        vigra::GreenAccessor<RGBValue> green;

        copyImage(srcImageRange(rgb, green), destImage(img));

        Image::ScanOrderIterator i = img.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i)
        {
            should(acc(i) == 2.0);
        }

        img = 42.0;
        col.setGreen(42.0);

        copyImage(srcImageRange(img), destImage(rgb, green));

        RGBImage::ScanOrderIterator ri = rgb.begin();
        RGBImage::Accessor racc = rgb.accessor();

        for(; ri != rgb.end(); ++ri)
        {
            should(racc(ri) == col);
        }
    }

    void copyBlueBandTest()
    {
        vigra::BlueAccessor<RGBValue> blue;

        copyImage(srcImageRange(rgb, blue), destImage(img));

        Image::ScanOrderIterator i = img.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i)
        {
            should(acc(i) == 3.0);
        }

        img = 42.0;
        col.setBlue(42.0);

        copyImage(srcImageRange(img), destImage(rgb, blue));

        RGBImage::ScanOrderIterator ri = rgb.begin();
        RGBImage::Accessor racc = rgb.accessor();

        for(; ri != rgb.end(); ++ri)
        {
            should(racc(ri) == col);
        }
    }

    void rgbToGrayTest()
    {
        vigra::RGBToGrayAccessor<RGBValue> gray;

        copyImage(srcImageRange(rgb, gray), destImage(img));

        Image::ScanOrderIterator i = img.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i)
        {
            shouldEqualTolerance(acc(i), col.luminance(), 1e-12);
        }

    }

    void initImageTest()
    {
        initImage(destImageRange(img), 42.0);

        Image::ScanOrderIterator i = img.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i)
        {
            should(acc(i) == 42.0);
        }

        initImage(View(img), 12.0);

        i = img.begin();
        for(; i != img.end(); ++i)
        {
            should(acc(i) == 12.0);
        }

        initImageWithFunctor(View(img), Counter());

        i = img.begin();
        for(int k=0; i != img.end(); ++i, ++k)
        {
            should(acc(i) == k);
        }
    }

    void initImageIfTest()
    {
        initImageIf(destImageRange(img), maskImage(mask), 42.0);

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator iend = img.end();
        Image::Accessor acc = img.accessor();

        should(acc(i) == 1.1);
        ++i;
        --iend;
        should(acc(iend) == 9.9);
        for(; i != iend; ++i)
        {
            should(acc(i) == 42.0);
        }

        initImageIf(View(img), mask, 12.0);
        i = img.begin();
        should(acc(i) == 1.1);
        should(acc(iend) == 9.9);
        ++i;
        for(; i != iend; ++i)
        {
            should(acc(i) == 12.0);
        }
    }

    void initImageBorderTest()
    {
        initImageBorder(destImageRange(img), 1, 42.0);

        Image::ScanOrderIterator i = img.begin();
        Image::Accessor acc = img.accessor();
        int k;

        for(k=0; i != img.end(); ++i, ++k)
        {
            if(k != 4)
                should(acc(i) == 42.0);
            else
                should(acc(i) == 5.5);
        }

        initImageBorder(View(img), 1, 12.0);

        i = img.begin();
        for(k=0; i != img.end(); ++i, ++k)
        {
            if(k != 4)
                should(acc(i) == 12.0);
            else
                should(acc(i) == 5.5);
        }
    }

    void findMinMaxTest()
    {
        vigra::FindMinMax<Image::value_type> minmax, minmax1;

        inspectImage(srcImageRange(img), minmax);

        should(minmax.count == 9);
        should(minmax.max == 9.9);
        should(minmax.min == 1.1);

        inspectImage(View(img), minmax1);

        should(minmax1.count == 9);
        should(minmax1.max == 9.9);
        should(minmax1.min == 1.1);
    }

    void findMinMaxIfTest()
    {
        vigra::FindMinMax<Image::value_type> minmax, minmax1;

        inspectImageIf(srcImageRange(img), maskImage(mask), minmax);

        should(minmax.count == 7);
        should(minmax.max == 8.8);
        should(minmax.min == 2.2);

        inspectImageIf(View(img), mask, minmax1);

        should(minmax1.count == 7);
        should(minmax1.max == 8.8);
        should(minmax1.min == 2.2);
    }

    void findAverageTest()
    {
        vigra::FindAverage<Image::value_type> average;

        inspectImage(srcImageRange(img), average);

        should(average.count() == 9);
        should(average() == 5.5);

        average(average); // combine

        should(average.count() == 18);
        should(average() == 5.5);
    }

    void findAverageIfTest()
    {
        vigra::FindAverage<Image::value_type> average;

        inspectImageIf(srcImageRange(img), maskImage(mask), average);

        should(average.count() == 7);
        should(average() == 5.5);
    }

    void findAverageWeightedTest()
    {
        vigra::FindAverage<Image::value_type> average;

        average(10);
        average(5,  2.0);
        average(42, 0.2);

        shouldEqualTolerance(average.count(), 3.2, 1e-12);
        shouldEqualTolerance(average(), (10+10+42*0.2)/(1+2+0.2), 1e-8);

        average(average); // combine

        shouldEqualTolerance(average.count(), 6.4, 1e-12);
        shouldEqualTolerance(average(), (10+10+42*0.2)/(1+2+0.2), 1e-8);
    }

    void findAverageAndVarianceTest()
    {
        vigra::FindAverageAndVariance<Image::value_type> averageAndVariance;

        inspectImage(srcImageRange(img), averageAndVariance);

        should(averageAndVariance.count() == 9);
        shouldEqualTolerance(averageAndVariance.average(), 5.5, 1e-14);

        // compute variance explicitly
        double sumOfSquares = 0.0;
        Image::ScanOrderIterator i = img.begin();
        for(; i != img.end(); ++i)
        {
            sumOfSquares += sq(*i - averageAndVariance.average());
        }

        shouldEqualTolerance(averageAndVariance.variance(), sumOfSquares / 9.0, 1e-14);
        shouldEqualTolerance(averageAndVariance(), sumOfSquares / 9.0, 1e-14);
        shouldEqualTolerance(averageAndVariance.variance(true), sumOfSquares / 8.0, 1e-14);

        // check merge of functors
        vigra::FindAverageAndVariance<Image::value_type> averageAndVariance1(averageAndVariance), averageAndVariance2;

        averageAndVariance1(5.31);
        averageAndVariance1(-0.3);
        averageAndVariance1(2.88);
        averageAndVariance1(10.521);

        averageAndVariance2(5.31);
        averageAndVariance2(-0.3);
        averageAndVariance2(2.88);
        averageAndVariance2(10.521);
        averageAndVariance2(averageAndVariance);

        should(averageAndVariance1.count() == 13);
        should(averageAndVariance2.count() == 13);
        shouldEqualTolerance(averageAndVariance1.average(), averageAndVariance2.average(), 1e-14);
        shouldEqualTolerance(averageAndVariance1.variance(), averageAndVariance2.variance(), 1e-14);
    }

    void reduceFunctorTest()
    {
        std::plus<double> p;
        ReduceFunctor<std::plus<double>, double> f(p, 0.0);

        inspectImage(srcImageRange(img), f);

        shouldEqual(f(), 49.5);
    }

    void meshGridTest()
    {
        // create an image whose values are equal to each pixel's distance from the image center
        int width = 5, height = 7;
        int xc = width/2, yc = height/2;

        DImage dist(width, height);
        Point2D upperLeft(-xc, -yc);

        using namespace vigra::functor;
        transformImage(meshGrid(upperLeft, upperLeft+dist.size()),
                       destImage(dist),
                       norm(Arg1()));
        for(int y=0; y<height; ++y)
            for(int x=0; x<width; ++x)
                shouldEqualTolerance(dist(x,y), vigra::hypot(x-xc, y-yc), 1e-14);
        dist.init(0.0);
        transformImage(meshGrid(Rect2D(upperLeft, upperLeft+dist.size())),
                       destImage(dist),
                       norm(Arg1()));
        for(int y=0; y<height; ++y)
            for(int x=0; x<width; ++x)
                shouldEqualTolerance(dist(x,y), vigra::hypot(x-xc, y-yc), 1e-14);
    }


    void findBoundingRectangleTest()
    {
        vigra::FindBoundingRectangle findRect;

        mask = 0;
        mask(1,1) = 1;

        inspectImageIf(srcIterRange<Diff2D>(vigra::Diff2D(0,0), img.size()),
                       maskImage(mask), findRect);

        Rect2D rect = findRect();
        should(rect.upperLeft().x == 1);
        should(rect.upperLeft().y == 1);
        should(rect.lowerRight().x == 2);
        should(rect.lowerRight().y == 2);

        mask(1,0) = 1;

        inspectImageIf(srcIterRange<Diff2D>(vigra::Diff2D(0,0), img.size()),
                       maskImage(mask), findRect);

        rect = findRect();
        should(rect.upperLeft().x == 1);
        should(rect.upperLeft().y == 0);
        should(rect.lowerRight().x == 2);
        should(rect.lowerRight().y == 2);

        mask(0,1) = 1;

        inspectImageIf(srcIterRange<Diff2D>(vigra::Diff2D(0,0), img.size()),
                       maskImage(mask), findRect);

        rect = findRect();
        should(rect.upperLeft().x == 0);
        should(rect.upperLeft().y == 0);
        should(rect.lowerRight().x == 2);
        should(rect.lowerRight().y == 2);

        mask(1,2) = 1;

        inspectImageIf(srcIterRange<Diff2D>(vigra::Diff2D(0,0), img.size()),
                       maskImage(mask), findRect);

        rect = findRect();
        should(rect.upperLeft().x == 0);
        should(rect.upperLeft().y == 0);
        should(rect.lowerRight().x == 2);
        should(rect.lowerRight().y == 3);

        mask(2,1) = 1;

        inspectImageIf(srcIterRange<Diff2D>(vigra::Diff2D(0,0), img.size()),
                       maskImage(mask), findRect);

        rect = findRect();
        should(rect.upperLeft().x == 0);
        should(rect.upperLeft().y == 0);
        should(rect.lowerRight().x == 3);
        should(rect.lowerRight().y == 3);

        vigra::FindBoundingRectangle findRect1;
        findRect1(vigra::Diff2D(4,4));

        findRect(findRect1);

        rect = findRect();
        should(rect.upperLeft().x == 0);
        should(rect.upperLeft().y == 0);
        should(rect.lowerRight().x == 5);
        should(rect.lowerRight().y == 5);
    }

    void lastValueFunctorTest()
    {
        typedef vigra::LastValueFunctor<int> LastValue;
        {
            LastValue a;
            a(10);
            shouldEqual(a(), 10);
            a(vigra::NumericTraits<int>::max());
            shouldEqual(a(), vigra::NumericTraits<int>::max());
            a(42);
            shouldEqual(a(), 42);
        }
        {
            // try to construct at same place on stack as a
            // (check whether ints are initialized to zero)
            LastValue b;
            shouldEqual(b(), 0);
        }
        LastValue c(23); // test init. with value
        shouldEqual(c(), 23);
    }

    void arrayOfRegionStatisticsTest()
    {
        vigra::MultiArray<2, unsigned char> labels(Shape2(3,3));
        labels = 1;
        labels[0] = 0;
        labels[8] = 0;

        vigra::ArrayOfRegionStatistics<vigra::FindMinMax<Image::value_type> > stats(1);

        inspectTwoImages(View(img), labels, stats);

        should(stats[0].count == 2);
        should(stats[0].min == 1.1);
        should(stats[0].max == 9.9);
        should(stats[1].count == 7);
        should(stats[1].min == 2.2);
        should(stats[1].max == 8.8);

        vigra::ArrayOfRegionStatistics<vigra::FindBoundingRectangle> stats1(1);
        labels.init(0);
        labels(1,1) = 1;

        inspectTwoImages(srcIterRange<Diff2D>(vigra::Diff2D(0,0), img.size()),
                         srcImage(labels), stats1);

        should(stats1[0].upperLeft.x == 0);
        should(stats1[0].upperLeft.y == 0);
        should(stats1[0].lowerRight.x == 3);
        should(stats1[0].lowerRight.y == 3);
        should(stats1[1].upperLeft.x == 1);
        should(stats1[1].upperLeft.y == 1);
        should(stats1[1].lowerRight.x == 2);
        should(stats1[1].lowerRight.y == 2);
    }

    void arrayOfRegionStatisticsIfTest()
    {
        vigra::MultiArray<2, unsigned char> labels(Shape2(3,3));
        labels = 1;
        labels[4] = 2;
        labels[5] = 2;

        vigra::ArrayOfRegionStatistics<vigra::FindMinMax<Image::value_type> > stats(2), stats1(2);

        inspectTwoImagesIf(srcImageRange(img), srcImage(labels), maskImage(mask), stats);

        should(stats[0].count == 0);
        should(stats[1].count == 5);
        should(stats[1].min == 2.2);
        should(stats[1].max == 8.8);
        should(stats[2].count == 2);
        should(stats[2].min == 5.5);
        should(stats[2].max == 6.6);

        inspectTwoImagesIf(View(img), labels, mask, stats1);

        should(stats1[0].count == 0);
        should(stats1[1].count == 5);
        should(stats1[1].min == 2.2);
        should(stats1[1].max == 8.8);
        should(stats1[2].count == 2);
        should(stats1[2].min == 5.5);
        should(stats1[2].max == 6.6);
    }

    void writeArrayOfRegionStatisticsTest()
    {
        vigra::BImage labels(3,3);
        labels = 1;
        labels.begin()[0] = 0;
        labels.begin()[8] = 0;

        vigra::ArrayOfRegionStatistics<vigra::FindAverage<Image::value_type> > stats(1);

        inspectTwoImages(srcImageRange(img), srcImage(labels), stats);

        should(stats[0].count() == 2);
        should(stats[0]() == 5.5);
        should(stats[1].count() == 7);
        should(stats[1]() == 5.5);

        transformImage(srcImageRange(labels), destImage(img), stats);

        Image::ScanOrderIterator i = img.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i)
        {
            should(acc(i) == 5.5);
        }

    }

    void linearIntensityTransformTest()
    {
        vigra::LinearIntensityTransform<Image::value_type> trans(2.0, -1.1);

        Image img1(3,3);

        transformImage(View(img), View(img1), trans);

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i, ++i1)
        {
            shouldEqual(2.0*(acc(i) - 1.1), acc(i1));
        }

        RGBImage img2(3,3);

        transformImage(srcImageRange(rgb), destImage(img2),
                        linearIntensityTransform(2.0, RGBValue(0.0, -1.0, -2.0)));

        shouldEqual(img2(0,0), RGBValue(2.0, 2.0, 2.0));

        transformImage(srcImageRange(rgb), destImage(img2),
                        linearIntensityTransform(RGBValue(2.0, 3.0, 4.0),
                                                 RGBValue(0.0, -1.0, -2.0)));

        shouldEqual(img2(0,0), RGBValue(2.0, 3.0, 4.0));
    }

    void scalarIntensityTransformTest()
    {
        Image img1(3,3);

        transformImage(srcImageRange(img), destImage(img1),
                       linearIntensityTransform<Image::value_type>(3.3));

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i, ++i1)
        {
            shouldEqual(3.3*acc(i), acc(i1));
        }

        RGBImage img2(3,3);

        transformImage(srcImageRange(rgb), destImage(img2),
                        linearIntensityTransform<RGBValue>(
                                                 RGBValue(1.0, 2.0, 3.0)));

        shouldEqual(img2(0,0), RGBValue(1.0, 4.0, 9.0));
    }

    void linearIntensityTransformIfTest()
    {
        vigra::LinearIntensityTransform<Image::value_type> trans(2.0, -1.1);

        Image img1(3,3);
        img1 = 42.0;

        transformImageIf(View(img), mask, View(img1), trans);

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::ScanOrderIterator i1end = img1.end();
        Image::Accessor acc = img.accessor();

        should(acc(i1) == 42.0);
        ++i;
        ++i1;
        --i1end;
        should(acc(i1end) == 42.0);

        for(; i1 != i1end; ++i, ++i1)
        {
            shouldEqual(2.0*(acc(i) - 1.1), acc(i1));
        }

    }

    void linearRangeMappingTest()
    {
        BImage img1(3,3);

        transformImage(srcImageRange(img), destImage(img1),
                       linearRangeMapping(1.0, 10.0, 0, 250));

        int res[] = {3, 33, 64, 94, 125, 156, 186, 217, 247 };

        shouldEqualSequence(img1.begin(), img1.end(), res);

        FindMinMax<double> minmax;
        inspectImage(srcImageRange(img), minmax);
        transformImage(srcImageRange(img), destImage(img1),
                       linearRangeMapping(minmax, 0, 250));
    
        int res1[] = {0, 31, 62, 94, 125, 156, 187, 219, 250 };
#if 0
        for(int i=0; i<9; ++i)
            std::cerr << (int)*(&img1(0,0)+i) << ", ";
        std::cerr << "\n";
#endif
        shouldEqualSequence(img1.begin(), img1.end(), res1);

        BRGBImage img2(3,3);

        transformImage(srcImageRange(rgb), destImage(img2),
                       linearRangeMapping(RGBValue(1.0, 1.0, 1.0),
                                          RGBValue(3.0, 3.0, 3.0),
                                          BRGBImage::value_type(0, 0, 0),
                                          BRGBImage::value_type(255, 255, 255)));

        shouldEqual(img2(0,0), BRGBImage::value_type(0, 128, 255));
    }

    void thresholdTest()
    {
        vigra::Threshold<Image::value_type, Image::value_type> trans(2.0, 9.0, 0.0, 1.0);

        Image img1(3,3);

        transformImage(srcImageRange(img), destImage(img1), trans);

        Image::ScanOrderIterator i1 = img1.begin();
        Image::ScanOrderIterator i1end = img1.end();
        Image::Accessor acc = img.accessor();

        should(acc(i1) == 0.0);
        ++i1;
        --i1end;
        should(acc(i1end) == 0.0);

        for(; i1 != i1end; ++i1)
        {
            should(1.0 == acc(i1));
        }

    }

    void brightnessContrastTest()
    {
        Image img1(3,3);

        transformImage(srcImageRange(img), destImage(img1),
              vigra::BrightnessContrastFunctor<double>(1.0, 1.0, 0.0, 255.0));

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();

        for(; i != img.end(); ++i, ++i1)
        {
            should(VIGRA_CSTD::fabs(*i - *i1) < 1.0e-10);
        }

        vigra::BrightnessContrastFunctor<unsigned char> charf(10.0, 1.0);

        for(int k=1; k < 255; ++k)
        {
            should(k < charf(k));
        }
        should(0 == charf(0));
        should(255 == charf(255));

        vigra::BrightnessContrastFunctor<RGBValue > rgbf(2.0, 1.0, RGBValue(0.0), RGBValue(255.0));

        should(col.red() < rgbf(col).red());
        should(col.green() < rgbf(col).green());
        should(col.blue() < rgbf(col).blue());
    }

    void gradientFunctionTest()
    {
        RGBImage in(3,3);
        int y;

        for(y=0; y<3; ++y)
        {
            for(int x=0; x<3; ++x)
            {
                in(x,y) = RGBValue(float(x+y));
            }
        }

        Image res(3,3);

        gradientBasedTransform(
          srcImageRange(in, vigra::RedAccessor<RGBValue >()),
          destImage(res), vigra::MagnitudeFunctor<double>());

        for(y=0; y<3; ++y)
        {
            for(int x=0; x<3; ++x)
            {
                should(VIGRA_CSTD::fabs(res(x,y) - VIGRA_CSTD::sqrt(2.0)) < 1e-6);
            }
        }

        gradientBasedTransform(RGBView(in), View(res), 
                               vigra::RGBGradientMagnitudeFunctor<double>());

        for(y=0; y<3; ++y)
        {
            for(int x=0; x<3; ++x)
            {
                should(VIGRA_CSTD::fabs(res(x,y) - VIGRA_CSTD::sqrt(6.0)) < 1e-6);
            }
        }
    }

    void additionTest()
    {
        Image img1(3,3), img2(3,3), img3(3,3);

        std::plus<Image::value_type> add;

        combineTwoImages(srcImageRange(img), srcImage(img), destImage(img1), add);
        combineTwoImages(View(img), View(img), View(img2), add);

        using namespace functor;
        combineThreeImages(View(img), View(img), View(img1), View(img3), Arg1() + Arg2() + Arg3());

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::ScanOrderIterator i2 = img2.begin();
        Image::ScanOrderIterator i3 = img3.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i, ++i1, ++i2, ++i3)
        {
            should(2.0*acc(i) == acc(i1));
            should(2.0*acc(i) == acc(i2));
            should(4.0*acc(i) == acc(i3));
        }
    }

    void additionIfTest()
    {
        Image img1(3,3), img2(3,3);
        img1 = 42.0;
        img2 = 42.0;

        std::plus<Image::value_type> add;

        combineTwoImagesIf(srcImageRange(img), srcImage(img),
                           maskImage(mask), destImage(img1), add);

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::ScanOrderIterator i1end = img1.end();
        Image::Accessor acc = img.accessor();

        should(acc(i1) == 42.0);
        ++i;
        ++i1;
        --i1end;
        should(acc(i1end) == 42.0);

        for(; i1 != i1end; ++i, ++i1)
        {
            should(2.0*acc(i) == acc(i1));
        }

        combineTwoImagesIf(View(img), View(img), mask, View(img2), add);
        should(View(img1) == View(img2));
    }

    void resizeNoInterpolationTest()
    {
        Image img1(5,5), img3(5,5);

        resizeImageNoInterpolation(srcImageRange(img), destImageRange(img1));

        Image::ScanOrderIterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        should(acc(i1) == 1.1);
        ++i1;
        should(acc(i1) == 2.2);
        ++i1;
        should(acc(i1) == 2.2);
        ++i1;
        should(acc(i1) == 3.3);
        ++i1;
        should(acc(i1) == 3.3);
        ++i1;
        should(acc(i1) == 4.4);
        ++i1;
        should(acc(i1) == 5.5);
        ++i1;
        should(acc(i1) == 5.5);
        ++i1;
        should(acc(i1) == 6.6);
        ++i1;
        should(acc(i1) == 6.6);
        ++i1;
        should(acc(i1) == 4.4);
        ++i1;
        should(acc(i1) == 5.5);
        ++i1;
        should(acc(i1) == 5.5);
        ++i1;
        should(acc(i1) == 6.6);
        ++i1;
        should(acc(i1) == 6.6);
        ++i1;
        should(acc(i1) == 7.7);
        ++i1;
        should(acc(i1) == 8.8);
        ++i1;
        should(acc(i1) == 8.8);
        ++i1;
        should(acc(i1) == 9.9);
        ++i1;
        should(acc(i1) == 9.9);
        ++i1;
        should(acc(i1) == 7.7);
        ++i1;
        should(acc(i1) == 8.8);
        ++i1;
        should(acc(i1) == 8.8);
        ++i1;
        should(acc(i1) == 9.9);
        ++i1;
        should(acc(i1) == 9.9);

        Image img2(3,3);
        img2 = 0;

        resizeImageNoInterpolation(View(img1), View(img2));

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i2 = img2.begin();

        for(; i != img.end(); ++i, ++i2)
        {
            should(acc(i) == acc(i2));
        }

        resizeImageSplineInterpolation(View(img), View(img3), BSpline<0>());
        shouldEqualSequence(img1.begin(), img1.end(), img3.begin());
    }

    void resizeLinearInterpolationTest()
    {
        Image img1(5,5);
        resizeImageLinearInterpolation(srcImageRange(img), destImageRange(img1));

        Image::ScanOrderIterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        shouldEqualTolerance(acc(i1), 1.1, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 1.65, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 2.2, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 2.75, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 3.3, 1e-14);
        ++i1;

        shouldEqualTolerance(acc(i1), 2.75, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 3.3, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 3.85, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 4.4, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 4.95, 1e-14);
        ++i1;

        shouldEqualTolerance(acc(i1), 4.4, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 4.95, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 5.5, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 6.05, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 6.6, 1e-14);
        ++i1;

        shouldEqualTolerance(acc(i1), 6.05, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 6.6, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 7.15, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 7.7, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 8.25, 1e-14);
        ++i1;

        shouldEqualTolerance(acc(i1), 7.7, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 8.25, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 8.8, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 9.35, 1e-14);
        ++i1;
        shouldEqualTolerance(acc(i1), 9.9, 1e-14);

        Image img2(3,3);
        img2 = 0;

        resizeImageNoInterpolation(View(img1), View(img2));

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i2 = img2.begin();

        for(; i != img.end(); ++i, ++i2)
        {
            should(acc(i) == acc(i2));
        }
    }


    Image img;
    vigra::MultiArray<2, unsigned char> mask;
    RGBImage rgb;
    RGBValue col;
};

struct ResizeImageTest
{
    typedef vigra::FImage Image;
    typedef vigra::FRGBImage RGBImage;
    typedef MultiArrayView<2, RGBValue<float> > View;

    ResizeImageTest()
    {
        ImageImportInfo ginfo("lenna128.xv");
        img.resize(ginfo.width(), ginfo.height());
        importImage(ginfo, destImage(img));

        ImageImportInfo cinfo("lenna128rgb.xv");
        rgb.resize(cinfo.width(), cinfo.height());
        importImage(cinfo, destImage(rgb));

    }

    void resizeLinearInterpolationReduceTest()
    {
        ImageImportInfo inforef("lenna42lin.xv");
        Image ref(inforef.size());
        importImage(inforef, destImage(ref));

        Image dest(inforef.size());

        resizeImageLinearInterpolation(srcImageRange(img), destImageRange(dest));

        shouldEqualSequenceTolerance(dest.begin(), dest.end(), ref.begin(), 1.0e-6f);

        ImageImportInfo inforgb("lenna42linrgb.xv");
        RGBImage rgbref(inforgb.size());
        importImage(inforgb, destImage(rgbref));

        RGBImage rgbdest(inforgb.size());

        resizeImageLinearInterpolation(View(rgb), View(rgbdest));

        shouldEqualSequenceTolerance(rgbdest.begin(), rgbdest.end(), rgbref.begin(),
                                     RGBImage::value_type(1.0e-6f));
    }

    void scalarExpand()
    {
        ImageImportInfo info("lenna288neu.xv");

        Image imgex(info.width(), info.height());
        importImage(info, destImage(imgex));

        Image img1(imgex.width(), imgex.height());
        resizeImageSplineInterpolation(srcImageRange(img), destImageRange(img1));

        shouldEqualSequenceTolerance(img1.begin(), img1.end(), imgex.begin(), 1e-4f);
    }

    void scalarReduce()
    {
        ImageImportInfo info("lenna42neu.xv");

        Image imgred(info.width(), info.height());
        importImage(info, destImage(imgred));

        Image img1(imgred.width(), imgred.height());

        resizeImageSplineInterpolation(srcImageRange(img), destImageRange(img1));

        shouldEqualSequenceTolerance(img1.begin(), img1.end(), imgred.begin(), 1e-4f);
   }

    void rgbExpand()
    {
        ImageImportInfo info("lenna288rgbneu.xv");

        RGBImage rgbex(info.width(), info.height());
        importImage(info, destImage(rgbex));

        RGBImage img1(rgbex.width(), rgbex.height());

        resizeImageSplineInterpolation(View(rgb), View(img1));

        RGBImage::ScanOrderIterator i1 = img1.begin();
        RGBImage::ScanOrderIterator iex = rgbex.begin();
        RGBImage::Accessor acc = img1.accessor();

        for(; i1 != img1.end(); ++i1, ++iex)
        {
            shouldEqualTolerance(acc.red(i1), acc.red(iex), 1e-4f);
            shouldEqualTolerance(acc.green(i1), acc.green(iex), 1e-4f);
            shouldEqualTolerance(acc.blue(i1), acc.blue(iex), 1e-4f);
        }
    }

    void rgbReduce()
    {
        ImageImportInfo info("lenna42rgbneu.xv");

        RGBImage rgbred(info.width(), info.height());
        importImage(info, destImage(rgbred));

        RGBImage img1(rgbred.width(), rgbred.height());

        resizeImageSplineInterpolation(View(rgb), View(img1));

        RGBImage::ScanOrderIterator i1 = img1.begin();
        RGBImage::ScanOrderIterator iref = rgbred.begin();
        RGBImage::Accessor acc = img1.accessor();

        for(; i1 != img1.end(); ++i1, ++iref)
        {
            shouldEqualTolerance(acc.red(i1), acc.red(iref), 1e-4f);
            shouldEqualTolerance(acc.green(i1), acc.green(iref), 1e-4f);
            shouldEqualTolerance(acc.blue(i1), acc.blue(iref), 1e-4f);
        }
    }

    /*Vergroesserungstest
    */
    void testCatmullRomInterpolationExtensionWithLena()
    {
        ImageImportInfo inforef("lenna367FIR.xv");
        Image ref(inforef.size());
        importImage(inforef, destImage(ref));

        Image dest(inforef.size());

        resizeImageCatmullRomInterpolation(srcImageRange(img), destImageRange(dest));

        shouldEqualSequenceTolerance(dest.begin(), dest.end(), ref.begin(), 1.0e-6f);
    }


    /*Verkleinerungstest
    */
    void testCatmullRomInterpolationReductionWithLena()
    {
        ImageImportInfo inforef("lenna42FIR.xv");
        Image ref(inforef.size());
        importImage(inforef, destImage(ref));

        Image dest(inforef.size());

        resizeImageCatmullRomInterpolation(srcImageRange(img), destImageRange(dest));

        shouldEqualSequenceTolerance(dest.begin(), dest.end(), ref.begin(), 1e-4f);

        ImageImportInfo inforgb("lennargb42FIR.xv");
        RGBImage rgbref(inforgb.size());
        importImage(inforgb, destImage(rgbref));

        RGBImage rgbdest(inforgb.size());

        resizeImageCatmullRomInterpolation(View(rgb), View(rgbdest));

        RGBImage::ScanOrderIterator i1 = rgbdest.begin();
        RGBImage::ScanOrderIterator iref = rgbref.begin();
        RGBImage::Accessor acc = rgbdest.accessor();

        for(; i1 != rgbdest.end(); ++i1, ++iref)
        {
            shouldEqualTolerance(acc.red(i1), acc.red(iref), 1e-4f);
            shouldEqualTolerance(acc.green(i1), acc.green(iref), 1e-4f);
            shouldEqualTolerance(acc.blue(i1), acc.blue(iref), 1e-4f);
        }
    }

    void testCubicInterpolationExtensionWithLena()
    {
        ImageImportInfo inforef("lenna367IIR.xv");
        Image ref(inforef.size());
        importImage(inforef, destImage(ref));

        Image dest(inforef.size());

        resizeImageSplineInterpolation(srcImageRange(img), destImageRange(dest));

        shouldEqualSequenceTolerance(dest.begin(), dest.end(), ref.begin(), 1e-4f);
    }

    /*Verkleinerungstest
    */
    void testCubicInterpolationReductionWithLena()
    {
        ImageImportInfo inforef("lenna42IIR.xv");
        Image ref(inforef.size());
        importImage(inforef, destImage(ref));

        Image dest(inforef.size());

        resizeImageSplineInterpolation(srcImageRange(img), destImageRange(dest));

        shouldEqualSequenceTolerance(dest.begin(), dest.end(), ref.begin(), 1e-4f);

        ImageImportInfo inforgb("lennargb42IIR.xv");
        RGBImage rgbref(inforgb.size());
        importImage(inforgb, destImage(rgbref));

        RGBImage rgbdest(inforgb.size());

        resizeImageSplineInterpolation(srcImageRange(rgb), destImageRange(rgbdest));

        RGBImage::ScanOrderIterator i1 = rgbdest.begin();
        RGBImage::ScanOrderIterator iref = rgbref.begin();
        RGBImage::Accessor acc = rgbdest.accessor();

        for(; i1 != rgbdest.end(); ++i1, ++iref)
        {
            shouldEqualTolerance(acc.red(i1), acc.red(iref), 1e-4f);
            shouldEqualTolerance(acc.green(i1), acc.green(iref), 1e-4f);
            shouldEqualTolerance(acc.blue(i1), acc.blue(iref), 1e-4f);
        }
    }

    void testCatmullRomInterpolationExtensionHandControled()
    {
        vigra::DImage src(6, 7), dest(10, 10, 145.346);

        for(int i = 0; i<src.width()*src.height(); i++)
        {
          src.begin()[i] = 0.25 + i;
        }

        static double refdata[100] =
        {           0.25,/**/    0.695816186557,   1.36111111111,       1.91666666667,      2.47222222222,    3.02777777778,    3.58333333333,      4.13888888889,      4.80418381344,/***/       5.25,/**/
           3.80555555556,        4.25137174211,    4.91666666667,       5.47222222222,      6.02777777778,    6.58333333333,    7.13888888889,      7.69444444444,      8.359739369,              8.80555555556,
                    8.25,        8.69581618656,    9.36111111111,       9.91666666667,      10.4722222222,    11.0277777778,    11.5833333333,      12.1388888889,      12.8041838134,            13.25,
                   12.25,/**/    12.6958161866,    13.3611111111,/***/  13.9166666667,      14.4722222222,    15.0277777778,    15.5833333333,      16.1388888889,      16.8041838134,            17.25,/**/
                   16.25,/***/   16.6958161866,    17.3611111111,       17.9166666667,/***/ 18.4722222222,    19.0277777778,    19.5833333333,      20.1388888889,      20.8041838134,            21.25,
                   20.25,        20.6958161866,    21.3611111111,       21.9166666667,      22.4722222222,    23.0277777778,    23.5833333333,      24.1388888889,      24.8041838134,            25.25,/***/
                   24.25,/**/    24.6958161866,    25.3611111111,       25.9166666667,      26.4722222222,    27.0277777778,    27.5833333333,      28.1388888889,      28.8041838134,            29.25,/**/
                   28.25,        28.6958161866,    29.3611111111,       29.9166666667,      30.4722222222,    31.0277777778,    31.5833333333,/***/ 32.1388888889,      32.8041838134,            33.25,
           32.6944444444,        33.140260631,     33.8055555556,       34.3611111111,      34.9166666667,    35.4722222222,    36.0277777778,      36.5833333333,/***/ 37.2486282579,            37.6944444444,
                   36.25,/**/    36.6958161866,    37.3611111111,/***/  37.9166666667,      38.4722222222,    39.0277777778,    39.5833333333,      40.1388888889,      40.8041838134,            41.25/**/
        };

        resizeImageCatmullRomInterpolation(srcImageRange(src), destImageRange(dest));

        shouldEqualSequenceTolerance(dest.begin(), dest.end(), refdata, 1e-11);
    }

    Image img;
    RGBImage rgb;
};

template <int N>
struct SplineImageViewTest
{
    typedef vigra::DImage Image;
    Image img;

    SplineImageViewTest()
    {
        ImageImportInfo ginfo("lenna128.xv");
        img.resize(ginfo.width(), ginfo.height());
        importImage(ginfo, destImage(img));
    }

    void testPSF()
    {
        int center = 10;
        Image img(2*center+1, 2*center+1);
        img.init(0.0);
        img(center, center) = 1.0;
        SplineImageView<N, double> view(srcImageRange(img), true);
        BSplineBase<N> spline;

        double d0 = center;
        double epsilon = 1.0e-10;
        shouldEqualTolerance(view(d0,d0), spline(0.0)*spline(0.0), epsilon);
        shouldEqualTolerance(view(d0,d0, 1, 0), spline(0.0, 1)*spline(0.0), epsilon);
        shouldEqualTolerance(view(d0,d0, 0, 1), spline(0.0, 1)*spline(0.0), epsilon);
        for(double d = 0.2; d < spline.radius(); d += 1.0)
        {
            double d1 = d + d0;
            shouldEqualTolerance(view(d1,d0), spline(d)*spline(0.0), epsilon);
            shouldEqualTolerance(view(d0,d1), spline(d)*spline(0.0), epsilon);
            shouldEqualTolerance(view(d1,d1), spline(d)*spline(d), epsilon);
            shouldEqualTolerance(view(d1, d0, 1, 0), spline(d, 1)*spline(0.0), epsilon);
            shouldEqualTolerance(view(d0, d1, 0, 1), spline(d, 1)*spline(0.0), epsilon);
            shouldEqualTolerance(view(d1, d1, 1, 0), spline(d, 1)*spline(d), epsilon);
            shouldEqualTolerance(view(d1, d1, 0, 1), spline(d, 1)*spline(d), epsilon);
            shouldEqualTolerance(view(d1, d1, 1, 1), spline(d, 1)*spline(d, 1), epsilon);
        }
    }

    void testCoefficientArray()
    {
        double x = 5.3, y = 7.85;
        double dx = (N % 2) ? x - VIGRA_CSTD::floor(x) : x - VIGRA_CSTD::floor(x + 0.5),
               dy = (N % 2) ? y - VIGRA_CSTD::floor(y) : y - VIGRA_CSTD::floor(y + 0.5);

        Image coefficients(N+1, N+1);
        SplineImageView<N, double> view(srcImageRange(img));
        // use the coefficients from a different point in the same facet -- should be the same!
        view.coefficientArray(x-0.1, y-0.1, coefficients);

        double f_x_y = 0.0;
        for(int ny = 0; ny < N + 1; ++ny)
            for(int nx = 0; nx < N + 1; ++nx)
                f_x_y += VIGRA_CSTD::pow(dx, nx) * VIGRA_CSTD::pow(dy, ny) * coefficients(nx, ny);

        shouldEqualTolerance(f_x_y, view(x, y), 1e-12);

        // test out-of-bounds cases
        x = -x;
        y += img.height();
        dx = (N % 2) ? x - VIGRA_CSTD::floor(x) : x - VIGRA_CSTD::floor(x + 0.5);
        dy = (N % 2) ? y - VIGRA_CSTD::floor(y) : y - VIGRA_CSTD::floor(y + 0.5);
        view.coefficientArray(x-0.1, y-0.1, coefficients);

        f_x_y = 0.0;
        for(int ny = 0; ny < N + 1; ++ny)
            for(int nx = 0; nx < N + 1; ++nx)
                f_x_y += VIGRA_CSTD::pow(dx, nx) * VIGRA_CSTD::pow(dy, ny) * coefficients(nx, ny);

        shouldEqualTolerance(f_x_y, view(x, y), 1e-12);

        x = 5.3 + img.width();
        y = -7.85;
        dx = (N % 2) ? x - VIGRA_CSTD::floor(x) : x - VIGRA_CSTD::floor(x + 0.5);
        dy = (N % 2) ? y - VIGRA_CSTD::floor(y) : y - VIGRA_CSTD::floor(y + 0.5);
        view.coefficientArray(x-0.1, y-0.1, coefficients);

        f_x_y = 0.0;
        for(int ny = 0; ny < N + 1; ++ny)
            for(int nx = 0; nx < N + 1; ++nx)
                f_x_y += VIGRA_CSTD::pow(dx, nx) * VIGRA_CSTD::pow(dy, ny) * coefficients(nx, ny);

        shouldEqualTolerance(f_x_y, view(x, y), 1e-12);
    }

    void testImageResize()
    {
        char name[200];
        sprintf(name, "splineimageview%d.xv", N);
        ImageImportInfo info(name);

        Image reference(info.width(), info.height());
        importImage(info, destImage(reference));

        SplineImageView<N, double> view(srcImageRange(img));

        for(int y=0; y<reference.height(); ++y)
        {
            for(int x=0; x<reference.width(); ++x)
            {
                double dx = (double)x / (reference.width() - 1) * (img.width() - 1);
                double dy = (double)y / (reference.height() - 1) * (img.height() - 1);
                shouldEqualTolerance(view(dx, dy), reference(x, y), 1e-4);
            }
        }
    }

    void testImageResize0()
    {
        Image reference(img.width()*2-1, img.height()*2-1);
        resizeImageNoInterpolation(srcImageRange(img), destImageRange(reference));

        SplineImageView<0, double> view(srcImageRange(img));
        SplineImageView<0, int> viewi(srcImageRange(img));

        // check that data have not been copied
        shouldEqual(view.image().width(), 0);
        shouldEqual(viewi.image().width(), img.width());

        for(int y=0; y<reference.height(); ++y)
        {
            for(int x=0; x<reference.width(); ++x)
            {
                double dx = (double)x / 2.0;
                double dy = (double)y / 2.0;
                shouldEqualTolerance(view(dx, dy), reference(x, y), 1e-4);
                shouldEqualTolerance(view.unchecked(dx, dy, 0, 0), reference(x, y), 1e-4);
                FixedPoint<10, 1> fx(dx), fy(dy);
                shouldEqual(viewi.unchecked(fx, fy, 0, 0), viewi.unchecked(dx, dy, 0, 0));
            }
        }
    }

    void testImageResize1()
    {
        Image reference(img.width()*2-1, img.height()*2-1);
        resizeImageLinearInterpolation(srcImageRange(img), destImageRange(reference));

        SplineImageView<1, double> view(srcImageRange(img));
        SplineImageView<1, short> viewi(srcImageRange(img));

        // check that data have not been copied unnecessarily
        shouldEqual(view.image().width(), 0);
        shouldEqual(viewi.image().width(), img.width());

        for(int y=0; y<reference.height(); ++y)
        {
            for(int x=0; x<reference.width(); ++x)
            {
                double dx = (double)x / 2.0;
                double dy = (double)y / 2.0;
                shouldEqualTolerance(view(dx, dy), reference(x, y), 1e-4);
                shouldEqualTolerance(view.unchecked(dx, dy, 0, 0), reference(x, y), 1e-4);
                FixedPoint<10, 1> fx(dx), fy(dy);
                shouldEqual(viewi.unchecked(fx, fy, 0, 0), viewi.unchecked(dx, dy, 0, 0));
                shouldEqual(viewi.unchecked(fx, fy), viewi.unchecked(dx, dy));
            }
        }
    }

    void testOutside()
    {
        int center = 10;
        Image img(2*center+1, 2*center+1);
        img.init(0.0);
        img(center, center) = 1.0;
        SplineImageView<N, double> view(srcImageRange(img), true);
        BSplineBase<N> spline;

        double epsilon = 1.0e-10;
        double cx = center-0.42;
        double cy = center-0.23;
        shouldEqualTolerance(view(-cx,cy),      view(cx,cy), epsilon);
        shouldEqualTolerance(view(-cx,cy,1,0), -view(cx,cy,1,0), epsilon);
        shouldEqualTolerance(view(-cx,cy,0,1),  view(cx,cy,0,1), epsilon);
        shouldEqualTolerance(view(-cx,cy,1,1), -view(cx,cy,1,1), epsilon);

        shouldEqualTolerance(view(cx,-cy),      view(cx,cy), epsilon);
        shouldEqualTolerance(view(cx,-cy,1,0),  view(cx,cy,1,0), epsilon);
        shouldEqualTolerance(view(cx,-cy,0,1), -view(cx,cy,0,1), epsilon);
        shouldEqualTolerance(view(cx,-cy,1,1), -view(cx,cy,1,1), epsilon);

        shouldEqualTolerance( view(view.width() + cx,cy),
                              view(view.width() - 2.0 - cx,cy), epsilon);
        shouldEqualTolerance( view(view.width() + cx,cy,1,0),
                             -view(view.width() - 2.0 - cx,cy,1,0), epsilon);
        shouldEqualTolerance( view(view.width() + cx,cy,0,1),
                              view(view.width() - 2.0 - cx,cy,0,1), epsilon);
        shouldEqualTolerance( view(view.width() + cx,cy,1,1),
                             -view(view.width() - 2.0 - cx,cy,1,1), epsilon);

        shouldEqualTolerance( view(cx, view.height() + cy),
                              view(cx, view.height() - 2.0 - cy), epsilon);
        shouldEqualTolerance( view(cx, view.height() + cy,1,0),
                              view(cx, view.height() - 2.0 - cy,1,0), epsilon);
        shouldEqualTolerance( view(cx, view.height() + cy,0,1),
                             -view(cx, view.height() - 2.0 - cy,0,1), epsilon);
        shouldEqualTolerance( view(cx, view.height() + cy,1,1),
                             -view(cx, view.height() - 2.0 - cy,1,1), epsilon);

        try
        {
            view(2*view.width(), 0);
            failTest("Out-of-range coordinate failed to throw exception");
        }
        catch(vigra::PreconditionViolation) {}
    }

    void testVectorSIV()
    {
        // (compile-time only test for now)
        BRGBImage rgb(20, 10);
        SplineImageView<N, TinyVector<float, 3> > view(srcImageRange(rgb));
        (void)view(4.5, 1.3);
    }

};

struct GeometricTransformsTest
{
    typedef vigra::DImage Image;
    typedef vigra::MultiArrayView<2, double> View;
    Image img;
    int w, h;

    GeometricTransformsTest()
    {
        ImageImportInfo ginfo("lenna128.xv");
        w = ginfo.width();
        h = ginfo.height();
        img.resize(w, h);
        importImage(ginfo, destImage(img));
    }

    void testSimpleGeometry()
    {
        Image res1(img.size()), res11(img.size()), res2(h, w), res21(h, w);

        rotateImage(srcImageRange(img), destImage(res1), 180);
        rotateImage(View(img), View(res11), 180);
        for(int y = 0; y < 10; ++y)
            for(int x = 0; x < 10; ++x)
            {
                shouldEqual(img(x,y), res1(w-x-1, h-y-1));
                shouldEqual(img(x,y), res11(w-x-1, h-y-1));
            }
        
        rotateImage(srcImageRange(img), destImage(res2), 90);
        for(int y = 0; y < 10; ++y)
            for(int x = 0; x < 10; ++x)
                shouldEqual(img(x,y), res2(y, w-x-1));
        rotateImage(srcImageRange(img), destImage(res2), -90);
        for(int y = 0; y < 10; ++y)
            for(int x = 0; x < 10; ++x)
                shouldEqual(img(x,y), res2(h-y-1, x));


        try
        {
            rotateImage(srcImageRange(img), destImage(res2), 22);
            failTest("rotateImage() failed to throw exception");
        }
        catch(vigra::PreconditionViolation) {}

        transposeImage(srcImageRange(img), destImage(res2), vigra::major);
        transposeImage(View(img), View(res21), vigra::major);
        View transposed(View(img).transpose());
        for(int y = 0; y < 10; ++y)
            for(int x = 0; x < 10; ++x)
            {
                shouldEqual(img(x,y), res2(y, x));
                shouldEqual(img(x,y), res21(y, x));
                shouldEqual(transposed(y,x), res21(y,x));
            }

        transposeImage(srcImageRange(img), destImage(res2), vigra::minor);
        for(int y = 0; y < 10; ++y)
            for(int x = 0; x < 10; ++x)
                shouldEqual(img(x,y), res2(h-y-1, w-x-1));

        reflectImage(srcImageRange(img), destImage(res1), vigra::horizontal);
        reflectImage(View(img), View(res11), vigra::horizontal);
        for(int y = 0; y < 10; ++y)
            for(int x = 0; x < 10; ++x)
            {
                shouldEqual(img(x,y), res1(x, h-y-1));
                shouldEqual(img(x,y), res11(x, h-y-1));
            }

        reflectImage(srcImageRange(img), destImage(res1), vigra::vertical);
        for(int y = 0; y < 10; ++y)
            for(int x = 0; x < 10; ++x)
                shouldEqual(img(x,y), res1(w-x-1, y));
        reflectImage(srcImageRange(img), destImage(res1), vigra::horizontal | vigra::vertical);
        for(int y = 0; y < 10; ++y)
            for(int x = 0; x < 10; ++x)
                shouldEqual(img(x,y), res1(w-x-1, h-y-1));

        double xfactor = 3.0;
        double yfactor = 2.0;
        Image res3((int)(w*xfactor), (int)(h*yfactor)), res31(res3.size());
        resampleImage(srcImageRange(img), destImage(res3), xfactor, yfactor);
        resampleImage(View(img), View(res31), xfactor, yfactor);
        for(int y = 0; y < 10; ++y)
            for(int x = 0; x < 10; ++x)
            {
                shouldEqual(img(x,y), res3(int(x*xfactor), int(y*yfactor)));
                shouldEqual(img(x,y), res31(int(x*xfactor), int(y*yfactor)));
            }
    }

    void testAffineMatrix()
    {
        typedef TinyVector<double, 2> Vector2;

        Matrix<double> point(3,1);
        point(0,0) = 1.6;
        point(1,0) = 2.9;
        point(2,0) = 1.0;

        Matrix<double> t = translationMatrix2D(Vector2(2.2, 4.1));
        Matrix<double> res = t * point;
        shouldEqualTolerance(res(0,0), 3.8, 1e-14);
        shouldEqualTolerance(res(1,0), 7.0, 1e-14);
        shouldEqual(res(2,0), 1.0);

        Matrix<double> r = rotationMatrix2DDegrees(-90.0);
        res = r * point;
        shouldEqualTolerance(res(0,0), 2.9, 1e-14);
        shouldEqualTolerance(res(1,0), -1.6, 1e-14);
        shouldEqual(res(2,0), 1.0);

        r = rotationMatrix2DDegrees(-90.0, Vector2(1.6, 2.9));
        res = r * point;
        shouldEqualTolerance(res(0,0), 1.6, 1e-14);
        shouldEqualTolerance(res(1,0), 2.9, 1e-14);
        shouldEqual(res(2,0), 1.0);

        Matrix<double> s = scalingMatrix2D(2.0);
        res = s * point;
        shouldEqualTolerance(res(0,0), 3.2, 1e-14);
        shouldEqualTolerance(res(1,0), 5.8, 1e-14);
        shouldEqual(res(2,0), 1.0);

        Matrix<double> sh = shearMatrix2D(2.0, 0.5);
        res = sh * point;
        shouldEqualTolerance(res(0,0), 7.4, 1e-14);
        shouldEqualTolerance(res(1,0), 3.7, 1e-14);
        shouldEqual(res(2,0), 1.0);
    }

    void testRotation()
    {
        Image res(img.size()), ref(img.size()), res1(img.size());
        importImage(vigra::ImageImportInfo("lenna_rotate.xv"), destImage(ref));

        SplineImageView<3, double> sp(srcImageRange(img));

        rotateImage(sp, destImage(res), 45.0);
        shouldEqualSequenceTolerance(res.begin(), res.end(), ref.begin(), 1e-12);

        rotateImage(sp, View(res1), 45.0);
        shouldEqualSequenceTolerance(res1.begin(), res1.end(), ref.begin(), 1e-12);

        TinyVector<double, 2> center((w-1.0)/2.0, (h-1.0)/2.0);
        affineWarpImage(sp, destImageRange(res), rotationMatrix2DDegrees(45.0, center));
        shouldEqualSequenceTolerance(res.begin(), res.end(), ref.begin(), 1e-12);

        affineWarpImage(sp, View(res1), rotationMatrix2DDegrees(45.0, center));
        shouldEqualSequenceTolerance(res1.begin(), res1.end(), ref.begin(), 1e-12);
    }

    void testScaling()
    {
        Image res(2*w-1, 2*h-1), ref(2*w-1, 2*h-1);
        resizeImageSplineInterpolation(srcImageRange(img), destImageRange(ref));

        SplineImageView<3, double> sp(srcImageRange(img));

        affineWarpImage(sp, destImageRange(res), scalingMatrix2D(0.5));
        shouldEqualSequenceTolerance(res.begin(), res.end(), ref.begin(), 1e-14);
    }
};

struct ImageFunctionsTestSuite
: public vigra::test_suite
{
    ImageFunctionsTestSuite()
    : vigra::test_suite("ImageFunctionsTestSuite")
    {
        add( testCase( &ImageFunctionsTest::copyImageTest));
        add( testCase( &ImageFunctionsTest::copyImageIfTest));
        add( testCase( &ImageFunctionsTest::copyRedBandTest));
        add( testCase( &ImageFunctionsTest::copyGreenBandTest));
        add( testCase( &ImageFunctionsTest::copyBlueBandTest));
        add( testCase( &ImageFunctionsTest::rgbToGrayTest));
        add( testCase( &ImageFunctionsTest::initImageTest));
        add( testCase( &ImageFunctionsTest::initImageIfTest));
        add( testCase( &ImageFunctionsTest::initImageBorderTest));
        add( testCase( &ImageFunctionsTest::findMinMaxTest));
        add( testCase( &ImageFunctionsTest::findMinMaxIfTest));
        add( testCase( &ImageFunctionsTest::findAverageTest));
        add( testCase( &ImageFunctionsTest::findAverageIfTest));
        add( testCase( &ImageFunctionsTest::findAverageWeightedTest));
        add( testCase( &ImageFunctionsTest::findAverageAndVarianceTest));
        add( testCase( &ImageFunctionsTest::reduceFunctorTest));
        add( testCase( &ImageFunctionsTest::meshGridTest));
        add( testCase( &ImageFunctionsTest::findBoundingRectangleTest));
        add( testCase( &ImageFunctionsTest::lastValueFunctorTest));
        add( testCase( &ImageFunctionsTest::arrayOfRegionStatisticsTest));
        add( testCase( &ImageFunctionsTest::arrayOfRegionStatisticsIfTest));
        add( testCase( &ImageFunctionsTest::writeArrayOfRegionStatisticsTest));
        add( testCase( &ImageFunctionsTest::linearIntensityTransformTest));
        add( testCase( &ImageFunctionsTest::scalarIntensityTransformTest));
        add( testCase( &ImageFunctionsTest::linearIntensityTransformIfTest));
        add( testCase( &ImageFunctionsTest::linearRangeMappingTest));
        add( testCase( &ImageFunctionsTest::thresholdTest));
        add( testCase( &ImageFunctionsTest::brightnessContrastTest));
        add( testCase( &ImageFunctionsTest::gradientFunctionTest));
        add( testCase( &ImageFunctionsTest::additionTest));
        add( testCase( &ImageFunctionsTest::additionIfTest));
        add( testCase( &ImageFunctionsTest::resizeNoInterpolationTest));
        add( testCase( &ImageFunctionsTest::resizeLinearInterpolationTest));

        add( testCase( &ResizeImageTest::resizeLinearInterpolationReduceTest));
        add( testCase( &ResizeImageTest::scalarExpand));
        add( testCase( &ResizeImageTest::scalarReduce));
        add( testCase( &ResizeImageTest::rgbExpand));
        add( testCase( &ResizeImageTest::rgbReduce));
        add( testCase( &ResizeImageTest::testCatmullRomInterpolationExtensionWithLena));
        add( testCase( &ResizeImageTest::testCatmullRomInterpolationReductionWithLena));
        add( testCase( &ResizeImageTest::testCubicInterpolationExtensionWithLena));
        add( testCase( &ResizeImageTest::testCubicInterpolationReductionWithLena));
        add( testCase( &ResizeImageTest::testCatmullRomInterpolationExtensionHandControled));
        add( testCase( &SplineImageViewTest<0>::testPSF));
        add( testCase( &SplineImageViewTest<0>::testCoefficientArray));
        add( testCase( &SplineImageViewTest<0>::testImageResize0));
        add( testCase( &SplineImageViewTest<0>::testOutside));
        add( testCase( &SplineImageViewTest<1>::testPSF));
        add( testCase( &SplineImageViewTest<1>::testCoefficientArray));
        add( testCase( &SplineImageViewTest<1>::testImageResize1));
        add( testCase( &SplineImageViewTest<1>::testOutside));
        add( testCase( &SplineImageViewTest<2>::testPSF));
        add( testCase( &SplineImageViewTest<2>::testCoefficientArray));
        add( testCase( &SplineImageViewTest<2>::testImageResize));
        add( testCase( &SplineImageViewTest<2>::testOutside));
        add( testCase( &SplineImageViewTest<3>::testPSF));
        add( testCase( &SplineImageViewTest<3>::testCoefficientArray));
        add( testCase( &SplineImageViewTest<3>::testImageResize));
        add( testCase( &SplineImageViewTest<3>::testOutside));
        add( testCase( &SplineImageViewTest<5>::testPSF));
        add( testCase( &SplineImageViewTest<5>::testCoefficientArray));
        add( testCase( &SplineImageViewTest<5>::testImageResize));
        add( testCase( &SplineImageViewTest<5>::testOutside));
        add( testCase( &SplineImageViewTest<5>::testVectorSIV));

        add( testCase( &GeometricTransformsTest::testSimpleGeometry));
        add( testCase( &GeometricTransformsTest::testAffineMatrix));
        add( testCase( &GeometricTransformsTest::testRotation));
        add( testCase( &GeometricTransformsTest::testScaling));
    }
};

int main(int argc, char ** argv)
{
    ImageFunctionsTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

