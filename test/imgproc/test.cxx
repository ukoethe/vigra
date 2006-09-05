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
#include <functional>
#include <cmath>
#include <time.h>
#include "unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/stdimagefunctions.hxx"
#include "vigra/splineimageview.hxx"
#include "vigra/impex.hxx"

using namespace vigra;

struct ImageFunctionsTest
{
    typedef vigra::DImage Image;
    typedef vigra::DRGBImage RGBImage;
    typedef Image::value_type GrayValue;
    typedef RGBImage::value_type RGBValue;
    
    ImageFunctionsTest()
    : img(3,3), mask(3,3), rgb(3,3), col(1.0, 2.0, 3.0)
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

    }

    void copyImageIfTest()
    {
        Image img1(3,3);
        img1 = 42.0;

        copyImageIf(srcImageRange(img), maskImage(mask), destImage(img1));

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
    }

    void copyRedBandTest()
    {
        vigra::RedAccessor<RGBValue> red;

        copyImage(srcImageRange(rgb, red), destImage(img));

        Image::ScanOrderIterator i = img.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i)
        {
            should(acc(i) == 1.0);
        }

        img = 42.0;
        col.setRed(42.0);

        copyImage(srcImageRange(img), destImage(rgb, red));

        RGBImage::ScanOrderIterator ri = rgb.begin();
        RGBImage::Accessor racc = rgb.accessor();

        for(; ri != rgb.end(); ++ri)
        {
            should(racc(ri) == col);
        }
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

    }

    void findMinMaxTest()
    {
        vigra::FindMinMax<Image::value_type> minmax;

        inspectImage(srcImageRange(img), minmax);

        should(minmax.count == 9);
        should(minmax.max == 9.9);
        should(minmax.min == 1.1);
    }

    void findMinMaxIfTest()
    {
        vigra::FindMinMax<Image::value_type> minmax;

        inspectImageIf(srcImageRange(img), maskImage(mask), minmax);

        should(minmax.count == 7);
        should(minmax.max == 8.8);
        should(minmax.min == 2.2);
    }

    void findAverageTest()
    {
        vigra::FindAverage<Image::value_type> average;

        inspectImage(srcImageRange(img), average);

        should(average.count == 9);
        should(average() == 5.5);
    }

    void findAverageIfTest()
    {
        vigra::FindAverage<Image::value_type> average;

        inspectImageIf(srcImageRange(img), maskImage(mask), average);

        should(average.count == 7);
        should(average() == 5.5);
    }

    void reduceFunctorTest()
    {
        std::plus<double> p;
        ReduceFunctor<std::plus<double>, double> f(p, 0.0);

        inspectImage(srcImageRange(img), f);

        shouldEqual(f(), 49.5);
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

    void arrayOfRegionStatisticsTest()
    {
        vigra::BImage labels(3,3);
        labels = 1;
        labels.begin()[0] = 0;
        labels.begin()[8] = 0;

        vigra::ArrayOfRegionStatistics<vigra::FindMinMax<Image::value_type> > stats(1);

        inspectTwoImages(srcImageRange(img), srcImage(labels), stats);

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
        vigra::BImage labels(3,3);
        labels = 1;
        labels.begin()[4] = 2;
        labels.begin()[5] = 2;

        vigra::ArrayOfRegionStatistics<vigra::FindMinMax<Image::value_type> > stats(2);

        inspectTwoImagesIf(srcImageRange(img), srcImage(labels), maskImage(mask), stats);

        should(stats[0].count == 0);
        should(stats[1].count == 5);
        should(stats[1].min == 2.2);
        should(stats[1].max == 8.8);
        should(stats[2].count == 2);
        should(stats[2].min == 5.5);
        should(stats[2].max == 6.6);
    }

    void writeArrayOfRegionStatisticsTest()
    {
        vigra::BImage labels(3,3);
        labels = 1;
        labels.begin()[0] = 0;
        labels.begin()[8] = 0;

        vigra::ArrayOfRegionStatistics<vigra::FindAverage<Image::value_type> > stats(1);

        inspectTwoImages(srcImageRange(img), srcImage(labels), stats);

        should(stats[0].count == 2);
        should(stats[0]() == 5.5);
        should(stats[1].count == 7);
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

        transformImage(srcImageRange(img), destImage(img1), trans);

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

        transformImageIf(srcImageRange(img), maskImage(mask), destImage(img1), trans);

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
					   linearRangeMapping(1.0, 10.0, 0, 255));

#if defined(__i386__)
        int res[] = {3, 34, 65, 96, 127, 159, 190, 221, 252 };
#else
        int res[] = {3, 34, 65, 96, 128, 159, 190, 221, 252 };
#endif
        shouldEqualSequence(img1.begin(), img1.end(), res);

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
          destImage(res), vigra::MagnitudeFunctor<float>());

        for(y=0; y<3; ++y)
        {
            for(int x=0; x<3; ++x)
            {
                should(VIGRA_CSTD::fabs(res(x,y) - VIGRA_CSTD::sqrt(2.0)) < 1e-6);
            }
        }

        gradientBasedTransform(
          srcImageRange(in),
          destImage(res), vigra::RGBGradientMagnitudeFunctor<float>());

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
        Image img1(3,3);

        std::plus<Image::value_type> add;

        combineTwoImages(srcImageRange(img), srcImage(img), destImage(img1), add);

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i, ++i1)
        {
            should(2.0*acc(i) == acc(i1));
        }

    }

    void additionIfTest()
    {
        Image img1(3,3);
        img1 = 42.0;

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

    }

    void resizeNoInterpolationTest()
    {
        Image img1(5,5);

        resizeImageNoInterpolation(srcImageRange(img), destImageRange(img1));

        Image::ScanOrderIterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        should(acc(i1) == 1.1);
        ++i1;
        should(acc(i1) == 1.1);
        ++i1;
        should(acc(i1) == 2.2);
        ++i1;
        should(acc(i1) == 2.2);
        ++i1;
        should(acc(i1) == 3.3);
        ++i1;
        should(acc(i1) == 1.1);
        ++i1;
        should(acc(i1) == 1.1);
        ++i1;
        should(acc(i1) == 2.2);
        ++i1;
        should(acc(i1) == 2.2);
        ++i1;
        should(acc(i1) == 3.3);
        ++i1;
        should(acc(i1) == 4.4);
        ++i1;
        should(acc(i1) == 4.4);
        ++i1;
        should(acc(i1) == 5.5);
        ++i1;
        should(acc(i1) == 5.5);
        ++i1;
        should(acc(i1) == 6.6);
        ++i1;
        should(acc(i1) == 4.4);
        ++i1;
        should(acc(i1) == 4.4);
        ++i1;
        should(acc(i1) == 5.5);
        ++i1;
        should(acc(i1) == 5.5);
        ++i1;
        should(acc(i1) == 6.6);
        ++i1;
        should(acc(i1) == 7.7);
        ++i1;
        should(acc(i1) == 7.7);
        ++i1;
        should(acc(i1) == 8.8);
        ++i1;
        should(acc(i1) == 8.8);
        ++i1;
        should(acc(i1) == 9.9);

        Image img2(3,3);
        img2 = 0;

        resizeImageNoInterpolation(srcImageRange(img1), destImageRange(img2));

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i2 = img2.begin();

        for(; i != img.end(); ++i, ++i2)
        {
            should(acc(i) == acc(i2));
        }
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

        resizeImageNoInterpolation(srcImageRange(img1), destImageRange(img2));

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i2 = img2.begin();

        for(; i != img.end(); ++i, ++i2)
        {
            should(acc(i) == acc(i2));
        }
    }


    Image img;
    vigra::BImage mask;
    RGBImage rgb;
    RGBValue col;
};

struct ResizeImageTest
{
    typedef vigra::FImage Image;
    typedef vigra::FRGBImage RGBImage;

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

        resizeImageLinearInterpolation(srcImageRange(rgb), destImageRange(rgbdest));

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

        resizeImageSplineInterpolation(srcImageRange(rgb), destImageRange(img1));

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

        resizeImageSplineInterpolation(srcImageRange(rgb), destImageRange(img1));

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

        resizeImageCatmullRomInterpolation(srcImageRange(rgb), destImageRange(rgbdest));

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

        Image coefficients;
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

        for(int y=0; y<reference.height()-1; ++y)
        {
            for(int x=0; x<reference.width()-1; ++x)
            {
                double dx = (double)x / 2.0;
                double dy = (double)y / 2.0;
                shouldEqualTolerance(view(dx, dy), reference(x+1, y+1), 1e-4);
                shouldEqualTolerance(view.unchecked(dx, dy, 0, 0), reference(x+1, y+1), 1e-4);
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
        add( testCase( &ImageFunctionsTest::reduceFunctorTest));
        add( testCase( &ImageFunctionsTest::findBoundingRectangleTest));
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
    }
};

int main()
{
    ImageFunctionsTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

