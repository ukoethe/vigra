#include <iostream>
#include <functional>
#include "unittest.h"
#include "vigra/stdimage.hxx"
#include "vigra/stdimagefunctions.hxx"
#include "vigra/viff.hxx"

struct ImageFunctionsTest
{
    typedef DImage Image;
    typedef FRGBImage RGBImage;

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

        mask = 1;
        mask.begin()[0] = 0;
        mask.begin()[8] = 0;

        rgb = col;
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
        RedAccessor<RGBImage::value_type> red;

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
        GreenAccessor<RGBImage::value_type> green;

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
        BlueAccessor<RGBImage::value_type> blue;

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
        RGBToGrayAccessor<RGBImage::value_type> gray;

        copyImage(srcImageRange(rgb, gray), destImage(img));

        Image::ScanOrderIterator i = img.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i)
        {
            should(acc(i) == col.luminance());
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
            if(k != 4) should(acc(i) == 42.0)
            else       should(acc(i) == 5.5);
        }

    }

    void findMinMaxTest()
    {
        FindMinMax<Image::value_type> minmax;

        inspectImage(srcImageRange(img), minmax);

        should(minmax.count == 9);
        should(minmax.max == 9.9);
        should(minmax.min == 1.1);
    }

    void findMinMaxIfTest()
    {
        FindMinMax<Image::value_type> minmax;

        inspectImageIf(srcImageRange(img), maskImage(mask), minmax);

        should(minmax.count == 7);
        should(minmax.max == 8.8);
        should(minmax.min == 2.2);
    }

    void findAverageTest()
    {
        FindAverage<Image::value_type> average;

        inspectImage(srcImageRange(img), average);

        should(average.count == 9);
        should(average() == 5.5);
    }

    void findAverageIfTest()
    {
        FindAverage<Image::value_type> average;

        inspectImageIf(srcImageRange(img), maskImage(mask), average);

        should(average.count == 7);
        should(average() == 5.5);
    }

    void arrayOfRegionStatisticsTest()
    {
        BImage labels(3,3);
        labels = 1;
        labels.begin()[0] = 0;
        labels.begin()[8] = 0;

        ArrayOfRegionStatistics<FindMinMax<Image::value_type> > stats(1);

        inspectTwoImages(srcImageRange(img), srcImage(labels), stats);

        should(stats[0].count == 2);
        should(stats[0].min == 1.1);
        should(stats[0].max == 9.9);
        should(stats[1].count == 7);
        should(stats[1].min == 2.2);
        should(stats[1].max == 8.8);
    }

    void arrayOfRegionStatisticsIfTest()
    {
        BImage labels(3,3);
        labels = 1;
        labels.begin()[4] = 2;
        labels.begin()[5] = 2;

        ArrayOfRegionStatistics<FindMinMax<Image::value_type> > stats(2);

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
        BImage labels(3,3);
        labels = 1;
        labels.begin()[0] = 0;
        labels.begin()[8] = 0;

        ArrayOfRegionStatistics<FindAverage<Image::value_type> > stats(1);

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
        LinearIntensityTransform<Image::value_type> trans(2.0, -1.1);

        Image img1(3,3);

        transformImage(srcImageRange(img), destImage(img1), trans);

        Image::ScanOrderIterator i = img.begin();
        Image::ScanOrderIterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        for(; i != img.end(); ++i, ++i1)
        {
            should(2.0*(acc(i) - 1.1) == acc(i1));
        }

    }

    void linearIntensityTransformIfTest()
    {
        LinearIntensityTransform<Image::value_type> trans(2.0, -1.1);

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
            should(2.0*(acc(i) - 1.1) == acc(i1));
        }

    }

    void thresholdTest()
    {
        Threshold<Image::value_type, Image::value_type> trans(2.0, 9.0, 0.0, 1.0);

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
        double epsilon = 0.00000001;

        resizeImageLinearInterpolation(srcImageRange(img), destImageRange(img1));

        Image::ScanOrderIterator i1 = img1.begin();
        Image::Accessor acc = img.accessor();

        should(abs(acc(i1) - 1.1) < epsilon);
        ++i1;
        should(abs(acc(i1) - 1.65) < epsilon);
        ++i1;
        should(abs(acc(i1) - 2.2) < epsilon);
        ++i1;
        should(abs(acc(i1) - 2.75) < epsilon);
        ++i1;
        should(abs(acc(i1) - 3.3) < epsilon);
        ++i1;

        should(abs(acc(i1) - 2.75) < epsilon);
        ++i1;
        should(abs(acc(i1) - 3.3) < epsilon);
        ++i1;
        should(abs(acc(i1) - 3.85) < epsilon);
        ++i1;
        should(abs(acc(i1) - 4.4) < epsilon);
        ++i1;
        should(abs(acc(i1) - 4.95) < epsilon);
        ++i1;

        should(abs(acc(i1) - 4.4) < epsilon);
        ++i1;
        should(abs(acc(i1) - 4.95) < epsilon);
        ++i1;
        should(abs(acc(i1) - 5.5) < epsilon);
        ++i1;
        should(abs(acc(i1) - 6.05) < epsilon);
        ++i1;
        should(abs(acc(i1) - 6.6) < epsilon);
        ++i1;

        should(abs(acc(i1) - 6.05) < epsilon);
        ++i1;
        should(abs(acc(i1) - 6.6) < epsilon);
        ++i1;
        should(abs(acc(i1) - 7.15) < epsilon);
        ++i1;
        should(abs(acc(i1) - 7.7) < epsilon);
        ++i1;
        should(abs(acc(i1) - 8.25) < epsilon);
        ++i1;

        should(abs(acc(i1) - 7.7) < epsilon);
        ++i1;
        should(abs(acc(i1) - 8.25) < epsilon);
        ++i1;
        should(abs(acc(i1) - 8.8) < epsilon);
        ++i1;
        should(abs(acc(i1) - 9.35) < epsilon);
        ++i1;
        should(abs(acc(i1) - 9.9) < epsilon);

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
    BImage mask;
    RGBImage rgb;
    RGBImage::value_type col;
};

struct ResizeImageSplineTest
{
    typedef FImage Image;
    typedef FRGBImage RGBImage;

    ResizeImageSplineTest()
    {
        ViffImage * viff = readViffImage("lenna128.xv");

        shouldMsg(viff != 0, "Unable to read test image \"lenna128.xv\"\n");

        img.resize(viff->row_size,viff->col_size);
        importViffImage(viff, destImage(img));
        freeViffImage(viff);

        viff = readViffImage("lenna128rgb.xv");

        shouldMsg(viff != 0, "Unable to read test image \"lenna128rgb.xv\"\n");

        rgb.resize(viff->row_size,viff->col_size);
        importViffImage(viff, destImage(rgb));
        freeViffImage(viff);

    }

    void scalarExpand()
    {
        ViffImage * viff = readViffImage("lenna288.xv");

        shouldMsg(viff != 0, "Unable to read test image \"lenna288.xv\"\n");

        Image imgex(viff->row_size,viff->col_size);
        importViffImage(viff, destImage(imgex));
        freeViffImage(viff);

        Image img1(imgex.width(), imgex.height());
        double epsilon = 0.000001;

        resizeImageSplineInterpolation(srcImageRange(img), destImageRange(img1));

        Image::ScanOrderIterator i1 = img1.begin();
        Image::ScanOrderIterator iex = imgex.begin();
        Image::Accessor acc = img1.accessor();

        for(; i1 != img1.end(); ++i1, ++iex)
        {
            should(abs(acc(i1) - acc(iex)) < epsilon);
        }
    }

    void scalarReduce()
    {
        ViffImage * viff = readViffImage("lenna42.xv");

        shouldMsg(viff != 0, "Unable to read test image \"lenna42.xv\"\n");

        Image imgred(viff->row_size,viff->col_size);
        importViffImage(viff, destImage(imgred));
        freeViffImage(viff);
        
        Image img1(imgred.width(), imgred.height());
        double epsilon = 0.000001;

        resizeImageSplineInterpolation(srcImageRange(img), destImageRange(img1));

        Image::ScanOrderIterator i1 = img1.begin();
        Image::ScanOrderIterator ired = imgred.begin();
        Image::Accessor acc = img1.accessor(); 

        for(; i1 != img1.end(); ++i1, ++ired)
        {
            should(abs(acc(i1) - acc(ired)) < epsilon);
        }
    }

    void rgbExpand()
    {
        ViffImage * viff = readViffImage("lenna288rgb.xv");

        shouldMsg(viff != 0, "Unable to read test image \"lenna288rgb.xv\"\n");

        RGBImage rgbex(viff->row_size,viff->col_size);
        importViffImage(viff, destImage(rgbex));
        freeViffImage(viff);

        RGBImage img1(rgbex.width(), rgbex.height());
        double epsilon = 0.0001;

        resizeImageSplineInterpolation(srcImageRange(rgb), destImageRange(img1));

        RGBImage::ScanOrderIterator i1 = img1.begin();
        RGBImage::ScanOrderIterator iex = rgbex.begin();
        RGBImage::Accessor acc = img1.accessor();

        for(; i1 != img1.end(); ++i1, ++iex)
        {   
            should(abs(acc.red(i1) - acc.red(iex)) < epsilon);
            should(abs(acc.green(i1) - acc.green(iex)) < epsilon);
            should(abs(acc.blue(i1) - acc.blue(iex)) < epsilon);
        }
    }

    void rgbReduce()
    {
        ViffImage * viff = readViffImage("lenna42rgb.xv");

        shouldMsg(viff != 0, "Unable to read test image \"lenna42rgb.xv\"\n");

        RGBImage rgbred(viff->row_size,viff->col_size);
        importViffImage(viff, destImage(rgbred));
        freeViffImage(viff);

        RGBImage img1(rgbred.width(), rgbred.height());
        double epsilon = 0.0001;

        resizeImageSplineInterpolation(srcImageRange(rgb), destImageRange(img1));

        RGBImage::ScanOrderIterator i1 = img1.begin();
        RGBImage::ScanOrderIterator ired = rgbred.begin();
        RGBImage::Accessor acc = img1.accessor();

        for(; i1 != img1.end(); ++i1, ++ired)
        {
            should(abs(acc.red(i1) - acc.red(ired)) < epsilon);
            should(abs(acc.green(i1) - acc.green(ired)) < epsilon);
            should(abs(acc.blue(i1) - acc.blue(ired)) < epsilon);
        }
    }

    Image img;
    RGBImage rgb;
};

        
        
struct ImageFunctionsTestSuite
: public TestSuite
{
    ImageFunctionsTestSuite()
    : TestSuite("ImageFunctionsTestSuite")
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
        add( testCase( &ImageFunctionsTest::arrayOfRegionStatisticsTest));
        add( testCase( &ImageFunctionsTest::arrayOfRegionStatisticsIfTest));
        add( testCase( &ImageFunctionsTest::writeArrayOfRegionStatisticsTest));
        add( testCase( &ImageFunctionsTest::linearIntensityTransformTest));
        add( testCase( &ImageFunctionsTest::linearIntensityTransformIfTest));
        add( testCase( &ImageFunctionsTest::thresholdTest));
        add( testCase( &ImageFunctionsTest::additionTest));
        add( testCase( &ImageFunctionsTest::additionIfTest));
        add( testCase( &ImageFunctionsTest::resizeNoInterpolationTest));
        add( testCase( &ImageFunctionsTest::resizeLinearInterpolationTest));
        add( testCase( &ResizeImageSplineTest::scalarExpand));
        add( testCase( &ResizeImageSplineTest::scalarReduce));
        add( testCase( &ResizeImageSplineTest::rgbExpand));
        add( testCase( &ResizeImageSplineTest::rgbReduce));
    }
};

int main()
{
    ImageFunctionsTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

