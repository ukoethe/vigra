#include <iostream>
#include <cmath>
#include "unittest.h"
#include "vigra/stdimage.hxx"
#include "vigra/convolution.hxx"
#include "vigra/impex.hxx"

struct ConvolutionTest
{
    typedef vigra::DImage Image;

    ConvolutionTest()
    : constimg(5,5), rampimg(5,1)
    {
        constimg = 1.0;
        
        ViffImage * viff = readViffImage("lenna128.xv");

        shouldMsg(viff != 0, "Unable to read test image \"lenna128.xv\"\n");

        lenna.resize(viff->row_size,viff->col_size);
        importViffImage(viff, destImage(lenna));
        freeViffImage(viff);

        Image::ScanOrderIterator i = rampimg.begin();
        Image::ScanOrderIterator end = rampimg.end();
        Image::Accessor acc = rampimg.accessor();
        
        for(int k=0; i != end; ++i, ++k)
        {
            acc.set(k, i);
        }
    }
    
    void separableConvolutionTest()
    {
        vigra::Kernel1D<double> binom;
        binom.initBinomial(1);
        
        vigra::Kernel1D<double>::Iterator center = binom.center();
        
        should(center[0] == 0.5);
        should(center[-1] == 0.25);
        should(center[1] == 0.25);
        
        Image tmp1(constimg);
        Image tmp2(constimg);
        tmp2 = 0.0;
        
        separableConvolveX(srcImageRange(constimg), destImage(tmp1), kernel1d(binom));
        separableConvolveY(srcImageRange(tmp1), destImage(tmp2), kernel1d(binom));
        
        Image::ScanOrderIterator i1 = constimg.begin();
        Image::ScanOrderIterator i1end = constimg.end();
        Image::ScanOrderIterator i2 = tmp2.begin();
        Image::Accessor acc = constimg.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(acc(i1) == acc(i2));
        }
    }
    
    void separableDerivativeRepeatTest()
    {
        vigra::Kernel1D<double> grad;
        grad.initSymmetricGradient();
        
        Image tmp1(rampimg);
        tmp1 = 0.0;
        Image tmp2(constimg);
        
        separableConvolveX(srcImageRange(rampimg), destImage(tmp1), kernel1d(grad));
        separableConvolveX(srcImageRange(constimg), destImage(tmp2), kernel1d(grad));
        
        Image::ScanOrderIterator i1 = tmp1.begin();
        Image::ScanOrderIterator i2 = tmp2.begin();
        Image::Accessor acc = tmp1.accessor();
        
        should(acc(i1) == 0.5);
        should(acc(i2) == 0.0);
        ++i1;
        ++i2;
        should(acc(i1) == 1.0);
        should(acc(i2) == 0.0);
        ++i1;
        ++i2;
        should(acc(i1) == 1.0);
        should(acc(i2) == 0.0);
        ++i1;
        ++i2;
        should(acc(i1) == 1.0);
        should(acc(i2) == 0.0);
        ++i1;
        ++i2;
        should(acc(i1) == 0.5);
        should(acc(i2) == 0.0);
    }
    
    void separableDerivativeReflectTest()
    {
        vigra::Kernel1D<double> grad;
        grad.initSymmetricGradient();
        grad.setBorderTreatment(vigra::BORDER_TREATMENT_REFLECT);
        
        Image tmp1(rampimg);
        tmp1 = 1000.0;
        
        separableConvolveX(srcImageRange(rampimg), destImage(tmp1), kernel1d(grad));
        
        Image::ScanOrderIterator i1 = tmp1.begin();
        Image::Accessor acc = tmp1.accessor();
        
        should(acc(i1) == 0.0);
        ++i1;
        should(acc(i1) == 1.0);
        ++i1;
        should(acc(i1) == 1.0);
        ++i1;
        should(acc(i1) == 1.0);
        ++i1;
        should(acc(i1) == 0.0);
    }
    
    void separableDerivativeAvoidTest()
    {
        vigra::Kernel1D<double> grad;
        grad.initSymmetricGradient();
        grad.setBorderTreatment(vigra::BORDER_TREATMENT_AVOID);
        
        Image tmp1(rampimg);
        tmp1 = 1000.0;
        
        separableConvolveX(srcImageRange(rampimg), destImage(tmp1), kernel1d(grad));
        
        Image::ScanOrderIterator i1 = tmp1.begin();
        Image::Accessor acc = tmp1.accessor();
        
        should(acc(i1) == 1000.0);
        ++i1;
        should(acc(i1) == 1.0);
        ++i1;
        should(acc(i1) == 1.0);
        ++i1;
        should(acc(i1) == 1.0);
        ++i1;
        should(acc(i1) == 1000.0);
    }
    
    void separableSmoothClipTest()
    {
        vigra::Kernel1D<double> binom;
        binom.initBinomial(1);
        binom.setBorderTreatment(vigra::BORDER_TREATMENT_CLIP);
        
        Image tmp1(rampimg);
        tmp1 = 1000.0;
        
        separableConvolveX(srcImageRange(rampimg), destImage(tmp1), kernel1d(binom));
        
        Image::ScanOrderIterator i1 = tmp1.begin();
        Image::Accessor acc = tmp1.accessor();
        
        should(acc(i1) == 1.0/3.0);
        ++i1;
        should(acc(i1) == 1.0);
        ++i1;
        should(acc(i1) == 2.0);
        ++i1;
        should(acc(i1) == 3.0);
        ++i1;
        should(acc(i1) == 11.0/3.0);
    }
    
    void separableSmoothWrapTest()
    {
        vigra::Kernel1D<double> binom;
        binom.initBinomial(1);
        binom.setBorderTreatment(vigra::BORDER_TREATMENT_WRAP);
        
        Image tmp1(rampimg);
        tmp1 = 1000.0;
        
        separableConvolveX(srcImageRange(rampimg), destImage(tmp1), kernel1d(binom));
        
        Image::ScanOrderIterator i1 = tmp1.begin();
        Image::Accessor acc = tmp1.accessor();
        
        should(acc(i1) == 1.25);
        ++i1;
        should(acc(i1) == 1.0);
        ++i1;
        should(acc(i1) == 2.0);
        ++i1;
        should(acc(i1) == 3.0);
        ++i1;
        should(acc(i1) == 2.75);
    }
    
    void gaussianSmoothingTest()
    {
        vigra::Kernel1D<double> gauss;
        gauss.initGaussian(1.0);

        Image tmp1(lenna.size());
        Image tmp2(lenna.size());

        separableConvolveX(srcImageRange(lenna), destImage(tmp1), kernel1d(gauss));
        separableConvolveY(srcImageRange(tmp1), destImage(tmp2), kernel1d(gauss));
        
        gaussianSmoothing(srcImageRange(lenna), destImage(tmp1), 1.0);
        
        Image::ScanOrderIterator i1 = tmp1.begin();
        Image::ScanOrderIterator i1end = tmp1.end();
        Image::ScanOrderIterator i2 = tmp2.begin();
        Image::Accessor acc = constimg.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(acc(i1) == acc(i2));
        }
    }
    
    void separableGradientTest()
    {
        Image sepgrad(lenna.size());
        importImage(vigra::ImageImportInfo("lenna128sepgrad.xv"), destImage(sepgrad));

        double epsilon = 0.00001;
        
        vigra::Kernel1D<double> gauss;
        gauss.initGaussian(1.0);
        vigra::Kernel1D<double> grad;
        grad.initGaussianDerivative(1.0, 1);

        Image tmp1(lenna);
        tmp1 = 0.0;
        Image tmp2(lenna);
        tmp2 = 0.0;
        Image tmp3(lenna);

        separableConvolveX(srcImageRange(lenna), destImage(tmp3), kernel1d(grad));
        separableConvolveY(srcImageRange(tmp3), destImage(tmp1), kernel1d(gauss));
        
        separableConvolveX(srcImageRange(lenna), destImage(tmp3), kernel1d(gauss));
        separableConvolveY(srcImageRange(tmp3), destImage(tmp2), kernel1d(grad));
        
        Image::ScanOrderIterator i1 = tmp1.begin();
        Image::ScanOrderIterator i1end = tmp1.end();
        Image::ScanOrderIterator i2 = tmp2.begin();
        Image::ScanOrderIterator i = sepgrad.begin();
        Image::Accessor acc = constimg.accessor();
        
        for(; i1 != i1end; ++i1, ++i2, ++i)
        {
            double grad = sqrt(acc(i1)*acc(i1)+acc(i2)*acc(i2));
            should(std::abs(grad - acc(i)) < epsilon);
        }
    }
    
    void gradientTest()
    {
        Image sepgrad(lenna.size());
        importImage(vigra::ImageImportInfo("lenna128sepgrad.xv"), destImage(sepgrad));

        double epsilon = 0.00001;
        
        Image tmpx(lenna.size());
        Image tmpy(lenna.size());

        gaussianGradient(srcImageRange(lenna), destImage(tmpx), destImage(tmpy), 1.0);

        Image::ScanOrderIterator i1 = tmpx.begin();
        Image::ScanOrderIterator i1end = tmpx.end();
        Image::ScanOrderIterator i2 = tmpy.begin();
        Image::ScanOrderIterator i = sepgrad.begin();
        Image::Accessor acc = constimg.accessor();
        
        for(; i1 != i1end; ++i1, ++i2, ++i)
        {
            double grad = sqrt(acc(i1)*acc(i1)+acc(i2)*acc(i2));
            should(std::abs(grad - acc(i)) < epsilon);
        }
    }
    
    void hessianTest()
    {
        Image resxx(lenna.size());
        Image resxy(lenna.size());
        Image resyy(lenna.size());
        Image refxx(lenna.size());
        Image refxy(lenna.size());
        Image refyy(lenna.size());

        hessianMatrixOfGaussian(srcImageRange(lenna), 
            destImage(resxx), destImage(resxy), destImage(resyy), 1.0);
            
        importImage(vigra::ImageImportInfo("lennahessxx.xv"), destImage(refxx));
        importImage(vigra::ImageImportInfo("lennahessyy.xv"), destImage(refyy));
        importImage(vigra::ImageImportInfo("lennahessxy.xv"), destImage(refxy));
        
        Image::ScanOrderIterator i1 = resxx.begin();
        Image::ScanOrderIterator i1end = resxx.end();
        Image::ScanOrderIterator i2 = resyy.begin();
        Image::ScanOrderIterator i3 = resxy.begin();
        Image::ScanOrderIterator r1 = refxx.begin();
        Image::ScanOrderIterator r2 = refyy.begin();
        Image::ScanOrderIterator r3 = refxy.begin();
        Image::Accessor acc = constimg.accessor();
        
        for(; i1 != i1end; ++i1, ++i2, ++i3, ++r1, ++r2, ++r3)
        {
            should(acc(i1) == acc(r1));
            should(acc(i2) == acc(r2));
            should(acc(i3) == acc(r3));
        }
    }
    
    void structureTensorTest()
    {
        Image resxx(lenna.size());
        Image resxy(lenna.size());
        Image resyy(lenna.size());
        Image refxx(lenna.size());
        Image refxy(lenna.size());
        Image refyy(lenna.size());

        structureTensor(srcImageRange(lenna), 
            destImage(resxx), destImage(resxy), destImage(resyy), 1.0, 2.0);
            
        importImage(vigra::ImageImportInfo("lennastxx.xv"), destImage(refxx));
        importImage(vigra::ImageImportInfo("lennastyy.xv"), destImage(refyy));
        importImage(vigra::ImageImportInfo("lennastxy.xv"), destImage(refxy));
        
        Image::ScanOrderIterator i1 = resxx.begin();
        Image::ScanOrderIterator i1end = resxx.end();
        Image::ScanOrderIterator i2 = resyy.begin();
        Image::ScanOrderIterator i3 = resxy.begin();
        Image::ScanOrderIterator r1 = refxx.begin();
        Image::ScanOrderIterator r2 = refyy.begin();
        Image::ScanOrderIterator r3 = refxy.begin();
        Image::Accessor acc = constimg.accessor();
        
        for(; i1 != i1end; ++i1, ++i2, ++i3, ++r1, ++r2, ++r3)
        {
            should(acc(i1) == acc(r1));
            should(acc(i2) == acc(r2));
            should(acc(i3) == acc(r3));
        }
    }
    
    void stdConvolutionTest()
    {
        vigra::Kernel1D<double> binom1;
        binom1.initBinomial(1);
        
        vigra::Kernel2D<double> binom2;
        binom2.initSeparable(binom1, binom1);
        
        Image tmp1(constimg);
        tmp1 = 0.0;

        convolveImage(srcImageRange(constimg), destImage(tmp1), kernel2d(binom2));
        
        Image::ScanOrderIterator i1 = constimg.begin();
        Image::ScanOrderIterator i1end = constimg.end();
        Image::ScanOrderIterator i2 = tmp1.begin();
        Image::Accessor acc = constimg.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(acc(i1) == acc(i2));
        }
    }
    
    void stdVersusSeparableConvolutionTest()
    {
        double epsilon = 0.00001;
        
        vigra::Kernel1D<double> gauss1;
        gauss1.initGaussian(2.0);
        
        vigra::Kernel2D<double> gauss2;
        gauss2.initSeparable(gauss1, gauss1);
        
        Image tmp1(lenna);
        tmp1 = 0.0;

        convolveImage(srcImageRange(lenna), destImage(tmp1), kernel2d(gauss2));
        
        Image tmp2(lenna);
        Image tmp3(lenna);
        tmp3 = 0.0;
        
        separableConvolveX(srcImageRange(lenna), destImage(tmp2), kernel1d(gauss1));
        separableConvolveY(srcImageRange(tmp2), destImage(tmp3), kernel1d(gauss1));
        
        Image::Iterator y1 = tmp1.upperLeft() - gauss2.upperLeft();
        Image::Iterator end = tmp1.lowerRight() - gauss2.lowerRight();
        Image::Iterator y2 = tmp3.upperLeft() - gauss2.upperLeft();
        Image::Accessor acc = tmp1.accessor();
        
        for(; y1.y != end.y; ++y1.y, ++y2.y)
        {
            Image::Iterator x1 = y1;
            Image::Iterator x2 = y2;
            for(; x1.x != end.x; ++x1.x, ++x2.x)
            {
                should(std::abs(acc(x1) - acc(x2)) < epsilon);
            }
        }
    }
    
    void recursiveSmoothTest()
    {
        Image tmp1(constimg);
        tmp1 = 0.0;

        recursiveSmoothX(srcImageRange(constimg), destImage(tmp1), 1.0);
        recursiveSmoothY(srcImageRange(tmp1), destImage(tmp1), 1.0);
        
        Image::ScanOrderIterator i1 = constimg.begin();
        Image::ScanOrderIterator i1end = constimg.end();
        Image::ScanOrderIterator i2 = tmp1.begin();
        Image::Accessor acc = constimg.accessor();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(acc(i1) == acc(i2));
        }
    }
    
    void recursiveGradientTest()
    {
        ViffImage * viff = readViffImage("lenna128recgrad.xv");

        shouldMsg(viff != 0, "Unable to read test image \"lenna128recgrad.xv\"\n");

        Image recgrad(viff->row_size,viff->col_size);
        importViffImage(viff, destImage(recgrad));
        freeViffImage(viff);

        double epsilon = 0.00001;
        
        Image tmp1(lenna);
        tmp1 = 0.0;
        Image tmp2(lenna);
        tmp2 = 0.0;

        recursiveFirstDerivativeX(srcImageRange(lenna), destImage(tmp1), 1.0);
        recursiveSmoothY(srcImageRange(tmp1), destImage(tmp1), 1.0);
        
        recursiveSmoothX(srcImageRange(lenna), destImage(tmp2), 1.0);
        recursiveFirstDerivativeY(srcImageRange(tmp2), destImage(tmp2), 1.0);
        
        Image::ScanOrderIterator i1 = tmp1.begin();
        Image::ScanOrderIterator i1end = tmp1.end();
        Image::ScanOrderIterator i2 = tmp2.begin();
        Image::ScanOrderIterator i = recgrad.begin();
        Image::Accessor acc = constimg.accessor();
        
        for(; i1 != i1end; ++i1, ++i2, ++i)
        {
            double grad = sqrt(acc(i1)*acc(i1)+acc(i2)*acc(i2));
            
            should(std::abs(grad - acc(i)) < epsilon);
        }
    }
    
    void recursiveSecondDerivativeTest()
    {
        double epsilon = 0.0001;
        double b = exp(-1.0);
        double factor = (1.0 - b) * (1.0 - b) / b;
        
        Image tmp1(rampimg);
        tmp1 = 0.0;
        Image tmp2(rampimg);
        tmp2 = 0.0;

        recursiveSmoothX(srcImageRange(rampimg), destImage(tmp1), 1.0);
        
        recursiveSecondDerivativeX(srcImageRange(rampimg), destImage(tmp2), 1.0);
        
        Image::ScanOrderIterator i1 = rampimg.begin();
        Image::ScanOrderIterator i1end = i1 + rampimg.width();
        Image::ScanOrderIterator i2 = tmp1.begin();
        Image::ScanOrderIterator i3 = tmp2.begin();
        Image::Accessor acc = constimg.accessor();
        
        for(; i1 != i1end; ++i1, ++i2, ++i3)
        {
            double diff = factor * (acc(i2) - acc(i1));
            should(std::abs(diff - acc(i3)) < epsilon);
        }
    }
    
    void nonlinearDiffusionTest()
    {
        double epsilon = 0.0000001;
        
        Image res(lenna.size());

        nonlinearDiffusion(srcImageRange(lenna), destImage(res),
                           vigra::DiffusivityFunctor<double>(4.0), 4.0);
        
        ViffImage * viff = readViffImage("lenna128nonlinear.xv");
        shouldMsg(viff != 0, "Unable to read test image \"lenna128nonlinear.xv\"\n");
        Image comp(lenna.size());
        importViffImage(viff, destImage(comp));        
        freeViffImage(viff);
        
        Image::ScanOrderIterator i1 = res.begin();
        Image::ScanOrderIterator i1end = res.end();
        Image::ScanOrderIterator i2 = comp.begin();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            should(std::abs(*i1 - *i2) < epsilon);
        }
    }
    
    Image constimg, lenna, rampimg;
};

struct ConvolutionTestSuite
: public vigra::test_suite
{
    ConvolutionTestSuite()
    : vigra::test_suite("ConvolutionTestSuite")
    {
        add( testCase( &ConvolutionTest::separableConvolutionTest));
        add( testCase( &ConvolutionTest::separableDerivativeRepeatTest));
        add( testCase( &ConvolutionTest::separableDerivativeReflectTest));
        add( testCase( &ConvolutionTest::separableDerivativeAvoidTest));
        add( testCase( &ConvolutionTest::separableSmoothClipTest));
        add( testCase( &ConvolutionTest::separableSmoothWrapTest));
        add( testCase( &ConvolutionTest::gaussianSmoothingTest));
        add( testCase( &ConvolutionTest::separableGradientTest));
        add( testCase( &ConvolutionTest::gradientTest));
        add( testCase( &ConvolutionTest::hessianTest));
        add( testCase( &ConvolutionTest::structureTensorTest));
        add( testCase( &ConvolutionTest::stdConvolutionTest));
        add( testCase( &ConvolutionTest::stdVersusSeparableConvolutionTest));
        add( testCase( &ConvolutionTest::recursiveSmoothTest));
        add( testCase( &ConvolutionTest::recursiveGradientTest));
        add( testCase( &ConvolutionTest::recursiveSecondDerivativeTest));
        add( testCase( &ConvolutionTest::nonlinearDiffusionTest));
    }
};

int main()
{
    ConvolutionTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

