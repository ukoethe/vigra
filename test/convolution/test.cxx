#include <iostream>
#include <cmath>
#include "unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/convolution.hxx"
#include "vigra/impex.hxx"
#include "vigra/combineimages.hxx"

using namespace vigra;

vigra::DImage getSymmetricImage(){
    vigra::DImage src(40, 1);

    vigra::DImage::Accessor acc_src = src.accessor();
    vigra::DImage::iterator iter_src = src.begin();

    int i = 0;
    for ( ; i < 20 ; i++, iter_src++){
      acc_src.set(i + 0.25, iter_src);
    }
    i--;
    for ( ; i >= 0 ; i--, iter_src++){
      acc_src.set(i + 0.25, iter_src);
    }
    return src;
}

struct ConvolutionTest
{
    typedef vigra::DImage Image;

    ConvolutionTest()
    : constimg(5,5), rampimg(5,1)
    {
        constimg.init(1.0);
        
        ImageImportInfo info("lenna128.xv");

        lenna.resize(info.width(), info.height());
        importImage(info, destImage(lenna));

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
        importImage(vigra::ImageImportInfo("lenna128sepgrad.tif"), destImage(sepgrad));
        
        vigra::Kernel1D<double> gauss;
        gauss.initGaussian(1.0);
        vigra::Kernel1D<double> grad;
        grad.initGaussianDerivative(1.0, 1);

        Image tmp1(lenna);
        tmp1.init(0.0);
        Image tmp2(lenna);
        tmp2.init(0.0);
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

            shouldEqualTolerance(grad, acc(i), 1e-7);
        }
    }
    
    void gradientTest()
    {
        Image sepgrad(lenna.size());
        importImage(vigra::ImageImportInfo("lenna128sepgrad.tif"), destImage(sepgrad));
        
        
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
            shouldEqualTolerance(grad, acc(i), 1e-7);
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
            
        importImage(vigra::ImageImportInfo("lennahessxx.tif"), destImage(refxx));
        importImage(vigra::ImageImportInfo("lennahessyy.tif"), destImage(refyy));
        importImage(vigra::ImageImportInfo("lennahessxy.tif"), destImage(refxy));
        
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
                shouldEqualTolerance(acc(i1), acc(r1), 1e-7);
                shouldEqualTolerance(acc(i2), acc(r2), 1e-7);
                shouldEqualTolerance(acc(i3), acc(r3), 1e-7);
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
            
        importImage(vigra::ImageImportInfo("lennastxx.tif"), destImage(refxx));
        importImage(vigra::ImageImportInfo("lennastyy.tif"), destImage(refyy));
        importImage(vigra::ImageImportInfo("lennastxy.tif"), destImage(refxy));
        
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
                shouldEqualTolerance(acc(i1), acc(r1), 1e-7);
                shouldEqualTolerance(acc(i2), acc(r2), 1e-7);
                shouldEqualTolerance(acc(i3), acc(r3), 1e-7);
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
                shouldEqualTolerance(acc(x1), acc(x2), 1e-7);
            }
        }
    }
    
    void recursiveFilterTestWithAvoid()
    {
        Image src_const(25, 25);
        Image dest(src_const);
        src_const.init(42.1);
        dest.init(1.12);

        recursiveFilterX(srcImageRange(src_const), destImage(dest), VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_AVOID);

        Image::Iterator is = src_const.upperLeft();
        Image::Iterator tmp_src = src_const.upperLeft();
        Image::Iterator isend = src_const.lowerRight();
        Image::Iterator id = dest.upperLeft();
        Image::Iterator tmp_dest = dest.upperLeft();
        Image::Accessor acc = src_const.accessor();

        for(int y = 0; is.y != isend.y; is.y++, id.y++, y++)
        {
            for (int x = 0 ; is.x != isend.x ; is.x++, id.x++, x++)
            {
                if (x < 11 || x > 13)
                {
                    should (acc(id) == 1.12);
                }
                else
                {
                    should (acc(id) == acc(is));
                }
            }
            is.x = tmp_src.x;
            id.x = tmp_dest.x;
        }

        recursiveFilterY(srcImageRange(src_const), destImage(dest), 
                         VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_AVOID);

        is = src_const.upperLeft();
        tmp_src = src_const.upperLeft();
        isend = src_const.lowerRight();
        id = dest.upperLeft();
        tmp_dest = dest.upperLeft();
        acc = src_const.accessor();

        for(int y = 0; is.y != isend.y; is.y++, id.y++, y++)
        {
            for (int x = 0 ; is.x != isend.x ; is.x++, id.x++, x++)
            {
                if ((x > 10 && x < 14) || (y > 10 && y < 14))
                {
                    should(acc(is) == acc(id));
                }
                else
                {
                    should(acc(id) == 1.12);
                }
            }
            is.x = tmp_src.x;
            id.x = tmp_dest.x;
        }

        // Hier wird an einem symmetrischen Bild /\ getestet 
        // ob die korrekten Daten eingehalten wurden.

        Image src(getSymmetricImage());
        dest = src;
        recursiveFilterX(srcImageRange(src), destImage(dest), 
               VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_AVOID);
        Image::value_type correct_data[40] = 
            {0.25, 1.25, 2.25, 3.25, 4.25, 5.25, 6.25, 7.25, 8.25, 9.25, 
             10.25, 11.249812, 12.249472, 13.248558, 14.246079, 15.239341, 
             16.221025, 17.171238, 18.035903, 18.668023, 
             18.668023, 18.035903, 17.171238, 16.221025, 15.239341, 
             14.246079, 13.248558, 12.249472, 11.249812, 10.25, 9.25, 
             8.25, 7.25, 6.25, 5.25, 4.25, 3.25, 2.25, 1.25, 0.25};

        Image::iterator dest_iter = dest.begin();

        for (int i = 0 ; i < 40; i++)
        {
            shouldEqualTolerance (correct_data[i], dest_iter[i], 1e-7);
        }

    }// end of recursiveFilterTestWithAvoid

    void recursiveFilterTestWithReflect()
    {
        Image src_const(25, 25);
        Image dest(src_const);
        src_const.init(42.1);
        dest.init(1.12);

        recursiveFilterX(srcImageRange(src_const), destImage(dest), 
                     VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_REFLECT);

        Image::Iterator is = src_const.upperLeft();
        Image::Iterator tmp_src = src_const.upperLeft();
        Image::Iterator isend = src_const.lowerRight();
        Image::Iterator id = dest.upperLeft();
        Image::Iterator tmp_dest = dest.upperLeft();
        Image::Accessor acc = src_const.accessor();

        for(int y = 0; is.y != isend.y; is.y++, id.y++, y++)
        {
            for (int x = 0 ; is.x != isend.x ; is.x++, id.x++, x++)
            {
                should (acc(id) == acc(is));
            }
            is.x = tmp_src.x;
            id.x = tmp_dest.x;
        }

        recursiveFilterY(srcImageRange(dest), destImage(dest), 
                   VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_REFLECT);

        is = src_const.upperLeft();
        tmp_src = src_const.upperLeft();
        isend = src_const.lowerRight();
        id = dest.upperLeft();
        tmp_dest = dest.upperLeft();
        acc = src_const.accessor();

        for(int y = 0; is.y != isend.y; is.y++, id.y++, y++)
        {
            for (int x = 0 ; is.x != isend.x ; is.x++, id.x++, x++)
            {
                should(acc(is) == acc(id));
            }
            is.x = tmp_src.x;
            id.x = tmp_dest.x;
        }

        // Hier wird an einem symmetrischen Bild /\ (Groesse 40x1)  getestet 
        // ob die korrekten Daten eingehalten wurden.

        Image src(getSymmetricImage());
        dest = src;
        recursiveFilterX(srcImageRange(src), destImage(dest), 
                VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_REFLECT);
        Image::value_type correct_data[40] = 
           {1.100911, 1.5630327, 2.3651583, 3.2923643, 4.2655848, 
           5.2557329, 6.2521079, 7.2507724, 8.2502757, 9.2500786, 
           10.249967, 11.249819, 12.249475, 13.248559, 14.246079, 
           15.239341, 16.221025, 17.171238, 18.035903, 18.668023, 
           18.668023, 18.035903, 17.171238, 16.221025, 15.239341, 
           14.246079, 13.248559, 12.249475, 11.249819, 10.249967, 
           9.2500786, 8.2502757, 7.2507724, 6.2521079, 5.2557329, 
           4.2655848, 3.2923643, 2.3651583, 1.5630327, 1.100911};

        Image::iterator dest_iter = dest.begin();
        for (int i = 0 ; i < 40; i++)
        {
            shouldEqualTolerance (correct_data[i], dest_iter[i], 1e-7);
        }
    }// end of recursiveFilterTestWithReflect


    void recursiveFilterTestWithClipOnConstImage()
    {
        Image src_const(25, 25);
        Image dest(src_const);
        src_const.init(42.1);
        dest.init(1.12);
        recursiveFilterX(srcImageRange(src_const), destImage(dest), 
                           VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_CLIP);

        Image::Iterator is = src_const.upperLeft();
        Image::Iterator tmp_src = src_const.upperLeft();
        Image::Iterator isend = src_const.lowerRight();
        Image::Iterator id = dest.upperLeft();
        Image::Iterator tmp_dest = dest.upperLeft();
        Image::Accessor acc = src_const.accessor();

        for(int y = 0; is.y != isend.y; is.y++, id.y++, y++)
        {
            for (int x = 0 ; is.x != isend.x ; is.x++, id.x++, x++)
            {
                shouldEqualTolerance (acc(id), acc(is), 0.000001);
            }
            is.x = tmp_src.x;
            id.x = tmp_dest.x;
        }

        recursiveFilterY(srcImageRange(src_const), destImage(dest), 
                VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_CLIP);

        is = src_const.upperLeft();
        tmp_src = src_const.upperLeft();
        isend = src_const.lowerRight();
        id = dest.upperLeft();
        tmp_dest = dest.upperLeft();
        acc = src_const.accessor();

        for(int y = 0; is.y != isend.y; is.y++, id.y++, y++)
        {
            for (int x = 0 ; is.x != isend.x ; is.x++, id.x++, x++)
            {
                shouldEqualTolerance (acc(id), acc(is), 0.000001);
            }
            is.x = tmp_src.x;
            id.x = tmp_dest.x;
        }
    }// end of recursiveFilterTestWithClipOnConstImage()

    void recursiveFilterTestWithClipOnNonConstImage(){
        Image src(40, 1);
        Image dest(src);
        Image::value_type correct_data[40] = 
            {0.831977, 1.53351, 2.3853, 3.31218, 4.27763, 5.26195, 
            6.25506, 7.25211, 8.25086, 9.25035, 10.2501, 11.2501, 
            12.25, 13.25, 14.25, 15.25, 16.25, 17.25, 18.25, 19.25, 
            20.25, 21.25, 22.25, 23.25, 24.25, 25.25, 26.25, 27.25, 
            28.2499, 29.2499, 30.2496, 31.2491, 32.2479, 33.2449, 
            34.2381, 35.2224, 36.1878, 37.1147, 37.9665, 38.668};

        Image::Accessor acc_src = src.accessor();
        Image::iterator iter_src = src.begin();

        for (int i = 0 ; i < 40 ; i++, iter_src++)
        {
            acc_src.set(i + 0.25, iter_src);
        }

        recursiveFilterX(srcImageRange(src), destImage(dest), 
                VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_CLIP);

        Image::iterator idest = dest.begin();
        Image::iterator idest_end = dest.end();
        Image::Accessor dest_acc = dest.accessor();

        for(int i = 0; idest != idest_end; idest++, i++)
        {
            shouldEqualTolerance (dest_acc(idest), correct_data[i], 0.00001);
        }

    }// end of recursiveFilterTestWithClipOnNonConstImage

    void recursiveFilterTestWithWrap()
    {
        Image src_const(25, 25);
        Image dest(src_const);
        src_const.init(42.1);
        dest.init(1.12);

        recursiveFilterX(srcImageRange(src_const), destImage(dest), 
                VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_WRAP);

        Image::Iterator is = src_const.upperLeft();
        Image::Iterator tmp_src = src_const.upperLeft();
        Image::Iterator isend = src_const.lowerRight();
        Image::Iterator id = dest.upperLeft();
        Image::Iterator tmp_dest = dest.upperLeft();
        Image::Accessor acc = src_const.accessor();

        for(int y = 0; is.y != isend.y; is.y++, id.y++, y++)
        {
            for (int x = 0 ; is.x != isend.x ; is.x++, id.x++, x++)
            {
                should (acc(id) == acc(is));
            }
            is.x = tmp_src.x;
            id.x = tmp_dest.x;
        }

        recursiveFilterY(srcImageRange(src_const), destImage(dest), 
                    VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_WRAP);

        is = src_const.upperLeft();
        tmp_src = src_const.upperLeft();
        isend = src_const.lowerRight();
        id = dest.upperLeft();
        tmp_dest = dest.upperLeft();
        acc = src_const.accessor();

        for(int y = 0; is.y != isend.y; is.y++, id.y++, y++)
        {
            for (int x = 0 ; is.x != isend.x ; is.x++, id.x++, x++)
            {
                should(acc(is) == acc(id));
            }
            is.x = tmp_src.x;
            id.x = tmp_dest.x;
        }


        // Hier wird an einem symmetrischen Bild /\ (Groesse 40x1)  getestet 
        // ob die korrekten Daten eingehalten wurden.

        Image src(getSymmetricImage());
        dest = src;
        recursiveFilterX(srcImageRange(src), destImage(dest), 
                    VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_WRAP);
        Image::value_type correct_data[40] = 
            {1.100911, 1.5630327, 2.3651583, 3.2923643, 4.2655848, 
            5.2557329, 6.2521079, 7.2507724, 8.2502757, 9.2500786, 
            10.249967, 11.249819, 12.249475, 13.248559, 14.246079, 
            15.239341, 16.221025, 17.171238, 18.035903, 18.668023, 
            18.668023, 18.035903, 17.171238, 16.221025, 15.239341, 
            14.246079, 13.248559, 12.249475, 11.249819, 10.249967, 
            9.2500786, 8.2502757, 7.2507724, 6.2521079, 5.2557329, 
            4.2655848, 3.2923643, 2.3651583, 1.5630327, 1.100911};

        Image::iterator dest_iter = dest.begin();

        for (int i = 0 ; i < 40; i++)
        {
            shouldEqualTolerance (correct_data[i], dest_iter[i], 1e-7);
        }

    }// end of recursiveFilterTestWithWrap

    void recursiveFilterTestWithRepeat()
    {
        Image src_const(25, 25);
        Image dest(src_const);
        src_const.init(42.1);
        dest.init(1.12);

        recursiveFilterX(srcImageRange(src_const), destImage(dest), 
                    VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_REPEAT);

        Image::Iterator is = src_const.upperLeft();
        Image::Iterator tmp_src = src_const.upperLeft();
        Image::Iterator isend = src_const.lowerRight();
        Image::Iterator id = dest.upperLeft();
        Image::Iterator tmp_dest = dest.upperLeft();
        Image::Accessor acc = src_const.accessor();

        for(int y = 0; is.y != isend.y; is.y++, id.y++, y++)
        {
            for (int x = 0 ; is.x != isend.x ; is.x++, id.x++, x++)
            {
                should (acc(id) == acc(is));
            }
            is.x = tmp_src.x;
            id.x = tmp_dest.x;
        }
        recursiveFilterY(srcImageRange(src_const), destImage(dest), 
                VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_REPEAT);

        is = src_const.upperLeft();
        tmp_src = src_const.upperLeft();
        isend = src_const.lowerRight();
        id = dest.upperLeft();
        tmp_dest = dest.upperLeft();
        acc = src_const.accessor();

        for(int y = 0; is.y != isend.y; is.y++, id.y++, y++)
        {
            for (int x = 0 ; is.x != isend.x ; is.x++, id.x++, x++)
            {
                should(acc(is) == acc(id));
            }
            is.x = tmp_src.x;
            id.x = tmp_dest.x;
        }

        // Hier wird an einem symmetrischen Bild /\ (Groesse 40x1)  getestet 
        // ob die korrekten Daten eingehalten wurden.

        Image src(getSymmetricImage());
        dest = src;
        recursiveFilterX(srcImageRange(src), destImage(dest), 
                VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_REPEAT);
        Image::value_type correct_data[40] = 
            {0.67545906, 1.4065176, 2.3075796, 3.2711823, 4.2577924, 
            5.2528662, 6.2510533, 7.2503844, 8.250133, 9.2500261, 
            10.249947, 11.249812, 12.249472, 13.248558, 14.246079, 
            15.239341, 16.221025, 17.171238, 18.035903, 18.668023, 
            18.668023, 18.035903, 17.171238, 16.221025, 15.239341, 
            14.246079, 13.248558, 12.249472, 11.249812, 10.249947, 
            9.2500261, 8.250133, 7.2503844, 6.2510533, 5.2528662, 
            4.2577924, 3.2711823, 2.3075796, 1.4065176, 0.67545906};

        Image::iterator dest_iter = dest.begin();

        for (int i = 0 ; i < 40; i++)
        {
            shouldEqualTolerance (correct_data[i], dest_iter[i], 1e-7);
        }

    }// end of recursiveFilterTestWithRepeat

    void recursiveFilterTestFromWrapWithReflect()
    {

        Image src_wrap(79, 1);
        Image src_reflect(40, 1);
        Image dest_wrap(src_wrap);
        Image dest_reflect(src_reflect);

        Image::Accessor acc_src_wrap = src_wrap.accessor();
        Image::iterator iter_src_wrap = src_wrap.begin();
        Image::Accessor acc_src_reflect = src_reflect.accessor();
        Image::iterator iter_src_reflect = src_reflect.begin();

        for (int i = 0 ; i < 40 ; i++, iter_src_wrap++, iter_src_reflect++)
        {
            acc_src_wrap.set(i + 0.25, iter_src_wrap);
            acc_src_reflect.set(i + 0.25, iter_src_reflect);
        }
        for (int j = 38 ; j >= 0 ; j--, iter_src_wrap++)
        {
            acc_src_wrap.set( j + 0.25, iter_src_wrap);
        }

        recursiveFilterX(srcImageRange(src_wrap), destImage(dest_wrap), 
                VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_WRAP);

        recursiveFilterX(srcImageRange(src_reflect), destImage(dest_reflect), 
                VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_REFLECT);

        Image::iterator iter_dest_wrap = dest_wrap.begin();
        Image::Accessor acc_dest_wrap = dest_wrap.accessor();
        Image::iterator iter_dest_reflect = dest_reflect.begin();
        Image::iterator end_dest_reflect = dest_reflect.end();
        Image::Accessor acc_dest_reflect = dest_reflect.accessor();

        while(iter_dest_reflect != end_dest_reflect)
        {
            shouldEqualTolerance(acc_dest_wrap(iter_dest_wrap), 
                           acc_dest_reflect(iter_dest_reflect), 1e-6);
            iter_dest_wrap++;
            iter_dest_reflect++;
        }
    }

    void recursiveFilterTestFromRepeatWithAvoid()
    {
        Image src_avoid(40, 1);
        src_avoid.init(11.47);
        Image src_repeat(18, 1);

        Image dest_repeat(src_repeat);
        Image dest_avoid(src_avoid);

        Image::Accessor acc_src_avoid = src_avoid.accessor();
        Image::iterator iter_src_avoid = src_avoid.begin();
        Image::Accessor acc_src_repeat = src_repeat.accessor();
        Image::iterator iter_src_repeat = src_repeat.begin();

        int i = 0;
        for ( ; i < 20 ; i++, iter_src_avoid++)
        {
            if(i > 10)
            {
                acc_src_repeat.set(i + 0.47, iter_src_repeat);
                acc_src_avoid.set(i + 0.47, iter_src_avoid);
                iter_src_repeat++;
            }
        }
        i--;

        for ( ; i >= 0 ; i--, iter_src_avoid++)
        {
            if(i > 10)
            {
                acc_src_repeat.set(i + 0.47, iter_src_repeat);
                acc_src_avoid.set(i + 0.47, iter_src_avoid);
                iter_src_repeat++;
            }
        }

        recursiveFilterX(srcImageRange(src_repeat), destImage(dest_repeat), VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_REPEAT);
        recursiveFilterX(srcImageRange(src_avoid), destImage(dest_avoid), VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_AVOID);

        Image::Accessor acc_dest_repeat = dest_repeat.accessor();
        Image::iterator dest_iter_repeat = dest_repeat.begin();

        Image::Accessor acc_dest_avoid = dest_avoid.accessor();
        Image::iterator dest_iter_avoid = dest_avoid.begin();

        for (int i = 0 ; i < 39 ; i++, dest_iter_avoid++)
        {
            if (i < 11 || i > 28)
            {
                should(acc_dest_avoid(dest_iter_avoid) == 11.47);
            }
            else
            {
                should(acc_dest_avoid(dest_iter_avoid) == acc_dest_repeat(dest_iter_repeat));
                dest_iter_repeat++;
            }
        }
    }

    /**
   * Es wird die Positionierung der einzelnen 
   * Punkte relativ zueinander getestet.
   */
    void recursiveFilterTestOfAllTreatmentsRelatively()
    {
        Image src(40, 1);

        Image::Accessor acc_src = src.accessor();
        Image::iterator iter_src = src.begin();

        int i = 0;
        for ( ; i < 20 ; i++, iter_src++)
        {
            acc_src.set(i + 0.25, iter_src);
        }
        i--;
        for ( ; i >= 0 ; i--, iter_src++)
        {
            acc_src.set(i + 0.25, iter_src);
        }

        Image dest_avoid(src);
        Image dest_repeat(src);
        Image dest_reflect(src);
        Image dest_wrap(src);
        Image dest_clip(src);

        recursiveFilterX(srcImageRange(src), destImage(dest_avoid), 
                VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_AVOID);
        recursiveFilterX(srcImageRange(src), destImage(dest_repeat), 
                VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_REPEAT);
        recursiveFilterX(srcImageRange(src), destImage(dest_reflect), 
                VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_REFLECT);
        recursiveFilterX(srcImageRange(src), destImage(dest_wrap), 
                VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_WRAP);
        recursiveFilterX(srcImageRange(src), destImage(dest_clip), 
                VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_CLIP);

        iter_src = src.begin();
        Image::iterator iter_dest_avoid = dest_avoid.begin();
        Image::iterator iter_dest_repeat = dest_repeat.begin();
        Image::iterator iter_dest_reflect = dest_reflect.begin();
        Image::iterator iter_dest_wrap = dest_wrap.begin();
        Image::iterator iter_dest_clip = dest_clip.begin();

        for (int x = 0 ;  x < 40 ; x++)
        {
            if(x > 9 && x < 30 )
            {
                shouldEqualTolerance(iter_dest_avoid[x], iter_dest_repeat[x], 1e-5);
                shouldEqualTolerance(iter_dest_avoid[x], iter_dest_reflect[x], 1e-5);
                shouldEqualTolerance(iter_dest_avoid[x], iter_dest_wrap[x], 1e-5);
                shouldEqualTolerance(iter_dest_avoid[x], iter_dest_clip[x], 1e-5);
            }
            else
            {
                should(iter_dest_avoid[x] == iter_src[x]);
                should(iter_dest_repeat[x] < iter_dest_reflect[x]);
                should(iter_dest_repeat[x] < iter_dest_clip[x]);
                if (x < 2 || x > 37)
                {
                    should(iter_dest_clip[x] < iter_dest_reflect[x]);
                }
                else
                {
                    should(iter_dest_clip[x] > iter_dest_reflect[x]);
                }
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
        ImageImportInfo info("lenna128recgrad.xv");

        Image recgrad(info.width(), info.height());
        importImage(info, destImage(recgrad));
        
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
            
            shouldEqualTolerance(grad, acc(i), 1e-7);
        }
    }
    
    void recursiveSecondDerivativeTest()
    {
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
            shouldEqualTolerance(diff, acc(i3), 1e-7);
        }
    }
    
    void nonlinearDiffusionTest()
    {
         
        Image res(lenna.size());

        nonlinearDiffusion(srcImageRange(lenna), destImage(res),
                           vigra::DiffusivityFunctor<double>(4.0), 4.0);
                
        Image comp(lenna.size());
        importImage(vigra::ImageImportInfo("lenna128nonlinear.tif"), destImage(comp));

        Image::ScanOrderIterator i1 = res.begin();
        Image::ScanOrderIterator i1end = res.end();
        Image::ScanOrderIterator i2 = comp.begin();
        
        for(; i1 != i1end; ++i1, ++i2)
        {
            shouldEqualTolerance((*i1), (*i2), 1e-7);
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
        add( testCase( &ConvolutionTest::recursiveFilterTestWithAvoid));
        add( testCase( &ConvolutionTest::recursiveFilterTestWithClipOnConstImage));
        add( testCase( &ConvolutionTest::recursiveFilterTestWithClipOnNonConstImage));
        add( testCase( &ConvolutionTest::recursiveFilterTestWithReflect));
        add( testCase( &ConvolutionTest::recursiveFilterTestWithWrap));
        add( testCase( &ConvolutionTest::recursiveFilterTestWithRepeat));
        add( testCase( &ConvolutionTest::recursiveFilterTestFromWrapWithReflect));
        add( testCase( &ConvolutionTest::recursiveFilterTestFromRepeatWithAvoid));
        add( testCase( &ConvolutionTest::recursiveFilterTestOfAllTreatmentsRelatively));
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

