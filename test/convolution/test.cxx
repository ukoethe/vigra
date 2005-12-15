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
#include <cmath>
#include "vigra/convolution.hxx"
#include "unittest.hxx"
#include "vigra/stdimage.hxx"
#include "vigra/impex.hxx"
#include "vigra/combineimages.hxx"
#include "vigra/resampling_convolution.hxx"

using namespace vigra;

vigra::DImage getSymmetricLine(){
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
vigra::DImage getSymmetricImage(){
    vigra::DImage src(10, 6);

    vigra::DImage::Accessor acc_src = src.accessor();
    vigra::DImage::iterator iter_src = src.begin();

    int i = 0;
    for(int y = 0; y < 6; y++){
        int x = 0;
        for ( ; x < 5 ; i++, x++, iter_src++){
            acc_src.set(i + 0.25, iter_src);
        }
        i--;
        for ( ; x < 10 ; i--, x++, iter_src++){
            acc_src.set(i + 0.25, iter_src);
        }
        (y%2 == 0)? i += 3 : i+= 2;
    }
    return src;
}

vigra::DImage getUnsymmetricImage(){
    vigra::DImage src(10, 6);

    vigra::DImage::Accessor acc_src = src.accessor();
    vigra::DImage::iterator iter_src = src.begin();

    int i = 0;
    for(int y = 0; y < 6; y++){
        for (int x = 0 ; x < 10 ; i++, x++, iter_src++){
            acc_src.set(i + 0.25, iter_src);
        }
        (y%2 == 0)? i++ : i+= 2;
    }
    return src;
}

struct ConvolutionTest
{
    typedef vigra::DImage Image;

    ConvolutionTest()
    : constimg(5,5), 
      rampimg(5,1), 
      sym_image(getSymmetricImage()), 
      unsym_image(getUnsymmetricImage())
    {
        constimg.init(1.0);
        
        vigra::Kernel1D<double> binom1;
        binom1.initBinomial(1);
        sym_kernel.initSeparable(binom1, binom1);

        unsym_kernel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) = 1,   2,   4,
                                                                  5,   11,  3,
                                                                  6,   8,   7 ;
        unsym_kernel.normalize(1);

        line_kernel.initExplicitly(Diff2D(-2,0), Diff2D(2,0)) = 1, 4,   12,   4, 1 ;
        line_kernel.normalize(1);
        
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

    void simpleSharpeningTest(){
         Image dest_lenna(lenna);
        simpleSharpening(srcImageRange(lenna), destImage(dest_lenna), 0.8);

        Image dest_correct;
        ImageImportInfo info_correct("lenna_simple_sharpening_orig.xv");
        dest_correct.resize(info_correct.width(), info_correct.height());
        importImage(info_correct, destImage(dest_correct));

        Image::iterator i_dest = dest_lenna.begin();
        Image::iterator end_dest = dest_lenna.end();
        Image::Accessor acc_dest = dest_lenna.accessor();
        Image::iterator i_dest_correct = dest_correct.begin();
        Image::Accessor acc_dest_correct = dest_correct.accessor();

        for (; i_dest != end_dest; i_dest++, i_dest_correct++ ){
            shouldEqualTolerance(acc_dest(i_dest) , acc_dest_correct(i_dest_correct), 1e-15);
        }

    }

    void gaussianSharpeningTest(){
         Image dest_lenna(lenna);
        gaussianSharpening(srcImageRange(lenna), destImage(dest_lenna), 3., 0.7);

        Image dest_correct;
        ImageImportInfo info_correct("lenna_gaussian_sharpening_orig.xv");
        dest_correct.resize(info_correct.width(), info_correct.height());
        importImage(info_correct, destImage(dest_correct));

        Image::iterator i_dest = dest_lenna.begin();
        Image::iterator end_dest = dest_lenna.end();
        Image::Accessor acc_dest = dest_lenna.accessor();
        Image::iterator i_dest_correct = dest_correct.begin();
        Image::Accessor acc_dest_correct = dest_correct.accessor();

        for (; i_dest != end_dest; i_dest++, i_dest_correct++ ){
            shouldEqualTolerance(acc_dest(i_dest) , acc_dest_correct(i_dest_correct), 1e-6);
        }

    }

    void stdConvolutionTestOnConstImage()
    {        
        Image tmp_clip(constimg);
        tmp_clip = 0.0;
        Image tmp_wrap(constimg);
        tmp_wrap = 0.0;
        Image tmp_avoid(constimg);
        tmp_avoid = 0.0;
        Image tmp_repeat(constimg);
        tmp_repeat = 0.0;
        Image tmp_reflect(constimg);
        tmp_reflect = 0.0;

        convolveImage(srcImageRange(constimg), destImage(tmp_clip), kernel2d(sym_kernel, BORDER_TREATMENT_CLIP));

        convolveImage(srcImageRange(constimg), destImage(tmp_avoid), kernel2d(sym_kernel, BORDER_TREATMENT_AVOID));
        
        convolveImage(srcImageRange(constimg), destImage(tmp_wrap), kernel2d(sym_kernel, BORDER_TREATMENT_WRAP));

        convolveImage(srcImageRange(constimg), destImage(tmp_repeat), kernel2d(sym_kernel, BORDER_TREATMENT_REPEAT));

        convolveImage(srcImageRange(constimg), destImage(tmp_reflect), kernel2d(sym_kernel, BORDER_TREATMENT_REFLECT));
        
        Image::ScanOrderIterator i_src = constimg.begin();
        Image::ScanOrderIterator i_src_end = constimg.end();
        Image::ScanOrderIterator i_clip = tmp_clip.begin();
        Image::ScanOrderIterator i_avoid = tmp_avoid.begin();
        Image::ScanOrderIterator i_wrap = tmp_wrap.begin();
        Image::ScanOrderIterator i_repeat = tmp_repeat.begin();
        Image::ScanOrderIterator i_reflect = tmp_reflect.begin();
        Image::Accessor acc = constimg.accessor();
        
        for(int y = 0; i_src != i_src_end; y++){
            for(int x = 0; x < constimg.size().x; x++, ++i_src, ++i_clip, ++i_wrap, ++i_repeat, ++i_reflect, ++i_avoid){
                should(acc(i_src) == acc(i_clip));
                should(acc(i_src) == acc(i_wrap));
                should(acc(i_src) == acc(i_repeat));
                should(acc(i_src) == acc(i_reflect));
                if(x != 0 && y != 0 && x != 4 && y != 4){
                    should(acc(i_src) == acc(i_avoid)); 
                }else{
                    should(acc(i_avoid) == 0);
                }
            }
        }


        //@todo Import<->Export funktioniert nicht wie frueher???
//         Image dest_lenna(lenna);
//         convolveImage(srcImageRange(lenna), destImage(dest_lenna), kernel2d(sym_kernel, BORDER_TREATMENT_REFLECT));
//         exportImage(srcImageRange(dest_lenna), ImageExportInfo("lenna_convolve_128x120.xv"));

    }
    

    void stdConvolutionTestWithAvoid(){
        Image dest(sym_image);
        dest.init(42.1);

        convolveImage(srcImageRange(sym_image), destImage(dest), kernel2d(sym_kernel, BORDER_TREATMENT_AVOID));

        Image::Iterator i_dest_2D = dest.upperLeft();
        Image::ScanOrderIterator i_dest = dest.begin();
        Image::ScanOrderIterator i_dest_end = dest.end();
        Image::Accessor acc = dest.accessor();

        //Kontrollierung der Randbehandlung und ein paar Pixel
        //aus der Mitte.
        should (acc(i_dest_2D + Diff2D(1,1)) == 3);
        should (acc(i_dest_2D + Diff2D(8,2)) == 4.5);
        for(int y = 0; i_dest != i_dest_end; y++){
            for(int x = 0; x < dest.size().x; x++, ++i_dest){
                if(x == 0 || y == 0 || x == 9 || y == 9){
                    should(acc(i_dest) == 42.1);
                }
            }
        }

    }
    
    void stdConvolutionTestWithClip(){
        Image dest(sym_image);
        dest.init(42.1);

        convolveImage(srcImageRange(sym_image), destImage(dest), kernel2d(unsym_kernel, BORDER_TREATMENT_CLIP));

        Image::Iterator i_dest_2D = dest.upperLeft();
        Image::ScanOrderIterator i_dest = dest.begin();
        Image::ScanOrderIterator i_dest_end = dest.end();
        Image::Accessor acc = dest.accessor();

        //Kontrollierung der Randbehandlung und ein paar Pixel
        //aus der Mitte.
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(0,0)), 0.88157,   1E-4);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(3,0)), 3.75,      1E-4);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(9,0)), 1.199999,  1E-4);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(9,2)), 3.5642913, 1E-4);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(9,5)), 7.5603463, 1E-4);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(3,5)), 10.225,    1E-4);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(0,5)), 7.6833271, 1E-4);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(0,2)), 3.3712039, 1E-4);

        shouldEqualTolerance (acc(i_dest_2D + Diff2D(3,2)), 6.05852,   1E-4);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(6,4)), 9.14363,   1E-4);
    }
    
    void stdConvolutionTestWithWrap(){
        Image dest(unsym_image);
        dest.init(42.1);

        convolveImage(srcImageRange(unsym_image), destImage(dest), kernel2d(unsym_kernel, BORDER_TREATMENT_WRAP));

        Image::Iterator i_dest_2D = dest.upperLeft();
        Image::ScanOrderIterator i_dest = dest.begin();
        Image::ScanOrderIterator i_dest_end = dest.end();
        Image::Accessor acc = dest.accessor();

        //Kontrollierung der Randbehandlung und ein paar Pixel
        //aus der Mitte.
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(0,0)), 30.2926,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(3,0)), 30.3138,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(9,0)), 33.7606,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(9,2)), 25.9309,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(9,5)), 50.25,     1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(3,5)), 46.8032,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(0,5)), 46.7819,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(0,2)), 22.4628,   1E-5);

        shouldEqualTolerance (acc(i_dest_2D + Diff2D(2,3)), 33.0798,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(6,4)), 48.4841,   1E-5);
    }
    
    void stdConvolutionTestWithReflect(){
        Image dest(unsym_image);
        dest.init(42.1);

        convolveImage(srcImageRange(unsym_image), destImage(dest), kernel2d(unsym_kernel, BORDER_TREATMENT_REFLECT));

        Image::Iterator i_dest_2D = dest.upperLeft();
        Image::ScanOrderIterator i_dest = dest.begin();
        Image::ScanOrderIterator i_dest_end = dest.end();
        Image::Accessor acc = dest.accessor();

        //Kontrollierung der Randbehandlung und ein paar Pixel
        //aus der Mitte.
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(0,0)), 7.35638,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(3,0)), 9.76064,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(9,0)), 15.25,     1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(9,2)), 27.9734,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(9,5)), 59.1436,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(3,5)), 53.6543,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(0,5)), 51.25,     1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(0,2)), 20.0798,   1E-5);

        shouldEqualTolerance (acc(i_dest_2D + Diff2D(2,3)), 33.0798,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(6,4)), 48.4841,   1E-5);
    }
    
    void stdConvolutionTestWithRepeat(){
        Image dest(unsym_image);
        dest.init(42.1);

        convolveImage(srcImageRange(unsym_image), destImage(dest), kernel2d(unsym_kernel, BORDER_TREATMENT_REPEAT));

        Image::Iterator i_dest_2D = dest.upperLeft();
        Image::ScanOrderIterator i_dest = dest.begin();
        Image::ScanOrderIterator i_dest_end = dest.end();
        Image::Accessor acc = dest.accessor();

        //Kontrollierung der Randbehandlung und ein paar Pixel
        //aus der Mitte.
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(0,0)), 2.14362,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(3,0)), 4.84574,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(9,0)), 10.5904,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(9,2)), 28.2287,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(9,5)), 61.0372,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(3,5)), 55.2926,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(0,5)), 52.5904,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(0,2)), 19.7819,   1E-5);

        shouldEqualTolerance (acc(i_dest_2D + Diff2D(2,3)), 33.0798,   1E-5);
        shouldEqualTolerance (acc(i_dest_2D + Diff2D(6,4)), 48.4841,   1E-5);
    }

    void stdConvolutionTestFromWrapWithReflect(){
    
        Image src_wrap(78, 1);
        Image src_reflect(40, 1);
        Image dest_wrap(src_wrap);
        Image dest_reflect(src_reflect);

        Image::Accessor acc_src_wrap = src_wrap.accessor();
        Image::iterator iter_src_wrap = src_wrap.begin();
        Image::Accessor acc_src_reflect = src_reflect.accessor();
        Image::iterator iter_src_reflect = src_reflect.begin();

        for (int i = 0 ; i < 40 ; i++, iter_src_wrap++, iter_src_reflect++){
            acc_src_wrap.set(i + 0.25, iter_src_wrap);
            acc_src_reflect.set(i + 0.25, iter_src_reflect);
        }
        for (int j = 38 ; j >= 1 ; j--, iter_src_wrap++){
            acc_src_wrap.set( j + 0.25, iter_src_wrap);
        }
    
        convolveImage(srcImageRange(src_wrap), destImage(dest_wrap), kernel2d(line_kernel, BORDER_TREATMENT_WRAP));
        convolveImage(srcImageRange(src_reflect), destImage(dest_reflect), kernel2d(line_kernel, BORDER_TREATMENT_REFLECT));

        Image::iterator iter_dest_wrap = dest_wrap.begin();
        Image::Accessor acc_dest_wrap = dest_wrap.accessor();
        Image::iterator iter_dest_reflect = dest_reflect.begin();
        Image::iterator end_dest_reflect = dest_reflect.end();
        Image::Accessor acc_dest_reflect = dest_reflect.accessor();

        while(iter_dest_reflect != end_dest_reflect){
            shouldEqualTolerance(acc_dest_wrap(iter_dest_wrap), acc_dest_reflect(iter_dest_reflect), 1e-6);
            iter_dest_wrap++;
            iter_dest_reflect++;
        }
    }
  
    void stdConvolutionTestFromRepeatWithAvoid(){
        Image src_avoid(40, 1);
        src_avoid.init(2.47);
        Image src_repeat(36, 1);
    
        Image dest_repeat(src_repeat);
        Image dest_avoid(src_avoid);

        Image::Accessor acc_src_avoid = src_avoid.accessor();
        Image::iterator iter_src_avoid = src_avoid.begin();
        Image::Accessor acc_src_repeat = src_repeat.accessor();
        Image::iterator iter_src_repeat = src_repeat.begin();

        int i = 0;
        for ( ; i < 20 ; i++, iter_src_avoid++){
            if(i > 1){
                acc_src_repeat.set(i + 0.47, iter_src_repeat);
                acc_src_avoid.set(i + 0.47, iter_src_avoid);
                iter_src_repeat++;
            }
        }
        i--;

        for ( ; i >= 0 ; i--, iter_src_avoid++){
            if(i > 1){
                acc_src_repeat.set(i + 0.47, iter_src_repeat);
                acc_src_avoid.set(i + 0.47, iter_src_avoid);
                iter_src_repeat++;
            }
        }

        convolveImage(srcImageRange(src_repeat), destImage(dest_repeat), kernel2d(line_kernel, BORDER_TREATMENT_REPEAT));
        convolveImage(srcImageRange(src_avoid), destImage(dest_avoid), kernel2d(line_kernel, BORDER_TREATMENT_AVOID));

        Image::Accessor acc_dest_repeat = dest_repeat.accessor();
        Image::iterator dest_iter_repeat = dest_repeat.begin();

        Image::Accessor acc_dest_avoid = dest_avoid.accessor();
        Image::iterator dest_iter_avoid = dest_avoid.begin();

        for (int i = 0 ; i < 39 ; i++, dest_iter_avoid++){
            if (i < 2 || i > 37){
                should(acc_dest_avoid(dest_iter_avoid) == 2.47);
            }else{
                shouldEqualTolerance(acc_dest_avoid(dest_iter_avoid) , acc_dest_repeat(dest_iter_repeat), 1e-15);
                dest_iter_repeat++;
            }
        }
    }
  
    /**
     * Es wird die Positionierung der einzelnen 
     * Punkte relativ zueinander getestet.
     */
    void stdConvolutionTestOfAllTreatmentsRelatively(){
        Image src(40, 1);

        Image::Accessor acc_src = src.accessor();
        Image::iterator iter_src = src.begin();

        int i = 0;
        for ( ; i < 20 ; i++, iter_src++){
            acc_src.set(i + 0.25, iter_src);
        }
        i--;
        for ( ; i >= 0 ; i--, iter_src++){
            acc_src.set(i + 0.25, iter_src);
        }

        Image dest_avoid(src);
        Image dest_repeat(src);
        Image dest_reflect(src);
        Image dest_wrap(src);
        Image dest_clip(src);

        convolveImage(srcImageRange(src), destImage(dest_avoid), kernel2d(line_kernel, BORDER_TREATMENT_AVOID));
        convolveImage(srcImageRange(src), destImage(dest_repeat), kernel2d(line_kernel, BORDER_TREATMENT_REPEAT));
        convolveImage(srcImageRange(src), destImage(dest_reflect), kernel2d(line_kernel, BORDER_TREATMENT_REFLECT));
        convolveImage(srcImageRange(src), destImage(dest_wrap), kernel2d(line_kernel, BORDER_TREATMENT_WRAP));
        convolveImage(srcImageRange(src), destImage(dest_clip), kernel2d(line_kernel, BORDER_TREATMENT_CLIP));

        iter_src = src.begin();
        Image::iterator iter_dest_avoid = dest_avoid.begin();
        Image::iterator iter_dest_repeat = dest_repeat.begin();
        Image::iterator iter_dest_reflect = dest_reflect.begin();
        Image::iterator iter_dest_wrap = dest_wrap.begin();
        Image::iterator iter_dest_clip = dest_clip.begin();

        for (int x = 0 ;  x < 40 ; x++){
            if(x > 1 && x < 38 ){
                shouldEqualTolerance(iter_dest_avoid[x], iter_dest_repeat[x], 1e-5);
                shouldEqualTolerance(iter_dest_avoid[x], iter_dest_reflect[x], 1e-5);
                shouldEqualTolerance(iter_dest_avoid[x], iter_dest_wrap[x], 1e-5);
                shouldEqualTolerance(iter_dest_avoid[x], iter_dest_clip[x], 1e-5);
            }else{
                should(iter_dest_avoid[x] == iter_src[x]);
                should(iter_dest_repeat[x] < iter_dest_reflect[x]);
                should(iter_dest_repeat[x] < iter_dest_clip[x]);
                if (x < 1 || x > 38){
                    should(iter_dest_clip[x] < iter_dest_reflect[x]);
                }
                else{
                    should(iter_dest_clip[x] > iter_dest_reflect[x]);
                }
            }
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
        gauss.setBorderTreatment(BORDER_TREATMENT_REFLECT);

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
                    shouldEqualTolerance (acc(id) , acc(is), 1e-15);
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
                    shouldEqualTolerance(acc(is) , acc(id), 1e-15);
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

        Image src(getSymmetricLine());
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
                shouldEqualTolerance (acc(id) , acc(is), 1e-15);
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
                shouldEqualTolerance(acc(is) , acc(id), 1e-15);
            }
            is.x = tmp_src.x;
            id.x = tmp_dest.x;
        }

        // Hier wird an einem symmetrischen Bild /\ (Groesse 40x1)  getestet 
        // ob die korrekten Daten eingehalten wurden.

        Image src(getSymmetricLine());
        dest = src;
        recursiveFilterX(srcImageRange(src), destImage(dest), 
                VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_REFLECT);
        Image::value_type correct_data[40] = 
           {1.10091101909, 1.56303266253, 2.36515826013, 3.29236429975, 4.26558480098, 
            5.25573290944, 6.25210788209, 7.25077235463, 8.25027572885, 9.25007858906, 
            10.2499668097, 11.2498189802, 12.2494745341, 13.2485593473, 14.2460793794, 
            15.2393409851, 16.2210251818, 17.171238053, 18.0359027479, 18.6680232997, 
            18.6680232997, 18.0359027479, 17.171238053, 16.2210251818, 15.2393409851, 
            14.2460793794, 13.2485593473, 12.2494745342, 11.2498189803, 10.24996681, 
            9.25007858994, 8.25027573123, 7.2507723611, 6.2521078997, 5.2557329573, 
            4.26558493107, 3.29236465337, 2.36515922136, 1.56303527544, 1.10091812172};

        Image::iterator dest_iter = dest.begin();
        for (int i = 0 ; i < 40; i++)
        {
            shouldEqualTolerance (correct_data[i], dest_iter[i], 1e-10);
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
                shouldEqualTolerance (acc(id) , acc(is), 1e-15);
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
                shouldEqualTolerance(acc(is) , acc(id), 1e-15);
            }
            is.x = tmp_src.x;
            id.x = tmp_dest.x;
        }


        // Hier wird an einem symmetrischen Bild /\ (Groesse 40x1)  getestet 
        // ob die korrekten Daten eingehalten wurden.

        Image src(getSymmetricLine());
        dest = src;
        recursiveFilterX(srcImageRange(src), destImage(dest), 
                    VIGRA_CSTD::exp(-1.0), BORDER_TREATMENT_WRAP);
        Image::value_type correct_data[40] = 
            {0.8319696, 1.4640946, 2.328761, 3.2789745, 4.260659,
            5.2539208, 6.2514412, 7.2505271, 8.2501855, 9.2500454,
            10.249955, 11.249814, 12.249473, 13.248559, 14.246079,
            15.239341, 16.221025, 17.171238, 18.035903, 18.668023,
            18.668023, 18.035903, 17.171238, 16.221025, 15.239341,
            14.246079, 13.248559, 12.249473, 11.249814, 10.249955,
            9.2500454, 8.2501855, 7.2505271, 6.2514412, 5.2539208,
            4.260659, 3.2789745, 2.328761, 1.4640946, 0.8319696};

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
                shouldEqualTolerance (acc(id) , acc(is), 1e-15);
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
                shouldEqualTolerance(acc(is) , acc(id), 1e-15);
            }
            is.x = tmp_src.x;
            id.x = tmp_dest.x;
        }

        // Hier wird an einem symmetrischen Bild /\ (Groesse 40x1)  getestet 
        // ob die korrekten Daten eingehalten wurden.

        Image src(getSymmetricLine());
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

        Image src_wrap(78, 1);
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
        for (int j = 38 ; j > 0 ; j--, iter_src_wrap++)
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
    
    Image constimg, lenna, rampimg, sym_image, unsym_image;
    vigra::Kernel2D<double> sym_kernel, unsym_kernel, line_kernel;
    
};

struct ResamplingConvolutionTest
{
    void testKernelsSpline()
    {
        BSpline<3, double> spline, dspline(1);
        ArrayVector<Kernel1D<double> > kernels(4);
        Rational<int> samplingRatio(4), offset(1,8);
        resampling_detail::MapTargetToSourceCoordinate 
                   mapCoordinate(samplingRatio, offset);
        createResamplingKernels(spline, mapCoordinate, kernels);
        
        for(unsigned int i = 0; i<kernels.size(); ++i)
        {
            double sum = 0.0;
            for(int k = kernels[i].left(); k <= kernels[i].right(); ++k)
            {
                double x = rational_cast<double>(k + i / samplingRatio + offset);
                shouldEqualTolerance(kernels[i][k], spline(x), 1e-14);
                sum += kernels[i][k];
            }
            shouldEqualTolerance(sum, 1.0, 1e-14);
        }

        createResamplingKernels(dspline, mapCoordinate, kernels);
        
        for(unsigned int i = 0; i<kernels.size(); ++i)
        {
            double sum = 0.0;
            for(int k = kernels[i].left(); k <= kernels[i].right(); ++k)
            {
                double x = rational_cast<double>(k + i / samplingRatio + offset);
                shouldEqualTolerance(kernels[i][k], dspline(x), 1e-14);
                sum += -x*kernels[i][k];
            }
            shouldEqualTolerance(sum, 1.0, 1e-14);
        }
    }
    
    void testKernelsGauss()
    {
        Gaussian<double> gauss(0.7), dgauss(0.7, 1);
        ArrayVector<Kernel1D<double> > kernels(4);
        Rational<int> samplingRatio(4), offset(1,8);
        resampling_detail::MapTargetToSourceCoordinate 
                   mapCoordinate(samplingRatio, offset);
        createResamplingKernels(gauss, mapCoordinate, kernels);
        
        for(unsigned int i = 0; i<kernels.size(); ++i)
        {
            double sum = 0.0;
            for(int k = kernels[i].left(); k <= kernels[i].right(); ++k)
            {
                double x = rational_cast<double>(k + i / samplingRatio + offset);
                shouldEqualTolerance(kernels[i][k], gauss(x), 1e-4);
                sum += kernels[i][k];
            }
            shouldEqualTolerance(sum, 1.0, 1e-14);
       }

        createResamplingKernels(dgauss, mapCoordinate, kernels);
        
        for(unsigned int i = 0; i<kernels.size(); ++i)
        {
            double sum = 0.0;
            double f = i == 0 || i == 3 ? 1.00218 : 0.99851;
            for(int k = kernels[i].left(); k <= kernels[i].right(); ++k)
            {
                double x = rational_cast<double>(k + i / samplingRatio + offset);
                shouldEqualTolerance(kernels[i][k], f*dgauss(x), 1e-5);
                sum += -x*kernels[i][k];
            }
            shouldEqualTolerance(sum, 1.0, 1e-14);
        }
    }
    
    void testOversamplingConstant()
    {
        BSpline<3, double> spline, dspline(1);
        Rational<int> samplingRatio(4,1), offset(1,8);
        
        FImage img(100, 100);
        img.init(1.0);
        
        int wnew = rational_cast<int>((img.width() - 1 - offset) * samplingRatio + 1);
        int hnew = rational_cast<int>((img.height() - 1 - offset) * samplingRatio + 1);
        
        FImage res(wnew, hnew);
        
        resamplingConvolveImage(srcImageRange(img), destImageRange(res),
             spline, samplingRatio, offset, spline, samplingRatio, offset);
        for(FImage::iterator i = res.begin(); i < res.end(); ++i)
            shouldEqual(*i, 1.0);
        
        resamplingConvolveImage(srcImageRange(img), destImageRange(res),
             dspline, samplingRatio, offset, spline, samplingRatio, offset);
        for(FImage::iterator i = res.begin(); i < res.end(); ++i)
            shouldEqual(*i, 0.0);
    }

    void testOversamplingReal()
    {
        Gaussian<double> gauss(0.7);
        Rational<int> samplingRatio(2,1), offset(1,4);
        
        ImageImportInfo info("lenna128.xv");
        FImage img(info.size());
        importImage(info, destImage(img));
        
        int wnew = rational_cast<int>((info.width() - 1 - offset) * samplingRatio + 1);
        int hnew = rational_cast<int>((info.height() - 1 - offset) * samplingRatio + 1);
        
        FImage res(wnew, hnew);        
        resamplingConvolveImage(srcImageRange(img), destImageRange(res),
             gauss, samplingRatio, offset, gauss, samplingRatio, offset);
             
        ImageImportInfo rinfo("resampling.xv");
        shouldEqual(rinfo.width(), wnew);
        shouldEqual(rinfo.height(), hnew);
        FImage ref(wnew, hnew);        
        importImage(rinfo, destImage(ref));

        for(FImage::iterator i = res.begin(), j = ref.begin(); i < res.end(); ++i, ++j)
            shouldEqualTolerance(*i, *j, 1e-6);
    }
};

struct ConvolutionTestSuite
: public vigra::test_suite
{
    ConvolutionTestSuite()
    : vigra::test_suite("ConvolutionTestSuite")
    {
        add( testCase( &ConvolutionTest::simpleSharpeningTest)); 
        add( testCase( &ConvolutionTest::gaussianSharpeningTest)); 
        add( testCase( &ConvolutionTest::stdConvolutionTestOnConstImage));
        add( testCase( &ConvolutionTest::stdConvolutionTestWithAvoid));
        add( testCase( &ConvolutionTest::stdConvolutionTestWithClip));
        add( testCase( &ConvolutionTest::stdConvolutionTestWithWrap));
        add( testCase( &ConvolutionTest::stdConvolutionTestWithReflect));
        add( testCase( &ConvolutionTest::stdConvolutionTestWithRepeat));
        add( testCase( &ConvolutionTest::stdConvolutionTestFromWrapWithReflect));
        add( testCase( &ConvolutionTest::stdConvolutionTestFromRepeatWithAvoid));
        add( testCase( &ConvolutionTest::stdConvolutionTestOfAllTreatmentsRelatively));

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

        add( testCase( &ResamplingConvolutionTest::testKernelsSpline));
        add( testCase( &ResamplingConvolutionTest::testKernelsGauss));
        add( testCase( &ResamplingConvolutionTest::testOversamplingConstant));
        add( testCase( &ResamplingConvolutionTest::testOversamplingReal));
    }
};

int main()
{
    ConvolutionTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}

