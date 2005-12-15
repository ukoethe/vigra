/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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
#include <stdio.h>
#include "vigra/stdimage.hxx"
#include "vigra/convolution.hxx"
#include "vigra/resizeimage.hxx"
#include "vigra/impex.hxx"

using namespace vigra; 

// Gaussian reduction to next pyramid level
template <class Image>
void reduceToNextLevel(Image & in, Image & out)
{    
    // image size at current level
    int width = in.width();
    int height = in.height();
    
    // image size at next smaller level
    int newwidth = (width + 1) / 2;
    int newheight = (height + 1) / 2;
    
    // resize result image to appropriate size
    out.resize(newwidth, newheight);
    
    // define a Gaussian kernel (size 5x1)
    vigra::Kernel1D<double> filter;
    filter.initExplicitly(-2, 2) = 0.05, 0.25, 0.4, 0.25, 0.05;
    
    vigra::BasicImage<typename Image::value_type> tmpimage1(width, height);
    vigra::BasicImage<typename Image::value_type> tmpimage2(width, height);
    
    // smooth (band limit) input image
    separableConvolveX(srcImageRange(in),
                       destImage(tmpimage1), kernel1d(filter));
    separableConvolveY(srcImageRange(tmpimage1),
                       destImage(tmpimage2), kernel1d(filter));
                       
    // downsample smoothed image
    resizeImageNoInterpolation(srcImageRange(tmpimage2), destImageRange(out));
    
}

int main(int argc, char ** argv)
{
    if(argc != 3)
    {
        std::cout << "Usage: " << argv[0] << " infile outfile" << std::endl;
        std::cout << "(supported formats: " << vigra::impexListFormats() << ")" << std::endl;
        
        return 1;
    }
    
    try
    {
        vigra::ImageImportInfo info(argv[1]);
        
        if(info.isGrayscale())
        {
            vigra::BImage levels[4];
        
            levels[0].resize(info.width(), info.height());
           
            importImage(info, destImage(levels[0]));
            
            for(int i=1; i<4; ++i)
            {
                // reduce gray image 3 times
                reduceToNextLevel(levels[i-1], levels[i]);

            }
            
            exportImage(srcImageRange(levels[3]), vigra::ImageExportInfo(argv[2]));
        }
        else
        {
            vigra::BRGBImage levels[4];
        
            levels[0].resize(info.width(), info.height());
           
            importImage(info, destImage(levels[0]));
            
            for(int i=1; i<4; ++i)
            {
                // reduce color image 3 times
                reduceToNextLevel(levels[i-1], levels[i]);
            }
            
            exportImage(srcImageRange(levels[3]), vigra::ImageExportInfo(argv[2]));
        }
    }
    catch (vigra::StdException & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
