/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/
 

#include <iostream>
#include <stdio.h>
#include "vigra/stdimage.hxx"
#include "vigra/convolution.hxx"
#include "vigra/resizeimage.hxx"
#include "vigra/impex.hxx"

using namespace vigra; // MSVC doesn't support Koenig lookup

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
