/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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
#include "vigra/stdimage.hxx"
#include "vigra/convolution.hxx"
#include "vigra/nonlineardiffusion.hxx"
#include "vigra/impex.hxx"

using namespace vigra; // MSVC doesn't support Koenig lookup


int main(int argc, char ** argv)
{
    if(argc != 3)
    {
        std::cout << "Usage: " << argv[0] << " infile outfile" << std::endl;
        std::cout << "(supported formats: " << vigra::impexListFormats() << ")" << std::endl;
        
        return 1;
    }
    
    // Type of smoothing: 
    int type;
    std::cout << "Type of smoothing (1 = Gauss, 2 = Exponential, 3 = nolinear) ? ";
    std::cin >> type;
    
    // input width of smoothing filter 
    double scale;
    std::cout << "Amount of smoothing (operator scale) ? ";
    std::cin >> scale;
    
    double edge_threshold;
    if(type == 3)
    {
        std::cout << "Edge threshold ? ";
        std::cin >> edge_threshold;
    }
    
    try
    {
        vigra::ImageImportInfo info(argv[1]);
        
        if(info.isGrayscale())
        {
            vigra::BImage in(info.width(), info.height());
            vigra::BImage out(info.width(), info.height());
           
            importImage(info, destImage(in));
            
            switch(type)
            {
              case 2:
              {
                // apply recursive filter (exponential filter) to gray image
                recursiveSmoothX(srcImageRange(in), destImage(out), scale);
                recursiveSmoothY(srcImageRange(out), destImage(out), scale);
                break;
              }
              case 3:
              {
                // apply nonlinear diffusion to gray image
                nonlinearDiffusion(srcImageRange(in), destImage(out),
                   vigra::DiffusivityFunctor<float>(edge_threshold), scale);
                break;
              }
              default:
              {
                vigra::FImage tmp(info.width(), info.height());

                // apply Gaussian filter to gray image
                vigra::Kernel1D<double> gauss;
                gauss.initGaussian(scale);
                separableConvolveX(srcImageRange(in), destImage(tmp), kernel1d(gauss));
                separableConvolveY(srcImageRange(tmp), destImage(out), kernel1d(gauss));
              }
            }
            
            exportImage(srcImageRange(out), vigra::ImageExportInfo(argv[2]));
        }
        else
        {
            vigra::BRGBImage in(info.width(), info.height());
            vigra::BRGBImage out(info.width(), info.height());
           
            importImage(info, destImage(in));
            
            switch(type)
            {
              case 2:
              {
                // apply recursive filter (exponential filter) to color image
                recursiveSmoothX(srcImageRange(in), destImage(out), scale);
                recursiveSmoothY(srcImageRange(out), destImage(out), scale);
                break;
              }
              case 3:
              {
                // apply nonlinear diffusion to color image
                nonlinearDiffusion(srcImageRange(in), destImage(out),
                   vigra::DiffusivityFunctor<float>(edge_threshold), scale);
                break;
              }
              default:
              {
                vigra::FRGBImage tmp(info.width(), info.height());

                // apply Gaussian filter to color image
                vigra::Kernel1D<double> gauss;
                gauss.initGaussian(scale);
                separableConvolveX(srcImageRange(in), destImage(tmp), kernel1d(gauss));
                separableConvolveY(srcImageRange(tmp), destImage(out), kernel1d(gauss));
              }
            }
            
            exportImage(srcImageRange(out), vigra::ImageExportInfo(argv[2]));
        }
    }
    catch (vigra::StdException & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
