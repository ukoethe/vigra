/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2000 by Ullrich Koethe                  */
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
#include "vigra/edgedetection.hxx"
#include "vigra/impex.hxx"


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
        
        precondition(info.isGrayscale(), "Sorry, cannot operate on color images");
        
        vigra::BImage in(info.width(), info.height());

        importImage(info, destImage(in));

        // input width of edge detection filter
        double scale;
        std::cout << "Operator scale ? ";
        std::cin >> scale;

        // input threshold for gradient magnitude
        double threshold;
        std::cout << "Gradient threshold ? ";
        std::cin >> threshold;
    
        // create output image of appropriate size
        vigra::BImage out(info.width(), info.height());
        
        // paint output image white
        out = 255;
        
        // call edge detection algorithm
        // edges will be marked black
        differenceOfExponentialEdgeImage(srcImageRange(in), destImage(out),
                       scale, threshold, 0);

        exportImage(srcImageRange(out), vigra::ImageExportInfo(argv[2]));
    }
    catch (vigra::StdException & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
