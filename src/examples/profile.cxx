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
#include "vigra/imageiteratoradapter.hxx"
#include "vigra/impex.hxx"

using namespace vigra; // MSVC doesn't support Koenig lookup


int main(int argc, char ** argv)
{
    if(argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " infile" << std::endl;
        std::cout << "(supported formats: " << vigra::impexListFormats() << ")" << std::endl;
        
        return 1;
    }
    
    try
    {
        vigra::ImageImportInfo info(argv[1]);
        
        vigra_precondition(info.isGrayscale(), "Sorry, cannot operate on color images");
        
        int w = info.width();
        int h = info.height();
            
        
        vigra::BImage in(w, h);
        importImage(info, destImage(in));
        
        int length = (w < h) ? h : w;

        // create output image of appropriate size
        vigra::BImage out(length, 256);
        
        
        // paint output image white
        out = 255;

        // create line iterator that iterates along the image diagonal
        vigra::LineIterator<vigra::BImage::Iterator> line(in.upperLeft(), in.lowerRight());
         
        // create line iterator that marks the end of iteration
        vigra::LineIterator<vigra::BImage::Iterator> end(in.lowerRight(), in.lowerRight());

        // create image iterator that points to the first pixel of the last
        // row of the destination image
        vigra::BImage::Iterator column = out.upperLeft() + vigra::Diff2D(0, 255);
        
        // iterate along the line and across the destination image
        for(; line != end; ++line, ++column.x)
        {
            vigra::BImage::Iterator row(column);
            // paint all pixels black whose coordinates are smaller than the
            // current gray value along the diagonal
            for(int y=0; y <= *line; ++y, --row.y)  *row = 0;
        }
        
        std::cout << "Writing profile.gif" << std::endl;
        exportImage(srcImageRange(out), vigra::ImageExportInfo("profile.gif"));
    }
    catch (vigra::StdException & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
