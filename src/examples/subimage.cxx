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
#include "vigra/stdimagefunctions.hxx"
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
    
    try
    {
        vigra::ImageImportInfo info(argv[1]);
        
        // define upper left and lower right corners of a 
        // subimage (region of interest)
        int sub_x0 = info.width() / 4;
        int sub_y0 = info.height() / 4;
        int sub_x1 = info.width() - sub_x0;
        int sub_y1 = info.height() - sub_y0;
            
        if(info.isGrayscale())
        {
            vigra::BImage in(info.width(), info.height());
            
            importImage(info, destImage(in));
            
            // create output image of appropriate size
            vigra::BImage out(sub_x1 - sub_x0, sub_y1 - sub_y0);
            
            // copy region of interest by moving the input 
            // iterators to the appropriate positions
            copyImage(srcIterRange(in.upperLeft() + vigra::Diff2D(sub_x0, sub_y0),
                                   in.upperLeft() + vigra::Diff2D(sub_x1, sub_y1)),
                      destImage(out));
            
            exportImage(srcImageRange(out), vigra::ImageExportInfo(argv[2]));
        }
        else
        {
            vigra::BRGBImage in(info.width(), info.height());
            
            importImage(info, destImage(in));
            
            // create output image of appropriate size
            vigra::BRGBImage out(sub_x1 - sub_x0, sub_y1 - sub_y0);
            
            // copy region of interest by moving the input 
            // iterators to the appropriate positions
            copyImage(srcIterRange(in.upperLeft() + vigra::Diff2D(sub_x0, sub_y0),
                                   in.upperLeft() + vigra::Diff2D(sub_x1, sub_y1)),
                      destImage(out));
            
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
