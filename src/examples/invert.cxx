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
#include "vigra/stdimagefunctions.hxx"
#include "vigra/impex.hxx"
#include <string.h>


int main(int argc, char ** argv)
{
    if(argc != 3)
    {
        std::cout << "Usage: " << argv[0] << " infile outfile" << std::endl;
        std::cout << "(supported fomats: " << vigraImpexListFormats() << ")" << std::endl;
        
        return 1;
    }
    
    try
    {
        ImageImportInfo info(argv[1]);
        
        if(info.isGrayscale())
        {
            BImage in(info.width(), info.height());
            BImage out(info.width(), info.height());
           
            importImage(info, destImage(in));
            
            // create a negative image by applying the expression
            //       newvalue = -1 * (oldvalue - 255)
            // to each pixel
            transformImage(srcImageRange(in), destImage(out),
                           linearIntensityTransform(-1, -255));
            
            if(strcmp(argv[2], "-") == 0)
            {
                // write stdout
                exportImage(srcImageRange(out), 
                 ImageExportInfo(argv[2]).setFileType(info.getFileType()));
            }
            else
            {
                exportImage(srcImageRange(out), ImageExportInfo(argv[2]));
            }
        }
        else
        {
            BRGBImage in(info.width(), info.height());
            BRGBImage out(info.width(), info.height());
           
            importImage(info, destImage(in));
            
            RGBValue<int> offset(-255, -255, -255);
            
            
            // create a negative image by applying the expression
            //       newvalue = -1 * (oldvalue + RGBValue<int>(-255, -255, -255))
            // to each pixel
            transformImage(srcImageRange(in), destImage(out),
                           linearIntensityTransform(-1, offset));
            
            if(strcmp(argv[2], "-") == 0)
            {
                // write stdout
                exportImage(srcImageRange(out), 
                 ImageExportInfo(argv[2]).setFileType(info.getFileType()));
            }
            else
            {
                exportImage(srcImageRange(out), ImageExportInfo(argv[2]));
            }
        }
    }
    catch (VigraStdException & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
