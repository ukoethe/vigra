/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
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
#include "vigra/stdimage.hxx"
#include "vigra/stdimagefunctions.hxx"
#include "vigra/impex.hxx"

using namespace vigra; 


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
            vigra::BImage in(info.width(), info.height());
            vigra::BImage out(info.width(), info.height());
           
            importImage(info, destImage(in));
            
            // create image iterator that points to upper left corner 
            // of source image
            vigra::BImage::Iterator sy = in.upperLeft();
            
            // create image iterator that points past the lower right corner of
            // source image (similarly to the past-the-end iterator in the STL)
            vigra::BImage::Iterator send = in.lowerRight();
            
            // create image iterator that points to upper left corner 
            // of destination image
            vigra::BImage::Iterator dy = out.upperLeft();
            
            // iterate down the first column of the images
            for(; sy.y != send.y; ++sy.y, ++dy.y)
            {
                // create image iterator that points to the first 
                // pixel of the current row of the source image
                vigra::BImage::Iterator sx = sy;

                // create image iterator that points to the first 
                // pixel of the current row of the destination image
                vigra::BImage::Iterator dx = dy;
                
                // iterate across current row
                for(; sx.x != send.x; ++sx.x, ++dx.x)
                {
                    // calculate negative gray value
                    *dx = 255 - *sx;
                }
            }
            
            exportImage(srcImageRange(out), vigra::ImageExportInfo(argv[2]));
        }
        else
        {
            vigra::BRGBImage in(info.width(), info.height());
            vigra::BRGBImage out(info.width(), info.height());
           
            importImage(info, destImage(in));
            
            vigra::RGBValue<int> offset(255, 255, 255);
            
            // create image iterator that points to upper left corner 
            // of source image
            vigra::BRGBImage::Iterator sy = in.upperLeft();
            
            // create image iterator that points past the lower right corner of
            // source image (similarly to the past-the-end iterator in the STL)
            vigra::BRGBImage::Iterator send = in.lowerRight();
            
            // create image iterator that points to upper left corner 
            // of destination image
            vigra::BRGBImage::Iterator dy = out.upperLeft();
            
            // iterate down the first column of the images
            for(; sy.y != send.y; ++sy.y, ++dy.y)
            {
                // create image iterator that points to the first 
                // pixel of the current row of the source image
                vigra::BRGBImage::Iterator sx = sy;

                // create image iterator that points to the first 
                // pixel of the current row of the destination image
                vigra::BRGBImage::Iterator dx = dy;
                
                // iterate across current row
                for(; sx.x != send.x; ++sx.x, ++dx.x)
                {
                    // calculate negative color
                    *dx = offset - *sx;
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
