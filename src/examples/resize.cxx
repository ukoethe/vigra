/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
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
#include "vigra/resizeimage.hxx"
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
        // read image given as first argument
        // file type is determined automatically
        vigra::ImageImportInfo info(argv[1]);

        double sizefactor;
        std::cerr << "Resize factor ? ";
        std::cin >> sizefactor;
        int method;
        std::cerr << "Method (0 - pixel repetition, 1 - linear, 2 - spline ? ";
        std::cin >> method;

        // calculate new image size
        int nw = (int)(sizefactor*(info.width()-1) + 1.5);
        int nh = (int)(sizefactor*(info.height()-1) + 1.5);

        if(info.isGrayscale())
        {
            // create a gray scale image of appropriate size
            vigra::BImage in(info.width(), info.height());
            vigra::BImage out(nw, nh);

            // import the image just read
            importImage(info, destImage(in));

            switch(method)
            {
              case 0:
                // resize the image, using a bi-cubic spline algorithms
                resizeImageNoInterpolation(srcImageRange(in),
                    destImageRange(out));
                break;
              case 1:
                // resize the image, using a bi-cubic spline algorithms
                resizeImageLinearInterpolation(srcImageRange(in),
                    destImageRange(out));
                break;
              default:
                // resize the image, using a bi-cubic spline algorithms
                resizeImageSplineInterpolation(srcImageRange(in),
                    destImageRange(out));
            }

            // write the image to the file given as second argument
            // the file type will be determined from the file name's extension
            exportImage(srcImageRange(out), vigra::ImageExportInfo(argv[2]));
        }
        else
        {
            // create a RGB image of appropriate size
            vigra::BRGBImage in(info.width(), info.height());
            vigra::BRGBImage out(nw, nh);

            // import the image just read
            importImage(info, destImage(in));

            switch(method)
            {
              case 0:
                // resize the image, using a bi-cubic spline algorithms
                resizeImageNoInterpolation(srcImageRange(in),
                    destImageRange(out));
                break;
              case 1:
                // resize the image, using a bi-cubic spline algorithms
                resizeImageLinearInterpolation(srcImageRange(in),
                    destImageRange(out));
                break;
              default:
                // resize the image, using a bi-cubic spline algorithms
                resizeImageSplineInterpolation(srcImageRange(in),
                    destImageRange(out));
            }

            // write the image to the file given as second argument
            // the file type will be determined from the file name's extension
            exportImage(srcImageRange(out), vigra::ImageExportInfo(argv[2]));
        }
    }
    catch (vigra::StdException & e)
    {
        // catch any errors that might have occurred and print their reason
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}
