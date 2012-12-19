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
#include "vigra/multi_array.hxx"
#include "vigra/rgbvalue.hxx"
#include "vigra/resizeimage.hxx"
#include "vigra/impex.hxx"


template<class ImageType>
bool resizeImageFile(const vigra::ImageImportInfo &info, const vigra::Shape2 &newSize,
                     int method, const char *outputFilename)
{
    // create a gray scale image of appropriate size
    ImageType in(info.shape());
    ImageType out(newSize);

    // import the image just read
    importImage(info, destImage(in));

    using vigra::BSpline;

    switch(method)
    {
      case 0:
        // equiv. to resizeImageSplineInterpolation with BSpline<0, double>:
        resizeImageNoInterpolation(srcImageRange(in), destImageRange(out));
        break;
      case 1:
        // equiv. to resizeImageSplineInterpolation with BSpline<1, double>:
        resizeImageLinearInterpolation(srcImageRange(in), destImageRange(out));
        break;
      case 2:
        resizeImageSplineInterpolation(srcImageRange(in), destImageRange(out),
                                       BSpline<2, double>());
        break;
      case 3:
        resizeImageSplineInterpolation(srcImageRange(in), destImageRange(out),
                                       BSpline<3, double>());
        break;
      case 4:
        resizeImageSplineInterpolation(srcImageRange(in), destImageRange(out),
                                       BSpline<4, double>());
        break;
      case 5:
        resizeImageSplineInterpolation(srcImageRange(in), destImageRange(out),
                                       BSpline<5, double>());
        break;
      case 6:
        resizeImageSplineInterpolation(srcImageRange(in), destImageRange(out),
                                       BSpline<6, double>());
        break;
      case 7:
        resizeImageSplineInterpolation(srcImageRange(in), destImageRange(out),
                                       BSpline<7, double>());
        break;
      default:
        std::cerr << "Invalid method " << method << " (must be 0..7)!\n";
        return false;
    }

    // write the image to the file given as second argument
    // the file type will be determined from the file name's extension
    exportImage(srcImageRange(out), vigra::ImageExportInfo(outputFilename));
    return true;
}


int main(int argc, char ** argv)
{
    using vigra::Shape2;

    if((argc < 3) || (argc > 5))
    {
        std::cout << "Usage: " << argv[0] << " infile outfile [factor] [method]" << std::endl;
        std::cout << "(supported formats: " << vigra::impexListFormats() << ")" << std::endl;
        std::cout << "If factor or method are not provided, you will be asked for\nthem on the command line." << std::endl;

        return 1;
    }

    try
    {
        // read image given as first argument
        // file type is determined automatically
        vigra::ImageImportInfo info(argv[1]);

        double sizefactor;
        if(argc > 3)
        {
            sizefactor = atof(argv[3]);
        }
        else
        {
            std::cerr << "Resize factor ? ";
            std::cin >> sizefactor;
        }

        int method;
        if(argc > 4)
        {
            method = atoi(argv[4]);
        }
        else
        {
            std::cerr << "Method (0: pixel repetition, 1: linear, 2-7: spline) ? ";
            std::cin >> method;
        }

        // calculate new image size
        Shape2 newSize((info.shape() - Shape2(1,1)) * sizefactor + Shape2(1,1));

        if(info.isGrayscale())
        {
            if(!resizeImageFile<vigra::MultiArray<2, unsigned char> >(info, newSize, method, argv[2]))
                return 1;
        }
        else
        {
            if(!resizeImageFile<vigra::MultiArray<2, vigra::RGBValue<unsigned char> > >(info, newSize, method, argv[2]))
                return 1;
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
