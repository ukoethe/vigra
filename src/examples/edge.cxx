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
#include <vigra/multi_array.hxx>
#include <vigra/edgedetection.hxx>
#include <vigra/impex.hxx>

using namespace vigra; 


int main(int argc, char ** argv)
{
    if(argc != 3)
    {
        std::cout << "Usage: " << argv[0] << " infile outfile" << std::endl;
        std::cout << "(supported formats: " << impexListFormats() << ")" << std::endl;
        
        return 1;
    }
    
    try
    {
        ImageImportInfo info(argv[1]);
        
        vigra_precondition(info.isGrayscale(), "Sorry, cannot operate on color images");
        
        MultiArray<2, UInt8> in(info.shape());

        importImage(info, in);

        // input width of edge detection filter
        int which;
        std::cout << "Use Canny or Shen-Castan detector (1 or 2) ? ";
        std::cin >> which;

        // input width of edge detection filter
        double scale;
        std::cout << "Operator scale ? ";
        std::cin >> scale;

        // input threshold for gradient magnitude
        double threshold;
        std::cout << "Gradient threshold ? ";
        std::cin >> threshold;
    
        // create output image of appropriate size
        MultiArray<2, UInt8> out(info.shape());
        
        // paint output image white
        out = 255;
        
        if(which == 2)
        {
            // call Shen-Castan edge detection algorithm
            // edges will be marked black
            differenceOfExponentialEdgeImage(in, out, scale, threshold, 0);
        }
        else
        {
            // call Canny edge detection algorithm
            // edges will be marked black
            cannyEdgeImage(in, out, scale, threshold, 0);
        }
        
        exportImage(out, ImageExportInfo(argv[2]));
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
