/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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
#include "vigra/edgedetection.hxx"
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
        
        vigra_precondition(info.isGrayscale(), "Sorry, cannot operate on color images");
        
        vigra::BImage in(info.width(), info.height());

        importImage(info, destImage(in));

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
        vigra::BImage out(info.width(), info.height());
        
        // paint output image white
        out = 255;
        
        if(which == 2)
        {
            // call edge detection algorithm
            // edges will be marked black
            differenceOfExponentialEdgeImage(srcImageRange(in), destImage(out),
                           scale, threshold, 0);
        }
        else
        {
            // call edge detection algorithm
            // edges will be marked black
            cannyEdgeImage(srcImageRange(in), destImage(out),
                           scale, threshold, 0);
        }
        
        exportImage(srcImageRange(out), vigra::ImageExportInfo(argv[2]));
    }
    catch (vigra::StdException & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
