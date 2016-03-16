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
 
#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/convolution.hxx>
#include <vigra/multi_blockwise.hxx>
#include <iostream>

using namespace vigra;

int main (int argc, char ** argv) 
{
    if(argc != 3) 
    {
        std::cout << "Usage: " << argv[0] << " infile outfile" << std::endl;
        std::cout << "(supported formats: " << impexListFormats() << ")" << std::endl;
        
        return 1;
    }
    try
    {
        // read image given as first argument
        ImageImportInfo info(argv[1]);
        
        // instantiate arrays for image data and for smoothed image of appropriate size
        if (info.isGrayscale()) 
        {
            MultiArray<2, float> imageArray(info.shape()),
                                            exportArray(info.shape());

            // copy image data into array
            importImage(info, imageArray);
            
            BlockwiseConvolutionOptions<2> opt = BlockwiseConvolutionOptions<2>();
            
            gaussianSmoothMultiArray(imageArray, exportArray, 
                                     BlockwiseConvolutionOptions<2>().stdDev(2.0));
            
            // write image data to the file given as second argument
            exportImage(exportArray, ImageExportInfo(argv[2]));
            
        }
        else
        {
            MultiArray<2, RGBValue<float> > imageArray(info.shape()),
                                            exportArray(info.shape());
        
            // copy image data into array
            importImage(info, imageArray);
            
            BlockwiseConvolutionOptions<2> opt = BlockwiseConvolutionOptions<2>();
            
            gaussianSmoothMultiArray(imageArray, exportArray, 
                                     BlockwiseConvolutionOptions<2>().stdDev(2.0));
            
            // write image data to the file given as second argument
            exportImage(exportArray, ImageExportInfo(argv[2]));
        }
        
        
    }
    catch (std::exception & e) 
    {
        // catch any errors that might have occurred and print their reason
        std::cout << e.what() << std::endl;
        return 1;
    }
    return 0;
}