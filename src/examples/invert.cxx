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
#include <vigra/stdimagefunctions.hxx>
#include <vigra/multi_math.hxx>
#include <vigra/impex.hxx>
#include <string.h>

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
        
        if(info.isGrayscale())
        {
            MultiArray<2, UInt8> in(info.width(), info.height()),
                                 out(info.width(), info.height());
           
            importImage(info, in);
            
            // create an inverted image by applying the expression
            //       newvalue = -1 * (oldvalue - 255)
            // to each pixel
            transformImage(in, out,
                           linearIntensityTransform(-1, -255));
            
            // the same can be achieved using array expressions
            using namespace multi_math;   // activate array expressions
            out = 255 - in;
            
            if(strcmp(argv[2], "-") == 0)
            {
                // write stdout
                exportImage(out, ImageExportInfo(argv[2]).setFileType(info.getFileType()));
            }
            else
            {
                exportImage(out, ImageExportInfo(argv[2]));
            }
        }
        else
        {
            MultiArray<2, RGBValue<UInt8> > in(info.width(), info.height()),
                                            out(info.width(), info.height());
           
            importImage(info, in);
            
            // create a negative image by applying the expression
            //       newvalue = -1 * (oldvalue + RGBValue<int>(-255, -255, -255))
            // to each pixel
            transformImage(in, out,
                           linearIntensityTransform(-1, RGBValue<int>(-255)));
            
            // the same can be achieved using array expressions
            using namespace multi_math;   // activate array expressions
            out = RGBValue<int>(255) - in;
            
            if(strcmp(argv[2], "-") == 0)
            {
                // write stdout
                exportImage(out, 
                 ImageExportInfo(argv[2]).setFileType(info.getFileType()));
            }
            else
            {
                exportImage(out, ImageExportInfo(argv[2]));
            }
        }
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
