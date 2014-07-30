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
#include <vigra/convolution.hxx>
#include <vigra/nonlineardiffusion.hxx>
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
    
    // Type of smoothing: 
    int type;
    std::cout << "Type of smoothing (1 = Gauss, 2 = Exponential, 3 = nonlinear) ? ";
    std::cin >> type;
    
    // input width of smoothing filter 
    double scale;
    std::cout << "Amount of smoothing (operator scale) ? ";
    std::cin >> scale;
    
    double edge_threshold;
    if(type == 3)
    {
        std::cout << "Edge threshold ? ";
        std::cin >> edge_threshold;
    }
    
    try
    {
        ImageImportInfo info(argv[1]);
        
        if(info.isGrayscale())
        {
            MultiArray<2, UInt8> in(info.width(), info.height());
            MultiArray<2, float> out(info.width(), info.height());
           
            importImage(info, destImage(in));
            
            switch(type)
            {
              case 2:
              {
                // apply recursive filter (exponential filter) to gray image
                recursiveSmoothX(in, out, scale);
                recursiveSmoothY(out, out, scale);
                break;
              }
              case 3:
              {
                // apply nonlinear diffusion to gray image
                nonlinearDiffusion(in, out,
                   DiffusivityFunctor<float>(edge_threshold), scale);
                break;
              }
              default:
              {
                gaussianSmoothing(in, out, scale);
              }
            }
            
            exportImage(out, ImageExportInfo(argv[2]));
        }
        else
        {
            MultiArray<2, RGBValue<UInt8> > in(info.shape());
            MultiArray<2, RGBValue<float> > out(info.shape());
           
            importImage(info, in);
            
            switch(type)
            {
              case 2:
              {
                // apply recursive filter (exponential filter) to color image
                recursiveSmoothX(in, out, scale);
                recursiveSmoothY(out, out, scale);
                break;
              }
              case 3:
              {
                // apply nonlinear diffusion to color image, one band at a time
                for(int band = 0; band<3; ++band)
                {
                    nonlinearDiffusion(in.bindElementChannel(band), out.bindElementChannel(band),
                                       DiffusivityFunctor<float>(edge_threshold), scale);
                }
                break;
              }
              default:
              {
                gaussianSmoothing(in, out, scale);
              }
            }
            
            exportImage(out, ImageExportInfo(argv[2]));
        }
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
