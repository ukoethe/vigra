/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2004 by Ullrich Koethe                  */
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
#include <functional>
#include "vigra/stdimage.hxx"
#include "vigra/combineimages.hxx"
#include "vigra/tensorutilities.hxx"
#include "vigra/boundarytensor.hxx"
#include "vigra/impex.hxx"

using namespace vigra; 


int main(int argc, char ** argv)
{
    if(argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " infile" << std::endl;
        std::cout << "(supported formats: " << impexListFormats() << ")" << std::endl;
        std::cout << "creates: boundarystrength.tif, cornerstrength.tif" << std::endl;
        
        return 1;
    }
    
    try
    {
        ImageImportInfo info(argv[1]);
        int w = info.width(), h = info.height();
        // create image of appropriate size for boundary tensor
        FVector3Image boundarytensor(w, h);
        
        // input scale of the bandpass to be used
        double scale;
        std::cout << "Operator scale ? ";
        std::cin >> scale;
        
        if(info.isGrayscale())
        {
            FImage in(w, h);
            importImage(info, destImage(in));

            boundaryTensor(srcImageRange(in), destImage(boundarytensor), scale);
        }
        else if(info.isColor())
        {
            FRGBImage in(w, h);
            importImage(info, destImage(in));
            
            // calculate the boundary tensor for every channel and add the results
            FVector3Image bandtensor(w, h);
            for(int b=0; b<3; ++b)
            {
                VectorElementAccessor<FRGBImage::Accessor> band(b);
                boundaryTensor(srcImageRange(in, band), destImage(bandtensor), scale);
                combineTwoImages(srcImageRange(boundarytensor), srcImage(bandtensor),
                                 destImage(boundarytensor), 
                                 std::plus<FVector3Image::value_type>());
            }
        }
        else
        {
            std::cerr << "Sorry, can only operate on gray and color images.\n";
            return 1;
        }
        
        FImage boundarystrength(w, h), cornerness(w, h);
        FVector2Image edgeness(w,h);
        
        tensorTrace(srcImageRange(boundarytensor), destImage(boundarystrength));
        tensorToEdgeCorner(srcImageRange(boundarytensor), destImage(edgeness), destImage(cornerness));

        exportImage(srcImageRange(boundarystrength), ImageExportInfo("boundarystrength.tif").setPixelType("UINT8"));
        exportImage(srcImageRange(cornerness), ImageExportInfo("cornerstrength.tif").setPixelType("UINT8"));
    }
    catch (std::exception & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
