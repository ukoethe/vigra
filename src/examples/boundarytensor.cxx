/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2004 by Ullrich Koethe                  */
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
