/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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
#include "vigra/stdimage.hxx"
#include "vigra/stdimagefunctions.hxx"
#include "vigra/localminmax.hxx"
#include "vigra/labelimage.hxx"
#include "vigra/seededregiongrowing.hxx"
#include "vigra/impex.hxx"

using namespace vigra; // MSVC doesn't support Koenig lookup

// define a functor that calsulates the squared magnitude of the gradient
// given the x- and y- components of the gradient
struct GradientSquaredMagnitudeFunctor
{
    float operator()(float const & g1, float const & g2) const 
    {
        return g1 * g1 + g2 * g2;
    }
    
    float operator()(vigra::RGBValue<float> const & rg1, vigra::RGBValue<float> const & rg2) const 
    {
        float g1 = rg1.squaredMagnitude();
        float g2 = rg2.squaredMagnitude();
        
        return g1 * g1 + g2 * g2;
    }
};

// generic implementation of the watershed algorithm
template <class InImage, class OutImage>
void watershedSegmentation(InImage & in, OutImage & out, double scale)
{
    typedef typename vigra::NumericTraits<typename InImage::value_type>::RealPromote
        TmpType;

    int w = in.width();
    int h = in.height();
            
    vigra::BasicImage<TmpType> gradientx(w, h);
    vigra::BasicImage<TmpType> gradienty(w, h);
    vigra::FImage gradientmag(w, h);

    // calculate the x- and y-components of the image gradient at given scale
    recursiveFirstDerivativeX(srcImageRange(in), destImage(gradientx), scale);
    recursiveSmoothY(srcImageRange(gradientx), destImage(gradientx), scale);

    recursiveFirstDerivativeY(srcImageRange(in), destImage(gradienty), scale);
    recursiveSmoothX(srcImageRange(gradienty), destImage(gradienty), scale);

    // transform components into gradient magnitude
    combineTwoImages(srcImageRange(gradientx), srcImage(gradienty), 
                     destImage(gradientmag), GradientSquaredMagnitudeFunctor());

    vigra::IImage labels(w, h);
    labels = 0;

    // find the local minima of the gradient magnitude
    // (might be larger than one pixel)
    extendedLocalMinima(srcImageRange(gradientmag), destImage(labels), 1); 

    // label the minima just found
    int max_region_label = 
        labelImageWithBackground(srcImageRange(labels), destImage(labels),
                                 false, 0);

    // create a statistics functor for region growing
    vigra::ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<float> > 
                                          gradstat(max_region_label);

    // perform region growing, starting from the minima of the gradient magnitude;
    // as the feature (first input) image contains the gradient magnitude,
    // this calculates the catchment basin of each minimum
    seededRegionGrowing(srcImageRange(gradientmag), srcImage(labels), 
                        destImage(labels), gradstat);

    // initialize a functor to determin the average gray-value or color
    // for each region (catchment basin) just found
    vigra::ArrayOfRegionStatistics<vigra::FindAverage<TmpType> > 
                                          averages(max_region_label);

    // calculate the averages
    inspectTwoImages(srcImageRange(in), srcImage(labels), averages);

    // write the averages into the destination image (the functor 'averages'
    // acts as a look-up table)
    transformImage(srcImageRange(labels), destImage(out), averages);
    
    // mark the watershaed (region boundaries) black
    regionImageToEdgeImage(srcImageRange(labels), destImage(out), 
                           vigra::NumericTraits<typename OutImage::value_type>::zero());

}


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
        
        // input width of gradient filter 
        double scale = 1.0;
        std::cout << "Scale for gradient calculation ? ";
        std::cin >> scale;
            
        if(info.isGrayscale())
        {
            int w = info.width();
            int h = info.height();
            
            vigra::BImage in(w, h);
            importImage(info, destImage(in));
            
            vigra::BImage out(w, h);

            // perform watershed segmentation on gray image
            // note that the watershed algorithm usually results in an 
            // oversegmentation (i.e., to many regions), but its boundary 
            // localization is quite good
            watershedSegmentation(in, out, scale);
           
            std::cout << "Writing " << argv[2] << std::endl;
            exportImage(srcImageRange(out), vigra::ImageExportInfo(argv[2]));
           
        }
        else
        {
            int w = info.width();
            int h = info.height();
            
            vigra::BRGBImage in(w, h);
            importImage(info, destImage(in));
            
            vigra::BRGBImage out(w, h);


            // perform watershed segmentation on color image
            // note that the watershed algorithm usually results in an 
            // oversegmentation (i.e., to many regions), but its boundary 
            // localization is quite good
            watershedSegmentation(in, out, scale);
           
            std::cout << "Writing " << argv[2] << std::endl;
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
