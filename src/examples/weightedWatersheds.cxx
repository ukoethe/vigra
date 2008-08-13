/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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
#include "vigra/stdimagefunctions.hxx"
#include "vigra/functorexpression.hxx"
#include "vigra/resizeimage.hxx"
#include "vigra/convolution.hxx"
#include "vigra/localminmax.hxx"
#include "vigra/labelimage.hxx"
#include "vigra/seededregiongrowing.hxx"
#include "vigra/impex.hxx"

using namespace vigra; 

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

        return g1 + g2;
    }
};

// generic implementation of the watershed algorithm
template <class InImage, class OutImage>
void weightedWatershedSegmentation(InImage & in, OutImage & out, double scale, unsigned int oversampling = 1)
{
    using namespace vigra::functor;
    
    typedef typename vigra::NumericTraits<typename InImage::value_type>::RealPromote
        TmpType;

    vigra_precondition(oversampling > 0,
       "weightedWatershedSegmentation(): oversampling must not be zero.");
    
    int w = oversampling*(in.width() - 1) + 1;
    int h = oversampling*(in.height() - 1) + 1;

    vigra::BasicImage<TmpType> gradientx(w, h);
    vigra::BasicImage<TmpType> gradienty(w, h);
    vigra::FImage gradientmag(w, h);

    // calculate the x- and y-components of the image gradient at given scale
    // optionally enlarge the image before
    if(oversampling > 1)
    {
        vigra::BasicImage<TmpType> big(w, h);
        resizeImageSplineInterpolation(srcImageRange(in), destImageRange(big));
        gaussianGradient(srcImageRange(big), destImage(gradientx), destImage(gradienty), scale*oversampling);
    }
    else
    {
        gaussianGradient(srcImageRange(in), destImage(gradientx), destImage(gradienty), scale);
    }
    
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
                        destImage(labels), gradstat, KeepContours);

    out.resize(w,h);
    
    // set boundary pixels to the corresponding gradient value, non-boundary pixels to zero
    combineTwoImages(srcImageRange(gradientmag), srcImage(labels),
                     destImage(out), 
                     ifThenElse(Arg2()==Param(0), sqrt(Arg1()), Param(0.0)));
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
        unsigned int oversampling = 1;
        std::cout << "Scale for gradient calculation ? ";
        std::cin >> scale;
        std::cout << "Oversampling ? ";
        std::cin >> oversampling;

        if(info.isGrayscale())
        {
            int w = info.width();
            int h = info.height();

            vigra::BImage in(w, h);
            importImage(info, destImage(in));

            vigra::FImage out;

            // perform watershed segmentation on gray image
            // note that the watershed algorithm usually results in an
            // oversegmentation (i.e., too many regions), but its boundary
            // localization is quite good
            weightedWatershedSegmentation(in, out, scale, oversampling);

            std::cout << "Writing " << argv[2] << std::endl;
            exportImage(srcImageRange(out), vigra::ImageExportInfo(argv[2]));
        }
        else
        {
            int w = info.width();
            int h = info.height();

            vigra::BRGBImage in(w, h);
            importImage(info, destImage(in));

            vigra::FImage out;


            // perform watershed segmentation on color image
            // note that the watershed algorithm usually results in an
            // oversegmentation (i.e., too many regions), but its boundary
            // localization is quite good
            weightedWatershedSegmentation(in, out, scale, oversampling);

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
