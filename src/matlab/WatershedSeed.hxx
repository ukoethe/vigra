#include <iostream>
#include <string>
#include <vigra/watersheds.hxx>
#include <vigra/watersheds3D.hxx>
#include <vigra/matlab.hxx>
#include <vigra/basicimageview.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/stdimagefunctions.hxx>
#include <vigra/localminmax.hxx>
#include <vigra/labelimage.hxx>
#include <vigra/seededregiongrowing.hxx>
#include <vigra/impex.hxx>

using namespace vigra;

struct GradientSquaredMagnitudeFunctor
{
//as
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



template <class InImage, class OutImage>
void watershedSeed(InImage & in, OutImage & out){
    typedef typename vigra::NumericTraits<typename InImage::value_type>::RealPromote
        TmpType;

    int w = in.width();
    int h = in.height();
	
	double scale = 1.0;
	
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
    extendedLocalMinima(srcImageRange(gradientmag), destImage(out), 1);

    // label the minima just found
    int max_region_label =
        labelImageWithBackground(srcImageRange(out), destImage(out),
                                 false, 0);

    // create a statistics functor for region growing
    vigra::ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<float> >
                                          gradstat(max_region_label);

    // perform region growing, starting from the minima of the gradient magnitude;
    // as the feature (first input) image contains the gradient magnitude,
    // this calculates the catchment basin of each minimum
    seededRegionGrowing(srcImageRange(gradientmag), srcImage(out),
                        destImage(out), gradstat);
/*
    // initialize a functor to determine the average gray-value or color
    // for each region (catchment basin) just found
    vigra::ArrayOfRegionStatistics<vigra::FindAverage<TmpType> >
                                          averages(max_region_label);

    // calculate the averages
    inspectTwoImages(srcImageRange(in), srcImage(labels), averages);

    // write the averages into the destination image (the functor 'averages'
    // acts as a look-up table)
    transformImage(srcImageRange(labels), destImage(out), averages);

    // mark the watersheds (region boundaries) black
    regionImageToEdgeImage(srcImageRange(labels), destImage(out),
                           vigra::NumericTraits<typename OutImage::value_type>::zero());
						   */
}

template <class InImage, class OutImage>
void watershedSeed4(InImage & in, OutImage & out){
    typedef typename vigra::NumericTraits<typename InImage::value_type>::RealPromote
        TmpType;

    int w = in.width();
    int h = in.height();
	
	double scale = 1.0;
	
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
    extendedLocalMinima(srcImageRange(gradientmag), destImage(out), 1, FourNeighborCode());
/*
    // label the minima just found
    int max_region_label =
        labelImageWithBackground(srcImageRange(out), destImage(out),
                                 false, 0);

    // create a statistics functor for region growing
    vigra::ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<float> >
                                          gradstat(max_region_label);

    // perform region growing, starting from the minima of the gradient magnitude;
    // as the feature (first input) image contains the gradient magnitude,
    // this calculates the catchment basin of each minimum
    seededRegionGrowing(srcImageRange(gradientmag), srcImage(out),
                        destImage(out), gradstat);

    // initialize a functor to determine the average gray-value or color
    // for each region (catchment basin) just found
    vigra::ArrayOfRegionStatistics<vigra::FindAverage<TmpType> >
                                          averages(max_region_label);

    // calculate the averages
    inspectTwoImages(srcImageRange(in), srcImage(out), averages);

    // write the averages into the destination image (the functor 'averages'
    // acts as a look-up table)
    transformImage(srcImageRange(out), destImage(out), averages);
	
    // mark the watersheds (region boundaries) black
    /*regionImageToEdgeImage(srcImageRange(out), destImage(out),
                           vigra::NumericTraits<typename OutImage::value_type>::zero());
/*						   */
}