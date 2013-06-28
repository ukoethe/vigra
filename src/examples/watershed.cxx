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
#include <vigra/multi_convolution.hxx>
#include <vigra/multi_watersheds.hxx>
#include <vigra/impex.hxx>

using namespace vigra; 

template <class InImage, class OutImage>
void watershedSegmentation(InImage & in, OutImage & out, double scale)
{
    // compute the gradient magnitude as a suitable boundary indicator
    MultiArray<2, float> gradient(in.shape());
    gaussianGradientMagnitude(in, gradient, scale);
    
    // Compute watershed segmentation using a region growing algorithm with 4-neighborhood.
    // Use option object to ask the watershed algorithm to compute seeds automatically at 
    // minima of the boundary indicator.
    MultiArray<2, unsigned int> labeling(in.shape());
    unsigned int max_region_label = 
        watershedsMultiArray(gradient, labeling, DirectNeighborhood,
                             WatershedOptions().seedOptions(SeedOptions().minima()));

    // Initialize a functor to determine the average gray-value or color
    // for each region (catchment basin) just found.
    // define a temporary type that can hold floating point values
    typedef typename NumericTraits<typename InImage::value_type>::RealPromote
        TmpType;
    ArrayOfRegionStatistics<FindAverage<TmpType> > averages(max_region_label);

    // calculate the averages
    inspectTwoImages(in, labeling, averages);

    // write the averages into the destination image (the functor 'averages'
    // acts as a look-up table)
    transformImage(labeling, out, averages);

    // mark the watersheds (region boundaries) black
    regionImageToEdgeImage(labeling, out,
                           NumericTraits<typename OutImage::value_type>::zero());
}


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

        // input width of gradient filter
        double scale = 1.0;
        std::cout << "Scale for gradient calculation ? ";
        std::cin >> scale;

        if(info.isGrayscale())
        {
            int w = info.width();
            int h = info.height();

            MultiArray<2, UInt8> in(w, h), out(w, h);
            
            importImage(info, in);

            // perform watershed segmentation on gray image
            // note that the watershed algorithm usually results in an
            // oversegmentation (i.e., too many regions), but its boundary
            // localization is quite good
            watershedSegmentation(in, out, scale);

            std::cout << "Writing " << argv[2] << std::endl;
            exportImage(out, ImageExportInfo(argv[2]));
        }
        else
        {
            int w = info.width();
            int h = info.height();

            MultiArray<2, RGBValue<UInt8> > in(w, h),  out(w, h);
            importImage(info, in);

            // perform watershed segmentation on color image
            // note that the watershed algorithm usually results in an
            // oversegmentation (i.e., too many regions), but its boundary
            // localization is quite good
            watershedSegmentation(in, out, scale);

            std::cout << "Writing " << argv[2] << std::endl;
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
