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
#include "vigra/stdimagefunctions.hxx"
#include "vigra/distancetransform.hxx"
#include "vigra/labelimage.hxx"
#include "vigra/seededregiongrowing.hxx"
#include "vigra/impex.hxx"

using namespace vigra; 

int main(int argc, char ** argv)
{
    try
    {
        int number_of_points = 25;
        int size = 512;
        
        // create input image
        vigra::FImage in(size, size);
        
        // paint it black
        in = 0;
        
        for(int i=1; i<=number_of_points; ++i)
        {
            // mark a number of points 
            int x = (int)((float)rand() / RAND_MAX * size);
            int y = (int)((float)rand() / RAND_MAX * size);
            
            // label each point with a unique number
            in(x,y) = i;
        }
        
        // create output image and paint it white
        vigra::IImage out(size, size);
        out = 255;
        
        // in the output image, paint the points black which were 
        // marked in the input image
        initImageIf(destImageRange(out), maskImage(in), 0);
        
        // create image to hold the distance transform
        vigra::FImage distances(size, size);
        
        // calculate Euclidean distance transform
        distanceTransform(srcImageRange(in), destImage(distances), 0, 2);
        
        exportImage(srcImageRange(distances), vigra::ImageExportInfo("distances.gif"));
        std::cout << "Wrote distance transform (distances.gif)" << std::endl;
        
        // initialize statistics functor for region growing
        vigra::ArrayOfRegionStatistics<vigra::SeedRgDirectValueFunctor<float> > 
            statistics(number_of_points);

        // perform region growing, starting at the points marked in the input
        // image; as the feature (first source) image contains the distances,
        // the will calculate the voronoi regions of these points
        seededRegionGrowing(srcImageRange(distances), srcImage(in),
                            destImage(in), statistics);
        
        // in the output image, mark the borders of the voronoi regions black 
        regionImageToEdgeImage(srcImageRange(in), destImage(out), 0);

        exportImage(srcImageRange(out), vigra::ImageExportInfo("voronoi.gif"));
        std::cout << "Wrote voronoi diagram (voronoi.gif)" << std::endl;
    }
    catch (vigra::StdException & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
