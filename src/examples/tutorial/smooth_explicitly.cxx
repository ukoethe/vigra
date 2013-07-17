#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>

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
        // Read the header information of the input image.
        ImageImportInfo info(argv[1]);

        if (info.isGrayscale())  // Process grayscale image.
        {   
            // Instantiate arrays of appropriate size for the original and smoothed image.
            MultiArray<2, float> inputImage(info.shape()),
                                 resultImage(info.shape());
            
            // Read the input data.
            importImage(info, inputImage);

            // Smooth image by calculating the mean of a 5x5 window around each pixel.
            Shape2 current; // the coordinate of the current pixel
            for(current[1] = 0; current[1] < inputImage.shape(1); ++current[1])
            {
                for (current[0] = 0; current[0] < inputImage.shape(0); ++current[0])
                {
                    // Compute the corners of the current smoothing window.
                    // Note:
                    //    * Near the image border, not all pixels in a 5x5 window are within
                    //      the image. The max() and min() operations make sure that the
                    //      window is appropriately truncated (for example, max(Shape2(0), ...) 
                    //      sets negative coordinates to zero).
                    //    * By convention, the lower right corner 'windowStop' must be outside the
                    //      window. Therefore, we have to add 'Shape2(3)' (a vector of 3s) to 'current'. 
                    //      In contrast, the upper left corner 'windowStart' must be inside the window, 
                    //      so we only subtract 'Shape2(2)' (a vector of 2s) from 'current'. 
                    Shape2 windowStart = max(Shape2(0),          current - Shape2(2));
                    Shape2 windowStop  = min(inputImage.shape(), current + Shape2(3));
                    
                    // Create a view to the current window.
                    MultiArrayView<2, float> window = inputImage.subarray(windowStart, windowStop);
                    
                    // Compute and store the average of the current window.
                    resultImage[current] = window.sum<float>() / window.size();
                }
            }

            // Write the result image to the file given as second argument.
            exportImage(resultImage, ImageExportInfo(argv[2]));
        }
        else // process color image
        {
            // Everything as before, but the pixel type is now RGBValue<float>
            MultiArray<2, RGBValue<float> > inputImage(info.shape());
            MultiArray<2, RGBValue<float> > resultImage(info.shape());
             
            importImage(info, inputImage);

            Shape2 current;
            for(current[1] = 0; current[1] < inputImage.shape(1); ++current[1])
            {
                for (current[0] = 0; current[0] < inputImage.shape(0); ++current[0])
                {
                    Shape2 windowStart = max(Shape2(0),          current - Shape2(2));
                    Shape2 windowStop  = min(inputImage.shape(), current + Shape2(3));
                    
                    MultiArrayView<2, RGBValue<float> > window = inputImage.subarray(windowStart, windowStop);
                    
                    resultImage[current] = window.sum<RGBValue<float> >() / window.size();
                }
            }

            exportImage(resultImage, ImageExportInfo(argv[2]));
        }
        return 0;
    }
    catch (std::exception & e) 
    {
        // Catch any errors that might have occurred and print their reason.
        std::cout << e.what() << std::endl;
        return 1;
    }
}
