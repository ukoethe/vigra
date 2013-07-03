#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/convolution.hxx>
#include <iostream>

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
        // read image given as first argument
        ImageImportInfo info(argv[1]);

        // choose convolving mode
        std::cout << "Which mode of convolution?\n";
        std::cout << "1 - disk average\n";
        std::cout << "2 - corners of a 3x3-box\n";
        std::cout << "3 - Gaussian filter\n";
        std::cout << "4 - separable Gaussian filter (x- and x-dimension separately)\n";
        int mode;
        std::cin >> mode;
        
        // instantiate convolution kernels
        Kernel2D<double> kernel2dim;
        Kernel1D<double> kernel1dim;
        
        // initialize kernel
        switch(mode) 
        {
            case 1: 
                // homogenous averaging within a disk of radius 2
                kernel2dim.initDisk(2); 
                break;
            case 2: 
                // strange custom-made filter: averaging the corners of a 3x3-box
                kernel2dim.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) = 
                        0.25, 0.0, 0.25, 
                        0.0,  0.0, 0.0,
                        0.25, 0.0, 0.25; 
                        break;
            case 3: 
                // gaussian convolution
                kernel2dim.initGaussian(1.5);
                break;
            case 4: 
                // separable gaussian convolution
                kernel1dim.initGaussian(1.5);
                break;
            default:
                vigra_precondition(false, "mode must be between 1 and 4.");
        }
           
        // process grayscale image
        if (info.isGrayscale()) 
        {
            // instantiate arrays for image data and for smoothed image of appropriate size
            MultiArray<2, float> imageArray(info.shape()),
                                 exportArray(info.shape());

            // copy image data into array
            importImage(info, imageArray);

            // convolve image
            if (mode == 4) 
            {
                 convolveImage(imageArray, exportArray, kernel1dim, kernel1dim);
            }
            else 
            {
                convolveImage(imageArray, exportArray, kernel2dim);
            }

            // write image data to the file given as second argument
            exportImage(exportArray, ImageExportInfo(argv[2]));
        }
        // process color image
        else 
        {
            // instantiate arrays for image data and for smoothed image of appropriate size
            MultiArray<2, RGBValue<float> > imageArray(info.shape());
            MultiArray<2, RGBValue<float> > exportArray(info.shape());

            // copy image data into array
            importImage(info, imageArray);

            // convolve image
            if (mode == 4) 
            {
                 convolveImage(imageArray, exportArray, kernel1dim, kernel1dim);
            }
            else 
            {
                convolveImage(imageArray, exportArray, kernel2dim);
            }

            // write image data to the file given as second argument
            exportImage(exportArray, ImageExportInfo(argv[2]));
        }
    }
    catch (std::exception & e) 
    {
        // catch any errors that might have occurred and print their reason
        std::cout << e.what() << std::endl;
        return 1;
    }
    return 0;
}
