#include <vigra/multi_array.hxx>
#include "vigra/impex.hxx"
#include <vigra/stdconvolution.hxx>
#include <iostream>

using namespace vigra;

int main (int argc, char ** argv) {

    // read image given as first argument
    vigra::ImageImportInfo info(argv[1]);

    // process grayscale image
    if (info.isGrayscale()) {
       
        // instantiate arrays for image data and for smoothed image of appropriate size
        vigra::MultiArray<2, float> imageArray(info.shape());
        vigra::MultiArray<2, float> exportArray(info.shape());

        // copy image data into array
        importImage(info, destImage(imageArray));

        // instantiate convolving kernel
        vigra::Kernel2D<double> kernel2dim;
        vigra::Kernel1D<double> kernel1dim;
        
        // choose of convolving mode
        std::cout << "Which mode of convolution?\n";
        std::cout << "1 - 3x3-average\n";
        std::cout << "2 - 3x3-box just corners\n";
        std::cout << "3 - gaussian\n";
        std::cout << "4 - X- and Y-dimension separately\n";
        int mode;
        std::cin >> mode;
        
        // initializing kernel       
        switch(mode) {
            // homogenous averaging of 3x3-box
            case 1: kernel2dim.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) = 1.0/9.0; break;
            // averaging corners of 3x3-box
            case 2: kernel2dim.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) = 
                        0.25, 0.0, 0.25, 
                        0.0,  0.0, 0.0,
                        0.25, 0.0, 0.25; break;
            // gaussian convolving
            case 3: kernel2dim.initGaussian(1.5); break;
            case 4: kernel1dim.initExplicitly(-1,1) = 1.0/3.0; break;
        }
       
        // convolve image
        if (mode == 4) {
             separableConvolveX(srcImageRange(imageArray), destImage(imageArray), kernel1d(kernel1dim));
             separableConvolveY(srcImageRange(imageArray), destImage(exportArray), kernel1d(kernel1dim));
        }
        else {
            convolveImage(srcImageRange(imageArray), destImage(exportArray), kernel2d(kernel2dim));
        }

        // write image data to the file given as second argument
        exportImage(srcImageRange(exportArray), vigra::ImageExportInfo(argv[2]));
    }

    // process color image
    else {
       
        // instantiate arrays for image data and for smoothed image of appropriate size
        vigra::MultiArray<2, vigra::RGBValue<float> > imageArray(info.shape());
        vigra::MultiArray<2, vigra::RGBValue<float> > exportArray(info.shape());

        // copy image data into array
        importImage(info, destImage(imageArray));

        // instantiate convolving kernel
        vigra::Kernel2D<double> kernel2dim;
        vigra::Kernel1D<double> kernel1dim;
        
        // choose of convolving mode
        std::cout << "Which mode of convolution?\n";
        std::cout << "1 - 3x3-average\n";
        std::cout << "2 - 3x3-box just corners\n";
        std::cout << "3 - gaussian\n";
        std::cout << "4 - X- and Y-dimension separately\n";
        int mode;
        std::cin >> mode;
        
        // initializing kernel       
        switch(mode) {
            // homogenous averaging of 3x3-box
            case 1: kernel2dim.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) = 1.0/9.0; break;
            // averaging corners of 3x3-box
            case 2: kernel2dim.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) = 
                        0.25, 0.0, 0.25, 
                        0.0,  0.0, 0.0,
                        0.25, 0.0, 0.25; break;
            // gaussian convolving
            case 3: kernel2dim.initGaussian(1.5); break;
            case 4: kernel1dim.initExplicitly(-1,1) = 1.0/3.0; break;
        }
       
        // convolve image
        if (mode == 4) {
             separableConvolveX(srcImageRange(imageArray), destImage(imageArray), kernel1d(kernel1dim));
             separableConvolveY(srcImageRange(imageArray), destImage(exportArray), kernel1d(kernel1dim));
        }
        else {
            convolveImage(srcImageRange(imageArray), destImage(exportArray), kernel2d(kernel2dim));
        }

        // write image data to the file given as second argument
        exportImage(srcImageRange(exportArray), vigra::ImageExportInfo(argv[2]));
        
    }
    return 0;
}
