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

        vigra::Kernel2D<float> kernel;
        std::cout << "Which mode of convolution?\n";
        std::cout << "1 - 3x3-average\n";
        std::cout << "2 - 3x3-box just corners\n";
        std::cout << "3 - gaussian\n";
        int mode;
        std::cin >> mode;
        switch(mode) {
            case 1: kernel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) = 1.0/9.0; break;
            case 2: kernel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) = 
                        0.25, 0.0, 0.25, 
                        0.0,  0.0, 0.0,
                        0.25, 0.0, 0.25; break;
            case 3: kernel.initGaussian(1.5); break;
        }
       
        convolveImage(srcImageRange(imageArray), destImage(exportArray), kernel2d(kernel));

        // write image data to the file given as second argument
        exportImage(srcImageRange(exportArray), vigra::ImageExportInfo(argv[2]));
    }

    // process color image
    else {
       
/*        // instantiate arrays for image data and for smoothed image of appropriate size
        vigra::MultiArray<2, vigra::RGBValue<float> > imageArray(info.shape());
        vigra::MultiArray<2, vigra::RGBValue<float> > exportArray(info.shape());

        // copy image data into array
        importImage(info, destImage(imageArray));

        // smooth image by calculating the mean of a 3x3-field around the pixel
        int x0, x1, y0, y1, elements;
        for (int i = 0; i<exportArray.size(0); i++) {

            // set top-left x-coordinate of the field
            if (i==0)   x0 = 0;         // pixel is on the left edge
            else        x0 = i-1;
            // set bottom-right x-coordinate of the field 
            // (remember: the bottom-right point is not included in subarray-method (that's why +2)
            if (i==(exportArray.size(0)-1)) x1 = i;         // pixel is on the right edge
            else                            x1 = i+2;

            for (int j = 0; j<exportArray.size(1); j++) {

                // set top-left y-coordinate of the field
                if (j==0)   y0 = 0;     // pixel is on the upper edge
                else        y0 = j-1;
                // set bottom-right x-coordinate of the field
                // (remember: the bottom-right point is not included in subarray-method (that's why +2)
                if (j==(exportArray.size(1)-1)) y1 = j;     // pixel is on the lower edge
                else                            y1 = j+2;

                // calculate sum of the field
                vigra::RGBValue<float> sum = 
                    imageArray.subarray(Shape2(x0, y0), Shape2(x1, y1)).sum<vigra::RGBValue<float> >();
                // determine number of elements in the field
                elements = imageArray.subarray(Shape2(x0, y0), Shape2(x1, y1)).size();
                // calculate the mean and set the value
                exportArray[Shape2(i,j)] = sum/elements;
            }
        }

        // write image data to the file given as second argument
        exportImage(srcImageRange(exportArray), vigra::ImageExportInfo(argv[2]));*/
    }
    return 0;
}
