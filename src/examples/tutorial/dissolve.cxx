#include <vigra/multi_array.hxx>
#include "vigra/impex.hxx"
#include <vigra/multi_math.hxx>
#include <iostream>

using namespace vigra;

int main (int argc, char ** argv) {

    // read the images given as first and second argument
    vigra::ImageImportInfo info1(argv[1]);
    vigra::ImageImportInfo info2(argv[2]);

    // process grayscale images
    if (info1.isGrayscale() && info2.isGrayscale()) {

        // instantiate array for image data
        vigra::MultiArray<2, unsigned char> imageArray1(info1.shape());
        vigra::MultiArray<2, unsigned char> imageArray2(info2.shape());

        // copy image data into array
        importImage(info1, destImage(imageArray1));
        importImage(info2, destImage(imageArray2));

        // calculate size of dissolved image
        int width = info1.width() < info2.width() ? info1.width() : info2.width();
        int height = info1.height() < info2.height() ? info1.height() : info2.height(); 
        
        // instantiate array of appropriate size for the dissolved image
        vigra::MultiArray<2, unsigned char> exportArray(Shape2(width, height));

        // create subviews for dissolving
        // the middle parts of the two images will be dissolved
        int x0, x1, y0, y1;
        x0 = 0.5*(info1.width()-width);
        y0 = 0.5*(info1.height()-height);
        x1 = 0.5*(info1.width()+width);
        y1 = 0.5*(info1.height()+height);
        vigra::MultiArray<2, unsigned char> subImageArray1 = 
            imageArray1.subarray(Shape2(x0,y0), Shape2(x1,y1));
  
        x0 = 0.5*(info2.width()-width);
        y0 = 0.5*(info2.height()-height);
        x1 = 0.5*(info2.width()+width);
        y1 = 0.5*(info2.height()+height);
        vigra::MultiArray<2, unsigned char> subImageArray2 = 
            imageArray2.subarray(Shape2(x0,y0), Shape2(x1,y1));

        // dissolve images
        using namespace vigra::multi_math;
        exportArray = 0.5*subImageArray1 + 0.5*subImageArray2;

        // write image data to the file given as third argument
        exportImage(srcImageRange(exportArray), vigra::ImageExportInfo(argv[3]));
    }

    // process color images
    else if (!info1.isGrayscale() && !info2.isGrayscale()) {
        // instantiate array for image data
        vigra::MultiArray<2, vigra::RGBValue<unsigned char> > imageArray1(info1.shape());
        vigra::MultiArray<2, vigra::RGBValue<unsigned char> > imageArray2(info2.shape());

        // copy image data into array
        importImage(info1, destImage(imageArray1));
        importImage(info2, destImage(imageArray2));

        // calculate size of dissolved image
        int width = info1.width() < info2.width() ? info1.width() : info2.width();
        int height = info1.height() < info2.height() ? info1.height() : info2.height(); 
        
        // instantiate array of appropriate size for the dissolved image
        vigra::MultiArray<2, vigra::RGBValue<unsigned char> > exportArray(Shape2(width, height));

        // create subviews for dissolving
        // the middle parts of the two images will be dissolved
        int x0, x1, y0, y1;
        x0 = 0.5*(info1.width()-width);
        y0 = 0.5*(info1.height()-height);
        x1 = 0.5*(info1.width()+width);
        y1 = 0.5*(info1.height()+height);
        vigra::MultiArray<2, vigra::RGBValue<unsigned char> > subImageArray1 = 
            imageArray1.subarray(Shape2(x0,y0), Shape2(x1,y1));
  
        x0 = 0.5*(info2.width()-width);
        y0 = 0.5*(info2.height()-height);
        x1 = 0.5*(info2.width()+width);
        y1 = 0.5*(info2.height()+height);
        vigra::MultiArray<2, vigra::RGBValue<unsigned char> > subImageArray2 = 
            imageArray2.subarray(Shape2(x0,y0), Shape2(x1,y1));

        // dissolve images
        using namespace vigra::multi_math;
        exportArray = 0.5*subImageArray1 + 0.5*subImageArray2;

        // write image data to the file given as third argument
        exportImage(srcImageRange(exportArray), vigra::ImageExportInfo(argv[3]));
    }

    else {
        std::cout << "Error: Cannot dissolve color and grayscale image!\n";
    }
    return 0;
}
