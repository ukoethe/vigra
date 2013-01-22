#include <vigra/multi_array.hxx>
#include "vigra/impex.hxx"

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
                float sum = imageArray.subarray(Shape2(x0, y0), Shape2(x1, y1)).sum<float>();
                // determine number of elements in the field
                elements = imageArray.subarray(Shape2(x0, y0), Shape2(x1, y1)).size();
                // calculate the mean and set the value
                exportArray[Shape2(i,j)] = sum/elements;
            }
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
        exportImage(srcImageRange(exportArray), vigra::ImageExportInfo(argv[2]));
    }
    return 0;
}
