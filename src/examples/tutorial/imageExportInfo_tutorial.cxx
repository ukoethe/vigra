#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/impex.hxx>

using namespace vigra; 

int main(int argc, char ** argv) 
{
    // instantiate array for image data of size 160x160 pixels
    MultiArray<2, UInt8> imageArray(Shape2(160,160));

    // create a black (0) and white (255) checker-board image
    for (int y = 0; y < imageArray.shape(1); y++) 
    {
        // note that the inner loop should go over the first (x-) axis
        // because elements along the first axis are consecutive in memory
        for (int x = 0; x < imageArray.shape(0); x++) 
        {
            if ((x%20)/10 == (y%20)/10)
                imageArray(x,y) = 0;
            else
                imageArray(x,y) = 255; 
        }
    }

    // write image data to "testimage.jpg" and set compression to 70%
    exportImage(imageArray, ImageExportInfo("testimage.jpg").setCompression("JPEG QUALITY=70"));
    
    // if you don't want to set any options in ImageExportInfo, you can simply write
    exportImage(imageArray, "testimage.gif");

    return 0;
}
