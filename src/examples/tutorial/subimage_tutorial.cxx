#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/impex.hxx>

using namespace vigra; 


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
        // read image given as first argument
        // file type is determined automatically
        ImageImportInfo info(argv[1]);
        
        if(info.isGrayscale()) 
        {
            // read image from file
			MultiArray<2, UInt8> imageArray(info.shape());
			importImage(info, imageArray);

            // we want to cut out the center of the original image, such that the 
            // size of the remaining image is half the original
            Shape2 upperLeft  = info.shape() / 4,
                   lowerRight = info.shape() - upperLeft;
        
            // create subimage around center for output
            MultiArrayView<2, UInt8> subimage = imageArray.subarray(upperLeft, lowerRight);
            
            // write the subimage to the file provided as second command line argument
            // the file type will be determined from the file name's extension
            exportImage(subimage, ImageExportInfo(argv[2]));
        }
        else 
        {
            // read image from file
			MultiArray<2, RGBValue<UInt8> > imageArray(Shape2(info.shape()));
			importImage(info, imageArray);

            // we want to cut out the center of the original image, such that the 
            // size of the remaining image is half the original
            Shape2 upperLeft  = info.shape() / 4,
                   lowerRight = info.shape() - upperLeft;
        
            // create subimage around center for output
            MultiArrayView<2, RGBValue<UInt8> > subimage = imageArray.subarray(upperLeft, lowerRight);
            
            // write the subimage to the file provided as second command line argument
            // the file type will be determined from the file name's extension
            exportImage(subimage, ImageExportInfo(argv[2]));
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
