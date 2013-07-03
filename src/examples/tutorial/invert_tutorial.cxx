#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/multi_math.hxx>

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
        // use multi_math functions
        using namespace multi_math;
        
        // read image given as first argument
        ImageImportInfo info(argv[1]);

        // process grayscale image
        if (info.isGrayscale()) 
        {
            // instantiate array for image data
            MultiArray<2, UInt8> imageArray(info.shape());
            
            // import image data into array
            importImage(info, imageArray);
            
            // invert image
            imageArray = 255-imageArray;
            
            // write image data to the file given as second argument
            exportImage(imageArray, ImageExportInfo(argv[2]));
        }
        // process color image
        else 
        {
            // instantiate array for image data
            MultiArray<2, RGBValue<UInt8> > imageArray(info.shape());
            
            // copy image data into array
            importImage(info, imageArray);
            
            // invert image
            RGBValue<UInt8> temp(255,255,255);
            imageArray = temp - imageArray;
            
            // write image data to the file given as second argument
            exportImage(imageArray, ImageExportInfo(argv[2]));
        }
    }
    catch (std::exception & e) 
    {
        // catch any errors that might have occurred (e.g. "file not found") and print their reason
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}
