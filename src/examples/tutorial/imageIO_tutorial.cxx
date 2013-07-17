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
        std::cout << "(grayscale only, supported formats: " << impexListFormats() << ")" << std::endl;
        
        return 1;
    }
    
    try 
    {
        char * in_filename  = argv[1];
        char * out_filename = argv[2];
        
        // read image header information from in_filename
        ImageImportInfo imageInfo(in_filename);

        // instantiate array for image data
        MultiArray<2, UInt8> imageArray(imageInfo.shape());

        // copy image data from file into array
        importImage(imageInfo, imageArray);
        
        // if you don't need the information from ImageImportInfo, you can also
        // simply pass the filename (this will resize imageArray internally!)
        importImage(in_filename, imageArray);

        // set every second horizontal line to black
        for (int i = 0; i<imageInfo.height(); i+=2) 
        {
            imageArray.bind<1>(i) = 0;   
        }

        // write image data to the file given as second argument
        exportImage(imageArray, ImageExportInfo(out_filename));
    }
    catch (std::exception & e) 
    {
        // catch any errors that might have occurred and print their reason
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}
