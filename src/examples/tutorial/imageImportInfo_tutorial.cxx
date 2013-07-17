#include <iostream>
#include <vigra/impex.hxx>

int main(int argc, char ** argv) 
{
   if(argc != 2)
    {
        std::cout << "Usage: " << argv[0] << " infile" << std::endl;
        std::cout << "(supported formats: " << vigra::impexListFormats() << ")" << std::endl;
        
        return 1;
    }
    
    try
    {
        // read image given as first command-line argument
        vigra::ImageImportInfo imageInfo(argv[1]);

        // print some information:
        std::cout << "Image information:\n";
        std::cout << "  file format: " << imageInfo.getFileType() << std::endl;
        std::cout << "  width:       " << imageInfo.width() << std::endl;
        std::cout << "  height:      " << imageInfo.height() << std::endl;
        std::cout << "  pixel type:  " << imageInfo.getPixelType() << std::endl;
        std::cout << "  color image: ";
        if (imageInfo.isColor())    std::cout << "yes (";
        else                        std::cout << "no  (";
        std::cout << "number of channels: " << imageInfo.numBands() << ")\n";
    }
    catch (std::exception & e) 
    {
        // catch any errors that might have occurred (e.g. "file not found") and print their reason
        std::cout << e.what() << std::endl;
        return 1;
    }
    return 0;
}
