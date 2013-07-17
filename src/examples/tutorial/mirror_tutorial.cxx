#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/basicgeometry.hxx>
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

    // choose axis to mirror about
    std::cout << "How to mirror image?\n";
    std::cout << "1 - top to bottom (= about horizontal axis)\n";
    std::cout << "2 - left to right (= about vertical axis)\n";
    int mode;
    std::cin >> mode;

    try 
    {
        ImageImportInfo info(argv[1]);
        
        MultiArray<2, UInt8> imageArray(info.shape());
        importImage(info, imageArray);
        
        MultiArray<2, UInt8> newImageArray(info.shape());

        if(mode == 1)
        {
            // mirror the image horizontally 
            // (for didactic reasons, we implement this variant explicitly,
            //  note that info.height()-1 is equal to the last y-index)
            for (int y=0; y<info.height(); y++) 
            {           
                newImageArray.bind<1>(y) = imageArray.bind<1>(info.height()-1-y);
            }
        }
        else
        {
            // mirror the image vertically
            reflectImage(imageArray, newImageArray, vertical);
        }
        
        exportImage(newImageArray, ImageExportInfo(argv[2]));
    }
    catch (std::exception & e) 
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

