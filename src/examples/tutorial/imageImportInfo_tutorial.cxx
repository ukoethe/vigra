#include <iostream>
#include "vigra/impex.hxx"

int main(int argc, char ** argv) {
    
    // read image given as first argument
    vigra::ImageImportInfo imageInfo(argv[1]);

    // Print some information:
    std::cout << "Image information:\n";
    std::cout << "Filetype: " << imageInfo.getFileType() << std::endl;
    std::cout << "Width: " << imageInfo.width() << std::endl;
    std::cout << "Height: " << imageInfo.height() << std::endl;
    std::cout << "Color Image: ";
    if (imageInfo.isColor())    std::cout << "yes\n";
    else                        std::cout << "no\n";

    return 0;
}
