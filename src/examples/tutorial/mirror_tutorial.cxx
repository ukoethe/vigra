#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/stdimage.hxx>
#include <vigra/impex.hxx>

using namespace vigra; 


int main(int argc, char ** argv) {
    if(argc != 3) {
        std::cout << "Usage: " << argv[0] << " infile outfile" << std::endl;
        std::cout << "(supported formats: " << vigra::impexListFormats() << ")" << std::endl;
        
        return 1;
    }

    try {
        vigra::ImageImportInfo info(argv[1]);
        
        if(info.isGrayscale()) {
			vigra::MultiArray<2, unsigned char> imageArray(Shape2(info.width(), info.height()));
			importImage(info, destImage(imageArray));
            vigra::MultiArray<2, unsigned char> newImageArray(Shape2(info.width(), info.height()));

            // mirror the image horizontal 
            // info.height() is equal to the index of the last row
			for (int i=0; i<info.height(); i++) {           
	        		newImageArray.bind<1>(i) = imageArray.bind<1>(info.height()-(i+1));
			}
            
            exportImage(srcImageRange(newImageArray), vigra::ImageExportInfo(argv[2]));
        }
        else {
            std::cout << "no gray scale image" << std::endl;
        }
    }
    catch (vigra::StdException & e) {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

