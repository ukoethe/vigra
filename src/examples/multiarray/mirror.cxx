
#include <iostream>
#include <vigra/multi_array.hxx>
#include "vigra/stdimage.hxx"
#include "vigra/impex.hxx"

using namespace vigra; 


int main(int argc, char ** argv)
{
    if(argc != 3)
    {
        std::cout << "Usage: " << argv[0] << " infile outfile" << std::endl;
        std::cout << "(supported formats: " << vigra::impexListFormats() << ")" << std::endl;
        
        return 1;
    }

    // ask which axis should be used
	int axis;
	std::cout << "Which axis? horizontal (1) or vertical (2): ";
	std::cin >> axis;
    
    try
    {
        // read image given as first argument
        // file type is determined automatically
        vigra::ImageImportInfo info(argv[1]);
        
        if(info.isGrayscale())
        {   
            // write image data to MultiArray
			vigra::MultiArray<2, unsigned char> tempImage(Shape2(info.width(), info.height()));
			importImage(info, destImage(tempImage));
            // create a gray scale image of appropriate size
            vigra::MultiArray<2, unsigned char> in(Shape2(info.width(), info.height()));

            // mirror horizontal
            if (axis == 1) {
				for (int i=0; i<info.height(); i++) {           
    	        		in.bind<1>(i) = tempImage.bind<1>(info.height()-(i+1));
				}
			}

            // mirror vertical
			if (axis == 2) {
				for (int i=0; i<info.width(); i++) {           
    	        		in.bind<0>(i) = tempImage.bind<0>(info.width()-(i+1));
				}
			}    
            
            // write the image to the file given as second argument
            // the file type will be determined from the file name's extension
            exportImage(srcImageRange(in), vigra::ImageExportInfo(argv[2]));
        }
        else
        {
            // write image data to MultiArray
           	vigra::MultiArray<2, vigra::RGBValue<unsigned char> > tempImage(Shape2(info.width(), info.height()));
			importImage(info, destImage(tempImage));
            // create a RGB-image of appropriate size
            vigra::MultiArray<2, vigra::RGBValue<unsigned char> > in(Shape2(info.width(), info.height()));
            
            // mirror horizontal
			if (axis == 1) {
				for (int i=0; i<info.height(); i++) {           
    	        		in.bind<1>(i) = tempImage.bind<1>(info.height()-(i+1));
				}
			}

            // mirror vertical
			if (axis == 2) {
				for (int i=0; i<info.width(); i++) {           
    	        		in.bind<0>(i) = tempImage.bind<0>(info.width()-(i+1));
				}
			}
            
            // write the image to the file given as second argument
            // the file type will be determined from the file name's extension
            exportImage(srcImageRange(in), vigra::ImageExportInfo(argv[2]));
        }
    }
    catch (vigra::StdException & e)
    {
        // catch any errors that might have occurred and print their reason
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

