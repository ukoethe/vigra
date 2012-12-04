
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

	int horFac, verFac;
	std::cout << "horizontal factor: ";
	std::cin >> horFac;
	std::cout << "vertical factor: ";
	std::cin >> verFac;
    
    try
    {
        // read image given as first argument
        // file type is determined automatically
        vigra::ImageImportInfo info(argv[1]);
        
        if(info.isGrayscale())
        {
			vigra::MultiArray<2, unsigned char> tempImage(Shape2(info.width(), info.height()));
			importImage(info, destImage(tempImage));
            // create a gray scale image of appropriate size
            vigra::MultiArray<2, unsigned char> in(Shape2(info.width()*horFac, info.height()*verFac));
	    	// import the image just read
            
			for (int j=0; j<verFac; j++) {           
				for (int i=0; i<horFac; i++) {
            		in.subarray(Shape2(info.width()*i, info.height()*j),Shape2(info.width()*(i+1), info.height()*(j+1))) = tempImage;
				}
			}
            
            // import the image just read
            //importImage(info, destImage(in));	    
            
            // write the image to the file given as second argument
            // the file type will be determined from the file name's extension
            exportImage(srcImageRange(in), vigra::ImageExportInfo(argv[2]));
        }
        else
        {
            
           	vigra::MultiArray<2, vigra::RGBValue<unsigned char> > tempImage(Shape2(info.width(), info.height()));
			importImage(info, destImage(tempImage));
            // create a gray scale image of appropriate size
            vigra::MultiArray<2, vigra::RGBValue<unsigned char> > in(Shape2(info.width()*horFac, info.height()*verFac));
	    	// import the image just read
            
			for (int j=0; j<verFac; j++) {           
				for (int i=0; i<horFac; i++) {
            		in.subarray(Shape2(info.width()*i, info.height()*j),Shape2(info.width()*(i+1), info.height()*(j+1))) = tempImage;
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

