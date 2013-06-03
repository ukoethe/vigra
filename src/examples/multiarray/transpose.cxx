
#include <iostream>
#include <vigra/multi_array.hxx>
#include "vigra/stdimage.hxx"
#include "vigra/impex.hxx"
#include <vigra/matrix.hxx>

using namespace vigra; 


int main(int argc, char ** argv)
{
    if(argc != 3)
    {
        std::cout << "Usage: " << argv[0] << " infile outfile" << std::endl;
        std::cout << "(supported formats: " << vigra::impexListFormats() << ")" << std::endl;
        
        return 1;
    }
    
    try
    {
        // read image given as first argument
        // file type is determined automatically
        vigra::ImageImportInfo info(argv[1]);
        
        if(info.isGrayscale())
        {
            // create a gray scale image of appropriate size
			vigra::MultiArray<2, unsigned char> temp(Shape2(info.width(), info.height()));
	    	// import the image just read
            importImage(info, destImage(temp));
            vigra::MultiArray<2, unsigned char> in = transpose(temp);

            // write the image to the file given as second argument
            // the file type will be determined from the file name's extension
            exportImage(srcImageRange(in), vigra::ImageExportInfo(argv[2]));
        }
        else
        {
            // create a RGB image of appropriate size
            vigra::MultiArray<2, vigra::RGBValue<unsigned char> > temp(Shape2(info.width(), info.height()));
            
            // import the image just read
            importImage(info, destImage(temp));
			vigra::MultiArray<2, vigra::RGBValue<unsigned char> > in = transpose(temp);
            
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
