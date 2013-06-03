
#include <iostream>
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
    
    try
    {
        // read image given as first argument
        // file type is determined automatically
        vigra::ImageImportInfo info(argv[1]);
        
        if(info.isGrayscale())
        {
            // create a gray scale image of appropriate size
            vigra::BImage in(info.width()*2, info.height()*2);
            
            // import the image just read
            importImage(info, destImage(in));
	    
            vigra::BImage::Iterator pos = in.upperLeft();
            vigra::BImage::Iterator end = in.lowerRight();
            int i = 0;	
            for (; pos.y != end.y; ++pos.y) {
                for(pos.x = in.width()/2; pos.x != end.x; ++pos.x) {
         		*pos = (i++)%255;
                }
            }
            
            // write the image to the file given as second argument
            // the file type will be determined from the file name's extension
            exportImage(srcImageRange(in), vigra::ImageExportInfo(argv[2]));
        }
        else
        {
            // create a RGB image of appropriate size
            vigra::BRGBImage in(info.width(), info.height());
            
            // import the image just read
            importImage(info, destImage(in));
            
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

