
#include <iostream>
#include <vigra/multi_array.hxx>
#include "vigra/stdimage.hxx"
#include "vigra/impex.hxx"
#include <vigra/rgbvalue.hxx>

using namespace vigra; 


int main(int argc, char ** argv)
{
    if(argc != 3)
    {
        std::cout << "Usage: " << argv[0] << " infile outfile" << std::endl;
        std::cout << "(supported formats: " << vigra::impexListFormats() << ")" << std::endl;
        
        return 1;
    }

	std::cout << "Which component? Red (1), Green (2) or Blue (3)? ";
	int comp;
	std::cin >> comp;

    try
    {
        // read image given as first argument
        // file type is determined automatically
        vigra::ImageImportInfo info(argv[1]);
        
        if(info.isGrayscale())
        {
			std::cout << "This is a gray scale image. Please use a color image!\n";
        }
        else
        {
            // create a gray scale image of appropriate size
            vigra::MultiArray<2, vigra::RGBValue<unsigned char> > in(Shape2(info.width(), info.height()));
	    	// import the image just read
			importImage(info, destImage(in));
            
			switch (comp) {
				case 1: 
					for (vigra::MultiArray<2, vigra::RGBValue<unsigned char> >::iterator i=in.begin(); i != in.end(); i++) {
						(*i).setBlue(0);
						(*i).setGreen(0);
					}
					break;
				case 2:
					for (vigra::MultiArray<2, vigra::RGBValue<unsigned char> >::iterator i=in.begin(); i != in.end(); i++) {
						(*i).setRed(0);
						(*i).setBlue(0);
					}
					break;
				case 3: 
					for (vigra::MultiArray<2, vigra::RGBValue<unsigned char> >::iterator i=in.begin(); i != in.end(); i++) {
						(*i).setRed(0);
						(*i).setGreen(0);
					}
					break;
				default: std::cout << "Wrong input. Please choose 1, 2 or 3!\n"; break;
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

