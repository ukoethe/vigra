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
        // read image given as first argument
        // file type is determined automatically
        vigra::ImageImportInfo info(argv[1]);
        
        if(info.isGrayscale()) {
            
            // write image data to MultiArray
			vigra::MultiArray<2, unsigned char> imageArray(Shape2(info.shape()));
			importImage(info, destImage(imageArray));

            // calculate upper-left (x0,y0) and lower-right (x1,y1) position of subimage
            int x0 = info.width()/4;
            int y0 = info.height()/4;
            int x1 = info.width()-x0;
            int y1 = info.height()-y0;

            // copy subimage into new array for output
            vigra::MultiArray<2, unsigned char> newImageArray = 
                    imageArray.subarray(Shape2(x0,y0), Shape2(x1, y1));
            
            // write the image to the file given as second argument
            // the file type will be determined from the file name's extension
            exportImage(srcImageRange(newImageArray), vigra::ImageExportInfo(argv[2]));
        }
        else {
            // write image data to MultiArray
			vigra::MultiArray<2, vigra::RGBValue<unsigned char> > imageArray(Shape2(info.shape()));
			importImage(info, destImage(imageArray));

            // calculate upper-left (x0,y0) and lower-right (x1,y1) position of subimage
            int x0 = info.width()/4;
            int y0 = info.height()/4;
            int x1 = info.width()-x0;
            int y1 = info.height()-y0;

            // copy subimage into new array for output
            vigra::MultiArray<2, vigra::RGBValue<unsigned char> > newImageArray = 
                    imageArray.subarray(Shape2(x0,y0), Shape2(x1, y1));
            
            // write the image to the file given as second argument
            // the file type will be determined from the file name's extension
            exportImage(srcImageRange(newImageArray), vigra::ImageExportInfo(argv[2]));
        }
    }
    catch (vigra::StdException & e) {
        // catch any errors that might have occurred and print their reason
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
