#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/multi_math.hxx>
#include <iostream>

using namespace vigra;

int main (int argc, char ** argv) 
{
    if(argc != 4)
    {
        std::cout << "Usage: " << argv[0] << " infile1 infile2 outfile" << std::endl;
        std::cout << "(supported formats: " << impexListFormats() << ")" << std::endl;
        
        return 1;
    }
    
    try 
    {
        // read the images given as first and second argument
        ImageImportInfo info1(argv[1]),
                        info2(argv[2]);
        
        // for simplicity, we always use RGB images, even when the inputs are grayscale
        // (in which case all three channels will be equal)
        MultiArray<2, RGBValue<UInt8> > imageArray1(info1.shape()),
                                        imageArray2(info2.shape());
        
        importImage(info1, imageArray1);
        importImage(info2, imageArray2);
        
        // calculate size of dissolved image (take the minimum along both axes)
        Shape2 resultShape = min(info1.shape(), info2.shape());
        
        // instantiate array of appropriate size for the dissolved image
        MultiArray<2, RGBValue<UInt8> > exportArray(resultShape);

        // create subviews for dissolving
        // the middle parts of the two images will be dissolved
        Shape2 start1 = (info1.shape() - resultShape) / 2,
               end1   = start1 + resultShape;
        MultiArray<2, RGBValue<UInt8> > subImageArray1 = imageArray1.subarray(start1, end1);

        Shape2 start2 = (info2.shape() - resultShape) / 2,
               end2   = start2 + resultShape;
        MultiArray<2, RGBValue<UInt8> > subImageArray2 = imageArray2.subarray(start2, end2);

        // dissolve images
        using namespace multi_math;
        exportArray = 0.5*subImageArray1 + 0.5*subImageArray2;

        // write image data to the file given as third argument
        exportImage(exportArray, ImageExportInfo(argv[3]));
    }
    catch (std::exception & e) 
    {
        // catch any errors that might have occurred and print their reason
        std::cout << e.what() << std::endl;
        return 1;
    }
    return 0;
}
