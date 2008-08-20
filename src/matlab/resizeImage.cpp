#include <vigra/matlab.hxx>
#include <vigra/basicimageview.hxx>
#include <vigra/resizeimage.hxx>

using namespace vigra;
 
/** MATLAB 
function resized = resizeImage(original, newShape)

resize a 2D array to a new shape using cubic spline interpolation.
    original - 2D matrix with original image
    newShape - int32-array of length 2 giving the new shape
*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs)
{    
    BasicImageView<double> in = matlab::getImage<double>(inputs[0]);    
    MultiArrayShape<2>::type newShape = matlab::getShape<2>(inputs[1]);
    BasicImageView<double> out = matlab::createImage<double>(newShape[0], newShape[0], outputs[0]);

    resizeImageSplineInterpolation(srcImageRange(in), destImageRange(out));
}
