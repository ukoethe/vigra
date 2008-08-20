#include <vigra/matlab.hxx>
#include <vigra/multi_resize.hxx>

using namespace vigra;

/** MATLAB 
function resized = resizeVolume(original, newShape)

resize a 3D array to a new shape using cubic spline interpolation.
    original - 3D matrix with original data
    newShape - int32-array of length 3 giving the new shape
*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs)
{
    MultiArrayView<3, double> in = matlab::getMultiArray<3, double>(inputs[0]);    
    MultiArrayShape<3>::type newShape = matlab::getShape<3>(inputs[1]);
    MultiArrayView<3, double> out = matlab::createMultiArray<3, double>(newShape, outputs[0]);

    resizeMultiArraySplineInterpolation(srcMultiArrayRange(in), destMultiArrayRange(out));
}
