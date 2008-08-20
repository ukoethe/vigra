#include <vigra/matlab.hxx>
#include <vigra/multi_resize.hxx>

using namespace vigra;
 
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs)
{
    MultiArrayView<3, double> in = matlab::getMultiArray<3, double>(inputs[0]);    
    MultiArrayShape<3>::type newShape = matlab::getShape<3>(inputs[1]);
    MultiArrayView<3, double> out = matlab::createMultiArray<3, double>(newShape, outputs[0]);

    resizeMultiArraySplineInterpolation(srcMultiArrayRange(in), destMultiArrayRange(out));
}
