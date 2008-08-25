#include <vigra/matlab.hxx>
#include <vigra/multi_resize.hxx>

using namespace vigra;

template <class T>
void vigraResizeVolumeImpl(matlab::OutputArray outputs, matlab::InputArray inputs)
{    
    MultiArrayView<4, T> in = matlab::getMultiArray<4, T>(inputs[0]);
    MultiArrayShape<3>::type newShape = matlab::getShape<3>(inputs[1]);
    MultiArrayShape<4>::type newShape4(newShape[0], newShape[1], newShape[2], in.shape(3));
    MultiArrayView<4, T> out = matlab::createMultiArray<4, T>(newShape4, outputs[0]);
    unsigned int splineOrder = inputs.isValid(2) 
                                   ? matlab::getScalar<unsigned int>(inputs[2])
                                   : 3;

    for(int k=0; k<in.shape(3); ++k)
    {
        MultiArrayView<3, T> ink = in.bindOuter(k);    
        MultiArrayView<3, T> outk = out.bindOuter(k);
        switch(splineOrder)
        {
          case 0:
            resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<0>());
            break;
          case 1:
            resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<1>());
            break;
          case 2:
            resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<2>());
            break;
          case 3:
            resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<3>());
            break;
          case 4:
            resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<4>());
            break;
          case 5:
            resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<5>());
            break;
          default:
            mexErrMsgTxt("Interpolation order must be between 0 and 5.");
        }
    }
}

 
/** MATLAB 
function resized = resizeVolume(original, newShape)
function resized = resizeVolume(original, newShape, splineOrder)

resize a 3D array to a new shape using spline interpolation of the given order 
(default: using cubic spline interpolation).
    original    - 3D matrix with original volume data 
                  (gray scale or multi-band/RGB, type 'uint8', 'single' or 'double')
    newShape    - int32-array of length 3 giving the new shape
    splineOrder - order of interpolation 
                  (0 <= splineOrder <= 5, default: 3, i.e. cubic splines)
*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs)
{    
    if(matlab::ValueType<UInt8>::check(inputs[0]))
        vigraResizeVolumeImpl<UInt8>(outputs, inputs);
    else if(matlab::ValueType<float>::check(inputs[0]))
        vigraResizeVolumeImpl<float>(outputs, inputs);
    else if(matlab::ValueType<double>::check(inputs[0]))
        vigraResizeVolumeImpl<double>(outputs, inputs);
    else
        mexErrMsgTxt("Input image must have type 'uint8', 'single' or 'double'.");
}
