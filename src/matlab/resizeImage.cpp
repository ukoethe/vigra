#include <vigra/matlab.hxx>
#include <vigra/basicimageview.hxx>
#include <vigra/resizeimage.hxx>

using namespace vigra;

template <class T>
void vigraResizeImageImpl(matlab::OutputArray outputs, matlab::InputArray inputs)
{    
    MultiArrayView<3, T> in = matlab::getMultiArray<3, T>(inputs[0]);
    MultiArrayShape<2>::type newShape = matlab::getShape<2>(inputs[1]);
    MultiArrayShape<3>::type newShape3(newShape[0], newShape[1], in.shape(2));
    MultiArrayView<3, T> out = matlab::createMultiArray<3, T>(newShape3, outputs[0]);
    unsigned int splineOrder = inputs.isValid(2) 
                                   ? matlab::getScalar<unsigned int>(inputs[2])
                                   : 3;

    for(int k=0; k<in.shape(2); ++k)
    {
        BasicImageView<T> ink = makeBasicImageView(in.bindOuter(k));    
        BasicImageView<T> outk = makeBasicImageView(out.bindOuter(k));    
        switch(splineOrder)
        {
          case 0:
            resizeImageNoInterpolation(srcImageRange(ink), destImageRange(outk));
            break;
          case 1:
            resizeImageLinearInterpolation(srcImageRange(ink), destImageRange(outk));
            break;
          case 2:
            resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), BSpline<2>());
            break;
          case 3:
            resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), BSpline<3>());
            break;
          case 4:
            resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), BSpline<4>());
            break;
          case 5:
            resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), BSpline<5>());
            break;
          default:
            mexErrMsgTxt("Interpolation order must be between 0 and 5.");
        }
    }
}

 
/** MATLAB 
function resized = resizeImage(original, newShape)
function resized = resizeImage(original, newShape, splineOrder)

resize a 2D array to a new shape using spline interpolation of the given order 
(default: using cubic spline interpolation).
    original - 2D matrix with original image 
               (gray scale or multi-band/RGB, type 'uint8', 'single' or 'double')
    newShape - int32-array of length 2 giving the new shape
    splineOrder - order of interpolation 
                  (0 <= splineOrder <= 5, default: 3, i.e. cubic splines)
*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs)
{    
    if(matlab::ValueType<UInt8>::check(inputs[0]))
        vigraResizeImageImpl<UInt8>(outputs, inputs);
    else if(matlab::ValueType<float>::check(inputs[0]))
        vigraResizeImageImpl<float>(outputs, inputs);
    else if(matlab::ValueType<double>::check(inputs[0]))
        vigraResizeImageImpl<double>(outputs, inputs);
    else
        mexErrMsgTxt("Input image must have type 'uint8', 'single' or 'double'.");
}
