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

    for(int k=0; k<in.shape(2); ++k)
    {
        BasicImageView<T> ink = makeBasicImageView(in.bindOuter(k));    
        BasicImageView<T> outk = makeBasicImageView(out.bindOuter(k));    
        resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk));
    }
}

 
/** MATLAB 
function resized = resizeImage(original, newShape)

resize a 2D array to a new shape using cubic spline interpolation.
    original - 2D matrix with original image 
               (gray scale or multi-band/RGB, type 'uint8', 'single' or 'double')
    newShape - int32-array of length 2 giving the new shape
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
