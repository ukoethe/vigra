/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <vigra/matlab_FLEXTYPE.hxx>
#include <string>
#include <vigra/resizeimage.hxx>
#include <iostream>


//this could be a typedef but if you want outType to be the same type as inType then you can just
//set outType to T



using namespace vigra;
using namespace matlab;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//This is not using base_data as base class. because input is 4D.



/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* This function does all the work
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


//#define RN_DEBUG
#define cP2_(a, b) cP<(int)a, b>::value
template <class T>
void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    //Options

    MultiArrayView<3,T>         in3D        = inputs.getMultiArray<3,T>(0, Required());

    LOAD_ENUM_OPTION3(method, BSpline,Catmull, Coscot);

    int                         splineOrder = (method == BSpline)?
                                          inputs.getScalarMinMax<int>("splineOrder", Optional(3),0, 5)
                                        : 0;

    TinyVector<double, 2>       defaultShape(2* (in3D.shape(0)), 2*(in3D.shape(1)));
    TinyVectorView<double, 2>   newShape  = inputs.getTinyVector<double, 2> ( 1, Optional(defaultShape));
    MultiArrayShape<3>::type    newShape3      (newShape[0], newShape[1], in3D.shape(2));


    MultiArrayView<3,T> out3D           = outputs.createMultiArray      <3,T>   (0, Required(), newShape3);
    // contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function)

    for(int k=0; k<in3D.shape(2); ++k)
    {

        BasicImageView<T> ink = makeBasicImageView(in3D.bindOuter(k));
        BasicImageView<T> outk = makeBasicImageView(out3D.bindOuter(k));

        switch(cantorPair(method, splineOrder))
        {
        case cP2_(BSpline, 0):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("BSpline 0");
            #endif
            resizeImageNoInterpolation(srcImageRange(ink), destImageRange(outk));
            break;
        case cP2_(BSpline, 1):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("BSpline 1");
            #endif
            resizeImageLinearInterpolation(srcImageRange(ink), destImageRange(outk));
            break;
        case cP2_(BSpline, 2):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("BSpline 2");
            #endif
            resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), vigra::BSpline<2>());
            break;
        case cP2_(BSpline, 3):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("BSpline 3");
            #endif
            resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), vigra::BSpline<3>());
            break;
        case cP2_(BSpline, 4):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("BSpline 4");
            #endif
            resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), vigra::BSpline<4>());
            break;
        case cP2_(BSpline, 5):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("BSpline 5");
            #endif
            resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), vigra::BSpline<5>());
            break;
        case cP2_(Catmull, 0):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("Catmull 0");
            #endif
            resizeImageCatmullRomInterpolation(srcImageRange(ink), destImageRange(outk));
            break;
        case cP2_(Coscot, 0):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("Coscot 0");
            #endif
            resizeImageCoscotInterpolation(srcImageRange(ink), destImageRange(outk));
            break;
        default:
            mexErrMsgTxt("Something went wrong");
        }
    }

}

void vigraMexFunction(vigra::matlab::OutputArray outputs, vigra::matlab::InputArray inputs)
{
    /*
    FLEXIBLE_TYPE_START(0, in);
        ALLOW_D;
    FLEXIBLE_TYPE_END;
    */

    mxClassID inClass;
    FLEX_TYPE(inClass, 0, in);
    switch(inClass)
    {
        ALLOW_D
        DEFAULT_ERROR;
    }
}


/** MATLAB
function resized = vigraResize2(original)
function resized = vigraResize2(original, newShape)
function resized = vigraResize2(original, options)
function resized = vigraResize2(original, newShape, options)

D = vigraResize2(inputImage)   # resizes original image data with default options.
D = vigraResize2(inputImage, [200 300], options)  # does the same with user options.

    original    - Array with original 2D image data
                    (gray scale or multi-band/RGB, numeric type)
    newShape    - int32-Array of length 2 that gives the new shape
                    (default: 2*size(original)-1 )
    options
        splineOrder - order of interpolation
            (0 <= splineOrder <= 5, default: 3, i.e. cubic splines)
            this option is only used for method 'BSpline'
        method - 'BSpline' (default), 'Coscot' or 'CatmullRom'
*/
