/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/resizeimage.hxx>



//this could be a typedef but if you want outType to be the same type as inType then you can just 
//set outType to T

#define vigraFunctor vigraResize3

using namespace vigra;
using namespace matlab;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//This is not using base_data as base class. because input is 4D.

template <class T>
struct data
: public base_data<T>
{
    declCharConstr3(method, BSpline, CatmullRom, Coscot);
    declScalarMinMax(int, splineOrder, 3, 0, 5);
    declOut(T);
    
    using base_data<T>::in3D;
    
    data(matlab::OutputArray outputs, matlab::InputArray inputs)
    :           base_data<T>(inputs),
                map(method),
                map(splineOrder)
    {
        typedef MultiArrayShape<2>::type Shape;
        
        // check if newShape was given
        bool newShapeGiven = inputs.size() == 3 ||
                             (inputs.size() == 2 && !this->options.isValid());                                
        Shape newShape = newShapeGiven
                             ? matlab::getShape<2>(inputs[1])
                             : Shape(in3D.shape(0)*2-1, in3D.shape(1)*2-1);
        mapOut_3D(T, newShape[0], newShape[1], in3D.shape(2));
    }
};

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* This function does all the work
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cP2_(a, b) cP<data<T>::a, b>::value
struct vigraFunctor
{
    template <class T>
    static void exec(matlab::OutputArray outputs, matlab::InputArray inputs){
        //Options
        data<T>  o(outputs, inputs);

        // contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function) 

        for(int k=0; k<o.in3D.shape(3); ++k)
        {
            BasicImageView<T> ink = makeBasicImageView(o.in3D.bindOuter(k));    
            BasicImageView<T> outk = makeBasicImageView(o.out3D.bindOuter(k));    
            switch(cantorPair(o.method, o.splineOrder))
            {
            case cP2_(BSpline, 0):
                resizeImageNoInterpolation(srcImageRange(ink), destImageRange(outk));
                break;
            case cP2_(BSpline, 1):
                resizeImageLinearInterpolation(srcImageRange(ink), destImageRange(outk));
                break;
            case cP2_(BSpline, 2):
                resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), BSpline<2>());
                break;
            case cP2_(BSpline, 3):
                resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), BSpline<3>());
                break;
            case cP2_(BSpline, 4):
                resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), BSpline<4>());
                break;
            case cP2_(BSpline, 5):
                resizeImageSplineInterpolation(srcImageRange(ink), destImageRange(outk), BSpline<5>());
                break;
            case cP2_(CatmullRom, 0):
                resizeImageCatmullRomInterpolation(srcImageRange(ink), destImageRange(outk));
                break;
            case cP2_(Coscot, 0):
                resizeImageCoscotInterpolation(srcImageRange(ink), destImageRange(outk));
                break;          
            default:
                mexErrMsgTxt("Something went wrong");
            }
        }
        
    }
};

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
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs){
    // 
    callMexFunctor<vigraFunctor>(outputs, inputs);
}
