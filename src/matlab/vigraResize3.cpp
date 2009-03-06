/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/multi_resize.hxx>



//this could be a typedef but if you want outType to be the same type as inType then you can just 
//set outType to T



using namespace vigra;
using namespace matlab;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//This is not using base_data as base class. because input is 4D.
template <class T>
struct data
: public base_data<T>
{
    MultiArrayView<4, T> in; 
    MultiArrayView<4, T> out;
    
    declCharConstr2(method, BSpline, CatmullRom);
    declScalarMinMax(int, splineOrder, 3, 0, 5);

    
    data(matlab::OutputArray outputs, matlab::InputArray inputs)
    :           base_data<T>(inputs),
                in(matlab::getMultiArray<4, T>(inputs[0])),
                initOption(method),
                initOption(splineOrder)
    {
        typedef MultiArrayShape<3>::type Shape;
        
        // check if newShape was given
        bool newShapeGiven = inputs.size() == 3 ||
                             (inputs.size() == 2 && !this->options.isValid());                                
        Shape newShape = newShapeGiven
                             ? matlab::getShape<3>(inputs[1])
                             : Shape(in.shape(0)*2-1, in.shape(1)*2-1, in.shape(2)*2-1);
        MultiArrayShape<4>::type newShape4(newShape[0], newShape[1], newShape[2], in.shape(3));
        
        out = matlab::createMultiArray<4,T>(newShape4,outputs[0]);
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

        for(int k=0; k<o.in.shape(3); ++k)
        {
            MultiArrayView<3, T> ink = o.in.bindOuter(k);    
            MultiArrayView<3, T> outk = o.out.bindOuter(k);
            switch(cantorPair(o.method, o.splineOrder))
            {
            case cP2_(BSpline, 0):
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<0>());
                break;
            case cP2_(BSpline, 1):
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<1>());
                break;
            case cP2_(BSpline, 2):
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<2>());
                break;
            case cP2_(BSpline, 3):
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<3>());
                break;
            case cP2_(BSpline, 4):
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<4>());
                break;
            case cP2_(BSpline, 5):
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), BSpline<5>());
                break;
            case cP2_(CatmullRom, 0):
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), CatmullRomSpline<double>());
                break;
             default:
                mexErrMsgTxt("Something went wrong");
            }
        }
        
    }
};


/*+++++++++++++++++++++++MexEntryFunc++++++++++++++++++++++++++++++++*/
/* Gatewayfunction - see matlab.hxx for details.
/* if a certain class is NOT supported - you will have to copy the 
/* body of the callMexFunctor function and edit it here.
/* Supports (u)int[8|16|32|64], float and double.
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/** MATLAB 
function resized = vigraResize3(original)
function resized = vigraResize3(original, newShape)
function resized = vigraResize3(original, options)
function resized = vigraResize3(original, newShape, options)

D = vigraResize3(inputVolume)   # resizes original volume data with default options.
D = vigraResize3(inputVolume, [200 300 100], options)  # does the same with user options.

    original    - Array with original 3D volume data 
                    (gray scale or multi-band/RGB, numeric type)
    newShape    - int32-Array of length 3 that gives the new shape
                    (default: 2*size(original)-1 )
    options
        splineOrder - order of interpolation 
            (0 <= splineOrder <= 5, default: 3, i.e. cubic splines)
            this option is only used for method 'BSpline'
        method - 'BSpline' (default) or 'CatmullRom'
*/

