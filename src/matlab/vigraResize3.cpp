/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by Rahul Nair and Ullrich Koethe          */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/



#include <vigra/matlab.hxx>
#include <vigra/matlab_FLEXTYPE.hxx>
#include <string>
#include <vigra/multi_resize.hxx>
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
    /***************************************************************************************************
    **              INIT PART                                                                         **
    ****************************************************************************************************/
    //Options

    MultiArrayView<4,T>         in4D        = inputs.getMultiArray<4,T>(0, v_required());

    VIGRA_CREATE_ENUM_AND_STD_MAP2(MapName, BSpline, Catmull)
    int method = inputs.getEnum("method", v_default(BSpline), MapName );

 //   Methods1 method2 = inputs.getEnum<Methods1>("method1", Optional(BSpline11));


    Int32                        splineOrder = (method == BSpline)?
                                          inputs.getScalarMinMax<int>("splineOrder", v_default(3),0, 5)
                                        : 0;

    TinyVector<double, 3>       defaultShape(2*in4D.shape(0) -1, 2*in4D.shape(1)-1, 2*in4D.shape(2)-1);
    TinyVectorView<double, 3>   newShape  = inputs.getTinyVector<double, 3> ( 1, v_default(defaultShape));
    MultiArrayShape<4>::type    newShape4      (Int32(newShape[0]), Int32(newShape[1]), Int32(newShape[2]), in4D.shape(3));


    MultiArrayView<4,T> out4D           = outputs.createMultiArray      <4,T>   (0, v_required(), newShape4);
    // contorPair maps 2 integers bijectively onto one dimension. (see Wikipedia Cantor pair Function)
    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/
        for(int k=0; k<in4D.shape(3); ++k)
        {
            MultiArrayView<3, T> ink = in4D.bindOuter(k);
            MultiArrayView<3, T> outk = out4D.bindOuter(k);
            switch(cantorPair(method, splineOrder))
            {
            case cP2_(BSpline, 0):
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk), vigra::BSpline<0>());
                break;
            case cP2_(BSpline, 1):
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk),  vigra::BSpline<1>());
                break;
            case cP2_(BSpline, 2):
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk),  vigra::BSpline<2>());
                break;
            case cP2_(BSpline, 3):
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk),  vigra::BSpline<3>());
                break;
            case cP2_(BSpline, 4):
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk),  vigra::BSpline<4>());
                break;
            case cP2_(BSpline, 5):
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk),  vigra::BSpline<5>());
                break;
            case cP2_(Catmull, 0):
                resizeMultiArraySplineInterpolation(srcMultiArrayRange(ink), destMultiArrayRange(outk),  vigra::CatmullRomSpline<double>());
                break;
             default:
                mexErrMsgTxt("Something went terribily wrong - Internal error");
            }
        }
}



/***************************************************************************************************
**           VIGRA GATEWAY                                                                        **
****************************************************************************************************/
void vigraMexFunction(vigra::matlab::OutputArray outputs, vigra::matlab::InputArray inputs)
{
    //Add classes as you feel
    switch(inputs.typeOf(0))
    {
        ALLOW_D
        default:
            mexErrMsgTxt("Type of input 0 not supported");
    }
}


/** MATLAB
function resized = vigraResize3(original)
function resized = vigraResize3(original, newShape)
function resized = vigraResize3(original, newShape, options)

D = vigraResize3(inputVolume)   # resizes original volume data with default options.
D = vigraResize3(inputVolume, [200 300 100], options)  # does the same with user options.

    original    - Array with original 3D volume data
                    (gray scale or multi-band/RGB, numeric type)
    newShape    - Array of length 3 that gives the new shape
                    (default: 2*size(original)-1 )
    options
        splineOrder - order of interpolation
            (0 <= splineOrder <= 5, default: 3, i.e. cubic splines)
            this option is only used for method 'BSpline'
        method - 'BSpline' (default) or 'Catmull'
*/
