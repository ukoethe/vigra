/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by Rahul Nair and Ullrich Koethe          */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
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

/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

//#define VIGRA_CHECK_BOUNDS
#include <vigra/matlab.hxx>
#include <vigra/labelimage.hxx>
#include <vigra/labelvolume.hxx>
#include <vigra/matlab_FLEXTYPE.hxx>


//this could be a typedef but if you want outType to be the same type as inType then you can just
//set outType to T


using namespace vigra;
using namespace matlab;



//#define RN_DEBUG
#define cP3_(a, b , c) cP3<a, b, c>::value
template <class T>
void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){

    /***************************************************************************************************
    **              INIT PART                                                                         **
    ***************************************************************************************************/
    typedef UInt32 OutputType;
    MultiArrayView<3,T> in3D        =       inputs.getMultiArray<3,T>(0, v_required());
    BasicImageView<T>   in          =       makeBasicImageView(in3D.bindOuter(0));
    Int32               numOfDim    =       inputs.getDimOfInput(0, v_required());

    Int32               v2Dconn[2]  = {8, 4};
    Int32               v3Dconn[2]  = {26, 6};
    Int32               connectivity= inputs.getScalarVals2D3D<Int32>("conn", 
                                                                      numOfDim == 2 ? v_default(8) : v_default(26),
                                                                      v2Dconn, v2Dconn+2,
                                                                      v3Dconn, v3Dconn+2,
                                                                      numOfDim);
    /*{
        if(numOfDim == 2 && connectivity!= 8 && connectivity != 4)
            mexErrMsgTxt("Connectivity for 2D data must be 8 or 4");
        else if(numOfDim == 3 && connectivity!= 26 && connectivity != 6)
            mexErrMsgTxt("Connectivity for 3D data must be 26 or 6");
    }*/


    bool                    hasBackground;
    T                       backgroundValue = inputs.getScalar<T>("backgroundValue", v_optional(hasBackground));

    MultiArrayView<3,OutputType>     out3D  = outputs.createMultiArray      <3,OutputType>   (0, v_required(), in3D.shape());
    BasicImageView<OutputType>       out(out3D.data(), in3D.shape(0), in3D.shape(1));

    Int32                     max_region_label = (hasBackground == true)? 1: 0;

    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/
    #ifdef RN_DEBUG
    mexPrintf("---%d---%d---%d---", max_region_label, numOfDim, connectivity);
    #endif
    switch(cantorPair(hasBackground , numOfDim, connectivity)){
        //cP is the templated version o f the cantorPair function first value is Dimension of Inputimage, second the connectivity setting
        //Code is basically the code on the VIGRA-reference page
        case cP3_(0, IMAGE, 8):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("0, IMAGE, 8");
            #endif
            max_region_label = labelImage(srcImageRange(in), destImage(out), true);
            break;
        case cP3_(0, IMAGE, 4):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("0, IMAGE, 4");
            #endif
            max_region_label = labelImage(srcImageRange(in), destImage(out), false);
            break;
        case cP3_(0, VOLUME, 6):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("0, VOLUME, 26");
            #endif
            max_region_label = labelVolumeSix(srcMultiArrayRange(in3D), destMultiArray(out3D));
            break;
        case cP3_(0, VOLUME, 26):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("0, VOLUME, 6");
            #endif
            max_region_label = labelVolume(srcMultiArrayRange(in3D), destMultiArray(out3D), NeighborCode3DTwentySix());
            break;
        case cP3_(1, IMAGE, 8):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("1, IMAGE, 8");
            #endif
            max_region_label = labelImageWithBackground(srcImageRange(in), destImage(out), true, backgroundValue);
            break;
        case cP3_(1, IMAGE, 4):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("1, IMAGE, 4");
            #endif
            max_region_label = labelImageWithBackground(srcImageRange(in), destImage(out), false, backgroundValue);
            break;
        case cP3_(1, VOLUME, 26):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("1, VOLUME, 26");
            #endif
            max_region_label = labelVolumeWithBackground(srcMultiArrayRange(in3D), destMultiArray(out3D),
                                                                                        NeighborCode3DTwentySix(), backgroundValue);
            break;
        case cP3_(1, VOLUME, 6):
            #ifdef RN_DEBUG
            mexWarnMsgTxt("1, VOLUME, 6");
            #endif
            max_region_label = labelVolumeWithBackground(srcMultiArrayRange(in3D), destMultiArray(out3D),
                                                                                        NeighborCode3DSix(), backgroundValue);
            #ifdef RN_DEBUG
            mexWarnMsgTxt("DONE");
            #endif
            break;
        default:
            mexErrMsgTxt("Something went wrong");
    }

    outputs.createScalar<int> (1, v_optional(), max_region_label);
}



/***************************************************************************************************
**           VIGRA GATEWAY                                                                        **
****************************************************************************************************/
void vigraMexFunction(vigra::matlab::OutputArray outputs, vigra::matlab::InputArray inputs)
{
    //Add classes as you feel
    switch(inputs.typeOf(0))
    {
        ALLOW_FD
        ALLOW_UINT_8_64
        ALLOW_INT_8_64
        default:
            mexErrMsgTxt("Type of input 0 not supported");
    }
}


/*+++++++++++++++++++++++MexEntryFunc++++++++++++++++++++++++++++++++*/
/* Gatewayfunction - see matlab.hxx for details.
/* if a certain class is NOT supported - you will have to copy the
/* body of the callMexFunctor function and edit it here.
/* Supports (u)int[8|16|32|64], float and double.
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/** MATLAB
function D = vigraConnectedComponents(inputArray)
function D = vigraConnectedComponents(inputArray, options);

D = vigraConnectedComponents(inputArray) performs connected components labeling
        using the default options.
D = vigraConnectedComponents(inputImage, options)  does the same with user options.

inputArray - a 2D or 3D array of numeric type

options    - is a struct with the following possible fields:
    'conn':            The neighborhood to be used
                       2D: 4 (default) or  8
                       3D: 6 (default) or 26
    'backgroundValue': Specify the value of a background region not to be labeled (will be labeled 0
                       in the result array, even when it doesn't form a single connected component).
                       If this option in not present, the entire image/volume will be labeled.

Usage:
    opt = struct('conn', 8);
    out = vigraConnectedComponents(in, opt);

*/
