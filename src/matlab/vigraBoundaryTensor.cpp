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

/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <vigra/matlab.hxx>
#include <string>
#include <vigra/boundarytensor.hxx>


using namespace vigra;
using namespace matlab;





template <class T>
void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    /***************************************************************************************************
    **              INIT PART                                                                         **
    ****************************************************************************************************/
    BasicImageView<T>   in      =   inputs.getImage<T>(0, v_required());
    double              scale   =   inputs.getScalarMinMax<double>(1, v_default(1.0), 0.0, "inf");

    MultiArrayView<3,T> res     =   outputs.createMultiArray<3,T>   (0, v_required(), 
                                                   MultiArrayShape<3>::type(3, in.width(), in.height()));
    
    vigra_precondition(sizeof(TinyVector<T, 3>) == sizeof(T)*res.stride(1),
           "vigraBoundaryTensor(): Internal error (unsuitable memory layout).");

    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/

    BasicImageView<TinyVector<T, 3> >  out(reinterpret_cast<TinyVector<T, 3> *>(res.data()),
                                            in.width(), in.height());
    boundaryTensor(srcImageRange(in), destImage(out), scale);
}



/***************************************************************************************************
**         VIGRA GATEWAY                                                                          **
****************************************************************************************************/
void vigraMexFunction(vigra::matlab::OutputArray outputs, vigra::matlab::InputArray inputs)
{
    //Add classes as you feel
    switch(inputs.typeOf(0))
    {
        ALLOW_FD
        default:
            mexErrMsgTxt("Type of input 0 not supported");
    }
}


/** MATLAB
function out = vigraBoundaryTensor(inputImage)
function out = vigraBoundaryTensor(inputImage, scale);

inputImage - 2D scalar input array
scale      - 1.0 (default), a positive floating point scale
             parameter for boundary tensor computation
             
out        - output boundary tensor image
             the first dimension holds the boundary tensor entries
                 B11(y,x) = out(1,y,x)
                 B21(y,x) = B12(y,x) = out(2,y,x)
                 B22(y,x) = out(3,y,x)

Usage:
    out = vigraBoundaryTensor(in, 2.0);
*/
