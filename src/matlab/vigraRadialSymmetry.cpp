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

/***************************************************************************************************
**         INCLUDES AND DEFS                                                                      **
****************************************************************************************************/

#include <vigra/matlab.hxx>
#include <vigra/symmetry.hxx>

using namespace vigra;
using namespace matlab;


template <class T>
void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    /***************************************************************************************************
    **              INIT PART                                                                         **
    ****************************************************************************************************/
    BasicImageView<T>   in      =   inputs.getImage<T>(0, v_required());
    double              scale   =   inputs.getScalarMinMax<double>(1, v_default(1.0), 0.0, "inf");

    BasicImageView<double> out  =   outputs.createImage<double>(0, v_required(), in.width(), in.height());

    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/
    radialSymmetryTransform(srcImageRange(in), destImage(out), scale);

}

/***************************************************************************************************
**         VIGRA GATEWAY                                                                          **
****************************************************************************************************/
void vigraMexFunction(vigra::matlab::OutputArray outputs, vigra::matlab::InputArray inputs)
{
    //Add classes as you feel
    switch(inputs.typeOf(0))
    {
        ALLOW_D
        default:
            mexErrMsgTxt("Type of input at position 0 not supported");
    }
}
/** MATLAB
function D = vigraRadialSymmetry(inputImage)
function D = vigraradialSymmetry(inputImage, options);

D = vigraRadialSymmetry(inputImage) computes the Fast Radial Symmetry Transform
            using default options, see vigra::RadialSymmetryTransform for more information.
D = vigraRadialSymmetry(inputImage, options)  does the same with user options.

inputImage - 2D input array
options    - a struct with following possible fields:
    'scale':    1.0 (default), any positive floating point value
                scale parameter for the vigraRadialSymmetry


Usage:
    opt = struct('method' ,value);
    out = vigraRadialSymmetry(in, opt);

*/
