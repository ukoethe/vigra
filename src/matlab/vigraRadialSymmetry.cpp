/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

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
    /*
    FLEXIBLE_TYPE_START(0, in);
        ALLOW_D;
    FLEXIBLE_TYPE_END;
    */
    //Add classes as you feel
    mxClassID inClass;
    FLEX_TYPE(inClass, 0, in);
    switch(inClass)
    {
        ALLOW_D
        DEFAULT_ERROR;
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
