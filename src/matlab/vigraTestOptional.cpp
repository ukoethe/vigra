/***************************************************************************************************
**         INCLUDES AND DEFS                                                                      **
****************************************************************************************************/

#include <vigra/matlab.hxx>
#include <vigra/multi_pointoperators.hxx>

using namespace vigra;
using namespace matlab;



void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    /***************************************************************************************************
    **              INIT PART                                                                         **
    ****************************************************************************************************/


    MultiArrayView<2,double>         in3D   = inputs.getMultiArray<2,double>(1, v_optional());
    bool agaIsSet;
    double                   aga            = inputs.getScalar<double>(0, v_optional(agaIsSet));

    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/

    MultiArrayShape<2>::type    newShape3 (3, 3);
    if(in3D.data() != 0)
    {
        MultiArrayView<2,double>   out3D       = outputs.createMultiArray      <2,double>   (1, v_optional(), in3D.shape());
        if(out3D.data() != 0)
            vigra::copyMultiArray(srcMultiArrayRange(in3D), destMultiArray(out3D));
    }
    else
    {
        MultiArrayView<2,double>   out3D       = outputs.createMultiArray      <2,double>   (1, v_optional(), newShape3);
    }


    if(agaIsSet)
        outputs.createScalar<double>(0, v_optional(), aga);
    else
        outputs.createScalar<double>(0, v_optional(), 2);
}

/***************************************************************************************************
**         VIGRA GATEWAY                                                                          **
****************************************************************************************************/
void vigraMexFunction(vigra::matlab::OutputArray outputs, vigra::matlab::InputArray inputs)
{
    vigraMain(outputs, inputs);
}
/** MATLAB
NODOC

*/

