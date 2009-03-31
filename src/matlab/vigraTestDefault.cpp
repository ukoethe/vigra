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
    MultiArrayShape<2>::type    newShape3 (3, 3);
    MultiArray<2,double> defMultArray(newShape3);
    MultiArrayView<2,double>         in3D   = inputs.getMultiArray<2,double>(2, v_default(MultiArrayView<2, double>(defMultArray)));
    double                   aga            = inputs.getScalar<double>(0, v_default(2));
    double                   sgs            = inputs.getScalar<double>(1, v_default(2, 3, int(aga) ));

    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/

    MultiArrayView<2,double>   out3D       = outputs.createMultiArray      <2,double>   (2, v_required(), in3D.shape());
    vigra::copyMultiArray(srcMultiArrayRange(in3D), destMultiArray(out3D));

    outputs.createScalar<double>(0, v_required(), aga);
    outputs.createScalar<double>(1, v_required(), sgs);
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

