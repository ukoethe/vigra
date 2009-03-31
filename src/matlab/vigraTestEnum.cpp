/***************************************************************************************************
**         INCLUDES AND DEFS                                                                      **
****************************************************************************************************/

#include <vigra/matlab.hxx>
using namespace vigra;
using namespace matlab;



void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    /***************************************************************************************************
    **              INIT PART                                                                         **
    ****************************************************************************************************/

    VIGRA_CREATE_ENUM_AND_STD_MAP3(Methods,MapName, first, second, third);
    double        method      = double(inputs.getEnum(0, v_required() , MapName ));

    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/

    outputs.createScalar<double>(0, v_required(), method);


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

