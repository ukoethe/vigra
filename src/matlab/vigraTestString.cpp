/***************************************************************************************************
**         INCLUDES AND DEFS                                                                      **
****************************************************************************************************/

#include <vigra/matlab.hxx>
#include <string>
using namespace vigra;
using namespace matlab;



void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    /***************************************************************************************************
    **              INIT PART                                                                         **
    ****************************************************************************************************/

    std::string   method      = inputs.getString(0, v_default("bla"));

    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/
    int out = 0;

    if(method == std::string("bla"))
        out = 1;
    else if (method == std::string("user"))
        out = 2;

    outputs.createScalar<double>(0, v_required(), out);


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

