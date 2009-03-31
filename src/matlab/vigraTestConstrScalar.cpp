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

    bool        Bool        =   inputs.getBool(0, v_required());
    double      MinMaxS     =   inputs.getScalarMinMax<double>(1, v_required(), 2.0, 3.0);
    double      vals1[]     =   {2, 4};
    double      vals2[]     =   {3, 5};
    double      ValsC       =   inputs.getScalarVals<double>(2, v_required(), vals1, vals1+2);
    double      Vals2D3DC   =   inputs.getScalarVals2D3D<double>(3, v_required(), vals2, vals2 +2,
                                                                                vals1, vals1 +2, MinMaxS);

    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/



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

