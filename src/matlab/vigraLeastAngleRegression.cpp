/***************************************************************************************************
**         INCLUDES AND DEFS                                                                      **
****************************************************************************************************/

#include <vigra/matlab.hxx>
#include <vigra/symmetry.hxx>
#include <vigra/regression.hxx>
//#include <vigra/multi_pointoperators.hxx>

using namespace vigra;
using namespace matlab;
using namespace linalg;

template <class T>
void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    /***************************************************************************************************
    **              INIT PART                                                                         **
    ****************************************************************************************************/
    typedef double OutputType;

    //A and b are references to input memory. changing them causes input changed.
    //MultiArrayView<2, T>   A_copy    =   inputs.getMultiArray<2, T>(0, v_required());
    //MultiArrayView<2, T>   b_copy    =   inputs.getMultiArray<2, T>(1, v_required());
    MultiArray<2, T>   A    =   inputs.getMultiArray<2, T>(0, v_required());
    MultiArray<2, T>   b    =   inputs.getMultiArray<2, T>(1, v_required());
    //vigra::copyMultiArray(srcMultiArrayRange(A), destMultiArray(A_copy));
    //vigra::copyMultiArray(srcMultiArrayRange(b), destMultiArray(b_copy));

    int max_solution_count      =   inputs.getScalarMinMax<int>("max_solution_count",v_default(0), 0, "inf");
    std::string mode            =   inputs.getString("mode", v_default(std::string("lasso")));


    ArrayVector<ArrayVector<int> > activeSets;
    ArrayVector<Matrix<OutputType> > lsq_solutions;
    ArrayVector<Matrix<OutputType> > lasso_solutions;
    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/

    // normalize the input

    int n = columnCount(A);
    Matrix<double> offset(1,n), scaling(1,n);
	prepareColumns(A, A, offset, scaling, linalg::DataPreparationGoals(ZeroMean|UnitVariance));
	prepareColumns(b, b, linalg::DataPreparationGoals(ZeroMean));

    unsigned int numSolutions = leastAngleRegression(A, b, activeSets, lasso_solutions, lsq_solutions,
                                                LeastAngleRegressionOptions()
                                                .maxSolutionCount(max_solution_count)
                                                .setMode(mode)
                                );

    MultiArrayView<2, OutputType> dense_lsq = outputs.createMultiArray<2, OutputType>(0,v_required(),
                                                    MultiArrayShape<2>::type(columnCount(A), numSolutions));
    MultiArrayView<2, OutputType> dense_lasso = outputs.createMultiArray<2, OutputType>(1,v_optional(),
                                                    MultiArrayShape<2>::type(columnCount(A), numSolutions));
    for (MultiArrayIndex k = 0; k < numSolutions; ++k)
    {
        for (unsigned int i = 0; i < activeSets[k].size(); ++i)
        {
            dense_lsq(activeSets[k][i], k) = lsq_solutions[k](i,0)*scaling(0, activeSets[k][i]);
        }
    }

    // Optional Output.
    if(dense_lasso.data() != 0)
    {
        for (MultiArrayIndex k = 0; k < numSolutions; ++k)
        {
            for (unsigned int i = 0; i < activeSets[k].size(); ++i)
            {
                dense_lasso(activeSets[k][i], k) = lasso_solutions[k](i,0)*scaling(0, activeSets[k][i]);
            }
        }
    }

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

    switch(inputs.typeOf(0))
    {
        ALLOW_D
        default:
	    mexErrMsgTxt("Type of input at position 0 not supported");
    }
}
/** MATLAB
function dense_lsq = vigraLeastAngleRegression(A, b)
function dense_lsq = vigraLeastAngleRegression(A, b, options);
function [dense_lsq dense_lasso] = vigraLeastAngleRegression(...);

Solves Equations of Type A*x = b using L1 regularisation
The columns of dense_lsq are the least squares Solutions in the
    (iterations?) of the Lasso
Columns of dense_lasso are...

options    - a struct with following possible fields:
    'mode':    default 'lasso', nnlasso (non negative lasso), lars (least angle regression)
    'max_solution_count':   default: Unused, Integral value > 0


*/
