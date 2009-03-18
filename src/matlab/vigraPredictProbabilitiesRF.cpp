#include <iostream>
#include <set>
#include <vigra/matlab.hxx>
#include "random_forest_impex.hxx"

using namespace vigra;
using namespace matlab;



void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    /* INPUT */
    if (inputs.size() != 2)
        mexErrMsgTxt("Two inputs required.");

    // get RF object
    std::auto_ptr<RandomForest<double> > rf =
                  matlab::importRandomForest<double>(matlab::getCellArray(inputs[0]));

    // get feature matrix
    MultiArrayView<2, double> features = inputs.getMultiArray<2, double> ( 1, v_required());
    if(rf->featureCount() != columnCount(features))
        mexErrMsgTxt("Feature array has wrong number of columns.");

    /* OUTPUT */
    MultiArrayView<2, double> probs = outputs.createMultiArray<2, double>(0, v_required(),
                                                                TinyVector<UInt32, 2>(rowCount(features), rf->labelCount()));

    rf->predictProbabilities(features, probs);

}




/***************************************************************************************************
**         VIGRA GATEWAY                                                                          **
****************************************************************************************************/
inline void vigraMexFunction(vigra::matlab::OutputArray outputs, vigra::matlab::InputArray inputs)
{
    vigraMain(outputs, inputs);
}
/** Matlab
function probs = predictProbabilitiesRF(RF, features)

Use a previously trained random forest classifier to predict labels for the given data
    RF        - MATLAB cell array representing the random forest classifier
    features  - M x N matrix, where M is the number of samples, N the number of features

    probs     - M x L matrix holding the predicted probabilities for each of
                the L possible labels
*/
