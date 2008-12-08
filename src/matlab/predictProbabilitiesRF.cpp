#include <iostream>
#include <set>
#include <vigra/matlab.hxx>
#include "random_forest_impex.hxx"

using namespace vigra;
using namespace matlab;

/** Matlab 
function probs = predictProbabilitiesRF(RF, features)

Use a previously trained random forest classifier to predict labels for the given data
    RF        - MATLAB cell array representing the random forest classifier
    features  - M x N matrix, where M is the number of samples, N the number of features

    probs     - M x L matrix holding the predicted probabilities for each of 
                the L possible labels
*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs)
{    
    /* INPUT */
    if (inputs.size() != 2)
        mexErrMsgTxt("Two inputs required.");
    if (outputs.size() > 1)
        mexErrMsgTxt("Too many output arguments");

    // get RF object
    std::auto_ptr<RandomForest<double> > rf =
                  matlab::importRandomForest<double>(matlab::getCellArray(inputs[0]));

    // get feature matrix
    MultiArrayView<2, double> features = matlab::getMatrix<double>(inputs[1]);
    if(rf->featureCount() != columnCount(features))
        mexErrMsgTxt("Feature array has wrong number of columns.");

    /* OUTPUT */
    MultiArrayView<2, double> probs = 
           matlab::createMatrix<double>(rowCount(features), rf->labelCount(), outputs[0]);
    
    rf->predictProbabilities(features, probs);
}
