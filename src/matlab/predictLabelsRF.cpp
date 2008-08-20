#include <iostream>
#include <set>
#include <vigra/matlab.hxx>
#include "create_random_forest.hxx"

using namespace vigra;

void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs)
{    
    /* INPUT */
    if (inputs.size() != 2)
        mexErrMsgTxt("Two inputs required.");
    if (outputs.size() > 1)
        mexErrMsgTxt("Too many output arguments");

    // get RF object
    std::auto_ptr<RandomForest<double> > rf =
                  matlab::createRandomForest<double>(matlab::getCellArray(inputs[0]));

    // get feature matrix
    MultiArrayView<2, double> features = matlab::getMatrix<double>(inputs[1]);
    if(rf->featureCount() != columnCount(features))
        mexErrMsgTxt("Feature array has wrong number of columns.");

    /* OUTPUT */
    MultiArrayView<2, double> labels = 
           matlab::createMatrix<double>(rowCount(features), 1, outputs[0]);
    
    rf->predictLabels(features, labels);
}
