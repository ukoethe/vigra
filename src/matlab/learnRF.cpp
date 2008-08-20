#include <iostream>
#include <set>
#include <vigra/matlab.hxx>
#include <vigra/random_forest.hxx>
#include "random_forest_impex.hxx"

using namespace vigra;

/** matlab 
function RF = learnRF(features, labels, treeCount)

Train a random forest classifier for the given data
    features  - M x N matrix, where M is the number of samples, N the number of features
    labels    - M x 1 matrix holding the true labels
    treeCount - number of trees to be used in the RF classifier
    
    RF        - MATLAB cell array representing the random forest classifier
*/
void vigraMexFunction(matlab::OutputArray outputs, matlab::InputArray inputs)
{    
    if (inputs.size() != 3)
        mexErrMsgTxt("Three inputs required.");
    if (outputs.size() > 1)
        mexErrMsgTxt("Too many output arguments");
                     
    MultiArrayView<2, double> features = matlab::getMatrix<double>(inputs[0]);
    MultiArrayView<2, double> labels   = matlab::getMatrix<double>(inputs[1]);
    if (features.shape(0) != labels.shape(0))
        mexErrMsgTxt("Sample number in labels and feature matrix must be equal");

    int numberOfSamples  = features.shape(0);
    int numberOfFeatures = features.shape(1);
    
    UInt32 Ntree = matlab::getScalar<UInt32>(inputs[2]);
    
    std::set<double> labelSet;
    for(int i=0; i<numberOfSamples; ++i)
        labelSet.insert(labels[i]);  
    std::size_t labelCount = labelSet.size();
    
    RandomForest<double> rf(labelSet.begin(), labelSet.end(), Ntree);
    
    rf.learn(features, labels);

    // OUTPUT
    matlab::exportRandomForest(rf, matlab::createCellArray(2*Ntree+2, outputs[0]));
}
