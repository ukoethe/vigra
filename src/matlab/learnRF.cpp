#include <iostream>
#include <set>
#include <vigra/matlab.hxx>
#include <vigra/random_forest.hxx>

using namespace vigra;

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

    /* OUTPUT */
    matlab::CellArray cells = matlab::createCellArray(2*Ntree+2, outputs[0]);

    // write RF parameters
    int parameterCount = 3;
    MultiArrayView<1, UInt32> parameters = matlab::createArray<UInt32>(parameterCount, cells[0]);
    parameters(0) = labelCount;
    parameters(1) = numberOfFeatures;
    parameters(2) = Ntree;
     
    // write array of possible class labels
    MultiArrayView<1, double> classLabels = matlab::createArray<double>(labelCount, cells[1]);
    for(unsigned int k =0; k<labelCount; ++k)
        classLabels(k) = rf.classes_[k];
    
    // for all decision trees
    for(unsigned int k=0; k<Ntree; ++k)
    {
        // write int tree array
        MultiArrayView<1, Int32> tree = 
            matlab::createArray<Int32>(rf.trees_[k].tree_.size(), cells[2*k+2]);
        for(unsigned int i =0; i<rf.trees_[k].tree_.size(); ++i)
            tree(i) = rf.trees_[k].tree_[i];
        
        // write double weight/threshold array
        MultiArrayView<1, double> weights = 
            matlab::createArray<double>(rf.trees_[k].terminalWeights_.size(), cells[2*k+3]);
        for(unsigned int i =0; i<rf.trees_[k].terminalWeights_.size(); ++i)
            weights(i) = rf.trees_[k].terminalWeights_[i];
    }
}
