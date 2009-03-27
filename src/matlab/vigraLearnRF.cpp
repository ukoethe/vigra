/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <iostream>
#include <set>
#include <string>
#include <vigra/matlab.hxx>
#include <vigra/random_forest.hxx>
#include "random_forest_impex.hxx"


using namespace vigra;
using namespace matlab;


void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    /***************************************************************************************************
    **              INIT PART                                                                         **
    ****************************************************************************************************/
    typedef double inputType;

    bool sample_with_replacement = inputs.getBool("sample_with_replacement", v_default(true));
    bool sample_classes_individually = inputs.getBool("sample_classes_individually", v_default(false));
    UInt32 min_split_node_size = UInt32(inputs.getScalarMinMax<double>("min_split_node_size",v_default(1.0), 0.0, "inf"));
    UInt32 treeCount = UInt32(inputs.getScalarMinMax<double>(2,v_default(255.0), 0.0, "inf"))/*treeCount*/;
    UInt32 mtry = UInt32(inputs.getScalar<double>("mtry",v_default(0)));
    UInt32 training_set_size = UInt32(inputs.getScalar<double>("training_set_size",v_default(0)));
    double training_set_proportion = inputs.getScalarMinMax<double>("training_set_proportion",v_default(1.0), 0.0, 1.0);

    MultiArrayView<2, inputType>  labels = inputs.getMultiArray<2, inputType>(1, v_required());
    MultiArrayView<2, inputType>  features = inputs.getMultiArray<2, inputType>(0, v_required());
    MultiArrayView<1, inputType>  weights = inputs.getMultiArray<1, inputType>("weights", v_optional());
    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/
    std::set<double> labelSet;
    for(int i=0; i<labels.size(); ++i)
        labelSet.insert(labels[i]);


    RandomForest<double> rf(labelSet.begin(), labelSet.end(),
                                            vigra::RandomForestOptions()\
                                            .trainingSetSizeProportional(training_set_proportion)\
                                            .featuresPerNode(mtry)\
                                            .minSplitNodeSize(min_split_node_size)\
                                            .trainingSetSizeAbsolute(training_set_size)\
                                            .sampleWithReplacement(sample_with_replacement)\
                                            .sampleClassesIndividually(sample_classes_individually)\
                                            .setTreeCount(treeCount)
                                            .weights(weights.data(), labelSet.size())
                            );
    double oobError = rf.learn(features, labels);

    matlab::exportRandomForest(rf, matlab::createCellArray(2*treeCount+2, outputs[0]));

    outputs.createScalar<double> (1, v_optional(), oobError);

}



/***************************************************************************************************
**         VIGRA GATEWAY                                                                          **
****************************************************************************************************/
inline void vigraMexFunction(vigra::matlab::OutputArray outputs, vigra::matlab::InputArray inputs)
{
    vigraMain(outputs, inputs);
}
