/*++++++++++++++++++++INCLUDES+and+Definitions++++++++++++++++++++++++*/

#include <iostream>
#include <set>
#include <string>
#include <vigra/matlab.hxx>
#include <vigra/random_forest.hxx>
#include "random_forest_impex.hxx"



//this could be a typedef but if you want outType to be the same type as inType then you can just
//set outType to T



using namespace vigra;
using namespace matlab;

/*+++++++++++++++++++User data structure+++++++++++++++++++++++++++++*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
//This is not using base_data as base class. because input is 4D.

template <class T>
struct data
: public base_data<T>
{
    declBool(sample_with_replacement,1 );
    declBool(sample_classes_individually,0);
    declScalarMinMax(unsigned int, min_split_node_size,1, 0, "inf");
    declScalarMinMax(unsigned int, treeCount, 255, 0, "inf");
    declScalar(unsigned int, mtry, 0);
    declScalar(unsigned int, training_set_size, 0);
    declScalarMinMax(double, training_set_proportion, 1.0, 0.0, 1.0);
    declMultiArrayDim(2, labels, double);
    declMultiArrayDim(2, features, double);

    int numberOfSamples, numberOfFeatures;


    data(matlab::OutputArray outputs, matlab::InputArray inputs)
    :           base_data<T>(inputs),
                initOption(sample_with_replacement),
                initOption(sample_classes_individually),
                initOption(min_split_node_size),
                initOption(treeCount),
                initOption(mtry),
                initOption(training_set_size),
                initOption(training_set_proportion),
                initOption(labels),
                initOption(features)
    {
        if (features.shape(0) != labels.shape(0))
            mexErrMsgTxt("Sample number in labels and feature matrix must be equal");
    }
};



struct vigraFunctor
{
    template <class T>
    static void exec(matlab::OutputArray outputs, matlab::InputArray inputs){
        //Options
        data<T>  o(outputs, inputs);

        std::set<double> labelSet;
        for(int i=0; i<o.labels.size(); ++i)
            labelSet.insert(o.labels[i]);


        RandomForest<double> rf(labelSet.begin(), labelSet.end(),
                                                vigra::RandomForestOptions()\
                                                .trainingSetSizeProportional(o.training_set_proportion)\
                                                .featuresPerNode(o.mtry)\
                                                .minSplitNodeSize(o.min_split_node_size)\
                                                .trainingSetSizeAbsolute(o.training_set_size)\
                                                .sampleWithReplacement(o.sample_with_replacement)\
                                                .sampleClassesIndividually(o.sample_classes_individually)\
                                                .setTreeCount(o.treeCount)
                                );
        double oobError = rf.learn(o.features, o.labels);

        matlab::exportRandomForest(rf, matlab::createCellArray(2*o.treeCount+2, outputs[0]));

        if(outputs.isValid(1))
            outputs[1] = matlab::createScalar<double>(oobError);

    }
};
