/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by Rahul Nair and Ullrich Koethe          */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/

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
    MultiArrayView<2, UInt8> oob_data =
                        outputs.createMultiArray<2, UInt8>(2, v_optional(),
                                                            MultiArrayShape<2>::type(features.shape(0), treeCount));
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
                                            .setTreeCount(treeCount)\
                                            .weights(weights.data(), labelSet.size())\
                                            .oobData(oob_data)
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


/** MATLAB
function RF = vigraLearnRF(features, labels) Trains a randomForest with Default TreeCount and options
function RF = vigraLearnRF(features, labels, treeCount)  does the same treeCount number of trees and default options.
function RF = vigraLearnRF(features, labels, treeCount, options)  does the same with user options.
function [RF oob] = vigraLearnRF(...)                Outputs the oob error estimate
function [RF oob oob_data] = vigraLearnRF(...)       Outputs additional oob data.

features    - A Nxp Matrix with N samples containing p features
labels      - A Nx1 Matrix with the corresponding Training labels
treeCount   - default: 255. An Integral Scalar Value > 0 - Number of Trees to be used in the Random Forest.
options     - a struct with the following possible fields (default will be used
              if field is not present)
    'sample_with_replacement'       logical, default : true
    'sample_classes_individually'   logical, default : false
    'min_split_node_size'           Scalar, default: 1.0 - controls size pruning of the tree while training.
    'mtry'                          Scalar, default: floor(sqrt(number of features))

    'training_set_size'             Scalar, default: Not used
    'training_set_proportion'       Scalar, default: 1.0
                                    The last two options exclude each other. if training_set_size always overrides
                                    training_set_proportional, if set.
                                    Controls the number of samples drawn to train an individual tree.
    'weights'
                                    Array containing training weights for each class. The size of the array is
                                    not checked so you may get wierd errors if you do not enforce the size constraints.

oob_data    - A NxNumberOfTrees Matrix. oob_data(i, j) = 0 if ith sample was not in the test set for the jth tree
                                                       = 1 if ith sample was correctly classified in jth tree
                                                       = 2 if ith sample was misclassified int jth tree
*/
