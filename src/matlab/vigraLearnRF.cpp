/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by Rahul Nair and Ullrich Koethe          */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
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
#include "random_forest_impex.hxx"
#include "RandomForestProgressVisitor.hxx"


using namespace vigra;
using namespace matlab;


void vigraMain(matlab::OutputArray outputs, matlab::InputArray inputs){
    /***************************************************************************************************
    **              INIT PART                                                                         **
    ****************************************************************************************************/
    typedef double inputType;
    typedef double inputLType;
    vigra::RandomForestOptions 	options;
    vigra::ProblemSpec<inputLType>	ext_param;
    options.sample_with_replacement(inputs.getBool("sample_with_replacement", 
                                                   v_default(true)));
    
    if(inputs.getBool("sample_classes_individually", v_default(false)))
        options.use_stratification(vigra::RF_EQUAL);
    options.min_split_node_size(inputs
            .getScalarMinMax<double>("min_split_node_size",
                                     v_default(1.0), 
                                     0.0, "inf"));
    options.tree_count(inputs
            .getScalarMinMax<double>(2,v_default(255.0), 0.0, "inf"));

    if(inputs.hasData("mtry"))
    {
        if(inputs.typeOf("mtry") == mxCHAR_CLASS)
        {
            std::map<std::string, int> map_;
                map_["RF_LOG"] 	= int(RF_LOG);
                map_["RF_SQRT"] = int(RF_SQRT);
                map_["RF_ALL"]	= int(RF_ALL);
            RF_OptionTag method  = RF_OptionTag(inputs.getEnum("mtry", v_default(RF_LOG), map_));
        }
        else
        {
            options
                .features_per_node(int(inputs.getScalar<double>("mtry",v_default(0))));
        }
    }

    double training_set_size 
        = inputs.getScalar<double>("training_set_size",v_default(0));
    if(training_set_size != 0)
    {
        options.samples_per_tree(int(training_set_size));
    }
    else
    {
        options.samples_per_tree(inputs
            .getScalarMinMax<double>("training_set_proportion",
                                     v_default(1.0), 0.0, 1.0));
    }

    MultiArrayView<2, inputLType>  labels 
        = inputs.getMultiArray<2, inputLType>(1, v_required());
    MultiArrayView<2, inputType>  features 
        = inputs.getMultiArray<2, inputType>(0, v_required());
    MultiArrayView<1, inputType>  weights 
        = inputs.getMultiArray<1, inputType>("weights", v_optional());

    if(weights.size() != 0)
        ext_param.class_weights(weights.data(), weights.data() + weights.size());

    double var_imp_rep
        = inputs.getScalar<double>("importance_repetition",v_default(10));
    
    VariableImportanceVisitor var_imp(var_imp_rep);
    RandomForestProgressVisitor progress;
    if(!outputs.isValid(2))
        var_imp.deactivate();


    /***************************************************************************************************
    **              CODE PART                                                                         **
    ****************************************************************************************************/

    RandomForest<inputLType> rf(options, ext_param);    
    double oobError = rf.learn(features, 
                               labels,
                               vigra::create_visitor(var_imp, progress));

    matlab::exportRandomForest(rf, matlab::createCellArray(2*options.tree_count_+2, outputs[0]));

    outputs.createScalar<double> (1, v_optional(), oobError);
    MultiArrayView<2, double> vari 
        = outputs.createMultiArray<2, double>(2, v_optional(), 
                            MultiArrayShape<2>::type(var_imp
                                                      .variable_importance_
                                                      .shape()));
    if(vari.size() != 0)
        vari = var_imp.variable_importance_;
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
function [RF oob var_imp] = vigraLearnRF(...)       Outputs variable importance.

features    - A Nxp Matrix with N samples containing p features
labels      - A Nx1 Matrix with the corresponding Training labels
treeCount   - default: 255. An Integral Scalar Value > 0 - Number of Trees to be used in the Random Forest.
options     - a struct with the following possible fields (default will be used
              if field is not present)
    'sample_with_replacement'       logical, default : true
    'sample_classes_individually'   logical, default : false
    'min_split_node_size'           Scalar, default: 1.0 - controls size pruning of the tree while training.
    'mtry'                          Scalar or String, 
                                    default: floor(sqrt(number of features)) ('RF_SQRT')
                                    if a Scalar value is specified it is taken as the 
                                    absolute value. Otherwise use one of the Tokens
                                    'RF_SQRT', 'RF_LOG' or 'RF_ALL'

    'training_set_size'             Scalar, default: Not used
    'training_set_proportion'       Scalar, default: 1.0
                                    The last two options exclude each other. if training_set_size always overrides
                                    training_set_proportional, if set.
                                    Controls the number of samples drawn to train an individual tree.
    'weights'
                                    Array containing training weights for each class. The size of the array is
                                    not checked so you may get wierd errors if you do not enforce the size constraints.
var_imp		- A FeatureCount x ClassCount +2 Matrix. 
                                    The last column is the variable importance based on mean decrease in impurity
                                    over all trees the end -1 column is the permutation based variable importance
                                    Columns 1 - ClassCount are the class wise permutation based variable importance
                                    scores.

//not yet supported
oob_data    - A NxNumberOfTrees Matrix. oob_data(i, j) = 0 if ith sample was not in the test set for the jth tree
                                                       = 1 if ith sample was correctly classified in jth tree
                                                       = 2 if ith sample was misclassified int jth tree
*/
