/************************************************************************/
/*                                                                      */
/*                 Copyright 2009 by Ullrich Koethe                     */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpylearning_PyArray_API
// #define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/random_forest.hxx>
#ifdef HasHDF5
# include <vigra/random_forest_hdf5_impex.hxx>
#endif
#include <set>
#include <cmath>
#include <memory>
#include <boost/python.hpp>

namespace python = boost::python;
namespace vigra
{

template<class FeatureType>
OnlinePredictionSet<FeatureType>* pythonConstructOnlinePredictioSet(NumpyArray<2,FeatureType> features,int num_sets)
{
    return new OnlinePredictionSet<FeatureType>(features,num_sets);
}

template<class LabelType, class FeatureType>
RandomForest<LabelType>*
pythonConstructRandomForest(int treeCount,
                            int mtry,
                            int min_split_node_size,
                            int training_set_size,
                            float training_set_proportions,
                            bool sample_with_replacement,
                            bool sample_classes_individually,
                            bool prepare_online)


{
    RandomForestOptions options;
    options .sample_with_replacement(sample_with_replacement)
            .tree_count(treeCount)
            .prepare_online_learning(prepare_online)
            .min_split_node_size(min_split_node_size);


    if(mtry  > 0)
        options.features_per_node(mtry);

    if(training_set_size != 0)
        options.samples_per_tree(training_set_size);
    else
        options.samples_per_tree(training_set_proportions);

    if(sample_classes_individually)
        options.use_stratification(RF_EQUAL);

    RandomForest<LabelType>* rf = new RandomForest<LabelType>(options);

    return rf;
}

#ifdef HasHDF5
template<class LabelType>
RandomForest<LabelType> * 
pythonImportRandomForestFromHDF5(std::string filename, 
					             std::string pathname = "")
{ 
    std::auto_ptr<RandomForest<LabelType> > rf(new RandomForest<LabelType>);
    
    vigra_precondition(rf_import_HDF5(*rf, filename, pathname),
           "RandomForest(): Unable to load from HDF5 file.");
           
    return rf.release();
}					
#endif // HasHDF5

template<class LabelType, class FeatureType>
python::tuple
pythonLearnRandomForestWithFeatureSelection(RandomForest<LabelType> & rf, 
                                            NumpyArray<2,FeatureType> trainData, 
                                            NumpyArray<2,LabelType> trainLabels)
{
    VariableImportanceVisitor var_imp;

    double oob;
	Py_BEGIN_ALLOW_THREADS
	oob = rf.learn(trainData, trainLabels, create_visitor(var_imp));
	Py_END_ALLOW_THREADS

    // std::cout << "out of bag: " << oob << std::endl;

    NumpyArray<2, double> varImp(MultiArrayShape<2>::type(var_imp.variable_importance_.shape(0),
                                                           var_imp.variable_importance_.shape(1))); 

    for (int x=0;x<varImp.shape(0);x++)
        for (int y=0;y<varImp.shape(1);y++)
            varImp(x,y)= var_imp.variable_importance_(x,y);

    return python::make_tuple(oob, varImp);
}

template<class LabelType, class FeatureType>
double
pythonLearnRandomForest(RandomForest<LabelType> & rf, 
                        NumpyArray<2,FeatureType> trainData, 
                        NumpyArray<2,LabelType> trainLabels)
{
  double oob;

  Py_BEGIN_ALLOW_THREADS
  oob = rf.learn(trainData, trainLabels);
  Py_END_ALLOW_THREADS
  
  //std::cout << "out of bag: " << oob << std::endl;
  return oob;
}

template<class LabelType,class FeatureType>
void 
pythonRFOnlineLearn(RandomForest<LabelType> & rf,
                    NumpyArray<2,FeatureType> trainData,
                    NumpyArray<2,LabelType> trainLabels,
                    int startIndex,
                    bool adjust_thresholds)
{
    Py_BEGIN_ALLOW_THREADS
    rf.onlineLearn(trainData, trainLabels, startIndex, adjust_thresholds);
    Py_END_ALLOW_THREADS
}

template<class LabelType,class FeatureType>
void 
pythonRFReLearnTree(RandomForest<LabelType> & rf,
                    NumpyArray<2,FeatureType> trainData,
                    NumpyArray<2,LabelType> trainLabels,
                    int treeId)
{
    Py_BEGIN_ALLOW_THREADS
    rf.reLearnTree(trainData,trainLabels,treeId);
    Py_END_ALLOW_THREADS
}

template<class LabelType,class FeatureType>
NumpyAnyArray 
pythonRFPredictLabels(RandomForest<LabelType> const & rf,
                      NumpyArray<2,FeatureType> testData,
                      NumpyArray<2,LabelType> res)
{
    //construct result
    res.reshapeIfEmpty(MultiArrayShape<2>::type(testData.shape(0),1),
                       "Output array has wrong dimensions.");
    rf.predictLabels(testData,res);
    return res;
}

template<class LabelType, class FeatureType>
NumpyAnyArray 
pythonRFPredictProbabilities(RandomForest<LabelType> & rf,
                             NumpyArray<2,FeatureType> testData, 
                             NumpyArray<2,float> res)
{
    //construct result
    res.reshapeIfEmpty(MultiArrayShape<2>::type(testData.shape(0), rf.ext_param_.class_count_),
                       "Output array has wrong dimensions.");
	Py_BEGIN_ALLOW_THREADS
    rf.predictProbabilities(testData,res);
	Py_END_ALLOW_THREADS
    return res;
}

template<class LabelType,class FeatureType>
NumpyAnyArray 
pythonRFPredictProbabilitiesOnlinePredSet(RandomForest<LabelType> & rf,
                                          OnlinePredictionSet<FeatureType> & predSet,
                                          NumpyArray<2,float> res)
{
    //construct result
    res.reshapeIfEmpty(MultiArrayShape<2>::type(predSet.features.shape(0),rf.ext_param_.class_count_),
                       "Output array has wrong dimenstions.");
    clock_t start=clock();
	Py_BEGIN_ALLOW_THREADS
    rf.predictProbabilities(predSet, res);
    Py_END_ALLOW_THREADS
    double duration=(clock()-start)/double(CLOCKS_PER_SEC);
    std::cerr<<"Prediction Time: "<<duration<<std::endl;
    return res;
}


void defineRandomForest()
{
    using namespace python;
    
    docstring_options doc_options(true, true, false);

    class_<OnlinePredictionSet<float> > pred_set_class("RF_OnlinePredictionSet",python::no_init);
    pred_set_class
        .def("__init__",
                python::make_constructor(registerConverters(&pythonConstructOnlinePredictioSet<float>),
                                         boost::python::default_call_policies(),
                                         (arg("features"))),
                                         "docu")

        .def("get_worsed_tree",&OnlinePredictionSet<float>::get_worsed_tree,
           "doku")
    
        .def("invalidateTree",&OnlinePredictionSet<float>::reset_tree,
         (arg("treeId")),
         "doku")
    ;

    enum_<RF_OptionTag>("RF_MTRY_SWITCH")
        .value("RF_MTRY_LOG",RF_LOG)
        .value("RF_MTRY_SQRT",RF_SQRT)
        .value("RF_MTRY_ALL",RF_ALL);

    class_<RandomForest<UInt32> > rfclass_new("RandomForest",python::no_init);

    rfclass_new
        .def("__init__",python::make_constructor(registerConverters(&pythonConstructRandomForest<UInt32,float>),
                                                 boost::python::default_call_policies(),
                                                 ( arg("treeCount")=255,
                                                   arg("mtry")= -1,
                                                   arg("min_split_node_size")=1,
                                                   arg("training_set_size")=0,
                                                   arg("training_set_proportions")=1.0,
                                                   arg("sample_with_replacement")=true,
                                                   arg("sample_classes_individually")=false,
                                                   arg("prepare_online_learning")=false)),
             "Constructor::\n\n"
             "  RandomForest(treeCount = 255, mtry=RF_SQRT, min_split_node_size=1,\n"
             "               training_set_size=0, training_set_proportions=1.0,\n"
             "               sample_with_replacement=True, sample_classes_individually=False,\n"
             "               prepare_online_learning=False)\n\n"
             "'treeCount' controls the number of trees that are created.\n\n"
             "See RandomForest_ and RandomForestOptions_ in the C++ documentation "
             "for the meaning of the other parameters.\n")
#ifdef HasHDF5
        .def("__init__",python::make_constructor(&pythonImportRandomForestFromHDF5<UInt32>,
                                                 boost::python::default_call_policies(),
                                                 ( arg("filename"),
                                                   arg("pathInFile")="")),
             "Load from HDF5 file::\n\n"
             "  RandomForest(filename, pathInFile)\n\n")
#endif // HasHDF5
        .def("featureCount",
            &RandomForest<UInt32>::column_count,
             "Returns the number of features the RandomForest works with.\n")
        .def("labelCount",
            &RandomForest<UInt32>::class_count,
             "Returns the number of labels, the RandomForest knows.\n")
        .def("treeCount",
             &RandomForest<UInt32>::tree_count,
             "Returns the 'treeCount', that was set when constructing the RandomForest.\n")
        .def("predictLabels",
             registerConverters(&pythonRFPredictLabels<UInt32,float>),
             (arg("testData"), arg("out")=object()),
             "Predict labels on 'testData'.\n\n"
             "The output is an array containing a labels for every test samples.\n")
        .def("predictProbabilities",
             registerConverters(&pythonRFPredictProbabilities<UInt32,float>),
             (arg("testData"), arg("out")=object()),
             "Predict probabilities for different classes on 'testData'.\n\n"
             "The output is an array containing a probability for every test sample and class.\n")
        .def("predictProbabilities",
             registerConverters(&pythonRFPredictProbabilitiesOnlinePredSet<UInt32,float>),
             (arg("testData"), arg("out")=object()),
             "The output is an array containing a probability for every test sample and class.\n")
        .def("learnRF",
             registerConverters(&pythonLearnRandomForest<UInt32,float>),
             (arg("trainData"), arg("trainLabels")),
             "Trains a random Forest using 'trainData' and 'trainLabels'.\n\n"
             "and returns the OOB. See the vigra documentation for the meaning af the rest of the paremeters.\n")
        .def("reLearnTree",
             registerConverters(&pythonRFReLearnTree<UInt32,float>),
            (arg("trainData"), arg("trainLabels"), arg("treeId")),
             "Re-learn one tree of the forest using 'trainData' and 'trainLabels'.\n\n"
             "and returns the OOB. This might be helpful in an online learning setup to improve the classifier.\n")
        .def("learnRFWithFeatureSelection",
             registerConverters(&pythonLearnRandomForestWithFeatureSelection<UInt32,float>),
             (arg("trainData"), arg("trainLabels")),
             "Train a random Forest using 'trainData' and 'trainLabels'.\n\n"
             "and returns the OOB and the Variable importance"
             "See the vigra documentation for the meaning af the rest of the paremeters.\n")
        .def("onlineLearn",
             registerConverters(&pythonRFOnlineLearn<UInt32,float>),
             (arg("trainData"),arg("trainLabels"),arg("startIndex")),
             "Learn online.\n\n"
             "Works only if forest has been created with prepare_online_learning=true. "
             "Needs the old training data and the new appened, starting at startIndex.\n\n")
#ifdef HasHDF5
        .def("writeHDF5", &rf_export_HDF5<UInt32>,
             (arg("filename"), arg("pathInFile")="", arg("overwriteflag")=false),
             "Store the random forest in the given HDF5 file 'filname' under the internal\n"
             "path 'pathInFile'. If a dataset already exists, 'overwriteflag' determines\n"
             "if the old data are overwritten.\n")
#endif // HasHDF5
        ;
}


void defineRandomForestOld();


} // namespace vigra


using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(learning)
{
    import_vigranumpy();
    defineRandomForest();
    defineRandomForestOld();
}


