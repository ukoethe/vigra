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
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/random_forest_3.hxx>

#ifdef HasHDF5
# include <vigra/random_forest_3_hdf5_impex.hxx>
#endif

#include <boost/python.hpp>

//#include <vigra/timing.hxx>
//#include <vigra/random.hxx>
//#include <set>
//#include <cmath>
//#include <memory>

namespace python = boost::python;

namespace vigra
{
namespace rf3
{

typedef float FeatureType;
typedef uint32_t LabelType;
typedef NumpyArray<2,FeatureType> FeatureArrayType;
typedef NumpyArray<1,LabelType> LabelArrayType; // FIXME TODO dunno if this needs to be 1 or 2d
typedef DefaultRF<FeatureArrayType,LabelArrayType>::type RandomForestType;

// random forest constructor
RandomForestType *
pythonConstructRandomForest3(
        FeatureArrayType features, 
        LabelArrayType labels,
        int treeCount,
        int mtry,
        int min_split_node_size,
        bool sample_with_replacement,
        bool sample_classes_individually,
        size_t resample_count,
        size_t max_depth,
        double tau,
        size_t n_threads
        )
{
    
    RandomForestOptions rf_opts;

    rf_opts.tree_count(treeCount);
    
    if(mtry > 0)
        rf_opts.features_per_node(mtry);
    
    // TODO training_set_size -> can't find the corresponding parameter, afaik this is the number of bootstrap samples used
    
    rf_opts.bootstrap_sampling(sample_with_replacement);
    rf_opts.min_num_instances(min_split_node_size);
    rf_opts.use_stratification(sample_classes_individually);
    
    // this is the number of instances that is resampled in each split / this disables bootstrap sampling and
    rf_opts.resample_count(resample_count);
    
    // expose the max tree depth, not in old rf afaik
    rf_opts.max_depth(max_depth);
    
    //expose node complexity, not in okd rf afaik
    rf_opts.node_complexity_tau(tau);
    
    // expose n_threads for multithreading
    rf_opts.n_threads(n_threads);
    
    // TODO expose class_weights
    // class_weights(std::vector<double> const & v)

    // TODO expose the split criterion (GINI, ENTTROPY, KSD), not in old rf afaik
    // split(RandomForestOptionTags p_split)
    
    // return pointer to the random forest
    PyAllowThreads _pythread;
    RandomForestType rf_tmp = random_forest(features, labels, rf_opts);
    RandomForestType * rf = new RandomForestType(rf_tmp);
    return rf;
}

// prediction
NumpyAnyArray 
pythonPredictProbabilities(const RandomForestType & rf,
                           FeatureArrayType features,
                           size_t n_threads,
                           NumpyArray<2,float> res) {

    res.reshapeIfEmpty(MultiArrayShape<2>::type(features.shape(0), rf.num_classes()),
                       "RandomForest.predictProbabilities(): Output array has wrong dimensions.");
    {
        PyAllowThreads _pythread;
        rf.predict_probabilities(features, res, n_threads);
    }
    return res;
}

NumpyAnyArray 
pythonPredictLabels(const RandomForestType & rf,
                           FeatureArrayType features,
                           size_t n_threads,
                           NumpyArray<1,LabelType> res) {

    res.reshapeIfEmpty(MultiArrayShape<1>::type(features.shape(0)),
                       "RandomForest.predictProbabilities(): Output array has wrong dimensions.");
    {
        PyAllowThreads _pythread;
        rf.predict(features, res, n_threads);
    }
    return res;
}

#ifdef HasHDF5
RandomForestType * 
pythonImportFromHDF5(const std::string & filename, const std::string & pathname )
{ 
    HDF5File h5ctx(filename);
    RandomForestType rf_tmp = random_forest_import_HDF5<FeatureArrayType,LabelArrayType>(h5ctx, pathname);
    RandomForestType * rf = new RandomForestType(rf_tmp);
    return rf;
}                   

void pythonExportHDF5(const RandomForestType & rf,
        const std::string & filename,
        const std::string & pathname)
{
    HDF5File h5ctx(filename, HDF5File::ReadWrite );
    random_forest_export_HDF5(rf, h5ctx, pathname);
}
#endif // HasHDF5
            
// expose the rf 3 to python
// expose the rf 3 as a python class + basic API
void exportRandomForest3()
{
    using namespace python;
    
    docstring_options doc_options(true, true, false);
    
    enum_<RandomForestOptionTags>("RF3_MTRY_SWITCH")
        .value("RF3_MTRY_LOG", vigra::rf3::RF_LOG)
        .value("RF3_MTRY_SQRT",vigra::rf3::RF_SQRT)
        .value("RF3_MTRY_ALL", vigra::rf3::RF_ALL);
    
    class_<RandomForestType> rfclass3("RandomForest3",python::no_init);
    //class_<RandomForestType> rfclass3("RandomForest3");

    rfclass3
       .def("__init__",
            python::make_constructor(registerConverters(&pythonConstructRandomForest3),
            boost::python::default_call_policies(),
            (arg("features"),
              arg("labels"),
              arg("treeCount")=255,
              arg("mtry")= -1,
              arg("min_split_node_size")=1,
              arg("sample_with_replacement")=true,
              arg("sample_classes_individually")=false,
              arg("resample_count")=0,
              arg("max_depth")=0,
              arg("tau")=-1,
              arg("n_threads")=1)
            ),
             "\nConstruct a new random forest::\n\n"
             "  RandomForest(features, labels, treeCount = 255, mtry=RF_SQRT, min_split_node_size=1,\n"
             "               training_set_size=0, training_set_proportions=1.0,\n"
             "               sample_with_replacement=True, sample_classes_individually=False,\n"
             "               )\n\n"
             "treeCount:\n"
             "     controls the number of trees that are created.\n"
             "See RandomForest_ and RandomForestOptions_ in the C++ documentation "
             "for the meaning of the other parameters.\n")
#ifdef HasHDF5
         .def("__init__",python::make_constructor(&pythonImportFromHDF5,
                                                 boost::python::default_call_policies(),
                                                 ( arg("filename"),
                                                   arg("pathInFile"))),
             "\nLoad from HDF5 file::\n\n"
             "  RandomForest(filename, pathInFile)\n\n")
        .def("writeHDF5", &pythonExportHDF5,
                    boost::python::default_call_policies(),
                    ( arg("filename"),
                      arg("pathInFile")),
                    "Store the random forest in the given HDF5 file 'filename' under the internal\n"
                    "path 'pathInFile'.\n")
#endif // HasHDF5
        .def("featureCount",
            &RandomForestType::num_features,
             "Returns the number of features the RandomForest works with.\n")
        .def("labelCount",
            &RandomForestType::num_classes,
             "Returns the number of labels, the RandomForest knows.\n")
        .def("treeCount",
             &RandomForestType::num_trees,
             "Returns the 'treeCount', that was set when constructing the RandomForest.\n")
        .def("predictProbabilities",
            registerConverters(&pythonPredictProbabilities),
            (
             arg("features"),
             arg("n_threads")=1,
             arg("out")=object()
            ),
             "Predict probabilities for different classes on 'testData'.\n\n"
             "The output is an array containing a probability for every test sample and class.\n")
        .def("predictLabels",
             registerConverters(&pythonPredictLabels),
             (
              arg("features"), 
             arg("n_threads")=1,
              arg("out")=object()),
             "Predict labels on 'features'.\n\n"
             "The output is an array containing a label for every test samples.\n")
    ;
}




} // end namespace rf3
} // end namepace vigra


