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
#include <vigra/random_forest_deprec.hxx>
#include <set>
#include <cmath>

#define RandomForest RandomForestDeprec
#define RandomForestOptions RandomForestOptionsDeprec
#define DecisionTree DecisionTreeDeprec

namespace python = boost::python;
namespace vigra
{

template<class LabelType,class FeatureType>
RandomForest<LabelType>*
pythonConstructRandomForest(NumpyArray<2, FeatureType> trainData,
                            NumpyArray<1, LabelType> trainLabels,
                            int treeCount,
                            int mtry,
                            int min_split_node_size,
                            int training_set_size,
                            float training_set_proportions,
                            bool sample_with_replacement,
                            bool sample_classes_individually)

{
    RandomForestOptions options;
    options
        .featuresPerNode(mtry)
        .sampleWithReplacement(sample_with_replacement)
        .setTreeCount(treeCount)
        .trainingSetSizeProportional(training_set_proportions)
        .trainingSetSizeAbsolute(training_set_size)
        .sampleClassesIndividually(sample_classes_individually)
        .minSplitNodeSize(min_split_node_size);
        
    std::set<LabelType> uniqueLabels(trainLabels.data(), trainLabels.data()+trainLabels.size());

    RandomForest<LabelType>* rf = 
        new RandomForest<LabelType>(uniqueLabels.begin(), uniqueLabels.end(),
                                    treeCount, options);
    double oob;

    {
        PyAllowThreads _pythread;
        oob = rf->learn(trainData, trainLabels);
    }
    
    std::cout << "Out-of-bag error " << oob << std::endl;
    return rf;
}

template<class LabelType,class FeatureType>
NumpyAnyArray 
pythonRFPredictLabels(RandomForest<LabelType> const & rf,
                      NumpyArray<2,FeatureType> testData,
                      NumpyArray<2,LabelType> res)
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(testData.shape(0), 1),
            "Output array has wrong dimensions.");
    
    PyAllowThreads _pythread;
    rf.predictLabels(testData, res);
    return res;
}

template<class LabelType, class FeatureType>
NumpyAnyArray 
pythonRFPredictProbabilities(RandomForest<LabelType> const & rf,
                             NumpyArray<2,FeatureType> testData,
                             NumpyArray<2,float> res)
{
    res.reshapeIfEmpty(MultiArrayShape<2>::type(testData.shape(0), rf.labelCount()),
            "Output array has wrong dimensions.");
    {
        PyAllowThreads _pythread;
        rf.predictProbabilities(testData, res);
    }
    return res;
}

void defineRandomForestOld()
{
    using namespace python;
    
    docstring_options doc_options(true, true, false);

    class_<RandomForest<UInt32> > rfclass("RandomForestOld",python::no_init);
    rfclass
        .def("__init__",python::make_constructor(registerConverters(&pythonConstructRandomForest<UInt32,float>),
                                                 boost::python::default_call_policies(),
                                                 ( arg("trainData"), arg("trainLabels"),
                                                   arg("treeCount")=255,
                                                   arg("mtry")=0,
                                                   arg("min_split_node_size")=1,
                                                   arg("training_set_size")=0,
                                                   arg("training_set_proportions")=1.0,
                                                   arg("sample_with_replacement")=true,
                                                   arg("sample_classes_individually")=false)),
             "Constructor::\n\n"
             "  RandomForestOld(trainData, trainLabels,\n"
             "                  treeCount = 255, mtry=0, min_split_node_size=1,\n"
             "                  training_set_size=0, training_set_proportions=1.0,\n"
             "                  sample_with_replacement=True, sample_classes_individually=False,)\n\n"
             "Construct and train a RandomForest using 'trainData' and 'trainLabels'. "
             "'treeCount' controls the number of trees that are created.\n\n"
             "See RandomForest_ and RandomForestOptions_ in the C++ documentation "
             "for the meaning of the other parameters.\n")
        .def("featureCount",
             &RandomForest<UInt32>::featureCount,
             "Returns the number of features the RandomForest works with.\n")
        .def("labelCount",
             &RandomForest<UInt32>::labelCount,
             "Returns the number of labels, the RanfomForest knows.\n")
        .def("treeCount",
             &RandomForest<UInt32>::treeCount,
             "Returns the 'treeCount', that was set when constructing the RandomForest.\n")
        .def("predictLabels",
             registerConverters(&pythonRFPredictLabels<UInt32,float>),
             (arg("testData"), arg("out")=object()),
             "Predict labels on 'testData'."
             "The output is an array containing a labels for every test samples.\n")
        .def("predictProbabilities",
             registerConverters(&pythonRFPredictProbabilities<UInt32,float>),
             (arg("testData"), arg("out")=object()),
             "Predict probabilities for different classes on 'testData'."
             "The output is an array containing a probability for every test sample and class.\n")
/*            .def("writeHDF5")
        .def("readHDF5")*/
        ;
}

} // namespace vigra

#undef RandomForest
#undef RandomForestOptions
#undef DecisionTree
