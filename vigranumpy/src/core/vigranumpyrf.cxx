/************************************************************************/
/*                                                                      */
/*                 Copyright 2009 by Ullrich Koethe                     */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpycmodule_PyArray_API
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/random_forest_deprec.hxx>
#include <set>
#include <cmath>

#ifdef HasHDF5
# include <vigra/random_forest_hdf5_impex.hxx>
#endif

namespace python = boost::python;
namespace vigra
{

template<class LabelType,class FeatureType>
RandomForest<LabelType>*
pythonConstructRandomForest(NumpyArray<2,FeatureType> trainData,NumpyArray<1,LabelType> trainLabels,
							int treeCount,
							int mtry,
							int min_split_node_size,
							int training_set_size,
							float training_set_proportions,
							bool sample_with_replacement,
							bool sample_classes_individually)

{
	RandomForestOptions options;
	options.featuresPerNode(mtry).sampleWithReplacement(sample_with_replacement).setTreeCount(treeCount)
		.trainingSetSizeProportional(training_set_proportions).trainingSetSizeAbsolute(training_set_size)
		.sampleClassesIndividually(sample_classes_individually).minSplitNodeSize(min_split_node_size);
	std::set<LabelType> uniqueLabels(trainLabels.data(),trainLabels.data()+trainLabels.size());

	RandomForest<LabelType>* rf=new RandomForest<LabelType>(uniqueLabels.begin(),uniqueLabels.end(),treeCount,options);
	rf->learn(trainData,trainLabels);

	return rf;
}

template<class LabelType,class FeatureType>
NumpyAnyArray pythonRFPredictLabels(RandomForest<LabelType>* rf,NumpyArray<2,FeatureType> testData,NumpyArray<2,LabelType> res)
{
	//construct result
	res.reshapeIfEmpty(MultiArrayShape<2>::type(testData.shape(0),1),"Output array has wrong dimensions.");
	rf->predictLabels(testData,res);
	return res;
}

template<class LabelType, class FeatureType>
NumpyAnyArray pythonRFPredictProbabilities(RandomForest<LabelType>* rf,NumpyArray<2,FeatureType> testData,NumpyArray<2,float> res)
{
	//construct result
	res.reshapeIfEmpty(MultiArrayShape<2>::type(testData.shape(0),rf->labelCount()),"Output array has wrong dimensions.");
	rf->predictProbabilities(testData,res);
	return res;
}

void defineRandomForest()
{
	using namespace python;

	class_<RandomForest<UInt32> > rfclass("RandomForest",python::no_init);
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
			 "Construct and train a random Forest using \"trainData\" and \"trainLabels\".\n"
			 "\"treeCount\" constrols the number of trees, that are created.\n"
			 "See the vigra documentation for the meaning af the rest of the paremeters.")
		.def("featureCount",
			 &RandomForest<UInt32>::featureCount,
			 "Returns the number of features the RandomForest works with.")
		.def("labelCount",
			 &RandomForest<UInt32>::labelCount,
			 "Returns the number of labels, the RanfomForest knows.")
		.def("treeCount",
			 &RandomForest<UInt32>::treeCount,
			 "Returns the \"treeCount\", that was set when constructing the RandomForest.")
		.def("predictLabels",
			 registerConverters(&pythonRFPredictLabels<UInt32,float>),
			 (arg("testData"), arg("out")=object()),
			 "Predict labels on \"testData\"."
			 "The output is an array containing a labels for every test samples.")
		.def("predictProbabilities",
			 registerConverters(&pythonRFPredictProbabilities<UInt32,float>),
			 (arg("testData"), arg("out")=object()),
			 "Predict probabilities for different classes on \"testData\"."
			 "The output is an array containing a probability for every test sample and class.")
/*			.def("writeHDF5")
		.def("readHDF5")*/
		;
}

} // namespace vigra
