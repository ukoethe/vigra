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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyclassification_PyArray_API
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/random_forest.hxx>
#include <set>
#include <cmath>
#include "tuples.hpp"


#ifdef HasHDF5
# include <vigra/random_forest_hdf5_impex.hxx>
#endif

namespace python = boost::python;
namespace vigra
{

template<class LabelType, class FeatureType>
RandomForest<LabelType>*
pythonConstructRandomForest(int treeCount,
							int mtry,
							int min_split_node_size,
							int training_set_size,
							float training_set_proportions,
							bool sample_with_replacement,
							bool sample_classes_individually)

{
	RandomForestOptions options;
	options.features_per_node(mtry).sample_with_replacement(sample_with_replacement).tree_count(treeCount)
	.min_split_node_size(min_split_node_size);


	if(training_set_size != 0)
		options.samples_per_tree(training_set_size);
	else
		options.samples_per_tree(training_set_proportions);

	if(sample_classes_individually)
		options.use_stratification(RF_EQUAL);

	RandomForest<LabelType>* rf = new RandomForest<LabelType>(options);

	return rf;
}

template<class LabelType, class FeatureType>
boost::tuples::tuple<double, NumpyArray<2,double> >
pythonLearnRandomForestWithFeatureSelection(RandomForest<LabelType>* rf, NumpyArray<2,FeatureType> trainData, NumpyArray<2,LabelType> trainLabels)
{
  VariableImportanceVisitor var_imp;
   
  double oob = rf->learn(trainData, trainLabels, create_visitor(var_imp));
  std::cout << "out of bag: " << oob << std::endl;
  
  NumpyArray<2, double> varImp(MultiArrayShape<2>::type(var_imp.variable_importance_.shape(0),var_imp.variable_importance_.shape(1))); 
  
  for (int x=0;x<varImp.shape(0);x++)
    for (int y=0;y<varImp.shape(1);y++)
      varImp(x,y)= var_imp.variable_importance_(x,y);

  return boost::tuples::tuple<double, NumpyArray<2, double> >(oob, varImp);
}

template<class LabelType, class FeatureType>
double
pythonLearnRandomForest(RandomForest<LabelType>* rf, NumpyArray<2,FeatureType> trainData, NumpyArray<2,LabelType> trainLabels)
{
  double oob = rf->learn(trainData, trainLabels);
  std::cout << "out of bag: " << oob << std::endl;

  return oob;
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
NumpyAnyArray pythonRFPredictProbabilities(RandomForest<LabelType>* rf,NumpyArray<2,FeatureType> testData, NumpyArray<2,float> res)
{
	//construct result
	res.reshapeIfEmpty(MultiArrayShape<2>::type(testData.shape(0),rf->ext_param_.class_count_),"Output array has wrong dimensions.");
	rf->predictProbabilities(testData,res);
	return res;
}

void defineRandomForest_new()
{
	using namespace python;

  typedef boost::tuples::tuple<double, NumpyArray<2,double> > trainReturn;
  boost::python::register_tuple< trainReturn >();

	class_<RandomForest<UInt32> > rfclass_new("RandomForest_new",python::no_init);
	rfclass_new
		.def("__init__",python::make_constructor(registerConverters(&pythonConstructRandomForest<UInt32,float>),
												 boost::python::default_call_policies(),
												 ( arg("treeCount")=255,
												   arg("mtry")=0,
												   arg("min_split_node_size")=1,
												   arg("training_set_size")=0,
												   arg("training_set_proportions")=1.0,
												   arg("sample_with_replacement")=true,
												   arg("sample_classes_individually")=false)),
			 "Constructs a random Forest \n"
			 "\"treeCount\" constrols the number of trees, that are created.\n"
			 "See the vigra documentation for the meaning af the rest of the paremeters.")
		.def("featureCount",
    &RandomForest<UInt32>::column_count,
			 "Returns the number of features the RandomForest works with.")
		.def("labelCount",
    &RandomForest<UInt32>::class_count,
			 "Returns the number of labels, the RanfomForest knows.")
		.def("treeCount",
			 &RandomForest<UInt32>::tree_count,
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
    .def("learnRF",
			 registerConverters(&pythonLearnRandomForest<UInt32,float>),
			 (arg("trainData"), arg("trainLabels")),
			 "Trains a random Forest using \"trainData\" and \"trainLabels\".\n"
			 "and returns the OOB. See the vigra documentation for the meaning af the rest of the paremeters.")
    .def("learnRFWithFeatureSelection",
			 registerConverters(&pythonLearnRandomForestWithFeatureSelection<UInt32,float>),
			 (arg("trainData"), arg("trainLabels")),
			 "Train a random Forest using \"trainData\" and \"trainLabels\".\n"
			 "and returns the OOB and the Variable importance"
			 "See the vigra documentation for the meaning af the rest of the paremeters.")

/*			.def("writeHDF5")
		.def("readHDF5")*/
		;
}

} // namespace vigra

