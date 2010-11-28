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
#define CLASSIFIER_TEST 1
#define HDF5 0
#define CROSSVAL 0

#ifdef _MSC_VER
# pragma warning (disable : 4244)
#endif

#include <iostream>
#include <fstream>
#include <functional>
#include <cmath>
#include <vigra/random_forest.hxx>
#include <vigra/random_forest_deprec.hxx>
#include <unittest.hxx>
#include <vector>
#include <limits>
//#include "data/RF_results.hxx"
#include "data/RF_data.hxx"
#include "test_visitors.hxx"

#if CROSSVAL
#include <vigra/crossvalidation.hxx>
#endif


#include <stdlib.h>

#ifdef HasHDF5
# include <vigra/hdf5impex.hxx>
# include <vigra/random_forest_hdf5_impex.hxx>
# if HDF5
#  include "vigra/dataset4.hxx"
# endif
#endif

using namespace vigra;


struct UnaryRandomFunctor
{
    double operator()(double)
    {
        return vigra::RandomTT800::global().uniform53();
    }
};


struct ClassifierTest
{
    RF_Test_Training_Data data;
    //RF_Test_Result_Data results;

    ClassifierTest():
        data()
    {}


/**     ClassifierTest::diffOnfiles()
    Helper method to check whether two files are binary equal.
**/
    void diffOnfiles(std::string oldName, std::string newName)
    {
        std::ifstream oldRF(oldName.c_str());
        if(oldRF.fail())
            vigra_fail("files with old values not found");
        std::ifstream newRF(newName.c_str());
        if(newRF.fail())
            vigra_fail("files with new values not found");
        std::string oldS, newS;
        int k = 0;
        while(std::getline(newRF, newS) && std::getline(oldRF, oldS))
        {

            std::ostringstream s1;
            s1 << "Error in line " << k << std::endl;
            //shouldEqualMessage(newS, oldS, s1.str().c_str() );
            vigra::detail::equal_impl(newS, oldS,s1.str().c_str(), __FILE__, __LINE__);
            ++k;
        }
        //clean up garbage.
        std::ofstream newRFclear(newName.c_str());
        newRFclear.close();
    }
    
/** Tests whether RF rejects Matrices containing Inf
 */	
	void RF_InfCheck()
    {
        //Create Test output by Random Forest
		typedef MultiArrayShape<2>::type Shp_t;
        {
            std::cerr << "RF_InfCheck(): checking contains_inf() and RF precondition on Inf\n";
			MultiArray<2, double> dfeatures(Shp_t(10, 5));
			MultiArray<2, float>  ffeatures(Shp_t(10, 5));
			
			for(int ii = 0; ii < 10; ++ii)
			{
				for(int jj = 0; jj < 5; ++jj)
				{
					dfeatures(ii, jj) = (ii*jj+jj)%97;
					ffeatures(ii, jj) = (ii*jj+jj)%97;
				}
			}
			shouldEqual(detail::contains_inf(dfeatures), false);
			shouldEqual(detail::contains_inf(ffeatures), false);
					
			dfeatures(2, 3) = std::numeric_limits<double>::infinity();
			ffeatures(2, 3) = std::numeric_limits<double>::infinity();
			
			shouldEqual(detail::contains_inf(dfeatures), true);
			shouldEqual(detail::contains_inf(ffeatures), true);
			
			for(int ii = 0; ii < 10; ++ii)
			{
				for(int jj = 0; jj < 5; ++jj)
				{
					dfeatures(ii, jj) = std::numeric_limits<double>::infinity();
					ffeatures(ii, jj) = std::numeric_limits<double>::infinity();
				}
			}
			
			shouldEqual(detail::contains_inf(dfeatures), true);
			shouldEqual(detail::contains_inf(ffeatures), true);
        }
		MultiArray<2, double> dfeatures(Shp_t(10, 5));
		MultiArray<2, double> dlabels(Shp_t(10, 1));
		for(int ii = 0; ii < 10; ++ii)
		{
			for(int jj = 0; jj < 5; ++jj)
			{
				dfeatures(ii, jj) = (ii*jj+jj)%97;
			}
			dlabels(ii, 0) = ii%2;
		}
		dlabels(9,0) = std::numeric_limits<double>::infinity();
		RandomForest<> rf;
		try
		{
			rf.learn(dfeatures, dlabels);
			shouldEqual(0, 1);
		}
		catch( ... )
		{
		}
		dfeatures(9,0) = dlabels(9,0);
		dlabels(9,0) = 9%2;
		try
		{
			rf.learn(dfeatures, dlabels);
			shouldEqual(0, 1);
		}
		catch( ... )
		{
		}

    }
	
	void RF_AlgorithmTest()
    {
		std::cerr << "WARNING: THIS TEST CURRENTLY ONLY CHECKS WHETHER ALGORITHMS COMPILE AND RUN WITHOUT ERROR...\n";
		int ii = data.size() - 3; // this is the pina_indians dataset
		std::cerr << "RF_AlgorithmTest()....";
		
		rf::algorithms::HClustering				linkage;
  		MultiArray<2, double>	distance;
		std::cerr << "cluster_permutation_importance()....\n";
  		cluster_permutation_importance(data.features(ii), data.labels(ii), linkage, distance);
		rf::algorithms::VariableSelectionResult  result1;
		std::cerr << "forward_selection()....\n";
		rf::algorithms::forward_selection(data.features(ii), data.labels(ii), result1);
		rf::algorithms::VariableSelectionResult  result;
		std::vector<int> r;
		r.push_back(0);
		r.push_back(1);
		r.push_back(2);
		r.push_back(3);
		r.push_back(4);
		r.push_back(5);
		r.push_back(6);
		r.push_back(7);
		result.init(data.features(ii), data.labels(ii), r.begin(), r.end());
		std::cerr << "rank_selection()....\n";
		rf::algorithms::rank_selection(data.features(ii), data.labels(ii), result);
		
		rf::algorithms::VariableSelectionResult  result2;
		std::cerr << "backward_elimination()....\n";
		rf::algorithms::backward_elimination(data.features(ii), data.labels(ii), result2);
		std::cerr << "DONE\n";
	}
	void RF_SpliceTest()
    {
		std::cerr << "RF_SpliceTest()....";

		MultiArray<2, double> a(MultiArrayShape<2>::type(20, 10), 1);
		for(int ii = 0; ii< 20; ii+=2)
		{
			for(int jj = 0; jj< 10; jj+=2)
			{
				a(ii, jj) = 3;
			}
		}
		MultiArray<2, double> b(_spl_shp(_spl(0,2,20),_spl(0,2,10)));
		copy_splice(_spl(0,2,20), _spl(0,2,10), a, b);

		MultiArray<2, double> c(b.shape(), 3);
		shouldEqual(b, c);
		std::cerr << "DONE\n";
	}

	
/** Check whether the ridge regression functor compiles and runs
 */
    void RFridgeRegressionTest()
    {
		std::cerr << "RF_ridgeRegressionTest()....";
		std::cerr << "WARNING: THIS TEST CURRENTLY ONLY CHECKS WHETHER RIDGE REGRESSION COMPILES AND RUNS WITHOUT ERROR...";
		int ii = data.size() - 3; // this is the pina_indians dataset
	    //the ridge split class
        vigra::RandomForest<> rf(vigra::RandomForestOptions().tree_count(1));
		vigra::GiniRidgeSplit ridgeSplit;
	    rf.learn(data.features(ii), data.labels(ii), rf_default(), ridgeSplit);
		std::cerr << "DONE\n";
    }
	
/** Tests whether RF rejects Matrices containing NaN
 */	
	void RF_NanCheck()
    {
        //Create Test output by Random Forest
		typedef MultiArrayShape<2>::type Shp_t;
        {
            std::cerr << "RF_NanCheck(): checking contains_nan() and RF precondition on NaN\n";
			MultiArray<2, double> dfeatures(Shp_t(10, 5));
			MultiArray<2, float>  ffeatures(Shp_t(10, 5));
			for(int ii = 0; ii < 10; ++ii)
			{
				for(int jj = 0; jj < 5; ++jj)
				{
					dfeatures(ii, jj) = (ii*jj+jj)%97;
					ffeatures(ii, jj) = (ii*jj+jj)%97;
				}
			}
			shouldEqual(detail::contains_nan(dfeatures), false);
			shouldEqual(detail::contains_nan(ffeatures), false);
					
			dfeatures(2, 3) = std::numeric_limits<double>::quiet_NaN();
			ffeatures(2, 3) = std::numeric_limits<float>::quiet_NaN();
			
			shouldEqual(detail::contains_nan(dfeatures), true);
			shouldEqual(detail::contains_nan(ffeatures), true);
			
			for(int ii = 0; ii < 10; ++ii)
			{
				for(int jj = 0; jj < 5; ++jj)
				{
					dfeatures(ii, jj) = std::numeric_limits<double>::quiet_NaN();
					ffeatures(ii, jj) = std::numeric_limits<double>::quiet_NaN();
				}
			}
			
			shouldEqual(detail::contains_nan(dfeatures), true);
			shouldEqual(detail::contains_nan(ffeatures), true);
			
			for(int ii = 0; ii < 10; ++ii)
			{
				for(int jj = 0; jj < 5; ++jj)
				{
					dfeatures(ii, jj) = (ii*jj+jj)%97;
					ffeatures(ii, jj) = (ii*jj+jj)%97;
				}
			}
			shouldEqual(detail::contains_nan(dfeatures), false);
			shouldEqual(detail::contains_nan(ffeatures), false);
					
			dfeatures(2, 3) = std::numeric_limits<double>::signaling_NaN();
			ffeatures(2, 3) = std::numeric_limits<double>::signaling_NaN();
			
			shouldEqual(detail::contains_nan(dfeatures), true);
			shouldEqual(detail::contains_nan(ffeatures), true);
			
			for(int ii = 0; ii < 10; ++ii)
			{
				for(int jj = 0; jj < 5; ++jj)
				{
					dfeatures(ii, jj) = std::numeric_limits<double>::signaling_NaN();
					ffeatures(ii, jj) = std::numeric_limits<double>::signaling_NaN();
				}
			}
			
			shouldEqual(detail::contains_nan(dfeatures), true);
			shouldEqual(detail::contains_nan(ffeatures), true);
        }
		MultiArray<2, double> dfeatures(Shp_t(10, 5));
		MultiArray<2, double> dlabels(Shp_t(10, 1));
		for(int ii = 0; ii < 10; ++ii)
		{
			for(int jj = 0; jj < 5; ++jj)
			{
				dfeatures(ii, jj) = (ii*jj+jj)%97;
			}
			dlabels(ii, 0) = ii%2;
		}
		dlabels(9,0) = std::numeric_limits<double>::quiet_NaN();
		RandomForest<> rf;
		try
		{
			rf.learn(dfeatures, dlabels);
			shouldEqual(0, 1);
		}
		catch( ... )
		{
		}
		dfeatures(9,0) = dlabels(9,0);
		dlabels(9,0) = 9%2;
		try
		{
			rf.learn(dfeatures, dlabels);
			shouldEqual(0, 1);
		}
		catch( ... )
		{
		}

    }
/**
        ClassifierTest::RFdefaultTest():
    Learns The Refactored Random Forest with a fixed Random Seed and default sampling Options on
    various UCI data sets. A print Visitor outputs Node and Split Information into a file
    (RandomForestNodeTest.log) This file is then compared with a similar file created with
    the old (functional) RandomForest (oldClassifier.log). As long as the order, in which the
    RND_ object is called is not changed the files should be binary identical and this Test
    should always succeed.
**/
    void RFdefaultTest()
    {
        //Create Test output by Random Forest
        {
            vigra::TestVisitor testVisitor;
            std::cerr << "RFdefaultTest(): Learning on Datasets\n";
            for(int ii = 0; ii < data.size() ; ii++)
            {
                vigra::RandomForest<>
					RF2(vigra::RandomForestOptions().tree_count(32));
                vigra::RandomForest<int>
                    RF3(vigra::RandomForestOptions().tree_count(5), vigra::ProblemSpec<int>().classes_(data.ClassIter(ii).begin(), data.ClassIter(ii).end()));

                RF3.learn(  data.features(ii),
                            data.labels(ii),
						   	rf_default(),
							rf_default(),
							rf_default(),
                            vigra::RandomMT19937(1));

                RF2.learn(  data.features(ii),
                            data.labels(ii),
						   	create_visitor(testVisitor),
							rf_default(),
							rf_default(),
                            vigra::RandomMT19937(1));

                testVisitor.fout <<  data.names(ii) << std::endl;
                std::cerr << "[";
                for(int ss = 0; ss < ii+1; ++ss)
                    std::cerr << "#";
                for(int ss = ii+1; ss < data.size(); ++ss)
                    std::cerr << " ";
                std::cerr << "] " << data.names(ii);
                std::cerr << "\n";
                
                // just to ensure that using other label classes does 
                // not destroy anything
                shouldEqual(RF2.tree(0).topology_, RF3.tree(0).topology_);
				
                shouldEqual(data.features(ii).shape(0),
						    RF2.ext_param_.row_count_);
				shouldEqual(data.features(ii).shape(1),
							RF2.ext_param_.column_count_);
				shouldEqual(data.ClassIter(ii).size(), 
							RF2.ext_param_.class_count_);
            }
        }
        std::cerr << std::endl;
        std::cerr << "RFdefaultTest(): Comparing with Working Version:";
        //Cheap diff on old and new Classifier.
        diffOnfiles("./data/oldClassifier.log", "RandomForestNodeTest.log");
        std::cerr << "DONE!\n\n";
    }

	/** checks whether the agglomeration of data and conversion from internal
	 * to external labels is working 
	 */
	void RFresponseTest()
	{
		int ii = 2; 
		// learn on glass data set and predict: 
		// this is interesting because there is no label with number 4
		// in this dataset.
		
		// check whether agglomeration of probabilities is done properlyy
	    std::cerr << "RFresponseTest(): Learning on Datasets\n";
		vigra::RandomForest<>
			RF(vigra::RandomForestOptions().tree_count(2)); 

		RF.learn( data.features(ii),
				  data.labels(ii),
				  rf_default(),
				  rf_default(),
				  rf_default(),
				  vigra::RandomMT19937(1));
        
		typedef MultiArrayShape<2>::type Shp;
		MultiArray<2, double> response(Shp(data.features(ii).shape(0),
									   data.ClassIter(ii).size()));
		RF.predictProbabilities(data.features(ii), response);
		for(int jj = 0; jj < response.shape(0); ++jj)
		{
			typedef ArrayVector<double>::iterator Iter;
			ArrayVector<double> tmp(data.ClassIter(ii).size(), 0.0);
			Iter a = RF.tree(0).predict(rowVector(data.features(ii),jj));
			double totalWeight = 0.0;
			std::transform(tmp.begin(), tmp.end(), 
						   a, tmp.begin(), std::plus<double>());
			totalWeight = std::accumulate(a, a + data.ClassIter(ii).size(),
										  totalWeight);
			a = RF.tree(1).predict(rowVector(data.features(ii),jj));
			std::transform(tmp.begin(), tmp.end(), 
						   a, tmp.begin(), std::plus<double>());
			totalWeight = std::accumulate(a, a + data.ClassIter(ii).size(),
										  totalWeight);
			std::transform(tmp.begin(), tmp.end(),tmp.begin(), 
						   std::bind2nd(std::divides<double>(), totalWeight));
			MultiArrayView<2, double> 
				should_resp(Shp(1, data.ClassIter(ii).size()), tmp.data());

			shouldEqual(rowVector(response, jj), should_resp);
		}
        
        StopAfterTree stopAfterTree(1);
        StopIfMargin  stopIfMargin(0.5);
        StopAfterVoteCount stopAfterVoteCount(0.5);
        StopIfConverging   stopIfConverging(0.5);
		RF.predictProbabilities(data.features(ii), response,stopAfterTree);
		RF.predictProbabilities(data.features(ii), response,stopIfMargin);
		RF.predictProbabilities(data.features(ii), response,stopAfterVoteCount);
		RF.predictProbabilities(data.features(ii), response,stopIfConverging);
		// to check whether labels are being currectly converted we use the
		// property of the random forest to almost surely have 0 prediction
		// error on the training data. with enough trees.
		vigra::RandomForest<>
			RF2(vigra::RandomForestOptions().tree_count(255)); 

		RF2.learn( data.features(ii),
				  data.labels(ii),
				  rf_default(),
				  rf_default(),
				  rf_default(),
				  vigra::RandomMT19937(1));
		MultiArray<2, double> dble_labels(Shp(data.features(ii).shape(0),
											  1));
		MultiArray<2, Int32>  int_labels(Shp(data.features(ii).shape(0),
											 1));

		RF2.predictLabels(data.features(ii), dble_labels);
	    RF2.predictLabels(data.features(ii), int_labels);
	
		for(int jj = 0; jj< data.features(ii).shape(0); ++jj)
		{
			shouldEqualTolerance(dble_labels[jj], data.labels(ii)[jj], 0.01);
			shouldEqualTolerance(int_labels[jj], data.labels(ii)[jj], 0.01);
		}
		std::cerr << "done \n";
	}

/** Learns The Refactored Random Forest with 100 trees 10 times and
 * 	calulates the mean oob error. The distribution of the oob error
 * 	is gaussian as a first approximation. The mean oob error should
 * 	thus be in a 3 sigma neighborhood of the oob error produced by the working
 * 	version of the Random Forest.
 */
    void RFoobTest()
    {

            std::cerr << "RFoobTest(): Learning each Datasets 10 times\n";
            typedef MultiArrayShape<2>::type _TTT;
            for(int ii = 0; ii <data.size() ; ii++)
            {
                double oob = 0.0;
                for(int jj = 0; jj < 10; ++jj)
                {
                    rf::visitors::OOB_PerTreeError oob_v;
                    vigra::RandomForest<> RF2(vigra::RandomForestOptions()
												 .tree_count(100));

                     RF2.learn(  data.features(ii), 
                                 data.labels(ii),
                                 rf::visitors::create_visitor(oob_v));
                     oob += oob_v.oobError;
                }
                oob = oob/10;
                oob = oob - data.oobError(ii);
                oob = oob<0? -oob : oob;
                std::ostringstream s1;
                s1 	<< "Error - mean oobError exceeds 3 sigma bound:  " << oob
				 	<< "<-->" <<  data.oobError(ii) << std::endl;
                vigra::detail::should_impl(oob < 3* data.oobSTD(ii),
										   s1.str().c_str(),
										   __FILE__,
										   __LINE__);
                std::cerr << "[";
                for(int ss = 0; ss < ii+1; ++ss)
                    std::cerr << "#";
                for(int ss = ii+1; ss < data.size(); ++ss)
                    std::cerr << " ";
                std::cerr << "] " << data.names(ii);
                std::cerr << "\n";
            }
            std::cerr << "done!\n";
    }


/**
        ClassifierTest::RFsetTest():
    Learns The Refactored Random Forest with 1200 Trees default options and random Seed for the
    XOR problem. There are 30 distinct Trees possible. This test checks whether all solutions were
    found.
    Since there are only 4 samples and 2 features This test should work in the majority of cases.
**/
    void RFsetTest()
    {
        //xor dataset.
        double features[] = {0, 0, 1, 1,
                         0, 1, 0, 1};
        int    labels[] = {1, 0, 0, 1};
        std::cerr << "RFsetTest(): Learning 1200 Trees on XOR problem. ";
        {
            vigra::SetTestVisitor testVisitor;
            vigra::RandomForest<> RF2(vigra::RandomForestOptions().tree_count(1200));
            RF2.learn(  MultiArrayView<2, double>(MultiArrayShape<2>::type(4,2), features),
                        MultiArrayView<2, int>(MultiArrayShape<2>::type(4,1), labels),
						create_visitor(testVisitor),
						rf_default(),
						rf_default(),
                        vigra::RandomTT800::global());

        }
        std::cerr << "DONE!\n";
        std::cerr << "RFsetTest(): Comparing with Working Version:";
            diffOnfiles("data/oldsetTest.log", "setTest.log");
        std::cerr << "DONE!\n\n";
    }

    void RFonlineTest()
    {
        //xor dataset.
        double features[] = {0, 0, 1, 1,
                         0, 1, 0, 1};
        int    labels[] = {1, 0, 0, 1};
        double features2[] = {0.4, 0.4, 0.6, 0.6,
                         0.4, 0.6, 0.4, 0.6};
        int    labels2[] = {1, 0, 0, 1};
        std::cerr << "RFsetTest(): Learning 1200 Trees on online problem. ";
        {
            vigra::SetTestVisitor testVisitor;
            vigra::RandomForest<> RF2(vigra::RandomForestOptions().tree_count(1200).prepare_online_learning(true));
            RF2.learn(  MultiArrayView<2, double>(MultiArrayShape<2>::type(4,2), features),
                        MultiArrayView<2, int>(MultiArrayShape<2>::type(4,1), labels),
						create_visitor(testVisitor),
						rf_default(),
						rf_default(),
                        vigra::RandomTT800::global());
            RF2.onlineLearn( MultiArrayView<2, double>(MultiArrayShape<2>::type(4,2), features),
                        MultiArrayView<2, int>(MultiArrayShape<2>::type(4,1), labels),  4);

        }
        std::cerr << "DONE!\n";
    }

/**
        ClassifierTest::RFnoiseTest():
    Learns The Refactored Random Forest with 100 Trees default options and random Seed for 64 dimensional
    unformily distributed samples
    oob error should be 0.5
**/
    void RFnoiseTest()
    {
        typedef vigra::MultiArrayShape<2>::type DiffT;
        vigra::MultiArray<2, double> features(DiffT(50000, 8));
        vigra::MultiArray<2, size_t> labels(DiffT(50000, 1), size_t(0));
        labels.subarray(DiffT(25000,0), DiffT(50000,1)).init(1);
        std::for_each(features.begin(), features.end(),UnaryRandomFunctor());
        std::cerr << "RFnoiseTest(): Learning 100 Trees on Noise.";
        vigra::RandomForest<> RF2(vigra::RandomForestOptions().tree_count(100));
        rf::visitors::OOB_PerTreeError oob_v;
        RF2.learn(features, labels,create_visitor(oob_v));
        std::cerr << "DONE!\n";
        std::cerr << "RFnoiseTest(): Comparing oob with 0.5:";
        shouldEqualTolerance(oob_v.oobError, 0.5, 0.01);
        std::cerr << "DONE!\n\n";
    }

	void RFvariableImportanceTest()
	{
		double pina_var_imp[] = 
		{
			0.017263, 0.040776, 0.003548, 0.003463, 0.005085, 0.015100, 0.005815, 0.019693, 
			0.000555, 0.034199, 0.000093, 0.001263, 0.000669, 0.014896, 0.002777, 0.007323, 
			0.017818, 0.074975, 0.003641, 0.004726, 0.005754, 0.029996, 0.008591, 0.027016, 
			13.743281, 48.682308, 15.098506, 10.868249, 11.145719, 29.414823, 22.270783, 23.060834 
		};

		vigra::MultiArrayView<2, double> p_imp(MultiArrayShape<2>::type(8, 4), pina_var_imp);
		vigra::MultiArray<2, double> zero(p_imp.shape(), 0.0);
        //Create Test output by Random Forest
        {
            std::cerr << "RFvariableImportanceTest(): "
							"Learning on Datasets\n";

			int ii = data.size() - 3; // this is the pina_indians dataset
            {
            	vigra::rf::visitors::VariableImportanceVisitor var_imp;

                vigra::RandomForest<>
					RF2(vigra::RandomForestOptions().tree_count(255));
                RF2.learn(  data.features(ii),
                            data.labels(ii),
						    create_visitor(var_imp),
							rf_default(),
							rf_default(),
                            vigra::RandomMT19937(1));
				
				var_imp.variable_importance_ -= p_imp;
				for(int jj = 0; jj < p_imp.shape(0);  ++jj)
					for(int gg = 0; gg < p_imp.shape(1); ++gg)
						shouldEqualTolerance(var_imp
											   .variable_importance_(jj, gg)
											 , 0.0,0.0001);
				std::cerr << std::endl;
                std::cerr << "[";
                for(int ss = 0; ss < ii+1; ++ss)
                    std::cerr << "#";
                for(int ss = ii+1; ss < data.size(); ++ss)
                    std::cerr << " ";
                std::cerr << "] " << data.names(ii);
                std::cerr << "\n";
            }
        }
        std::cerr << std::endl;
        std::cerr << "DONE!\n\n";
	}

    void RFwrongLabelTest()
    {
        double rawfeatures [] = 
        {
            0, 1, 1, 1, 1,
            0, 1, 1, 1, 1, 
        };
        double rawlabels [] =
        {
            0, 1, 1, 0, 0
        };
        typedef MultiArrayShape<2>::type Shp;
        MultiArrayView<2, double> features(Shp(5, 2), rawfeatures);
        MultiArrayView<2, double> labels(Shp(5, 1), rawlabels);
       
        RandomForest<> rf(RandomForestOptions().tree_count(10).sample_with_replacement(false));
        rf.learn(features, labels);

        MultiArray<2, double> prob(Shp(1, 2));
        rf.predictProbabilities(rowVector(features, 2), prob);

        shouldEqual(prob[0], prob[1]);
        shouldEqual(prob[1], 0.5);

    }
#ifdef HasHDF5
	/** checks whether hdf5 import export is working
	 */
    void HDF5ImpexTest()
    {

            for(int ii = 0; ii <data.size() ; ii++)
            {
				std::cerr << "Running HDF5 Impex Test\n";
				std::string filename = data.names(ii) + "_rf.hdf5";
				std::remove(filename.c_str());
				vigra::RandomForest<> RF(vigra::RandomForestOptions()
											 .tree_count(100));

				 RF.learn(data.features(ii), data.labels(ii));
				 rf_export_HDF5(RF, data.names(ii) + "_rf.hdf5");

				 vigra::RandomForest<> RF2;
				 rf_import_HDF5(RF2,data.names(ii) + "_rf.hdf5");

				 should(RF.ext_param_== RF2.ext_param_);
				 should(RF.options_ ==  RF2.options_);
				 shouldEqual(RF.trees_.size(), RF2.trees_.size());
				 for(int jj = 0; jj < int(RF.trees_.size()); ++jj)
				 {
					 should(RF.trees_[jj].topology_ == 
							RF2.trees_[jj].topology_);

					 should(RF.trees_[jj].parameters_ ==
							RF2.trees_[jj].parameters_);

				 }
                 std::cerr << "[";
                 for(int ss = 0; ss < ii+1; ++ss)
                     std::cerr << "#";
                 for(int ss = ii+1; ss < data.size(); ++ss)
                     std::cerr << " ";
                 std::cerr << "] " << data.names(ii);
                 std::cerr << "\n";
            }
            std::cerr << "done!\n";
    }
#endif
};


struct ClassifierTestSuite
: public vigra::test_suite
{
    ClassifierTestSuite()
    : vigra::test_suite("ClassifierTestSuite")
    {
	
        add( testCase( &ClassifierTest::RFdefaultTest));
#ifndef FAST
        add( testCase( &ClassifierTest::RFsetTest));
        add( testCase( &ClassifierTest::RFonlineTest));
        add( testCase( &ClassifierTest::RFoobTest));
        add( testCase( &ClassifierTest::RFnoiseTest));
        add( testCase( &ClassifierTest::RFvariableImportanceTest));
        add( testCase( &ClassifierTest::RF_NanCheck));
        add( testCase( &ClassifierTest::RF_InfCheck));
        add( testCase( &ClassifierTest::RF_SpliceTest));
        add( testCase( &ClassifierTest::RF_AlgorithmTest));
#endif
        add( testCase( &ClassifierTest::RFresponseTest));
		
        add( testCase( &ClassifierTest::RFridgeRegressionTest));
#ifdef HasHDF5
		add( testCase( &ClassifierTest::HDF5ImpexTest));
#endif
    }
};


#if CLASSIFIER_TEST
int main(int argc, char ** argv)
{
    ClassifierTestSuite test;

    int failed = test.run(vigra::testsToBeExecuted(argc, argv));

    std::cout << test.report() << std::endl;
    return (failed != 0);
}
#endif

#if CROSSVAL
int main(int argc, char ** argv)
{
    RF_Test_Training_Data data;

    std::vector<size_t> indices;
    for(int ii = 1; ii < argc; ++ii)
    {
        indices.push_back(atoi(argv[ii]));
    }
           //int ii = 1;
    for(int ii = 0; ii < data.size(); ii++)
    {
        if(indices.size() != 0)
            if(find(indices.begin(), indices.end(), ii) == indices.end())
                continue;
        MultiArray<3, double> err(MultiArrayShape<3>::type(51, 4, 5));
        for(int ss = 0; ss <= 50; ++ss)
        {
            for(int jj = 0; jj < 5; ++jj)
            {
                vigra::RandomForest<int> RF2(vigra::RandomForestOptions().setTreeCount(255).trainingSetSizeProportional(double(ss)*0.02));
                vigra::RandomForest<int> RF(vigra::RandomForestOptions().setTreeCount(255).trainingSetSizeProportional(double(ss)*0.02).sampleWithReplacement(false));
                TinyVector<double, 4> vc2 = crossvalidate( data.features(ii),
                                                         data.labels(ii),
                                                         RF2);
                TinyVector<double, 4> vc = crossvalidate( data.features(ii),
                                                         data.labels(ii),
                                                         RF);
                err(ss,0 ,jj) = vc[0];
                err(ss,1 ,jj) = vc[3];
                err(ss,2 ,jj) = vc2[0];
                err(ss,3 ,jj) = vc2[3];
                std::cerr << jj << " " << ss << " " << ii << std::endl;
            }
        }
        std::string name = "gini_bootstrap_" + data.names(ii) + ".hdf5";
        writeToHDF5File(name.c_str(), "results", err);
    }

}
#endif

#if HDF5
int main(int argc, char ** argv)
{
    RF_Test_Training_Data data;

    for(int ii = 0; ii < data.size(); ii++)
    {
		H5_dataset set_t(data.names(ii)+std::string(".h5"));
		set_t.set_source("1","features", data.features(ii));
		set_t.set_source("1","labels", data.labels(ii));
	}

    for(int ii = 0; ii < data.size(); ++ii)
    {
		H5_dataset set_t(data.names(ii)+std::string(".h5"));
		std::set<std::string>   grpses = set_t.get_groups();
		grpses.erase("labels");
		std::string a = "labels";
		std::set<std::string>  srces = set_t.get_sources();
		MultiArray<2, double> feats(set_t.shape(srces.begin(), srces.end(), grpses.begin(), grpses.end()));
		MultiArray<2, int> labels(set_t.shape(srces.begin(), srces.end(), &a, &a +1));
		set_t.get_source("1","features", feats);
		set_t.get_source("1","labels", labels);
		std::cerr << (feats == data.features(ii)) <<	std::endl;
		std::cerr << (labels == data.labels(ii)) <<	std::endl;
	}
}
#endif
