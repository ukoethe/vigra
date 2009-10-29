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
#define CLASSIFIER_TEST 1
#define HDF5 0
#define CROSSVAL 0


#include <iostream>
#include <fstream>
#include <functional>
#include <cmath>
#include <vigra/random_forest.hxx>
#include <unittest.hxx>
#include <vector>
//#include "data/RF_results.hxx"
#include "data/RF_data.hxx"
#include "test_visitors.hxx"

#if CROSSVAL
#include <vigra/crossvalidation.hxx>
#endif


#include <stdlib.h>

#if HDF5
#include <vigra/dataset4.hxx>
#include "vigra/hdf5impex.hxx"
#endif
using namespace vigra;


struct UnaryRandomFunctor
{
    double operator()(double in)
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
            vigra::StopVisiting stop;
            vigra::TestVisitor<StopVisiting> testVisitor(stop);
            std::cerr << "RFdefaultTest(): Learning on Datasets\n";
            for(int ii = 0; ii < data.size() ; ii++)
            {

                vigra::RandomForest<>
					RF2(vigra::RandomForestOptions().tree_count(32));

                RF2.learn(  data.features(ii),
                            data.labels(ii),
							rf_default(),
							rf_default(),
						   	testVisitor,
                            vigra::RandomMT19937(1));

                testVisitor.fout <<  data.names(ii) << std::endl;
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
        std::cerr << "RFdefaultTest(): Comparing with Working Version:";
        //Cheap diff on old and new Classifier.
        diffOnfiles("./data/oldClassifier.log", "RandomForestNodeTest.log");
        std::cerr << "DONE!\n\n";
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
                    vigra::RandomForest<> RF2(vigra::RandomForestOptions()
												 .tree_count(100));

                     oob  += RF2.learn(  data.features(ii),
                                            data.labels(ii));
                    //std::cerr << oob <<" " <<  data.oobError(ii) <<std::endl;
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
            vigra::StopVisiting stop;
            vigra::SetTestVisitor<StopVisiting> testVisitor(stop);
            vigra::RandomForest<> RF2(vigra::RandomForestOptions().tree_count(1200));
            RF2.learn(  MultiArrayView<2, double>(MultiArrayShape<2>::type(4,2), features),
                        MultiArrayView<2, int>(MultiArrayShape<2>::type(4,1), labels),
						rf_default(),
						rf_default(),
						testVisitor,
                        vigra::RandomTT800::global());

        }
        std::cerr << "DONE!\n";
		system("ls");
        std::cerr << "RFsetTest(): Comparing with Working Version:";
            diffOnfiles("data/oldsetTest.log", "setTest.log");
        std::cerr << "DONE!\n\n";
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
        double oob = RF2.learn(features, labels);
        std::cerr << "DONE!\n";
        std::cerr << "RFnoiseTest(): Comparing oob with 0.5:";
        shouldEqualTolerance(oob, 0.5, 0.01);
        std::cerr << "DONE!\n\n";
    }
};


struct ClassifierTestSuite
: public vigra::test_suite
{
    ClassifierTestSuite()
    : vigra::test_suite("ClassifierTestSuite")
    {
        add( testCase( &ClassifierTest::RFdefaultTest));
        add( testCase( &ClassifierTest::RFsetTest));
        add( testCase( &ClassifierTest::RFoobTest));
        add( testCase( &ClassifierTest::RFnoiseTest));
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
