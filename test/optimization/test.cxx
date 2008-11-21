/************************************************************************/
/*                                                                      */
/*                  Copyright 2008 by Ullrich Koethe                    */
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

#include "unittest.hxx"

#include <iostream>
#include <iterator>
#include <algorithm>

#include "vigra/matrix.hxx"
#include "vigra/regression.hxx"

#include "larsdata.hxx"

using namespace vigra;
using namespace vigra::linalg;

struct OptimizationTest {
  	typedef vigra::linalg::Matrix<double>::difference_type Shape;
  	
  	ArrayVector<Matrix<double> > x, y, lars, larslsq, lasso, lassolsq, nnlasso, nnlassolsq;
  	int size;
  	static double w[100];

    OptimizationTest()
    {
        loadData(x, y, lars, larslsq, lasso, lassolsq, nnlasso, nnlassolsq);
        size = x.size();
    }

    void testLarsImpl(LeastAngleRegressionOptions const & options, 
                       Matrix<double> X, Matrix<double> y, 
                       Matrix<double> const & ref, const char * message)
    {
        double epsilon = 1e-10;
        
	    ArrayVector<Matrix<double> > results;
	    ArrayVector<ArrayVector<int> > activeSets;
    	
	    Matrix<double> offset(1,50);
	    Matrix<double> scaling(1,50, 1.0); // trivial init

	    prepareColumns(y, y, DataPreparationGoals(ZeroMean));
	    prepareColumns(X, X, offset, scaling, DataPreparationGoals(ZeroMean|UnitVariance));

	    int numSolutions = leastAngleRegression(X, y, activeSets, results, options);
	    shouldMsg(numSolutions == ref.columnCount(), (std::string("wrong number of solutions in ") + message).c_str());
	    
	    for (MultiArrayIndex j = 0; j < numSolutions; ++j) 
	    {
		    Matrix<double> B(50, 1);
		    for (unsigned int i = 0; i < activeSets[j].size(); ++i) 
		    {
			    // activeSets[j][i] is the true index of the i-th result 
			    B(activeSets[j][i],0) = results[j](i,0)*scaling(0, activeSets[j][i]);
		    }
		    std::ostringstream s;
		    s << "solution " << j << " differs in " << message;
		    shouldMsg((B - columnVector(ref, j)).norm(0) < epsilon, s.str().c_str());
		}
    }

#define VIGRA_TEST_LARS(options, k, name, msg) \
    testLarsImpl(options, x[k], y[k], name[k], msg);

    void testLSQ()
    {
        double epsilon = 1e-10;
        char * methods[3] = { "QR", "SVD", "NE" };
        
        for(int m=0; m<3; ++m)
        {
		    for(int k=0; k<size; ++k)
		    {
	            Matrix<double> result(50, 1);
                leastSquares(x[k], y[k], result, methods[m]);
    	        
	            // check the KKT conditions
	            Matrix<double> r = transpose(x[k])*(y[k] - x[k]*result);
	            r /= x[k].norm(2)*y[k].norm(2);
	            std::ostringstream s;
	            s << "failure in problem " << k << " of LSQ test with " << methods[m] << " solver";
	            shouldMsg(r.norm(0) < epsilon, s.str().c_str());
		    }
		}
    }

    void testWeightedLSQ()
    {
        double epsilon = 1e-10;
        Matrix<double> weights(100, 1, w);
                                
        char * methods[3] = { "QR", "SVD", "NE" };
        
        for(int m=0; m<3; ++m)
        {
		    for(int k=0; k<size; ++k)
		    {
	            Matrix<double> result(50, 1);
                weightedLeastSquares(x[k], y[k], weights, result, methods[m]);
    	        
	            // check the KKT conditions
	            Matrix<double> r = transpose(x[k])*(weights*pointWise(y[k] - x[k]*result));
	            r /= x[k].norm(2)*y[k].norm(2);
	            std::ostringstream s;
	            s << "failure in problem " << k << " of weighted LSQ test with " << methods[m] << " solver";
	            shouldMsg(r.norm(0) < epsilon, s.str().c_str());
		    }
		}
    }

    void testRidgeRegression()
    {
        double epsilon = 1e-10;
        double la[4] = {0.01, 1.0, 100.0, 10000.0};
        ArrayVector<double> lambdas(la, la+4);
        Matrix<double> weights(100, 1, w);
        
	    for(int k=0; k<size; ++k)
	    {
            Matrix<double> results(50, 4);
            ridgeRegressionSeries(x[k], y[k], results, lambdas);
            for(int m=0; m<4; ++m)
            {
	            Matrix<double> result(50, 1);
                ridgeRegression(x[k], y[k], result, lambdas[m]);
    	        
	            // check the KKT conditions
	            Matrix<double> r = transpose(x[k])*(y[k] - x[k]*result) - lambdas[m]*result;
	            r /= x[k].norm(2)*y[k].norm(2);
	            std::ostringstream s;
	            s << "failure in problem " << k << " of ridge regression test";
	            shouldMsg(r.norm(0) < epsilon, s.str().c_str());
	            shouldMsg((result - columnVector(results, m)).norm(0) < epsilon, s.str().c_str());
		    }
            for(int m=0; m<4; ++m)
            {
	            Matrix<double> result(50, 1);
                weightedRidgeRegression(x[k], y[k], weights, result, lambdas[m]);
    	        
	            // check the KKT conditions
	            Matrix<double> r = transpose(x[k])*(weights*pointWise(y[k] - x[k]*result)) - lambdas[m]*result;
	            r /= x[k].norm(2)*y[k].norm(2);
	            std::ostringstream s;
	            s << "failure in problem " << k << " of weighted ridge regression test";
	            shouldMsg(r.norm(0) < epsilon, s.str().c_str());
		    }
		}
    }

    void testLars()
    {
        LeastAngleRegressionOptions larsOptions;
		larsOptions = larsOptions.lars().leastSquaresSolutions(false);
		
		for(int k=0; k<size; ++k)
		{
		    std::ostringstream s;
		    s << "lars " << k;
    		VIGRA_TEST_LARS(larsOptions, k, lars, s.str().c_str());
		}
    }

    void testLarsLSQ()
    {
        LeastAngleRegressionOptions larsOptions;
		larsOptions = larsOptions.lars().leastSquaresSolutions(true);
		
		for(int k=0; k<size; ++k)
		{
		    std::ostringstream s;
		    s << "larslsq " << k;
    		VIGRA_TEST_LARS(larsOptions, k, larslsq, s.str().c_str());
		}
    }
    
    void testLasso()
    {
        LeastAngleRegressionOptions larsOptions;
		larsOptions = larsOptions.lasso().leastSquaresSolutions(false);
		
		for(int k=0; k<size; ++k)
		{
		    std::ostringstream s;
		    s << "lasso " << k;
    		VIGRA_TEST_LARS(larsOptions, k, lasso, s.str().c_str());
		}
    }

    void testLassoLSQ()
    {
        LeastAngleRegressionOptions larsOptions;
		larsOptions = larsOptions.lasso().leastSquaresSolutions(true);
		
		for(int k=0; k<size; ++k)
		{
		    std::ostringstream s;
		    s << "lassolsq " << k;
    		VIGRA_TEST_LARS(larsOptions, k, lassolsq, s.str().c_str());
		}
    }
    
    void testNNLasso()
    {
        LeastAngleRegressionOptions larsOptions;
		larsOptions = larsOptions.nnlasso().leastSquaresSolutions(false);
		
		for(int k=0; k<size; ++k)
		{
		    std::ostringstream s;
		    s << "nnlasso " << k;
    		VIGRA_TEST_LARS(larsOptions, k, nnlasso, s.str().c_str());
		}
    }

    void testNNLassoLSQ()
    {
        LeastAngleRegressionOptions larsOptions;
		larsOptions = larsOptions.nnlasso().leastSquaresSolutions(true);
		
		for(int k=0; k<size; ++k)
		{
		    std::ostringstream s;
		    s << "nnlassolsq " << k;
    		VIGRA_TEST_LARS(larsOptions, k, nnlassolsq, s.str().c_str());
		}
    }

    void testNNLSQ()
    {
        double epsilon = 1e-10;
		for(int k=0; k<size; ++k)
		{
	        Matrix<double> result(50, 1);
            nonnegativeLeastSquares(x[k], y[k], result);
	        
	        // check the KKT conditions
	        Matrix<double> r = transpose(x[k])*(y[k] - x[k]*result);
	        r /= x[k].norm(2)*y[k].norm(2);
	        
	        std::ostringstream s;
	        s << "failure in problem " << k << " of NNLSQ test";
	        for(int l=0; l<50; ++l)
    	        shouldMsg((r(l,0) <= 0.0 && result(l,0) == 0.0) || (abs(r(l,0)) < epsilon && result(l,0) > 0.0), s.str().c_str());
		}
    }
};

double OptimizationTest::w[100] =
               {3.099060, 2.241629, 0.551814, 9.882278, 9.069642, 6.263997, 14.999678, 7.299786,
                6.164468, 0.876546, 14.982520, 1.723430, 0.083692, 0.054570, 3.551352, 0.602337,
                3.629894, 4.573182, 1.967879, 1.218648, 0.517033, 15.419422, 1.523039, 11.806089, 
                1.568337, 0.590645, 0.180897, 0.228799, 13.380690, 1.392058, 8.199124, 0.036260, 
                12.429014, 0.588538, 4.029622, 1.356536, 11.249162, 9.584928, 1.193385, 0.564644, 
                14.343903, 0.772126, 2.372715, 0.865896, 3.259700, 6.267629, 2.852698, 9.231133, 
                0.293141, 14.956796, 7.291935, 11.719283, 0.020900, 0.698194, 0.641376, 7.234362, 
                2.282460, 0.124404, 9.981096, 0.000032, 3.045198, 2.892267, 1.918650, 4.422195, 
                0.012394, 15.208476, 8.253837, 15.333965, 11.425123, 12.661551, 14.815347, 
                0.009620, 0.024826, 0.036806, 10.438017, 0.148247, 4.152140, 0.835896, 1.313798, 
                0.062187, 10.299775, 7.222942, 5.173045, 14.041757, 0.666464, 7.385935, 3.626625, 
                15.613323, 3.012300, 14.823977, 10.285900, 1.518461, 6.132636, 0.381827,
                10.990995, 0.087974, 10.590009, 0.481253, 5.965124, 5.592337};



struct OptimizationTestSuite : public vigra::test_suite {
	OptimizationTestSuite() : vigra::test_suite("Optimization Tests") {
		add( testCase(&OptimizationTest::testLSQ));
		add( testCase(&OptimizationTest::testWeightedLSQ));
		add( testCase(&OptimizationTest::testRidgeRegression));
		add( testCase(&OptimizationTest::testLars));
		add( testCase(&OptimizationTest::testLarsLSQ));
		add( testCase(&OptimizationTest::testLasso));
		add( testCase(&OptimizationTest::testLassoLSQ));
		add( testCase(&OptimizationTest::testNNLasso));
		add( testCase(&OptimizationTest::testNNLassoLSQ));
		add( testCase(&OptimizationTest::testNNLSQ));
	}
};

int main() {
	OptimizationTestSuite test;
	int failed = test.run();
	std::cout << test.report() << std::endl;
	return failed;
}


