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

struct LarsTest {
  	typedef vigra::linalg::Matrix<double>::difference_type Shape;
  	
  	ArrayVector<Matrix<double> > x, y, lars, larslsq, lasso, lassolsq, nnlasso, nnlassolsq;
  	int size;

    LarsTest()
    {
        loadData(x, y, lars, larslsq, lasso, lassolsq, nnlasso, nnlassolsq);
        size = x.size();
    }

    void testLarsImpl(LeastAngleRegressionOptions const & options, 
                       Matrix<double> X, Matrix<double> y, 
                       Matrix<double> const & ref, const char * message)
    {
        double espilon = 1e-10;
        
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
		    shouldMsg((B - columnVector(ref, j)).norm(0) < espilon, s.str().c_str());
		}
    }

#define VIGRA_TEST_LARS(options, k, name, msg) \
    testLarsImpl(options, x[k], y[k], name[k], msg);

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
};


struct OptimizationTestSuite : public vigra::test_suite {
	OptimizationTestSuite() : vigra::test_suite("Optimization Tests") {
		add( testCase(&LarsTest::testLars));
		add( testCase(&LarsTest::testLarsLSQ));
		add( testCase(&LarsTest::testLasso));
		add( testCase(&LarsTest::testLassoLSQ));
		add( testCase(&LarsTest::testNNLasso));
		add( testCase(&LarsTest::testNNLassoLSQ));
	}
};

int main() {
	OptimizationTestSuite test;
	int failed = test.run();
	std::cout << test.report() << std::endl;
	return failed;
}


