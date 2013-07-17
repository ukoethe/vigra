/************************************************************************/
/*                                                                      */
/*                  Copyright 2008 by Ullrich Koethe                    */
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

#include "unittest.hxx"

#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>

#include "vigra/matrix.hxx"
#include "vigra/regression.hxx"
#include "vigra/quadprog.hxx"

#include "larsdata.hxx"

using namespace vigra;
using namespace vigra::linalg;

struct GaussianFitModel
{
    // pass in scalar data
    template <class T>
    T operator()(double data, TinyVector<T, 3> const & p) const
    {
        return p[1] * exp(-0.5 * sq(data / p[0])) + p[2];
    }
    
    // pass in vector data
    template <class T>
    T operator()(MultiArrayView<1, double> const & data, TinyVector<T, 3> const & p) const
    {
        return p[1] * exp(-0.5 * sq(data(0) / p[0])) + p[2];
    }
};


struct OptimizationTest 
{
    typedef vigra::linalg::Matrix<double>::difference_type Shape;

    ArrayVector<Matrix<double> > x, y, lars, larslsq, lasso, lassolsq, nnlasso, nnlassolsq;
    int size;

    static double w[100];
    static double gaussianFitData[400];

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
        ArrayVector<ArrayVector<MultiArrayIndex> > activeSets;

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
        const char * methods[3] = { "QR", "SVD", "NE" };

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

        const char * methods[3] = { "QR", "SVD", "NE" };

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

    void testNonlinearLSQ()
    {
        int size = 200;

        // data are generated from a * exp(-0.5*sq(x/b)) + c + random.normal(0, 0.05)
        // with (a, b, c) = (1, 1, 0.5)
        Matrix<double> points(Shape2(size, 2), gaussianFitData);
        MultiArrayView<1, double> features(points.bindOuter(0)),
                                  response(points.bindOuter(1));
        {
            TinyVector<double, 3> p(2.0, 0.5, 0.0);
            // test 1D feature vector
            nonlinearLeastSquares(features, response, p, GaussianFitModel());
            shouldEqualTolerance(p[0], 0.99661454591763432, 1e-14);
            shouldEqualTolerance(p[1], 1.0070468569580588, 1e-14);
            shouldEqualTolerance(p[2], 0.49323974179926072, 1e-14);
        }
        {
            TinyVector<double, 3> p(2.0, 0.5, 0.0);
            // test 2D feature matrix (which happens to have just one column)
            nonlinearLeastSquares(points.subarray(Shape2(), Shape2(size,1)), response, p, GaussianFitModel());
            shouldEqualTolerance(p[0], 0.99661454591763432, 1e-14);
            shouldEqualTolerance(p[1], 1.0070468569580588, 1e-14);
            shouldEqualTolerance(p[2], 0.49323974179926072, 1e-14);
        }
    }

    void testQuadProg()
    {
        {
            double Gdata[] = { 2.1, 0.0, 1.0,
                               0.0, 2.2, 0.0,
                               1.0, 0.0, 3.1};

            double gdata[] = {6.0, 1.0, 1.0};

            double CEdata[] = {1.0, 2.0, -1.0};

            double cedata[] = {4.0};

            double CIdata[] = { 1.0,  0.0, 0.0,
                                0.0,  1.0, 0.0,
                                0.0,  0.0, 1.0,
                               -1.0, -1.0, 0.0};

            double cidata[] = {0.0, 0.0, 0.0, -10.0};

            double xrefdata[] = {0.0, 2.0, 0.0};

            Matrix<double> G(3,3, Gdata),
                           g(3,1, gdata),
                           CE(1,3, CEdata),
                           ce(1,1, cedata),
                           CI(4,3, CIdata),
                           ci(4,1, cidata),
                           x(3,1),
                           xref(3,1, xrefdata);

            // Solution: x = [0 2 0]'
            shouldEqualTolerance(quadraticProgramming(G, g, CE, ce, CI, ci, x), 6.4, 1e-14);
            shouldEqualSequenceTolerance(x.data(), x.data()+3, xrefdata, 1e-14);
        }
        {
            double Gdata[] = {13.0, 12.0, -2.0,
                              12.0, 17.0,  6.0,
                              -2.0,  6.0, 12.0};

            double gdata[] = {-22.0, -14.5, 13.0};

            double CIdata[] = { 1.0,  0.0,  0.0,
                                0.0,  1.0,  0.0,
                                0.0,  0.0,  1.0,
                                0.0,  0.0,  1.0,  // add a redundant constraint
                               -1.0,  0.0,  0.0,
                                0.0, -1.0,  0.0,
                                0.0,  0.0, -1.0};

            double cidata[] = {-1.0, -1.0, -1.0, -2.0, -1.0, -1.0, -1.0};

            double xrefdata[] = {1.0, 0.5, -1.0};

            Matrix<double> G(3,3, Gdata),
                           g(3,1, gdata),
                           CE,
                           ce,
                           CI(7,3, CIdata),
                           ci(7,1, cidata),
                           x(3,1),
                           xref(3,1, xrefdata);

            // Solution: x = [1 0.5 -1]'
            shouldEqualTolerance(quadraticProgramming(G, g, CE, ce, CI, ci, x), -22.625, 1e-14);
            shouldEqualSequenceTolerance(x.data(), x.data()+3, xrefdata, 1e-14);
        }
        {
            double Gdata[] = {11.9557487988864, 2.37931476086954, -0.376766571133756, 0.223144794097961, -2.05504905104678, -3.64568978396075, -5.59430271562319, -1.69330364188253,
                               2.37931476086954, 17.1587043263327, 1.74450492549782, 6.72064470789053, 2.58420114935295, 3.98719003404737, -4.68188629728773, -0.895269022607405,
                              -0.376766571133756, 1.74450492549782, 15.2278954316710, -4.29485935734649, -8.44276270167430, -1.41513468539932, 0.888903479635589, -2.05250235130216,
                               0.223144794097961, 6.72064470789053, -4.29485935734649, 14.3556037317358, 3.81380696097931, -1.80980093568364, -2.37884660765149, 0.837452262645372,
                              -2.05504905104678, 2.58420114935295, -8.44276270167430, 3.81380696097931, 11.6645675398446, 4.09567324453584, -1.30279988467841, -0.235847705894946,
                              -3.64568978396075, 3.98719003404737, -1.41513468539932, -1.80980093568364, 4.09567324453584, 19.1274437292212, -1.29389377140431, -5.21988129175308,
                              -5.59430271562319, -4.68188629728773, 0.888903479635589, -2.37884660765149, -1.30279988467841, -1.29389377140431, 15.2120567578534, 1.62881095102307,
                              -1.69330364188253, -0.895269022607405, -2.05250235130216, 0.837452262645372, -0.235847705894946, -5.21988129175308, 1.62881095102307, 9.35509425897139};

            double gdata[] = {1.69056816359673, -5.99644841532114, -2.18039280131329, 7.36068310570474, -1.56285526687905, 0.556686337262202, -1.39556850584431, -4.21389206565463};

            double CEdata[] = {-0.941710232692403, 1.62075129330797, -0.635514086344284, -2.04726855027881, -0.0879725045073330, -0.0439787231114959, -0.214829528926237, -0.746695232460116,
                                0.384997223212006, -3.05182514896839, -1.02873559051701, -1.12930451304791, 1.07377722353829, -0.799868398167378, 0.00731478600822687, 0.349275518947309,
                               -0.278886986977952, -0.0484535232017660, 1.64138030542441, -2.35558595029120, -0.311909083814529, -0.865157920346768, -1.03947227889259, 0.484013190408825,
                               -0.982943684246934, 0.318202298763980, 0.0194950686266762, -0.561248752100947, -1.47877377059869, -0.119006984268877, 0.832835771260238, -1.00785901222435};

            double cedata[] = {1.00346920730319, -2.67608891438430, 0.0168223926638857, -1.44324544697726};

            double xrefdata[] = {0.11925229791508, 1.15958300340908, 0.14739510825795, 0.00000000000000, 0.72868950990454, 0.09811723636230, 0.18105311441297, 0.75350474556599};

            Matrix<double> G(8,8, Gdata),
                           g(8,1, gdata),
                           CE(4,8, CEdata),
                           ce(4,1, cedata),
                           CI(identityMatrix<double>(8)),
                           ci(8,1), x(8,1),
                           xref(8,1, xrefdata);

            shouldEqualTolerance(quadraticProgramming(G, g, CE, ce, CI, ci, x), 5.96700441471631, 1e-10);
            shouldEqualSequenceTolerance(x.data(), x.data()+8, xrefdata, 1e-10);
        }
        double epsilon = 1e-10;
        for(int k=0; k<size; ++k)
        {
            Matrix<double> result(50, 1), ref(50,1);
            Matrix<double> G = transpose(x[k])*x[k],
                           g = -transpose(x[k])*y[k];
            double m = quadraticProgramming(G, g, Matrix<double>(), Matrix<double>(),
                                 identityMatrix<double>(50), Matrix<double>(50, 1), result);
            nonnegativeLeastSquares(x[k], y[k], ref);
            shouldEqualTolerance(2.0*m, squaredNorm(x[k]*ref-y[k]) - squaredNorm(y[k]), epsilon);
            shouldEqualSequenceTolerance(ref.data(), ref.data()+50, result.data(), epsilon);
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

double OptimizationTest::gaussianFitData[400] = {
-4.000000, 0.444541,
-3.959799, 0.538727,
-3.919598, 0.518276,
-3.879397, 0.412113,
-3.839196, 0.518404,
-3.798995, 0.541461,
-3.758794, 0.503801,
-3.718593, 0.491741,
-3.678392, 0.460771,
-3.638191, 0.429009,
-3.597990, 0.541560,
-3.557789, 0.486328,
-3.517588, 0.490383,
-3.477387, 0.589003,
-3.437186, 0.536945,
-3.396985, 0.521662,
-3.356784, 0.510677,
-3.316583, 0.580087,
-3.276382, 0.590646,
-3.236181, 0.551794,
-3.195980, 0.535165,
-3.155779, 0.402148,
-3.115578, 0.513988,
-3.075377, 0.502330,
-3.035176, 0.514688,
-2.994975, 0.558430,
-2.954774, 0.375726,
-2.914573, 0.485836,
-2.874372, 0.529562,
-2.834171, 0.494678,
-2.793970, 0.449333,
-2.753769, 0.566007,
-2.713568, 0.539023,
-2.673367, 0.479504,
-2.633166, 0.546958,
-2.592965, 0.575754,
-2.552764, 0.538718,
-2.512563, 0.582603,
-2.472362, 0.550975,
-2.432161, 0.532177,
-2.391960, 0.499256,
-2.351759, 0.558655,
-2.311558, 0.578853,
-2.271357, 0.619601,
-2.231156, 0.577235,
-2.190955, 0.613576,
-2.150754, 0.550746,
-2.110553, 0.568697,
-2.070352, 0.611762,
-2.030151, 0.574626,
-1.989950, 0.679089,
-1.949749, 0.672611,
-1.909548, 0.675465,
-1.869347, 0.691201,
-1.829146, 0.788757,
-1.788945, 0.678421,
-1.748744, 0.606669,
-1.708543, 0.742304,
-1.668342, 0.746127,
-1.628141, 0.739815,
-1.587940, 0.734492,
-1.547739, 0.779914,
-1.507538, 0.830060,
-1.467337, 0.815632,
-1.427136, 0.981810,
-1.386935, 0.834181,
-1.346734, 0.864140,
-1.306533, 0.811485,
-1.266332, 0.961096,
-1.226131, 0.870745,
-1.185930, 0.968019,
-1.145729, 1.004959,
-1.105528, 1.007269,
-1.065327, 1.153907,
-1.025126, 1.141014,
-0.984925, 1.181630,
-0.944724, 1.095902,
-0.904523, 1.220690,
-0.864322, 1.213103,
-0.824121, 1.250634,
-0.783920, 1.286928,
-0.743719, 1.212951,
-0.703518, 1.259559,
-0.663317, 1.345653,
-0.623116, 1.190765,
-0.582915, 1.419420,
-0.542714, 1.390718,
-0.502513, 1.379101,
-0.462312, 1.409671,
-0.422111, 1.363268,
-0.381910, 1.412171,
-0.341709, 1.498303,
-0.301508, 1.520465,
-0.261307, 1.601247,
-0.221106, 1.472156,
-0.180905, 1.450842,
-0.140704, 1.464438,
-0.100503, 1.444060,
-0.060302, 1.494291,
-0.020101, 1.518935,
0.020101, 1.498086,
0.060302, 1.553001,
0.100503, 1.483252,
0.140704, 1.472778,
0.180905, 1.454707,
0.221106, 1.394221,
0.261307, 1.388047,
0.301508, 1.396606,
0.341709, 1.508361,
0.381910, 1.474431,
0.422111, 1.483513,
0.462312, 1.332036,
0.502513, 1.282955,
0.542714, 1.330059,
0.582915, 1.352545,
0.623116, 1.348480,
0.663317, 1.354922,
0.703518, 1.294989,
0.743719, 1.345523,
0.783920, 1.224325,
0.824121, 1.166410,
0.864322, 1.104242,
0.904523, 1.119812,
0.944724, 1.152129,
0.984925, 1.071238,
1.025126, 1.138132,
1.065327, 1.137580,
1.105528, 0.924277,
1.145729, 1.061945,
1.185930, 0.883011,
1.226131, 0.991640,
1.266332, 1.009765,
1.306533, 0.929159,
1.346734, 0.839811,
1.386935, 0.852935,
1.427136, 0.848105,
1.467337, 0.831661,
1.507538, 0.810848,
1.547739, 0.796379,
1.587940, 0.794108,
1.628141, 0.705262,
1.668342, 0.736556,
1.708543, 0.808252,
1.748744, 0.697509,
1.788945, 0.679672,
1.829146, 0.741614,
1.869347, 0.546297,
1.909548, 0.720579,
1.949749, 0.617860,
1.989950, 0.646273,
2.030151, 0.632174,
2.070352, 0.664405,
2.110553, 0.594448,
2.150754, 0.565075,
2.190955, 0.655597,
2.231156, 0.464782,
2.271357, 0.576826,
2.311558, 0.501742,
2.351759, 0.524873,
2.391960, 0.657789,
2.432161, 0.549708,
2.472362, 0.556816,
2.512563, 0.453497,
2.552764, 0.502001,
2.592965, 0.544503,
2.633166, 0.548955,
2.673367, 0.558903,
2.713568, 0.525611,
2.753769, 0.548909,
2.793970, 0.542868,
2.834171, 0.426533,
2.874372, 0.517917,
2.914573, 0.552697,
2.954774, 0.542204,
2.994975, 0.493085,
3.035176, 0.469709,
3.075377, 0.452920,
3.115578, 0.501249,
3.155779, 0.563532,
3.195980, 0.408463,
3.236181, 0.472325,
3.276382, 0.447676,
3.316583, 0.543335,
3.356784, 0.475859,
3.396985, 0.479589,
3.437186, 0.491873,
3.477387, 0.524637,
3.517588, 0.482437,
3.557789, 0.349477,
3.597990, 0.528711,
3.638191, 0.523288,
3.678392, 0.490176,
3.718593, 0.446792,
3.758794, 0.518444,
3.798995, 0.519696,
3.839196, 0.477128,
3.879397, 0.489703,
3.919598, 0.453953,
3.959799, 0.491464,
4.000000, 0.422814
};

struct OptimizationTestSuite : public vigra::test_suite 
{
    OptimizationTestSuite() : vigra::test_suite("Optimization Tests") 
    {
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
        add( testCase(&OptimizationTest::testNonlinearLSQ));
        add( testCase(&OptimizationTest::testQuadProg));
    }
};

int main() 
{
    OptimizationTestSuite test;
    int failed = test.run();
    std::cout << test.report() << std::endl;
    return failed;
}


