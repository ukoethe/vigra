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


#ifndef VIGRA_REGRESSION_HXX
#define VIGRA_REGRESSION_HXX

#include "matrix.hxx"
#include "linear_solve.hxx"
#include "singular_value_decomposition.hxx"
#include "numerictraits.hxx"
#include "functorexpression.hxx"


namespace vigra
{

namespace linalg
{

/** \addtogroup Optimization Optimization and Regression
 */
//@{
   /** Ordinary Least Squares Regression.

       Given a matrix \a A with <tt>m</tt> rows and <tt>n</tt> columns (with <tt>m \>= n</tt>),
       and a column vector \a b of length <tt>m</tt> rows, this function computes
       the column vector \a x of length <tt>n</tt> rows that solves the optimization problem

        \f[ \tilde \textrm{\bf x} = \textrm{argmin}
            \left|\left|\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right|\right|_2^2
        \f]

       When \a b is a matrix with <tt>k</tt> columns, \a x must also have
       <tt>k</tt> columns, which will contain the solutions for the corresponding columns of
       \a b. Note that all matrices must already have the correct shape.

       This function is just another name for \ref linearSolve(), perhaps
       leading to more readable code when \a A is a rectangular matrix. It returns
       <tt>false</tt> when the rank of \a A is less than <tt>n</tt>.
       See \ref linearSolve() for more documentation.

    <b>\#include</b> \<vigra/regression.hxx\>
        Namespaces: vigra and vigra::linalg
   */
template <class T, class C1, class C2, class C3>
inline bool
leastSquares(MultiArrayView<2, T, C1> const & A,
             MultiArrayView<2, T, C2> const &b, MultiArrayView<2, T, C3> &x,
             std::string method = "QR")
{
    return linearSolve(A, b, x, method);
}

   /** Weighted Least Squares Regression.

       Given a matrix \a A with <tt>m</tt> rows and <tt>n</tt> columns (with <tt>m \>= n</tt>),
       a vector \a b of length <tt>m</tt>, and a weight vector \a weights of length <tt>m</tt>
       with non-negative entries, this function computes the vector \a x of length <tt>n</tt>
       that solves the optimization problem

        \f[ \tilde \textrm{\bf x} = \textrm{argmin}
            \left(\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right)^T
             \textrm{diag}(\textrm{\bf weights})
             \left(\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right)
        \f]

       where <tt>diag(weights)</tt> creates a diagonal matrix from \a weights.
       The algorithm calls \ref leastSquares() on the equivalent problem

        \f[ \tilde \textrm{\bf x} = \textrm{argmin}
             \left|\left|\textrm{diag}(\textrm{\bf weights})^{1/2}\textrm{\bf A} \textrm{\bf x} -
                  \textrm{diag}(\textrm{\bf weights})^{1/2} \textrm{\bf b}\right|\right|_2^2
        \f]

       where the square root of \a weights is just taken element-wise.

       When \a b is a matrix with <tt>k</tt> columns, \a x must also have
       <tt>k</tt> columns, which will contain the solutions for the corresponding columns of
       \a b. Note that all matrices must already have the correct shape.

       The function returns
       <tt>false</tt> when the rank of the weighted matrix \a A is less than <tt>n</tt>.

    <b>\#include</b> \<vigra/regression.hxx\>
        Namespaces: vigra and vigra::linalg
   */
template <class T, class C1, class C2, class C3, class C4>
bool
weightedLeastSquares(MultiArrayView<2, T, C1> const & A,
             MultiArrayView<2, T, C2> const &b, MultiArrayView<2, T, C3> const &weights,
             MultiArrayView<2, T, C4> &x, std::string method = "QR")
{
    typedef T Real;

    const unsigned int rows = rowCount(A);
    const unsigned int cols = columnCount(A);
    const unsigned int rhsCount = columnCount(b);
    vigra_precondition(rows >= cols,
       "weightedLeastSquares(): Input matrix A must be rectangular with rowCount >= columnCount.");
    vigra_precondition(rowCount(b) == rows,
       "weightedLeastSquares(): Shape mismatch between matrices A and b.");
    vigra_precondition(rowCount(b) == rowCount(weights) && columnCount(weights) == 1,
       "weightedLeastSquares(): Weight matrix has wrong shape.");
    vigra_precondition(rowCount(x) == cols && columnCount(x) == rhsCount,
       "weightedLeastSquares(): Result matrix x has wrong shape.");

    Matrix<T> wa(A.shape()), wb(b.shape());

    for(unsigned int k=0; k<rows; ++k)
    {
        vigra_precondition(weights(k,0) >= 0,
           "weightedLeastSquares(): Weights must be positive.");
        T w = std::sqrt(weights(k,0));
        for(unsigned int l=0; l<cols; ++l)
            wa(k,l) = w * A(k,l);
        for(unsigned int l=0; l<rhsCount; ++l)
            wb(k,l) = w * b(k,l);
    }

    return leastSquares(wa, wb, x, method);
}

   /** Ridge Regression.

       Given a matrix \a A with <tt>m</tt> rows and <tt>n</tt> columns (with <tt>m \>= n</tt>),
       a vector \a b of length <tt>m</tt>, and a regularization parameter <tt>lambda \>= 0.0</tt>,
       this function computes the vector \a x of length <tt>n</tt>
       that solves the optimization problem

        \f[ \tilde \textrm{\bf x} = \textrm{argmin}
            \left|\left|\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right|\right|_2^2 +
            \lambda \textrm{\bf x}^T\textrm{\bf x}
        \f]

       This is implemented by means of \ref singularValueDecomposition().

       When \a b is a matrix with <tt>k</tt> columns, \a x must also have
       <tt>k</tt> columns, which will contain the solutions for the corresponding columns of
       \a b. Note that all matrices must already have the correct shape.

       The function returns <tt>false</tt> if the rank of \a A is less than <tt>n</tt>
       and <tt>lambda == 0.0</tt>.

    <b>\#include</b> \<vigra/regression.hxx\>
        Namespaces: vigra and vigra::linalg
   */
template <class T, class C1, class C2, class C3>
bool
ridgeRegression(MultiArrayView<2, T, C1> const & A,
                MultiArrayView<2, T, C2> const &b, MultiArrayView<2, T, C3> &x, double lambda)
{
    typedef T Real;

    const unsigned int rows = rowCount(A);
    const unsigned int cols = columnCount(A);
    const unsigned int rhsCount = columnCount(b);
    vigra_precondition(rows >= cols,
       "ridgeRegression(): Input matrix A must be rectangular with rowCount >= columnCount.");
    vigra_precondition(rowCount(b) == rows,
       "ridgeRegression(): Shape mismatch between matrices A and b.");
    vigra_precondition(rowCount(x) == cols && columnCount(x) == rhsCount,
       "ridgeRegression(): Result matrix x has wrong shape.");
    vigra_precondition(lambda >= 0.0,
       "ridgeRegression(): lambda >= 0.0 required.");

    unsigned int m = rows;
    unsigned int n = cols;

    Matrix<T> u(m, n), s(n, 1), v(n, n);

    unsigned int rank = singularValueDecomposition(A, u, s, v);
    if(rank < n && lambda == 0.0)
        return false;

    Matrix<T> t = transpose(u)*b;
    for(unsigned int k=0; k<cols; ++k)
        for(unsigned int l=0; l<rhsCount; ++l)
            t(k,l) *= s(k,0) / (sq(s(k,0)) + lambda);
    x = v*t;
    return true;
}

   /** Weighted ridge Regression.

       Given a matrix \a A with <tt>m</tt> rows and <tt>n</tt> columns (with <tt>m \>= n</tt>),
       a vector \a b of length <tt>m</tt>, a weight vector \a weights of length <tt>m</tt>
       with non-negative entries, and a regularization parameter <tt>lambda >= 0.0</tt>
       this function computes the vector \a x of length <tt>n</tt>
       that solves the optimization problem

        \f[ \tilde \textrm{\bf x} = \textrm{argmin}
            \left(\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right)^T
             \textrm{diag}(\textrm{\bf weights})
             \left(\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right) +
             \lambda \textrm{\bf x}^T\textrm{\bf x}
        \f]

       where <tt>diag(weights)</tt> creates a diagonal matrix from \a weights.
       The algorithm calls \ref ridgeRegression() on the equivalent problem

        \f[ \tilde \textrm{\bf x} = \textrm{argmin}
            \left|\left|\textrm{diag}(\textrm{\bf weights})^{1/2}\textrm{\bf A} \textrm{\bf x} -
                  \textrm{diag}(\textrm{\bf weights})^{1/2} \textrm{\bf b}\right|\right|_2^2 +
             \lambda \textrm{\bf x}^T\textrm{\bf x}
        \f]

       where the square root of \a weights is just taken element-wise.  This solution is
       computed by means of \ref singularValueDecomposition().

       When \a b is a matrix with <tt>k</tt> columns, \a x must also have
       <tt>k</tt> columns, which will contain the solutions for the corresponding columns of
       \a b. Note that all matrices must already have the correct shape.

       The function returns <tt>false</tt> if the rank of \a A is less than <tt>n</tt>
       and <tt>lambda == 0.0</tt>.

    <b>\#include</b> \<vigra/regression.hxx\>
        Namespaces: vigra and vigra::linalg
   */
template <class T, class C1, class C2, class C3, class C4>
bool
weightedRidgeRegression(MultiArrayView<2, T, C1> const & A,
             MultiArrayView<2, T, C2> const &b, MultiArrayView<2, T, C3> const &weights,
             MultiArrayView<2, T, C4> &x, double lambda)
{
    typedef T Real;

    const unsigned int rows = rowCount(A);
    const unsigned int cols = columnCount(A);
    const unsigned int rhsCount = columnCount(b);
    vigra_precondition(rows >= cols,
       "weightedRidgeRegression(): Input matrix A must be rectangular with rowCount >= columnCount.");
    vigra_precondition(rowCount(b) == rows,
       "weightedRidgeRegression(): Shape mismatch between matrices A and b.");
    vigra_precondition(rowCount(b) == rowCount(weights) && columnCount(weights) == 1,
       "weightedRidgeRegression(): Weight matrix has wrong shape.");
    vigra_precondition(rowCount(x) == cols && columnCount(x) == rhsCount,
       "weightedRidgeRegression(): Result matrix x has wrong shape.");
    vigra_precondition(lambda >= 0.0,
       "weightedRidgeRegression(): lambda >= 0.0 required.");

    Matrix<T> wa(A.shape()), wb(b.shape());

    for(unsigned int k=0; k<rows; ++k)
    {
        vigra_precondition(weights(k,0) >= 0,
           "weightedRidgeRegression(): Weights must be positive.");
        T w = std::sqrt(weights(k,0));
        for(unsigned int l=0; l<cols; ++l)
            wa(k,l) = w * A(k,l);
        for(unsigned int l=0; l<rhsCount; ++l)
            wb(k,l) = w * b(k,l);
    }

    return ridgeRegression(wa, wb, x, lambda);
}

   /** Ridge Regression with many lambdas.

       This executes \ref ridgeRegression() for a sequence of regularization parameters. This
       is implemented so that the \ref singularValueDecomposition() has to be executed only once.
       \a lambda must be an array conforming to the <tt>std::vector</tt> interface, i.e. must
       support <tt>lambda.size()</tt> and <tt>lambda[k]</tt>. The columns of the matrix \a x
       will contain the solutions for the corresponding lambda, so the  number of columns of
       the matrix \a x must be equal to <tt>lambda.size()</tt>, and \a b must be a columns vector,
       i.e. cannot contain several right hand sides at once.

       The function returns <tt>false</tt> when the matrix \a A is rank deficient. If this
       happens, and one of the lambdas is zero, the corresponding column of \a x will be skipped.

    <b>\#include</b> \<vigra/regression.hxx\>
        Namespaces: vigra and vigra::linalg
   */
template <class T, class C1, class C2, class C3, class Array>
bool
ridgeRegressionSeries(MultiArrayView<2, T, C1> const & A,
          MultiArrayView<2, T, C2> const &b, MultiArrayView<2, T, C3> &x, Array const & lambda)
{
    typedef T Real;

    const unsigned int rows = rowCount(A);
    const unsigned int cols = columnCount(A);
    const unsigned int lambdaCount = lambda.size();
    vigra_precondition(rows >= cols,
       "ridgeRegressionSeries(): Input matrix A must be rectangular with rowCount >= columnCount.");
    vigra_precondition(rowCount(b) == rows && columnCount(b) == 1,
       "ridgeRegressionSeries(): Shape mismatch between matrices A and b.");
    vigra_precondition(rowCount(x) == cols && columnCount(x) == lambdaCount,
       "ridgeRegressionSeries(): Result matrix x has wrong shape.");

    unsigned int m = rows;
    unsigned int n = cols;

    Matrix<T> u(m, n), s(n, 1), v(n, n);

    unsigned int rank = singularValueDecomposition(A, u, s, v);

    Matrix<T> xl = transpose(u)*b;
    Matrix<T> xt(cols,1);
    for(unsigned int i=0; i<lambdaCount; ++i)
    {
        vigra_precondition(lambda[i] >= 0.0,
           "ridgeRegressionSeries(): lambda >= 0.0 required.");
        if(lambda[i] == 0.0 && rank < rows)
            continue;
        for(unsigned int k=0; k<cols; ++k)
            xt(k,0) = xl(k,0) * s(k,0) / (sq(s(k,0)) + lambda[i]);
        columnVector(x, i) = v*xt;
    }
    return (rank == n);
}

/** \brief Pass options to leastAngleRegression().

    <b>\#include</b> \<vigra/regression.hxx\>
        Namespaces: vigra and vigra::linalg
*/
class LeastAngleRegressionOptions
{
  public:
    enum Mode { LARS, LASSO, NNLASSO };

        /** Initialize all options with default values.
        */
    LeastAngleRegressionOptions()
    : max_solution_count(0),
      unconstrained_dimension_count(0),
      mode(LASSO),
      least_squares_solutions(true)
    {}

        /** Maximum number of solutions to be computed.

            If \a n is 0 (the default), the number of solutions is determined by the length
            of the solution array. Otherwise, the minimum of maxSolutionCount() and that
            length is taken.<br>
            Default: 0 (use length of solution array)
        */
    LeastAngleRegressionOptions & maxSolutionCount(unsigned int n)
    {
        max_solution_count = (int)n;
        return *this;
    }

        /** Set the mode of the algorithm.

            Mode must be one of "lars", "lasso", "nnlasso". The function just calls
            the member function of the corresponding name to set the mode.

            Default: "lasso"
        */
    LeastAngleRegressionOptions & setMode(std::string mode)
    {
        for(unsigned int k=0; k<mode.size(); ++k)
            mode[k] = (std::string::value_type)tolower(mode[k]);
        if(mode == "lars")
            this->lars();
        else if(mode == "lasso")
            this->lasso();
        else if(mode == "nnlasso")
            this->nnlasso();
        else
            vigra_fail("LeastAngleRegressionOptions.setMode(): Invalid mode.");
        return *this;
    }


        /** Use the plain LARS algorithm.

            Default: inactive
        */
    LeastAngleRegressionOptions & lars()
    {
        mode = LARS;
        return *this;
    }

        /** Use the LASSO modification of the LARS algorithm.

            This allows features to be removed from the active set under certain conditions.<br>
            Default: active
        */
    LeastAngleRegressionOptions & lasso()
    {
        mode = LASSO;
        return *this;
    }

        /** Use the non-negative LASSO modification of the LARS algorithm.

            This enforces all non-zero entries in the solution to be positive.<br>
            Default: inactive
        */
    LeastAngleRegressionOptions & nnlasso()
    {
        mode = NNLASSO;
        return *this;
    }

        /** Compute least squares solutions.

            Use least angle regression to determine active sets, but
            return least squares solutions for the features in each active set,
            instead of constrained solutions.<br>
            Default: <tt>true</tt>
        */
    LeastAngleRegressionOptions & leastSquaresSolutions(bool select = true)
    {
        least_squares_solutions = select;
        return *this;
    }

    int max_solution_count, unconstrained_dimension_count;
    Mode mode;
    bool least_squares_solutions;
};

namespace detail {

template <class T, class C1, class C2>
struct LarsData
{
    typedef typename MultiArrayShape<2>::type Shape;

    int activeSetSize;
    MultiArrayView<2, T, C1> A;
    MultiArrayView<2, T, C2> b;
    Matrix<T> R, qtb, lars_solution, lars_prediction, next_lsq_solution, next_lsq_prediction, searchVector;
    ArrayVector<MultiArrayIndex> columnPermutation;

    // init data for a new run
    LarsData(MultiArrayView<2, T, C1> const & Ai, MultiArrayView<2, T, C2> const & bi)
    : activeSetSize(1),
      A(Ai), b(bi), R(A), qtb(b),
      lars_solution(A.shape(1), 1), lars_prediction(A.shape(0), 1),
      next_lsq_solution(A.shape(1), 1), next_lsq_prediction(A.shape(0), 1), searchVector(A.shape(0), 1),
      columnPermutation(A.shape(1))
    {
        for(unsigned int k=0; k<columnPermutation.size(); ++k)
            columnPermutation[k] = k;
    }

    // copy data for the recursive call in nnlassolsq
    LarsData(LarsData const & d, int asetSize)
    : activeSetSize(asetSize),
      A(d.R.subarray(Shape(0,0), Shape(d.A.shape(0), activeSetSize))), b(d.qtb), R(A), qtb(b),
      lars_solution(d.lars_solution.subarray(Shape(0,0), Shape(activeSetSize, 1))), lars_prediction(d.lars_prediction),
      next_lsq_solution(d.next_lsq_solution.subarray(Shape(0,0), Shape(activeSetSize, 1))), 
      next_lsq_prediction(d.next_lsq_prediction), searchVector(d.searchVector),
      columnPermutation(A.shape(1))
    {
        for(unsigned int k=0; k<columnPermutation.size(); ++k)
            columnPermutation[k] = k;
    }
};

template <class T, class C1, class C2, class Array1, class Array2, class Array3>
unsigned int 
leastAngleRegressionMainLoop(LarsData<T, C1, C2> & d,
                             Array1 & activeSets, 
                             Array2 * lars_solutions, Array3 * lsq_solutions,
                             LeastAngleRegressionOptions const & options)
{
    using namespace vigra::functor;

    typedef typename MultiArrayShape<2>::type Shape;
    typedef typename Matrix<T>::view_type Subarray;
    typedef ArrayVector<MultiArrayIndex> Permutation;
    typedef typename Permutation::view_type ColumnSet;

    vigra_precondition(d.activeSetSize > 0,
       "leastAngleRegressionMainLoop() must not be called with empty active set.");

    bool enforce_positive = (options.mode == LeastAngleRegressionOptions::NNLASSO);
    bool lasso_modification = (options.mode != LeastAngleRegressionOptions::LARS);

    const MultiArrayIndex rows = rowCount(d.R);
    const MultiArrayIndex cols = columnCount(d.R);
    const MultiArrayIndex maxRank = std::min(rows, cols);

    MultiArrayIndex maxSolutionCount = options.max_solution_count;
    if(maxSolutionCount == 0)
        maxSolutionCount = lasso_modification
                                ? 10*maxRank
                                : maxRank;

    bool needToRemoveColumn = false;
    MultiArrayIndex columnToBeAdded = 0, columnToBeRemoved = 0;
    MultiArrayIndex currentSolutionCount = 0;
    while(currentSolutionCount < maxSolutionCount)
    {
        ColumnSet activeSet = d.columnPermutation.subarray(0, (unsigned int)d.activeSetSize);
        ColumnSet inactiveSet = d.columnPermutation.subarray((unsigned int)d.activeSetSize, (unsigned int)cols);

        // find next dimension to be activated
        Matrix<T> cLARS = transpose(d.A) * (d.b - d.lars_prediction),      // correlation with LARS residual
                  cLSQ  = transpose(d.A) * (d.b - d.next_lsq_prediction);  // correlation with LSQ residual

        // In theory, all vectors in the active set should have the same correlation C, and
        // the correlation of all others should not exceed this. In practice, we may find the
        // maximum correlation in any variable due to tiny numerical inaccuracies. Therefore, we
        // determine C from the entire set of variables.
        MultiArrayIndex cmaxIndex = enforce_positive
                                      ? argMax(cLARS)
                                      : argMax(abs(cLARS));
        T C = abs(cLARS(cmaxIndex, 0));

        Matrix<T> ac(cols - d.activeSetSize, 1);
        for(MultiArrayIndex k = 0; k<cols-d.activeSetSize; ++k)
        {
            T rho = cLSQ(inactiveSet[k], 0), 
              cc  = C - sign(rho)*cLARS(inactiveSet[k], 0);

            if(rho == 0.0)  // make sure that 0/0 cannot happen in the other cases
                ac(k,0) = 1.0; // variable k is linearly dependent on the active set
            else if(rho > 0.0)
                ac(k,0) = cc / (cc + rho); // variable k would enter the active set with positive sign
            else if(enforce_positive)
                ac(k,0) = 1.0; // variable k cannot enter the active set because it would be negative
            else
                ac(k,0) = cc / (cc - rho); // variable k would enter the active set with negative sign
        }

        // in the non-negative case: make sure that a column just removed cannot re-enter right away
        // (in standard LASSO, this is allowed, because the variable may re-enter with opposite sign)
        if(enforce_positive && needToRemoveColumn)
                ac(columnToBeRemoved-d.activeSetSize,0) = 1.0;

        // find candidate
        // Note: R uses Arg1() > epsilon, but this is only possible because it allows several variables to
        //       join the active set simultaneously, so that gamma = 0 cannot occur.
        columnToBeAdded = argMin(ac);

        // if no new column can be added, we do a full step gamma = 1.0 and then stop, unless a column is removed below
        T gamma = (d.activeSetSize == maxRank)
                     ? 1.0
                     : ac(columnToBeAdded, 0);

        // adjust columnToBeAdded: we skipped the active set
        if(columnToBeAdded >= 0)
            columnToBeAdded += d.activeSetSize;

        // check whether we have to remove a column from the active set
        needToRemoveColumn = false;
        if(lasso_modification)
        {
            // find dimensions whose weight changes sign below gamma*searchDirection
            Matrix<T> s(Shape(d.activeSetSize, 1), NumericTraits<T>::max());
            for(MultiArrayIndex k=0; k<d.activeSetSize; ++k)
            {
                if(( enforce_positive && d.next_lsq_solution(k,0) < 0.0) ||
                   (!enforce_positive && sign(d.lars_solution(k,0))*sign(d.next_lsq_solution(k,0)) == -1.0))
                        s(k,0) = d.lars_solution(k,0) / (d.lars_solution(k,0) - d.next_lsq_solution(k,0));
            }

            columnToBeRemoved = argMinIf(s, Arg1() <= Param(gamma));
            if(columnToBeRemoved >= 0)
            {
                needToRemoveColumn = true; // remove takes precedence over add
                gamma = s(columnToBeRemoved, 0);
            }
        }

        // compute the current solutions
        d.lars_prediction  = gamma * d.next_lsq_prediction + (1.0 - gamma) * d.lars_prediction;
        d.lars_solution    = gamma * d.next_lsq_solution   + (1.0 - gamma) * d.lars_solution;
        if(needToRemoveColumn)
            d.lars_solution(columnToBeRemoved, 0) = 0.0;  // turn possible epsilon into an exact zero

        // write the current solution
        ++currentSolutionCount;
        activeSets.push_back(typename Array1::value_type(d.columnPermutation.begin(), d.columnPermutation.begin()+d.activeSetSize));

        if(lsq_solutions != 0)
        {
            if(enforce_positive)
            {
                ArrayVector<Matrix<T> > nnresults;
                ArrayVector<ArrayVector<MultiArrayIndex> > nnactiveSets;
                LarsData<T, C1, C2> nnd(d, d.activeSetSize);

                leastAngleRegressionMainLoop(nnd, nnactiveSets, &nnresults, (Array3*)0,
                                             LeastAngleRegressionOptions().leastSquaresSolutions(false).nnlasso());
                //Matrix<T> nnlsq_solution(d.activeSetSize, 1);
                typename Array2::value_type nnlsq_solution(Shape(d.activeSetSize, 1));
                for(unsigned int k=0; k<nnactiveSets.back().size(); ++k)
                {
                    nnlsq_solution(nnactiveSets.back()[k],0) = nnresults.back()[k];
                }
                //lsq_solutions->push_back(nnlsq_solution);
                lsq_solutions->push_back(typename Array3::value_type());
                lsq_solutions->back() = nnlsq_solution;
            }
            else
            {
                //lsq_solutions->push_back(d.next_lsq_solution.subarray(Shape(0,0), Shape(d.activeSetSize, 1)));
                lsq_solutions->push_back(typename Array3::value_type());
                lsq_solutions->back() = d.next_lsq_solution.subarray(Shape(0,0), Shape(d.activeSetSize, 1));
            }
        }
        if(lars_solutions != 0)
        {
            //lars_solutions->push_back(d.lars_solution.subarray(Shape(0,0), Shape(d.activeSetSize, 1)));
            lars_solutions->push_back(typename Array2::value_type());
            lars_solutions->back() = d.lars_solution.subarray(Shape(0,0), Shape(d.activeSetSize, 1));
        }

        // no further solutions possible
        if(gamma == 1.0)
            break;

        if(needToRemoveColumn)
        {
            --d.activeSetSize;
            if(columnToBeRemoved != d.activeSetSize)
            {
                // remove column 'columnToBeRemoved' and restore triangular form of R
                // note: columnPermutation is automatically swapped here
                detail::upperTriangularSwapColumns(columnToBeRemoved, d.activeSetSize, d.R, d.qtb, d.columnPermutation);

                // swap solution entries
                std::swap(d.lars_solution(columnToBeRemoved, 0), d.lars_solution(d.activeSetSize,0));
                std::swap(d.next_lsq_solution(columnToBeRemoved, 0), d.next_lsq_solution(d.activeSetSize,0));
                columnToBeRemoved = d.activeSetSize; // keep track of removed column
            }
            d.lars_solution(d.activeSetSize,0) = 0.0;
            d.next_lsq_solution(d.activeSetSize,0) = 0.0;
        }
        else
        {
            vigra_invariant(columnToBeAdded >= 0,
                "leastAngleRegression(): internal error (columnToBeAdded < 0)");
            // add column 'columnToBeAdded'
            if(d.activeSetSize != columnToBeAdded)
            {
                std::swap(d.columnPermutation[d.activeSetSize], d.columnPermutation[columnToBeAdded]);
                columnVector(d.R, d.activeSetSize).swapData(columnVector(d.R, columnToBeAdded));
                columnToBeAdded = d.activeSetSize; // keep track of added column
            }

            // zero the corresponding entries of the solutions
            d.next_lsq_solution(d.activeSetSize,0) = 0.0;
            d.lars_solution(d.activeSetSize,0) = 0.0;

            // reduce R (i.e. its newly added column) to triangular form
            detail::qrColumnHouseholderStep(d.activeSetSize, d.R, d.qtb);
            ++d.activeSetSize;
        }

        // compute the LSQ solution of the new active set
        Subarray Ractive = d.R.subarray(Shape(0,0), Shape(d.activeSetSize, d.activeSetSize));
        Subarray qtbactive = d.qtb.subarray(Shape(0,0), Shape(d.activeSetSize, 1));
        Subarray next_lsq_solution_view = d.next_lsq_solution.subarray(Shape(0,0), Shape(d.activeSetSize, 1));
        linearSolveUpperTriangular(Ractive, qtbactive, next_lsq_solution_view);

        // compute the LSQ prediction of the new active set
        d.next_lsq_prediction.init(0.0);
        for(MultiArrayIndex k=0; k<d.activeSetSize; ++k)
            d.next_lsq_prediction += next_lsq_solution_view(k,0)*columnVector(d.A, d.columnPermutation[k]);
    }

    return (unsigned int)currentSolutionCount;
}

template <class T, class C1, class C2, class Array1, class Array2>
unsigned int
leastAngleRegressionImpl(MultiArrayView<2, T, C1> const & A, MultiArrayView<2, T, C2> const &b,
                         Array1 & activeSets, Array2 * lasso_solutions, Array2 * lsq_solutions,
                         LeastAngleRegressionOptions const & options)
{
    using namespace vigra::functor;

    const MultiArrayIndex rows = rowCount(A);

    vigra_precondition(rowCount(b) == rows && columnCount(b) == 1,
       "leastAngleRegression(): Shape mismatch between matrices A and b.");

    bool enforce_positive = (options.mode == LeastAngleRegressionOptions::NNLASSO);

    detail::LarsData<T, C1, C2> d(A, b);

    // find dimension with largest correlation
    Matrix<T> c = transpose(A)*b;
    MultiArrayIndex initialColumn = enforce_positive
                                       ? argMaxIf(c, Arg1() > Param(0.0))
                                       : argMax(abs(c));
    if(initialColumn == -1)
        return 0; // no solution found

    // prepare initial active set and search direction etc.
    std::swap(d.columnPermutation[0], d.columnPermutation[initialColumn]);
    columnVector(d.R, 0).swapData(columnVector(d.R, initialColumn));
    detail::qrColumnHouseholderStep(0, d.R, d.qtb);
    d.next_lsq_solution(0,0) = d.qtb(0,0) / d.R(0,0);
    d.next_lsq_prediction = d.next_lsq_solution(0,0) * columnVector(A, d.columnPermutation[0]);
    d.searchVector = d.next_lsq_solution(0,0) * columnVector(A, d.columnPermutation[0]);

    return leastAngleRegressionMainLoop(d, activeSets, lasso_solutions, lsq_solutions, options);
}

} // namespace detail

   /** Least Angle Regression.

    <b>\#include</b> \<vigra/regression.hxx\>
        Namespaces: vigra and vigra::linalg

   <b> Declarations:</b>

    \code
    namespace vigra {
      namespace linalg {
        // compute either LASSO or least squares solutions
        template <class T, class C1, class C2, class Array1, class Array2>
        unsigned int
        leastAngleRegression(MultiArrayView<2, T, C1> const & A, MultiArrayView<2, T, C2> const &b,
                             Array1 & activeSets, Array2 & solutions,
                             LeastAngleRegressionOptions const & options = LeastAngleRegressionOptions());

        // compute LASSO and least squares solutions
        template <class T, class C1, class C2, class Array1, class Array2>
        unsigned int
        leastAngleRegression(MultiArrayView<2, T, C1> const & A, MultiArrayView<2, T, C2> const &b,
                             Array1 & activeSets, Array2 & lasso_solutions, Array2 & lsq_solutions,
                             LeastAngleRegressionOptions const & options = LeastAngleRegressionOptions());
      }
      using linalg::leastAngleRegression;
    }
    \endcode

       This function implements Least Angle Regression (LARS) as described in

       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
       B.Efron, T.Hastie, I.Johnstone, and R.Tibshirani: <i>"Least Angle Regression"</i>,
       Annals of Statistics 32(2):407-499, 2004.

       It is an efficient algorithm to solve the L1-regularized least squares (LASSO) problem

        \f[ \tilde \textrm{\bf x} = \textrm{argmin}
             \left|\left|\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right|\right|_2^2
           \textrm{ subject to } \left|\left|\textrm{\bf x}\right|\right|_1\le s
        \f]

       and the L1-regularized non-negative least squares (NN-LASSO) problem

        \f[ \tilde \textrm{\bf x} = \textrm{argmin} \left|\left|\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right|\right|_2^2
           \textrm{ subject to } \left|\left|\textrm{\bf x}\right|\right|_1\le s \textrm{ and } \textrm{\bf x}\ge \textrm{\bf 0}
        \f]

       where \a A is a matrix with <tt>m</tt> rows and <tt>n</tt> columns (often with <tt>m \< n</tt>),
       \a b a vector of length <tt>m</tt>, and a regularization parameter s \>= 0.0.
       L1-regularization has the desirable effect that it causes the solution \a x to be sparse, i.e. only
       the most important variables (called the <em>active set</em>) have non-zero values. The
       key insight of the LARS algorithm is the following: When the solution vector is considered
       as a function of the regularization parameter s, then <b>x</b>(s) is a piecewise
       linear function, i.e. a polyline in n-dimensional space. The knots of the polyline
       occur precisely at those values of s where one variable enters or leaves the active set,
       and can be efficiently computed.

       Therefore, leastAngleRegression() returns the entire solution path as a sequence of knot points, starting
       at \f$\textrm{\bf x}(s=0)\f$ (where the only feasible solution is obviously <b>x</b> = 0) and ending at
       \f$\textrm{\bf x}(s=\infty)\f$ (where the solution becomes the ordinary least squares solution). Actually,
       the initial null solution is not explicitly returned, i.e. the sequence starts at the first non-zero
       solution with one variable in the active set. The function leastAngleRegression() returns the number
       of solutions( i.e. knot points) computed.

       The sequences of active sets and corresponding variable weights are returned in \a activeSets and
       \a solutions respectively. That is, <tt>activeSets[i]</tt> is an \ref vigra::ArrayVector "ArrayVector\<int\>"
       containing the indices of the variables that are active at the i-th knot, and <tt>solutions</tt> is a
       \ref vigra::linalg::Matrix "Matrix\<T\>" containing the weights of those variables, in the same order (see
       example below). Variables not contained in <tt>activeSets[i]</tt> are zero at this solution.

       The behavior of the algorithm can be adapted by \ref vigra::linalg::LeastAngleRegressionOptions
       "LeastAngleRegressionOptions":
        <DL>
        <DT><b>options.lasso()</b> (active by default)
                          <DD> Compute the LASSO solution as described above.
        <DT><b>options.nnlasso()</b> (inactive by default)
                          <DD> Compute non-negative LASSO solutions, i.e. use the additional constraint that
                               <b>x</b> \>= 0 in all solutions.
        <DT><b>options.lars()</b> (inactive by default)
                          <DD> Compute a solution path according to the plain LARS rule, i.e. never remove
                               a variable from the active set once it entered.
        <DT><b>options.leastSquaresSolutions(bool)</b> (default: true)
                          <DD> Use the algorithm mode selected above
                               to determine the sequence of active sets, but then compute and return an
                               ordinary (unconstrained) least squares solution for every active set.<br>
                               <b>Note:</b> The second form of leastAngleRegression() ignores this option and
                               does always compute both constrained and unconstrained solutions (returned in
                               \a lasso_solutions and \a lsq_solutions respectively).
        <DT><b>maxSolutionCount(unsigned int n)</b> (default: n = 0, i.e. compute all solutions)
                          <DD> Compute at most <tt>n</tt> solutions.
        </DL>

        <b>Usage:</b>

        \code
        int m = ..., n = ...;
        Matrix<double> A(m, n), b(m, 1);
        ... // fill A and b

        // normalize the input
        Matrix<double> offset(1,n), scaling(1,n);
        prepareColumns(A, A, offset, scaling, DataPreparationGoals(ZeroMean|UnitVariance));
        prepareColumns(b, b, DataPreparationGoals(ZeroMean));

        // arrays to hold the output
        ArrayVector<ArrayVector<int> > activeSets;
        ArrayVector<Matrix<double> > solutions;

        // run leastAngleRegression() in non-negative LASSO mode
        int numSolutions = leastAngleRegression(A, b, activeSets, solutions,
                                    LeastAngleRegressionOptions().nnlasso());

        // print results
        Matrix<double> denseSolution(1, n);
        for (MultiArrayIndex k = 0; k < numSolutions; ++k)
        {
            // transform the sparse solution into a dense vector
            denseSolution.init(0.0); // ensure that inactive variables are zero
            for (unsigned int i = 0; i < activeSets[k].size(); ++i)
            {
                // set the values of the active variables;
                // activeSets[k][i] is the true index of the i-th variable in the active set
                denseSolution(0, activeSets[k][i]) = solutions[k](i,0);
            }

            // invert the input normalization
            denseSolution = denseSolution * pointWise(scaling);

            // output the solution
            std::cout << "solution " << k << ":\n" << denseSolution << std::endl;
        }
        \endcode

        <b>Required Interface:</b>

        <ul>
        <li> <tt>T</tt> must be numeric type (compatible to double)
        <li> <tt>Array1 a1;</tt><br>
             <tt>a1.push_back(ArrayVector\<int\>());</tt>
        <li> <tt>Array2 a2;</tt><br>
             <tt>a2.push_back(Matrix\<T\>());</tt>
        </ul>
   */
doxygen_overloaded_function(template <...> unsigned int leastAngleRegression)

template <class T, class C1, class C2, class Array1, class Array2>
inline unsigned int
leastAngleRegression(MultiArrayView<2, T, C1> const & A, MultiArrayView<2, T, C2> const &b,
                     Array1 & activeSets, Array2 & solutions,
                     LeastAngleRegressionOptions const & options = LeastAngleRegressionOptions())
{
    if(options.least_squares_solutions)
        return detail::leastAngleRegressionImpl(A, b, activeSets, (Array2*)0, &solutions, options);
    else
        return detail::leastAngleRegressionImpl(A, b, activeSets, &solutions, (Array2*)0, options);
}

template <class T, class C1, class C2, class Array1, class Array2>
inline unsigned int
leastAngleRegression(MultiArrayView<2, T, C1> const & A, MultiArrayView<2, T, C2> const &b,
                     Array1 & activeSets, Array2 & lasso_solutions, Array2 & lsq_solutions,
                     LeastAngleRegressionOptions const & options = LeastAngleRegressionOptions())
{
    return detail::leastAngleRegressionImpl(A, b, activeSets, &lasso_solutions, &lsq_solutions, options);
}

   /** Non-negative Least Squares Regression.

       Given a matrix \a A with <tt>m</tt> rows and <tt>n</tt> columns (with <tt>m \>= n</tt>),
       and a column vector \a b of length <tt>m</tt> rows, this function computes
       a column vector \a x of length <tt>n</tt> with <b>non-negative entries</b> that solves the optimization problem

        \f[ \tilde \textrm{\bf x} = \textrm{argmin}
            \left|\left|\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right|\right|_2^2
            \textrm{ subject to } \textrm{\bf x} \ge \textrm{\bf 0}
        \f]

       Both \a b and \a x must be column vectors (i.e. matrices with <tt>1</tt> column).
       Note that all matrices must already have the correct shape. The solution is computed by means
       of \ref leastAngleRegression() with non-negativity constraint.

    <b>\#include</b> \<vigra/regression.hxx\>
        Namespaces: vigra and vigra::linalg
   */
template <class T, class C1, class C2, class C3>
inline void
nonnegativeLeastSquares(MultiArrayView<2, T, C1> const & A,
                        MultiArrayView<2, T, C2> const &b, MultiArrayView<2, T, C3> &x)
{
    vigra_precondition(columnCount(A) == rowCount(x) && rowCount(A) == rowCount(b),
        "nonnegativeLeastSquares(): Matrix shape mismatch.");
    vigra_precondition(columnCount(b) == 1 && columnCount(x) == 1,
        "nonnegativeLeastSquares(): RHS and solution must be vectors (i.e. columnCount == 1).");

    ArrayVector<ArrayVector<MultiArrayIndex> > activeSets;
    ArrayVector<Matrix<T> > results;

    leastAngleRegression(A, b, activeSets, results,
                         LeastAngleRegressionOptions().leastSquaresSolutions(false).nnlasso());
    x.init(NumericTraits<T>::zero());
    if(activeSets.size() > 0)
        for(unsigned int k=0; k<activeSets.back().size(); ++k)
            x(activeSets.back()[k],0) = results.back()[k];
}


//@}

} // namespace linalg

using linalg::leastSquares;
using linalg::weightedLeastSquares;
using linalg::ridgeRegression;
using linalg::weightedRidgeRegression;
using linalg::ridgeRegressionSeries;
using linalg::nonnegativeLeastSquares;
using linalg::leastAngleRegression;
using linalg::LeastAngleRegressionOptions;

} // namespace vigra

#endif // VIGRA_REGRESSION_HXX
