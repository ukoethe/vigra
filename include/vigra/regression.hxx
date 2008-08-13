/************************************************************************/
/*                                                                      */
/*                  Copyright 2007 by Ullrich Koethe                    */
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

/** \addtogroup MatrixAlgebra
 */
//@{
   /** Ordinary Least Squares Regression.

       Given a matrix \a A with <tt>m</tt> rows and <tt>n</tt> columns (with <tt>m \>= n</tt>),
       and a column vector \a b of length <tt>m</tt> rows, this function computes 
       a column vector \a x of length <tt>n</tt> rows such that the residual

        \f[ \left|\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right|^2
        \f]

       is minimized. When \a b is a matrix with <tt>k</tt> columns, \a x must also have 
       <tt>k</tt> columns, which will contain the solutions for the corresponding columns of 
       \a b. Note that all matrices must already have the correct shape.
       
       This function is just another name for \ref linearSolve(), perhaps 
       leading to more readable code when \a A is a rectangular matrix. It returns
       <tt>false</tt> when the rank of \a A is less than <tt>n</tt>.
       See \ref linearSolve() for more documentation.

    <b>\#include</b> \<<a href="regression_8hxx-source.html">vigra/regression_8hxx.hxx</a>\> or<br>
    <b>\#include</b> \<<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>\><br>
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
       with non-negative entries, this function computes a vector \a x of length <tt>n</tt> 
       such that the weighted residual

        \f[  \left(\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right)^T 
             \textrm{diag}(\textrm{\bf weights}) 
             \left(\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right)
        \f]

       is minimized, where <tt>diag(weights)</tt> creates a diagonal matrix from \a weights.
       The algorithm calls \ref leastSquares() on the equivalent problem 

        \f[ \left|\textrm{diag}(\textrm{\bf weights})^{1/2}\textrm{\bf A} \textrm{\bf x} - 
                  \textrm{diag}(\textrm{\bf weights})^{1/2} \textrm{\bf b}\right|^2
        \f]
        
       where the square root of \a weights is just taken element-wise. 
       
       When \a b is a matrix with <tt>k</tt> columns, \a x must also have 
       <tt>k</tt> columns, which will contain the solutions for the corresponding columns of 
       \a b. Note that all matrices must already have the correct shape.

       The function returns
       <tt>false</tt> when the rank of the weighted matrix \a A is less than <tt>n</tt>.

    <b>\#include</b> \<<a href="regression_8hxx-source.html">vigra/regression.hxx</a>\> or<br>
    <b>\#include</b> \<<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>\><br>
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
       this function computes a vector \a x of length <tt>n</tt> such that the residual

        \f[ \left|\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right|^2 + 
            \lambda \textrm{\bf x}^T\textrm{\bf x}
        \f]

       is minimized. This is implemented by means of \ref singularValueDecomposition().
       
       When \a b is a matrix with <tt>k</tt> columns, \a x must also have 
       <tt>k</tt> columns, which will contain the solutions for the corresponding columns of 
       \a b. Note that all matrices must already have the correct shape.
       
       The function returns <tt>false</tt> if the rank of \a A is less than <tt>n</tt>
       and <tt>lambda == 0.0</tt>.

    <b>\#include</b> \<<a href="regression_8hxx-source.html">vigra/regression.hxx</a>\> or<br>
    <b>\#include</b> \<<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>\><br>
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
       this function computes a vector \a x of length <tt>n</tt> such that the weighted residual

        \f[  \left(\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right)^T 
             \textrm{diag}(\textrm{\bf weights}) 
             \left(\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right) +
             \lambda \textrm{\bf x}^T\textrm{\bf x}
        \f]

       is minimized, where <tt>diag(weights)</tt> creates a diagonal matrix from \a weights.
       The algorithm calls \ref ridgeRegression() on the equivalent problem 

        \f[ \left|\textrm{diag}(\textrm{\bf weights})^{1/2}\textrm{\bf A} \textrm{\bf x} - 
                  \textrm{diag}(\textrm{\bf weights})^{1/2} \textrm{\bf b}\right|^2 + 
             \lambda \textrm{\bf x}^T\textrm{\bf x}
        \f]
        
       where the square root of \a weights is just taken element-wise.  This solution is 
       computed by means of \ref singularValueDecomposition().     
       
       When \a b is a matrix with <tt>k</tt> columns, \a x must also have 
       <tt>k</tt> columns, which will contain the solutions for the corresponding columns of 
       \a b. Note that all matrices must already have the correct shape.
 
       The function returns <tt>false</tt> if the rank of \a A is less than <tt>n</tt>
       and <tt>lambda == 0.0</tt>.

    <b>\#include</b> \<<a href="regression_8hxx-source.html">vigra/regression.hxx</a>\> or<br>
    <b>\#include</b> \<<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>\><br>
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

    <b>\#include</b> \<<a href="regression_8hxx-source.html">vigra/regression.hxx</a>\> or<br>
    <b>\#include</b> \<<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>\><br>
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
        if(lambda == 0.0 && rank < rows)
            continue;
        for(unsigned int k=0; k<cols; ++k)
            xt(k,0) = xl(k,0) * s(k,0) / (sq(s(k,0)) + lambda[i]);
        columnVector(x, i) = v*xt;
    }
    return (rank < n);
}

/** \brief Pass options to leastAngleRegression().

*/
class LeastAngleRegressionOptions
{
  public:
        /** Initialize all options with default values.
        */
    LeastAngleRegressionOptions()
    : bic_variance(0.0),
      max_solution_count(0), 
      unconstrained_dimension_count(0),
      bic_delay(5),
      lasso_modification(true), 
      enforce_positive(false),
      least_squares_solutions(true)
    {}

        /** Minimum number of solutions to be computed.
        
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

#if 0 // option currently disabled
        /** Number of unconstrained dimensions in the feature matrix.
        
            The first \a n columns in the feature matrix will be considered as unconstrained,
            i.e. they will always be in the active set, and there is no restriction on the size
            or sign of the corresponding solution coefficients.<br>
            Default: 0 
        */
    LeastAngleRegressionOptions & unconstrainedDimensionCount(unsigned int n)
    {
        unconstrained_dimension_count = (int)n;
        return *this;
    }
#endif

        /** Use the LASSO modification of the LARS algorithm.
        
            This allows features to be removed from the active set under certain conditions.<br>
            Default: <tt>true</tt>
        */
    LeastAngleRegressionOptions & useLassoModification(bool select)
    {
        lasso_modification = select;
        return *this;
    }

        /** Enforce all solution coefficients in the active set to be positive.
        
            This implies the LASSO extension. The constraint does not apply to
            unconstrained dimensions.<br>
            Default: <tt>false</tt>
        */
    LeastAngleRegressionOptions & enforcePositiveSolutions(bool select)
    {
        if(select)
            lasso_modification = true;
        enforce_positive = select;
        return *this;
    }

        /** Compute least squares solutions.
        
            Use least angle regression to determine active sets, but
            return least squares solutions for the features in each active set,
            instead of constrained solutions.<br>
            Default: <tt>true</tt>
        */
    LeastAngleRegressionOptions & leastSquaresSolutions(bool select)
    {
        least_squares_solutions = select;
        return *this;
    }

        /** Stop iterations at minimum of BIC (Bayesian Information Criterium).
        
            The BIC is calculated according to 
            
            \f[
                \textrm{BIC}_k = \frac{|\textrm{\bf b} - \textrm{\bf A} \textrm{\bf x}_k|^2}{variance} + K_k \log N
            \f]
            
            where \f$x_k\f$ is the k-th solution, \f$K_k\f$ is the number of non-zero coefficients
            in the k-th solution, are N is the length of b. Passing a non-positive value for
            \a variance indicates that the BIC is not to be used. \a delay > 0 specifies the number 
            of subsequent steps where the BIC must be larger than the current minimal BIC before
            it is actually accepted as a minimum (in order to avoid spurious minima due to noise).
            
            Defaults: <br>
            <tt>variance = 0.0</tt> meaning that BIC is not used.<br>
            <tt>delay = 5</tt>
        */
    LeastAngleRegressionOptions & stopAtMinimumOfBIC(double variance, unsigned int delay = 5)
    {
        bic_variance = variance;
        bic_delay = (int)delay;
        return *this;
    }

    double bic_variance;
    int max_solution_count, unconstrained_dimension_count, bic_delay;
    bool lasso_modification, enforce_positive, least_squares_solutions;
};

template <class T, class C1, class C2, class Array1, class Array2>
unsigned int
leastAngleRegression(MultiArrayView<2, T, C1> const & A, MultiArrayView<2, T, C2> const &b, 
                     Array1 & solutions, Array2 & activeSets,
                     LeastAngleRegressionOptions const & options = LeastAngleRegressionOptions())
{
    using namespace vigra::functor;
    using namespace vigra::linalg;
    
    typedef typename MultiArrayShape<2>::type Shape;
    typedef typename Matrix<T>::view_type Subarray;

    typedef typename Array2::value_type Permutation;
    typedef typename Permutation::view_type ColumnSet;
    
    if(options.enforce_positive && !options.lasso_modification)
        vigra_precondition(false,
              "leastAngleRegression(): Positive solutions can only be enforced whan LASSO modification is active.");

    const MultiArrayIndex rows = rowCount(A);
    const MultiArrayIndex cols = columnCount(A);
    const MultiArrayIndex maxRank = std::min(rows, cols);
    const unsigned int ucols = (unsigned int)cols;

    vigra_precondition(rowCount(b) == rows && columnCount(b) == 1,
       "leastAngleRegression(): Shape mismatch between matrices A and b.");
       
    MultiArrayIndex maxSolutionCount = options.max_solution_count;
    if(maxSolutionCount == 0)
        maxSolutionCount = options.lasso_modification
                                ? 10*maxRank
                                : maxRank;
    
    Matrix<T> R(A), qtb(b);

    // the first activeSetSize entries will hold the active set indices,
    // the other entries are the inactive set, all permuted in the same way as the
    // columns of the matrix R
    Permutation columnPermutation(ucols);
    for(int k=0; k<cols; ++k)
        columnPermutation[k] = k;
        
    // find dimension with largest correlation
    Matrix<T> c = transpose(A)*b;
    MultiArrayIndex initialColumn;
    if(options.enforce_positive)
        initialColumn = argMaxIf(c, Arg1() > Param(0.0));
    else
        initialColumn = argMax(abs(c));
    if(initialColumn == -1)
        return 0; // no solution found
    
    // prepare initial active set and search direction etc.
    MultiArrayIndex activeSetSize = 1;
    std::swap(columnPermutation[0], columnPermutation[initialColumn]);
    columnVector(R, 0).swapData(columnVector(R, initialColumn));
    detail::qrColumnHouseholderStep(0, R, qtb);

    Matrix<T> lsq_solution(cols, 1), lars_solution(cols,1), lsq_prediction(rows,1), lars_prediction(rows,1); // initially zero
    
    Matrix<T> next_lsq_solution(cols, 1);
    next_lsq_solution(0,0) = qtb(0,0) / R(0,0);
    Matrix<T> searchVector = 
         next_lsq_solution(0,0) * columnVector(A, (MultiArrayIndex)columnPermutation[0]);
    
    double minimal_bic = NumericTraits<double>::max();
    int minimal_bic_solution = -1;
    
    Permutation columnsToBeRemoved;
    MultiArrayIndex currentSolutionCount = 0;
    while(currentSolutionCount < maxSolutionCount)
    {
        ColumnSet activeSet = columnPermutation.subarray(0, (unsigned int)activeSetSize);
        ColumnSet inactiveSet = columnPermutation.subarray((unsigned int)activeSetSize, ucols);
        
        // find next dimension to be activated
        Matrix<T> c(cols - activeSetSize, 1), ac(cols - activeSetSize, 1);
        Matrix<T> lars_residual = b - lars_prediction;

        T C = abs(dot(columnVector(A, activeSet[0]), lars_residual));
        
        for(MultiArrayIndex k = 0; k<cols-activeSetSize; ++k)
        {
            // perform permutation on A explicitly, so that we need not store a permuted copy of A
            c(k, 0) = dot(columnVector(A, inactiveSet[k]), lars_residual);
            T a  = dot(columnVector(A, inactiveSet[k]), searchVector),
              am = (C - c(k, 0)) / (C - a),
              ap = (C + c(k, 0)) / (C + a);

            if(!options.enforce_positive && ap > 0.0 && ap < am)
                ac(k, 0) = ap;
            else
                ac(k, 0) = am;
        }

        MultiArrayIndex columnToBeAdded = argMinIf(ac, Arg1() > Param(0.0));
        if(columnToBeAdded == -1)
            break;  // no further solution possible
        
        T gamma = ac(columnToBeAdded, 0);

        // adjust columnToBeAdded: we skipped the active set
        columnToBeAdded += activeSetSize; 
        
        // check whether we have to remove columns from the active set
        bool needToRemoveColumns = false;
        if(options.lasso_modification)
        {
            // find dimensions whose weight changes sign below gamma*searchDirection
            Matrix<T> d(Shape(activeSetSize, 1), NumericTraits<T>::max());
            for(MultiArrayIndex k=0; k<activeSetSize; ++k)
                if(sign(lsq_solution(k,0))*sign(next_lsq_solution(k,0)) == -1.0)
                    d(k,0) = lsq_solution(k,0) / (lsq_solution(k,0) - next_lsq_solution(k,0));
            int changesSign = argMinIf(d, Arg1() < Param(gamma));
            if(changesSign >= 0)
            {
                needToRemoveColumns = true;
                gamma = d(changesSign, 0);
            }
        }
        
//        gamma = std::min(gamma, NumericTraits<T>::one()); // is this ever necessary ??

        // compute the current solution
        Subarray current_lsq_solution;
        Subarray current_lars_solution;
        if(needToRemoveColumns)
        {
            T tolerance = NumericTraits<T>::epsilon();  // FIXME: adjust tolerance to problem
            lars_solution = gamma * next_lsq_solution + (1.0 - gamma) * lars_solution;
            
            columnsToBeRemoved.clear();
            for(MultiArrayIndex k=0; k<activeSetSize; ++k)
                if((options.enforce_positive && lars_solution(k,0) <= tolerance) ||
                   abs(lars_solution(k,0)) <= tolerance)
                    columnsToBeRemoved.push_back(k);
            
            for(MultiArrayIndex k=0; k<(MultiArrayIndex)columnsToBeRemoved.size(); ++k)
            {
                // remove column 'columnsToBeRemoved[k]' and restore triangular from of R
                detail::upperTriangularSwapColumns(columnsToBeRemoved[k], activeSetSize, R, qtb, columnPermutation);

                // remove entry 'columnsToBeRemoved[k]' from the LARS solution
                std::swap(lars_solution(activeSetSize,0), lars_solution(columnsToBeRemoved[k], 0));
                --activeSetSize;
            }
            
            // remember the active subvector of the solutions
            current_lars_solution = lars_solution.subarray(Shape(0,0), Shape(activeSetSize, 1));
            current_lsq_solution = lsq_solution.subarray(Shape(0,0), Shape(activeSetSize, 1));

            // compute the LSQ solution of the reduced active set
            Subarray Ractive = R.subarray(Shape(0,0), Shape(activeSetSize, activeSetSize));
            Subarray qtbactive = qtb.subarray(Shape(0,0), Shape(activeSetSize, 1));
            linearSolveUpperTriangular(Ractive, qtbactive, current_lsq_solution);

            // compute the predictions of the reduced active set
            lsq_prediction.init(NumericTraits<T>::zero()); 
            lars_prediction.init(NumericTraits<T>::zero()); 
            for(MultiArrayIndex k=0; k<activeSetSize; ++k)
            {
               lsq_prediction += current_lsq_solution(k,0)*columnVector(A, columnPermutation[k]);
               lars_prediction += current_lars_solution(k,0)*columnVector(A, columnPermutation[k]);
            }
        }
        else
        {
            lars_solution = gamma * next_lsq_solution + (1.0 - gamma) * lars_solution;
            current_lars_solution = lars_solution.subarray(Shape(0,0), Shape(activeSetSize, 1));
            current_lsq_solution = next_lsq_solution.subarray(Shape(0,0), Shape(activeSetSize, 1));

            lsq_prediction = lars_prediction + searchVector;
            lars_prediction += gamma*searchVector;
        }
            
        ++currentSolutionCount;
        activeSets.push_back(Permutation(columnPermutation.subarray(0, (unsigned int)activeSetSize)));
        
        double residual;
        if(options.least_squares_solutions)
        {
            solutions.push_back(current_lsq_solution);
            residual = squaredNorm(b - lsq_prediction);
        }
        else
        {
            solutions.push_back(current_lars_solution);
            residual = squaredNorm(b - lars_prediction);
        }
        
        // determine BIC and possibly stop iteration
        if(options.bic_variance > 0.0)
        {
            double bic = residual / options.bic_variance + activeSetSize*std::log((double)rows);
            if(bic < minimal_bic)
            {
                minimal_bic = bic;
                minimal_bic_solution = currentSolutionCount;
            }
            if(currentSolutionCount - minimal_bic_solution >= options.bic_delay)
                return (unsigned int)minimal_bic_solution;
        }

        if(needToRemoveColumns)
        {
            searchVector = lsq_prediction - lars_prediction;
        }
        else
        {
            // add column 'columnToBeAdded'
            std::swap(columnPermutation[activeSetSize], columnPermutation[columnToBeAdded]);
            columnVector(R, activeSetSize).swapData(columnVector(R, columnToBeAdded));

            // zero the corresponding entry of the solutions
            lsq_solution(activeSetSize,0) = 0.0;
            lars_solution(activeSetSize,0) = 0.0;            
            
            // reduce R (i.e. its newly added column) to triangular form
            bool singular = !detail::qrColumnHouseholderStep(activeSetSize, R, qtb);
            if(singular || closeAtTolerance(qtb(activeSetSize,0) / R(activeSetSize, activeSetSize), 0.0)) // FIXME: use tolerance???
                break; // no further solutions possible
            ++activeSetSize;
 
            // compute LSQ solution of new active set
            Subarray Ractive = R.subarray(Shape(0,0), Shape(activeSetSize, activeSetSize));
            Subarray qtbactive = qtb.subarray(Shape(0,0), Shape(activeSetSize, 1));
            Subarray next_lsq_solution_k = next_lsq_solution.subarray(Shape(0,0), Shape(activeSetSize, 1));
            linearSolveUpperTriangular(Ractive, qtbactive, next_lsq_solution_k);
            
            if(activeSetSize == maxRank)
            {
                // if all columns are active, LARS solution and LSQ solution are identical, and no further solution is possible
                ++currentSolutionCount;
                activeSets.push_back(Permutation(columnPermutation.subarray(0, (unsigned int)activeSetSize)));
                solutions.push_back(next_lsq_solution_k);
                break;
            }

            // compute new search direction
            searchVector = -lars_prediction;
            for(MultiArrayIndex k=0; k<activeSetSize; ++k)
                searchVector += next_lsq_solution_k(k,0)*columnVector(A, columnPermutation[k]);
        }
    }
    
    if(options.bic_variance > 0.0 && minimal_bic_solution != (int)currentSolutionCount)
        return (unsigned int)minimal_bic_solution;
    else
        return (unsigned int)currentSolutionCount;
}

} // namespace linalg

using linalg::leastSquares;
using linalg::weightedLeastSquares;
using linalg::ridgeRegression;
using linalg::weightedRidgeRegression;
using linalg::ridgeRegressionSeries;
using linalg::leastAngleRegression;
using linalg::LeastAngleRegressionOptions;

} // namespace vigra

#endif // VIGRA_REGRESSION_HXX
