/************************************************************************/
/*                                                                      */
/*                  Copyright 2007 by Ullrich Koethe                    */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

/** \addtogroup LinearAlgebraFunctions Matrix functions
 */
//@{
   /** Ordinary Least Squares Regression.

       Given a matrix \a A with <tt>m</tt> rows and tt>n</tt> columns (with <tt>m >= n</tt>),
       and a column vector \a b of length <tt>m</tt> rows, this function computes 
       a column vector \a x of length <tt>n</tt> rows such that the residual

        \f[ \left|\textrm{\bf A} \textrm{\bf x} - \textrm{\bf b}\right|^2
        \f]

       is minimized. When \a b is a matrix with <tt>k</tt> columns, \a x must also have 
       <tt>k</tt> columns, which will contain the solutions for the corresponding columns of 
       \a b. Note that all matrices must already have the correct shape.
       
       This function is just another name for \ref linearSolve(), perhaps 
       leading to more readable code when \A is a rectangular matrix. It returns
       <tt>false</tt> when the rank of \a A is less than <tt>n</tt>.
       See \ref linearSolve() for more documentation.

    <b>\#include</b> "<a href="regression_8hxx-source.html">vigra/regression_8hxx.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
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

       Given a matrix \a A with <tt>m</tt> rows and tt>n</tt> columns (with <tt>m >= n</tt>),
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

    <b>\#include</b> "<a href="regression_8hxx-source.html">vigra/regression_8hxx.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
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

       Given a matrix \a A with <tt>m</tt> rows and tt>n</tt> columns (with <tt>m >= n</tt>),
       a vector \a b of length <tt>m</tt>, and a regularization parameter <tt>lambda >= 0.0</tt>, 
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

    <b>\#include</b> "<a href="regression_8hxx-source.html">vigra/regression_8hxx.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
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

       Given a matrix \a A with <tt>m</tt> rows and tt>n</tt> columns (with <tt>m >= n</tt>),
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

    <b>\#include</b> "<a href="regression_8hxx-source.html">vigra/regression_8hxx.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
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

    <b>\#include</b> "<a href="regression_8hxx-source.html">vigra/regression_8hxx.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
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
        max_solution_count = n;
        return *this;
    }

#if 0 // option currently disabled
        /** Number of unconstrained dimensions in the feature matrix.
        
            The first \a n columns in the feature matrix will be considered as unconstrained,
            i.e. they will always be in the active set, and there is no restriction on the size
            or sign of the coirresponding solution coefficients.<br>
            Default: 0 
        */
    LeastAngleRegressionOptions & unconstrainedDimensionCount(unsigned int n)
    {
        unconstrained_dimension_count = n;
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
            
            where \f$x_k\f$ is the k-th solution, \f$K_k\f$ is the number of non-zero coeffocients
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
        vigra_precondition(delay > 0,
           "LeastAngleRegressionOptions::stopAtMinimumOfBIC(): delay must be > 0.");
        
        bic_variance = variance;
        bic_delay = delay;
        return *this;
    }

    double bic_variance;
    unsigned int max_solution_count, unconstrained_dimension_count, bic_delay;
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
    
    typedef typename MultiArrayView<2, T, C1>::difference_type Shape;
    typedef typename Matrix<T>::view_type Subarray;
    
    if(options.enforce_positive && !options.lasso_modification)
        vigra_precondition(false,
              "leastAngleRegression(): Positive solutions can only be enforced whan LASSO modification is active.");

    const unsigned int rows = rowCount(A);
    const unsigned int cols = columnCount(A);

    vigra_precondition(rowCount(b) == rows && columnCount(b) == 1,
       "leastAngleRegression(): Shape mismatch between matrices A and b.");
       
    unsigned int maxSolutionCount;
    if(options.max_solution_count == 0)
        maxSolutionCount = solutions.size();
    else
        maxSolutionCount = std::min(solutions.size(), options.max_solution_count);
    vigra_precondition(maxSolutionCount <= activeSets.size(),
       "leastAngleRegression(): Active sets array too small.");
    
    Matrix<T> R(A), qtb(b);

    // set first activeSetSize entries will hold the active set indices,
    // the other entries are the inactive set, all order in the same way as the
    // columns of the matrix R
    unsigned int k;
    ArrayVector<int> activeSet(cols);
    for(k=0; k<cols; ++k)
        activeSet[k] = k;
        
    // find dimension with largest correlation
    Matrix<T> c = transpose(A)*b;
    int initialDimension;
    if(options.enforce_positive)
        initialDimension = argMaxIf(c, Arg1() > Param(0.0));
    else
        initialDimension = argMax(abs(c));
    if(initialDimension == -1)
        return 0; // no solution found
    T C = abs(c(initialDimension, 0));
    
    // prepare initial active set and search direction etc.
    int activeSetSize = 1;
    std::swap(activeSet[0], activeSet[initialDimension]);
    columnVector(R, 0).swapData(columnVector(R, initialDimension));
    detail::qrLinearSolveOneStep(0, R, qtb);

    Matrix<T> lsq_solution(cols, 1), lars_solution(cols,1), mu(b.shape()); // initially zero
    Matrix<T> next_lsq_solution(cols, 1);
    next_lsq_solution(0,0) = qtb(0,0) / R(0,0);
    Matrix<T> searchVector = next_lsq_solution(0,0) * columnVector(A, activeSet[0]);
    
    double minimal_bic = NumericTraits<double>::max();
    int minimal_bic_solution = -1;
    
    for(k=0; k < maxSolutionCount; ++k)
    {
        if(activeSetSize == std::min(rows, cols))
        {
            // cannot have more solutions than the size of the matrix A
            // last solution is then always the LSQ solution
            solutions[k] = next_lsq_solution.subarray(Shape(0,0), Shape(activeSetSize, 1));
            activeSets[k].resize(activeSetSize);
            std::copy(activeSet.begin(), activeSet.begin()+activeSetSize, activeSets[k].begin());
            ++k;
            break;
        }
        
        // find next dimension to be activated
        Matrix<T> c(cols - activeSetSize, 1), a(cols - activeSetSize, 1);
        Matrix<T> bmu = b - mu;
        for(unsigned int l = activeSetSize; l<cols; ++l)
        {
            // perform permutation on A explicitly, so that we need not store a permuted copy of A
            c(l-activeSetSize, 0) = dot(columnVector(A, activeSet[l]), bmu);
            a(l-activeSetSize, 0) = dot(columnVector(A, activeSet[l]), searchVector);
        }
        Matrix<T> ac = (C - c) / pointWise(C - a);
        if(!options.enforce_positive)
            ac = joinVertically(ac, (C + c) / pointWise(C + a));
        
        int limitingDimension = argMinIf(ac, Arg1() > Param(0.0));
        if(limitingDimension == -1)
            break;  // no further solution found
        
        T gamma = ac(limitingDimension, 0);
        
        // adjust limitingDimension: we possibly joined two ac vectors
        limitingDimension %= (cols - activeSetSize);
        C = abs(c(limitingDimension, 0));

        // adjust limitingDimension: we skipped the active set
        limitingDimension += activeSetSize; 
        
        // check whether we have to remove a dimension from the active set
        Subarray lsq_solution_k = lsq_solution.subarray(Shape(0,0), Shape(activeSetSize, 1));
        Subarray next_lsq_solution_k = next_lsq_solution.subarray(Shape(0,0), Shape(activeSetSize, 1));
        if(options.lasso_modification)
        {
            // find dimensions whose weight changes sign below gamma*searchDirection
            Matrix<T> d(Shape(activeSetSize, 1), NumericTraits<T>::max());
            for(int l=0; l<activeSetSize; ++l)
                if(sign(lsq_solution_k(l,0))*sign(next_lsq_solution_k(l,0)) == -1.0)
                    d(l,0) = lsq_solution_k(l,0) / (lsq_solution_k(l,0) - next_lsq_solution_k(l,0));
            int changesSign = argMinIf(d, Arg1() < Param(gamma));
            if(changesSign >= 0)
            {
                limitingDimension = changesSign;
                gamma = d(changesSign, 0);
            }
        }

        // compute and write the current solution
        lsq_solution_k = next_lsq_solution_k;
        Subarray lars_solution_k = lars_solution.subarray(Shape(0,0), Shape(activeSetSize, 1));
        lars_solution_k = gamma * lsq_solution_k + (1.0 - gamma) * lars_solution_k;
        
        double residual;
        if(options.least_squares_solutions)
        {
            solutions[k] = lsq_solution_k;
            residual = squaredNorm(bmu - searchVector);
        }
        else
        {
            solutions[k] = lars_solution_k;
            residual = squaredNorm(bmu - gamma*searchVector);
        }
        
        activeSets[k].resize(activeSetSize);
        std::copy(activeSet.begin(), activeSet.begin()+activeSetSize, activeSets[k].begin());
        
        if(options.bic_variance > 0.0)
        {
            double bic = residual / options.bic_variance + activeSetSize*std::log(rows);
            if(bic < minimal_bic)
            {
                minimal_bic = bic;
                minimal_bic_solution = k;
            }
            if(k - minimal_bic_solution >= options.bic_delay)
                break;
        }

        // update the active set and its QR factorization
        std::swap(activeSet[activeSetSize], activeSet[limitingDimension]);
        if(limitingDimension < activeSetSize)
        {
            std::swap(lsq_solution(activeSetSize,0), lsq_solution(limitingDimension, 0));
            std::swap(lars_solution(activeSetSize,0), lars_solution(limitingDimension, 0));
            detail::qrLinearSolveSwap(limitingDimension, activeSetSize, R, qtb);
            --activeSetSize;
        }
        else
        {
            lsq_solution(activeSetSize,0) = 0.0;
            lars_solution(activeSetSize,0) = 0.0;
            columnVector(R, activeSetSize).swapData(columnVector(R, limitingDimension));
            bool singular = !detail::qrLinearSolveOneStep(activeSetSize, R, qtb);
            if(singular || closeAtTolerance(qtb(activeSetSize,0) / R(activeSetSize, activeSetSize), 0.0))
            {
                ++k;
                break; // no further solutions possible
            }
            ++activeSetSize;
        }
        
        // compute LSQ solution of new active set
        Subarray Ractive = R.subarray(Shape(0,0), Shape(activeSetSize, activeSetSize));
        Subarray qtbactive = qtb.subarray(Shape(0,0), Shape(activeSetSize, 1));
        Subarray next_lsq_solution_k1 = next_lsq_solution.subarray(Shape(0,0), Shape(activeSetSize, 1));
        reverseElimination(Ractive, qtbactive, next_lsq_solution_k1);
        
        // compute new search direction
        mu += gamma*searchVector;
        
        searchVector = -mu;
        for(unsigned int l=0; l<activeSetSize; ++l)
            searchVector += next_lsq_solution_k1(l,0)*columnVector(A, activeSet[l]);
    }
    
    if(options.bic_variance > 0.0 && minimal_bic_solution + 1 != k)
        return minimal_bic_solution + 1;
    else
        return k;
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
