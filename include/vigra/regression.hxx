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

template <class T, class C1, class C2, class Array1, class Array2>
unsigned int
leastAngleRegression(MultiArrayView<2, T, C1> const & A, MultiArrayView<2, T, C2> const &b, 
                     Array1 & solutions, Array2 & activeSets
                     /* LARSOptions const & options */)
{
    using namespace vigra::functor;
    using namespace vigra::linalg;
    
    typedef typename MultiArrayView<2, T, C1>::difference_type Shape;
    typedef typename Matrix<T>::view_type Subarray;
    
    const unsigned int rows = rowCount(A);
    const unsigned int cols = columnCount(A);
    const unsigned int maxSolutionCount = std::min(solutions.size(), std::min(rows, cols));
    vigra_precondition(rowCount(b) == rows && columnCount(b) == 1,
       "leastAngleRegression(): Shape mismatch between matrices A and b.");
    vigra_precondition(maxSolutionCount <= activeSets.size(),
       "leastAngleRegression(): Active sets array too small.");
       
    Matrix<T> X(A);
    Matrix<T> mu(b.shape());

    unsigned int k = 0;
    ArrayVector<int> activeSet(cols);
    for(k=0; k<cols; ++k)
        activeSet[k] = k;
        
    T C = 0.0;
    int best = -1,
        activeSetSize = 0;
    for(k=0; k < maxSolutionCount; ++k)
    {
        Subarray Xinactive = X.subarray(Shape(0, activeSetSize), Shape(rows, cols));
        Matrix<T> c = transpose(Xinactive)*(b - mu);
        if(activeSetSize == 0)
        {
            // find initial active column
            if(false) // FIXME: positive LASSO restriction
                best = argMaxIf(c, Arg1() > Param(0.0));
            else
                best = argMax(abs(c));
            if(best == -1)
                break; // no solution found
            C = abs(c(best, 0));
        }
        else
        {
            Subarray Xactive = X.subarray(Shape(0,0), Shape(rows, activeSetSize));
            Matrix<T> u = Xactive * solutions[k-1] - mu;
            Matrix<T> a = transpose(Xinactive)*u;
            Matrix<T> ac = (C - c) / pointWise(C - a);
            if(true) // FIXME: not positive LASSO restriction
                ac = joinColumns(ac, (C + c) / pointWise(C + a));
            best = argMinIf(ac, Arg1() > Param(0.0));
            if(best == -1)
                break; // no solution found
            T gamma = ac(best, 0);
            mu += gamma*u;
            
            // adjust best: we possibly joined two ac vectors
            best %= (cols - activeSetSize);
            C = abs(c(best, 0));

            // adjust best: we skipped the active set
            best += activeSetSize; 
        }
        
        columnVector(X, k).swapData(columnVector(X, best));
        std::swap(activeSet[k], activeSet[best]);
        ++activeSetSize;

        Subarray Xactive = X.subarray(Shape(0,0), Shape(rows, activeSetSize));
        solutions[k].reshape(Shape(activeSetSize, 1));
        leastSquares(Xactive, b, solutions[k]);

        activeSets[k].resize(activeSetSize);
        std::copy(activeSet.begin(), activeSet.begin()+activeSetSize, activeSets[k].begin());
    }
    
    return k;
}

//@}

} // namespace linalg

using linalg::leastSquares;
using linalg::weightedLeastSquares;
using linalg::ridgeRegression;
using linalg::weightedRidgeRegression;
using linalg::ridgeRegressionSeries;
using linalg::leastAngleRegression;

} // namespace vigra

#endif // VIGRA_REGRESSION_HXX
