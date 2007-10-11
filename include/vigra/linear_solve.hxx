/************************************************************************/
/*                                                                      */
/*        Copyright 2004 by Gunnar Kedenburg and Ullrich Koethe         */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
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


#ifndef VIGRA_LINEAR_SOLVE_HXX
#define VIGRA_LINEAR_SOLVE_HXX

#include <ctype.h>
#include <string>
#include "mathutil.hxx"
#include "matrix.hxx"
#include "singular_value_decomposition.hxx"


namespace vigra
{

namespace linalg
{

/** \addtogroup LinearAlgebraFunctions Matrix functions
 */
//@{
    /** invert matrix \a v.
        If the matrix \a v is square, \a res must have the same shape and will contain the
        inverse of \a v. If \a v is rectangular, it must have more rows than columns, and \a res
        must have the transposed shape of \a v. The inverse is then computed in the least-squares 
        sense, i.e. \a res will be the pseudo-inverse (Moore-Penrose inverse).
        The function returns <tt>true</tt> upon success, and <tt>false</tt> if \a v 
        is not invertible (has not full rank). The inverse is computed by means of QR 
        decomposition. 
        
    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear_solve.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2>
bool inverse(const MultiArrayView<2, T, C1> &v, MultiArrayView<2, T, C2> &res)
{
    const unsigned int n = columnCount(v);
    vigra_precondition(n <= rowCount(v),
       "inverse(): input matrix must have at least as many rows as columns.");
    vigra_precondition(n == rowCount(res) && rowCount(v) == columnCount(res),
       "inverse(): shape of output matrix must be the transpose of the input matrix' shape.");

    Matrix<T> q(v.shape()), r(n, n);
    if(!qrDecomposition(v, q, r))
        return false; // a didn't have full rank
    reverseElimination(r, transpose(q), res); 
    return true;
}

    /** create the inverse of matrix \a v.

        The result is returned as a temporary matrix. If the matrix \a v is square, 
        the result will have the same shape and containa the inverse of \a v. 
        If \a v is rectangular, it must have more rows than columns, and the result will
        have the transposed shape of \a v. The inverse is then computed in the least-squares 
        sense, i.e. \a res will be the pseudo-inverse (Moore-Penrose inverse).
        The inverse is computed by means of QR decomposition. If \a v
        is not invertible, <tt>vigra::PreconditionViolation</tt> exception is thrown.
        Usage:

        \code
        vigra::Matrix<double> v(n, n);
        v = ...;

        vigra::Matrix<double> m = inverse(v);
        \endcode

    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear_solve.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C>
TemporaryMatrix<T> inverse(const MultiArrayView<2, T, C> &v)
{
    TemporaryMatrix<T> ret(columnCount(v), rowCount(v));  // transpose shape
    vigra_precondition(inverse(v, ret),
        "inverse(): matrix is not invertible.");
    return ret;
}

template <class T>
TemporaryMatrix<T> inverse(const TemporaryMatrix<T> &v)
{
    if(columnCount(v) == rowCount(v))
    {
        vigra_precondition(inverse(v, const_cast<TemporaryMatrix<T> &>(v)),
            "inverse(): matrix is not invertible.");
        return v;      
    }
    else
    {
        TemporaryMatrix<T> ret(columnCount(v), rowCount(v));  // transpose shape
        vigra_precondition(inverse(v, ret),
            "inverse(): matrix is not invertible.");
        return ret;
    }
}

namespace detail {

template <class T, class C1>
T
determinantByLUDecomposition(MultiArrayView<2, T, C1> const & a)
{
    unsigned int m = rowCount(a), n = columnCount(a);
    vigra_precondition(n == m,
       "determinant(): square matrix required.");
       
    Matrix<T> LU(a);
    double permutationSign = 1.0;
    T det = 1.0;

    for (unsigned int j = 0; j < n; ++j) 
    {
        // Apply previous transformations.
        for (unsigned int i = 0; i < m; ++i) 
        {
            // Most of the time is spent in the following dot product.
            unsigned int kmax = i < j ? i : j;
            T s = 0.0;
            for (unsigned int k = 0; k < kmax; ++k) 
            {
                s += LU(i,k)*LU(k,j);
            }

            LU(i,j) = LU(i,j) -= s;
        }

        // Find pivot and exchange if necessary.
        unsigned int p = j;
        for (unsigned int i = j+1; i < m; ++i) 
            if (abs(LU(i,j)) > abs(LU(p,j))) 
                p = i;

        if (p != j) 
        {
            Matrix<T> t = rowVector(LU, p);
            rowVector(LU, p) = rowVector(LU, j);
            rowVector(LU, j) = t;
            permutationSign = -permutationSign;
        }
        
        det *= LU(j,j);

        // Compute multipliers.

        if ((j < m) && (LU(j,j) != 0.0)) 
        {
            for (unsigned int i = j+1; i < m; ++i) 
            {
               LU(i,j) /= LU(j,j);
            }
        }
    }
    return det * permutationSign;
}

} // namespace detail



    /** Compute the determinant of a square matrix.

        \a method must be one of the following:
        <DL>
        <DT>"Cholesky"<DD> Compute the solution by means of Cholesky decomposition. This
                           method is faster than "LU", but requres the matrix \a a 
                           to be symmetric positive definite. If this is 
                           not the case, a <tt>ContractViolation</tt> exception is thrown.
                           
        <DT>"LU"<DD> (default) Compute the solution by means of LU decomposition.
        </DL>

    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear_solve.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1>
T
determinant(MultiArrayView<2, T, C1> const & a, std::string method = "LU")
{
    unsigned int n = columnCount(a);
    vigra_precondition(rowCount(a) == n,
               "determinant(): Square matrix required.");    

    if(n == 1)
        return a(0,0);
    if(n == 2)
        return a(0,0)*a(1,1) - a(0,1)*a(1,0);
    if(method == "LU")
    {
        return detail::determinantByLUDecomposition(a);
    }
    else if(method == "Cholesky")
    {
        Matrix<T> L(a.shape());
        vigra_precondition(choleskyDecomposition(a, L),
           "determinant(): Cholesky method requires symmetric positive definite matrix.");
        T det = L(0,0);
        for(unsigned int k=1; k<n; ++k)
            det *= L(k,k);
        return sq(det);
    }
    else
    {
        vigra_precondition(false, "determinant(): Unknown solution method.");
    }
    return T();
}


    /** Cholesky decomposition.

        \a A must be a symmetric positive definite matrix, and \a L will be a lower
        triangular matrix, such that (up to round-off errors):

        \code
        A == L * transpose(L);
        \endcode

        This implementation cannot be applied in-place, i.e. <tt>&L == &A</tt> is an error.

    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear_solve.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2>
bool choleskyDecomposition(MultiArrayView<2, T, C1> const & A,
                           MultiArrayView<2, T, C2> &L)
{
    typedef T Real;
    
	unsigned int n = columnCount(A);
	
    vigra_precondition(rowCount(A) == n,
                       "choleskyDecomposition(): Input matrix must be square.");
    vigra_precondition(n == columnCount(L) && n == rowCount(L),
                       "choleskyDecomposition(): Output matrix must have same shape as input matrix.");

     for (unsigned int j = 0; j < n; ++j) 
	 {
        Real d(0.0);
        for (unsigned int k = 0; k < j; ++k) 
		{
            Real s(0.0);
            for (unsigned int i = 0; i < k; ++i) 
			{
               s += L(k, i)*L(j, i);
            }
            L(j, k) = s = (A(j, k) - s)/L(k, k);
            d = d + s*s;
            if(A(k, j) != A(j, k))
                return false;  // A is not symmetric 
         }
         d = A(j, j) - d;
         if(d <= 0.0)
            return false;  // A is not positive definite
         L(j, j) = std::sqrt(d);
         for (unsigned int k = j+1; k < n; ++k) 
		 {
            L(j, k) = 0.0;
         }
	}
    return true;
}

    /** QR decomposition.

        \a a contains the original matrix, results are returned in \a q and \a r, where
        \a q is a orthogonal matrix, and \a r is an upper triangular matrix, such that 
        (up to round-off errors):

        \code
        a == q * r;
        \endcode

        If \a a dosn't have full rank, the function returns <tt>false</tt>. 
        The decomposition is computed by householder transformations. It can be applied in-place,
        i.e. <tt>&a == &q</tt> or <tt>&a == &r</tt> are allowed.

    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear_solve.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
bool qrDecomposition(MultiArrayView<2, T, C1> const & a,
                     MultiArrayView<2, T, C2> &q, MultiArrayView<2, T, C3> &r)
{
    typedef T Real;
    
    const unsigned int m = rowCount(a);
    const unsigned int n = columnCount(a);
    vigra_precondition(m >= n &&
                       n == columnCount(r) && n == rowCount(r) &&
                       n == columnCount(q) && m == rowCount(q),
                       "qrDecomposition(): Matrix shape mismatch.");

    Matrix<T> qr = a;
    
    bool fullRank = true;

    // Main loop.
    for (unsigned int k = 0; k < n; ++k) 
    {
        // Compute 2-norm of k-th column without under/overflow.
        Real nrm = 0.0;
        for (unsigned int i = k; i < m; ++i) 
        {
            nrm = hypot(nrm, qr(i, k));
        }

        if (nrm != 0.0) 
        {
            // Form k-th Householder vector.
            if (qr(k, k) < 0.0) 
            {
                nrm = -nrm;
            }
            for (unsigned int i = k; i < m; ++i) 
            {
                qr(i, k) /= nrm;
            }
            qr(k, k) += 1.0;

            // Apply transformation to remaining columns.
            for (unsigned int j = k+1; j < n; ++j) 
            {
                Real s = 0.0; 
                for (unsigned int i = k; i < m; ++i) 
                {
                    s += qr(i,k)*qr(i,j);
                }
                s = -s/qr(k,k);
                for (unsigned int i = k; i < m; ++i) 
                {
                    qr(i,j) += s*qr(i,k);
                }
            }
        }
        r(k,k) = -nrm;
        if(nrm == 0.0)
            fullRank = false;
    }
    for (unsigned int i = 0; i < n; ++i) 
    {
        for (unsigned int j = i+1; j < n; ++j) 
        {
            r(i,j) = qr(i,j);
            r(j,i) = 0.0;
        }
    }
    for (int k = n-1; k >= 0; --k) 
    {
        for (unsigned int i = 0; i < m; ++i) 
        {
            q(i,k) = 0.0;
        }
        q(k,k) = 1.0;
        for (unsigned int j = k; j < n; ++j) 
        {
            if (qr(k,k) != 0.0) 
            {
                Real s = 0.0;
                for (unsigned int i = k; i < m; ++i) 
                {
                    s += qr(i,k)*q(i,j);
                }
                s = -s/qr(k,k);
                for (unsigned int i = k; i < m; ++i) 
                {
                    q(i,j) += s*qr(i,k);
                }
            }
        }
    }
    return fullRank;
}

    /** Solve a linear system with right-triangular coefficient matrix.

        The square matrix \a r must be a right-triangular coefficient matrix as can,
        for example, be obtained by means of QR decomposition. If \a r doesn't have full rank
        the function fails and returns <tt>false</tt>, otherwise it returns <tt>true</tt>.
        
        The column vectors in \a b are the right-hand sides of the equation (so, several equations
        with the same coefficients can be solved in one go). The result is returned
        int \a x, whose columns contain the solutions for the correspoinding
        columns of \a b. The number of columns of \a a must equal the number of rows of
        both \a b and \a x, and the number of columns of \a b and \a x must be
        equal. This implementation can be applied in-place, i.e. <tt>&b == &x</tt> is allowed.

    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear_solve.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
bool reverseElimination(const MultiArrayView<2, T, C1> &r, const MultiArrayView<2, T, C2> &b,
                        MultiArrayView<2, T, C3> & x)
{
    unsigned int m = columnCount(r);
    unsigned int n = columnCount(b);
    vigra_precondition(m == rowCount(r),
        "reverseElimination(): square coefficient matrix required.");
    vigra_precondition(m == rowCount(b) && m == rowCount(x) && n == columnCount(x),
        "reverseElimination(): matrix shape mismatch.");

    for(unsigned int k = 0; k < n; ++k)
    {
        for(int i=m-1; i>=0; --i)
        {
            if(r(i,i) == NumericTraits<T>::zero())
                return false;  // r doesn' have full rank
            T sum = b(i, k);
            for(unsigned int j=i+1; j<m; ++j)
                 sum -= r(i, j) * x(j, k);
            x(i, k) = sum / r(i, i);
        }
    }
    return true;
}

    /** Solve a linear system with left-triangular coefficient matrix.

        The square matrix \a l must be a left-triangular coefficient matrix. If \a l 
        doesn't have full rank the function fails and returns <tt>false</tt>, 
        otherwise it returns <tt>true</tt>.
        
        The column vectors in \a b are the right-hand sides of the equation (so, several equations
        with the same coefficients can be solved in one go). The result is returned
        int \a x, whose columns contain the solutions for the correspoinding
        columns of \a b. The number of columns of \a a must equal the number of rows of
        both \a b and \a x, and the number of columns of \a b and \a x must be
        equal. This implementation can be applied in-place, i.e. <tt>&b == &x</tt> is allowed.

    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear_solve.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
bool leftReverseElimination(const MultiArrayView<2, T, C1> &l, const MultiArrayView<2, T, C2> &b,
                            MultiArrayView<2, T, C3> & x)
{
    unsigned int m = columnCount(l);
    unsigned int n = columnCount(b);
    vigra_precondition(m == rowCount(l),
        "leftReverseElimination(): square coefficient matrix required.");
    vigra_precondition(m == rowCount(b) && m == rowCount(x) && n == columnCount(x),
        "leftReverseElimination(): matrix shape mismatch.");

    for(unsigned int k = 0; k < n; ++k)
    {
        for(unsigned int i=0; i<m; ++i)
        {
            if(l(i,i) == NumericTraits<T>::zero())
                return false;  // l doesn' have full rank
            T sum = b(i, k);
            for(unsigned int j=0; j<i; ++j)
                 sum -= l(i, j) * x(j, k);
            x(i, k) = sum / l(i, i);
        }
    }
    return true;
}

    /** Solve a linear system.

        The \a a is the coefficient matrix, and the column vectors
        in \a b are the right-hand sides of the equation (so, several equations
        with the same coefficients can be solved in one go). The result is returned 
        in \a res, whose columns contain the solutions for the corresponding
        columns of \a b. The number of columns of \a a must equal the number of rows of
        both \a b and \a res, and the number of columns of \a b and \a res must be
        equal. 
        
        \a method must be one of the following:
        <DL>
        <DT>"Cholesky"<DD> Compute the solution by means of Cholesky decomposition. The 
                           coefficient matrix \a a must by symmetric positive definite. If
                           this is not the case, the function returns <tt>false</tt>.
                           
        <DT>"QR"<DD> (default) Compute the solution by means of QR decomposition.  The 
                           coefficient matrix \a a can be square or rectangular. In the latter case,
                           it must have more rows than columns, and the solution will be computed in the 
                           least squares sense. If \a a doesn't have full rank, the function 
                           returns <tt>false</tt>.

        <DT>"SVD"<DD> Compute the solution by means of singular value decomposition.  The 
                           coefficient matrix \a a can be square or rectangular. In the latter case,
                           it must have more rows than columns, and the solution will be computed in the 
                           least squares sense. If \a a doesn't have full rank, the function 
                           returns <tt>false</tt>.
        </DL>
        
        This function can be applied in-place, i.e. <tt>&b == &res</tt> or <tt>&a == &res</tt> are allowed
        (provided they have the required shapes).

    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear_solve.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
bool linearSolve(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b,
                 MultiArrayView<2, T, C3> & res, std::string method = "QR")
{
    vigra_precondition(columnCount(a) <= rowCount(a),
        "linearSolve(): Coefficient matrix a must have at least as many rows as columns.");
    vigra_precondition(columnCount(a) == rowCount(res) && 
                       rowCount(a) == rowCount(b) && columnCount(b) == columnCount(res),
        "linearSolve(): matrix shape mismatch.");

    for(unsigned int k=0; k<method.size(); ++k)
        method[k] = tolower(method[k]);
    
    if(method == "cholesky")
    {
        vigra_precondition(columnCount(a) == rowCount(a),
            "linearSolve(): Cholesky method requires square coefficient matrix.");
        Matrix<T> L(a.shape());
        if(!choleskyDecomposition(a, L))
            return false; // false if a wasn't symmetric positive definite
        leftReverseElimination(L, b, res);
        reverseElimination(transpose(L), res, res);
    }
    else if(method == "qr")
    {
        Matrix<T> q(a.shape()), r(columnCount(a), columnCount(a));
        if(!qrDecomposition(a, q, r))
            return false; // a didn't have full rank
        reverseElimination(r, transpose(q) * b, res);
    }
    else if(method == "svd")
    {
        unsigned int n = rowCount(b);
        unsigned int m = columnCount(b);
	    Matrix<T> u(a.shape()), s(n, 1), v(n, n);

        unsigned int rank = singularValueDecomposition(a, u, s, v);
        if(rank < n)
            return false; // a didn't have full rank

        Matrix<T> t = transpose(u)*b;
        for(unsigned int k=0; k<n; ++k)
            for(unsigned int l=0; l<m; ++l)
                t(k,l) /= s(k,0);
        res = v*t;
    }
    else
    {
        vigra_precondition(false, "linearSolve(): Unknown solution method.");
    }
    return true;
}

//@}

} // namespace linalg

using linalg::inverse;
using linalg::determinant;
using linalg::linearSolve;
using linalg::choleskyDecomposition;
using linalg::qrDecomposition;
using linalg::reverseElimination;
using linalg::leftReverseElimination;

} // namespace vigra


#endif // VIGRA_LINEAR_SOLVE_HXX
