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

#include "matrix.hxx"


namespace vigra
{

namespace linalg
{

/** \addtogroup LinearAlgebraFunctions Matrix functions
 */
//@{
    /** invert square matrix \a v.
        The result is written into \a r which must have the same shape.
        The inverse is calculated by means of QR decomposition. If \a v
        is not invertible, <tt>vigra::PreconditionViolation</tt> exception is thrown.

    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear_solve.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2>
void inverse(const MultiArrayView<2, T, C1> &v, MultiArrayView<2, T, C2> &r)
{
    const unsigned int n = rowCount(r);
    vigra_precondition(n == columnCount(v) && n == rowCount(v) && n == columnCount(r),
       "inverse(): matrices must be square.");
    vigra_precondition(linearSolve(v, identityMatrix<T>(n), r),
        "inverse(): matrix is not invertible.");
}

    /** create the inverse of square matrix \a v.
        The result is returned as a temporary matrix.
        The inverse is calculated by means of QR decomposition. If \a v
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
    const unsigned int n = rowCount(v);
    vigra_precondition(n == columnCount(v),
       "inverse(): matrix must be square.");
    TemporaryMatrix<T> ret = identityMatrix<T>(n);
    vigra_precondition(linearSolve(v, ret, ret),
        "inverse(): matrix is not invertible.");
    return ret;
}

template <class T>
TemporaryMatrix<T> inverse(const TemporaryMatrix<T> &v)
{
    const unsigned int n = v.rowCount();
    vigra_precondition(n == v.columnCount(),
       "inverse(): matrix must be square.");
    vigra_precondition(linearSolve(v, identityMatrix<T>(n), const_cast<TemporaryMatrix<T> &>(v)),
        "inverse(): matrix is not invertible.");
    return v;
}

    /** QR decomposition.

        \a a contains the original matrix, results are returned in \a q and \a r, where
        \a q is a orthogonal matrix, and \a r is an upper triangular matrix, and
        the following relation holds (up to round-off errors):

        \code
        assert(a == q * r);
        \endcode

        This implementation uses householder transformations. It can be applied in-place,
        i.e. <tt>&a == &q</tt> or <tt>&a == &r</tt> are allowed.

    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear_solve.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
void qrDecomposition(MultiArrayView<2, T, C1> const & a,
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
}

    /** Solve a linear system with right-triangular defining matrix.

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

    for(unsigned int k=0; k<m; ++k)
        if(r(k,k) == NumericTraits<T>::zero())
            return false; // r doesn't have full rank.
    for(unsigned int k = 0; k < n; ++k)
    {
        x(m-1, k) = b(m-1, k) / r(m-1, m-1);
        if(m >= 2)
        {
            for(int i = m-2; i >= 0; --i)
            {
                // compute the i'th inner product, excluding the diagonal entry.
                T sum = NumericTraits<T>::zero();
                for(unsigned int j = i+1; j < m; ++j)
                    sum += r(i, j) * x(j, k);
                if(r(i, i) != NumericTraits<T>::zero())
                    x(i, k) = (b(i, k) - sum) / r(i, i);
                else
                    x(i, k) = NumericTraits<T>::zero();
            }
        }
    }
    return true;
}

    /** Solve a linear system.

        The \a a is the coefficient matrix, and the column vectors
        in \a b are the right-hand sides of the equation (so, several equations
        with the same coefficients can be solved in one go). When \a s is rectangular,
        it must have more rows than columns, and the solution is computed in the least-squares sense.
        
        The result is returned int \a res, whose columns contain the solutions for the corresponding
        columns of \a b. The number of columns of \a a must equal the number of rows of
        both \a b and \a res, and the number of columns of \a b and \a res must be
        equal. The algorithm uses QR decomposition of \a a. It fails and returns
        <tt>false</tt> if \a a doesn't have full rank. This implementation can be
        applied in-place, i.e. <tt>&b == &res</tt> or <tt>&a == &res</tt> are allowed
        (as long as the shapes match).

    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear_solve.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
bool linearSolve(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b,
                 MultiArrayView<2, T, C3> & res)
{
    vigra_precondition(columnCount(a) <= rowCount(a),
        "linearSolve(): Coefficient matrix a must have at least as many rows as columns.");
    vigra_precondition(columnCount(a) == rowCount(res) && 
                       rowCount(a) == rowCount(b) && columnCount(b) == columnCount(res),
        "linearSolve(): matrix shape mismatch.");

    Matrix<T> q(a.shape()), r(columnCount(a), columnCount(a));
    qrDecomposition(a, q, r);
    return reverseElimination(r, transpose(q) * b, res); // false if a didn't have full rank
}

//@}

} // namespace linalg

using linalg::inverse;
using linalg::linearSolve;
using linalg::qrDecomposition;
using linalg::reverseElimination;

} // namespace vigra


#endif // VIGRA_LINEAR_SOLVE_HXX
