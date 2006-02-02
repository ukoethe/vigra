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

#include "vigra/matrix.hxx"


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
    
    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear__solve.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear__solve.hxx</a>" or<br>
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
    vigra_precondition(linearSolve(v, identityMatrix<T>(n), v),
        "inverse(): matrix is not invertible.");
    return v;
}

    /** QR decomposition.
     
        \a a contains the original matrix, results are returned in \a q and \a r, where
        \a q is a orthogonal matrix, and \a r is an upper triangular matrix, and
        the following relation holds:
        
        \code
        assert(a == q * r);
        \endcode
        
        This implementation uses householder transformations. It can be applied in-place,
        i.e. <tt>&a == &r</tt> is allowed.
    
    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear__solve.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
void qrDecomposition(MultiArrayView<2, T, C1> const & a,
                     MultiArrayView<2, T, C2> &q, MultiArrayView<2, T, C3> &r)
{
    typedef typename MultiArrayView<2, T, C2>::difference_type MatrixShape;
    typedef typename MultiArray<1, T>::difference_type VectorShape;
    
    // the orthogonal matrix q will have as many rows and columns as
    // the original matrix has columns.
    const unsigned int rows = rowCount(a);
    const unsigned int cols = columnCount(a);
    vigra_precondition(cols == columnCount(r) && cols == rowCount(r) &&
                       cols == columnCount(q) && cols == rowCount(q),
                       "qrDecomposition(): Matrix shape mismatch.");

    identityMatrix(q);
    r.copy(a);   // does nothing if &r == &a
    
    for(unsigned int k = 0; (k < cols) && (k < rows - 1); ++k) {

        const unsigned int rows_left = rows - k;
        const unsigned int cols_left = cols - k;

        // create a view on the remaining part of r
        MatrixShape rul(k, k);
        MultiArrayView<2, T, C2> rsub = r.subarray(rul, r.shape());

        // decompose the first row
        MultiArrayView <1, T, C2 > vec = rsub.bindOuter(0);

        // defining householder vector
        VectorShape ushape(rows_left);
        MultiArray<1, T> u(ushape);
        for(unsigned int i = 0; i < rows_left; ++i)
            u(i) = vec(i);
        u(0) += norm(vec);

        const T divisor = squaredNorm(u);
        const T scal = (divisor == 0) ? 0.0 : 2.0 / divisor;

        // apply householder elimination on rsub
        for(unsigned int i = 0; i < cols_left; ++i) {

            // compute the inner product of the i'th column of rsub with u
            T sum = dot(u, rsub.bindOuter(i));

            // add rsub*(uu')/(u'u)
            sum *= scal;
            for(unsigned int j = 0; j < rows_left; ++j)
                rsub(j, i) -= sum * u(j); 
        }
        
        MatrixShape qul(0, k);
        MultiArrayView <2, T, C3 > qsub = q.subarray(qul, q.shape());

        // apply the (self-inverse) householder matrix on q
        for(unsigned int i = 0; i < cols; ++i) {

            // compute the inner product of the i'th row of q with u
            T sum = dot(qsub.bindInner(i), u);
            
            // add q*(uu')/(u'u)
            sum *= scal;
            for(unsigned int j = 0; j < rows_left; ++j)
                qsub(i, j) -= sum * u(j);
        }
    }
}

    /** Solve a linear system with right-triangular defining matrix.
     
        The square matrix \a a must be a right-triangular coefficient matrix as can,
        for example, be obtained by means of QR decomposition. The column vectors
        in \a b are the right-hand sides of the equation (so, several equations
        with the same coefficients can be solved in one go). The result is returned
        int \a x, whose columns contain the solutions for the correspoinding 
        columns of \a b. The number of columns of \a a must equal the number of rows of
        both \a b and \a x, and the number of columns of \a b and \a x must be
        equal. This implementation can be applied in-place, i.e. <tt>&b == &x</tt> is allowed.
    
    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear__solve.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
void reverseElimination(const MultiArrayView<2, T, C1> &r, const MultiArrayView<2, T, C2> &b,
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
}

    /** Solve a linear system.
     
        The square matrix \a a is the coefficient matrix, and the column vectors
        in \a b are the right-hand sides of the equation (so, several equations
        with the same coefficients can be solved in one go). The result is returned
        int \a res, whose columns contain the solutions for the correspoinding 
        columns of \a b. The number of columns of \a a must equal the number of rows of
        both \a b and \a res, and the number of columns of \a b and \a res must be
        equal. The algorithm uses QR decomposition of \a a. The algorithm returns
        <tt>false</tt> if \a a doesn't have full rank. This implementation  can be 
        applied in-place, i.e. <tt>&b == &res</tt> or <tt>&a == &res</tt> are allowed.
    
    <b>\#include</b> "<a href="linear__solve_8hxx-source.html">vigra/linear__solve.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
bool linearSolve(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b,
                 MultiArrayView<2, T, C3> & res)
{
    unsigned int acols = columnCount(a);
    unsigned int bcols = columnCount(b);
    vigra_precondition(acols == rowCount(a),
        "linearSolve(): square coefficient matrix required.");
    vigra_precondition(acols == rowCount(b) && acols == rowCount(res) && bcols == columnCount(res),
        "linearSolve(): matrix shape mismatch.");
    
    Matrix<T> q(acols, acols), r(a);
    qrDecomposition(r, q, r);
    for(unsigned int k=0; k<acols; ++k)
        if(r(k,k) == NumericTraits<T>::zero())
            return false; // a didn't have full rank.
    q.transpose();
    reverseElimination(r, q * b, res);
    return true;
}

//@}

} // namespace linalg

using linalg::inverse;
using linalg::linearSolve;
using linalg::qrDecomposition;
using linalg::reverseElimination;

} // namespace vigra


#endif // VIGRA_LINEAR_SOLVE_HXX
