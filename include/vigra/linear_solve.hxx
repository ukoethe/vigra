/************************************************************************/
/*                                                                      */
/*     Copyright 2003-2008 by Gunnar Kedenburg and Ullrich Koethe       */
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

namespace detail {

template <class T, class C1>
T determinantByLUDecomposition(MultiArrayView<2, T, C1> const & a)
{
    typedef MultiArrayShape<2>::type Shape;

    MultiArrayIndex m = rowCount(a), n = columnCount(a);
    vigra_precondition(n == m,
       "determinant(): square matrix required.");
       
    Matrix<T> LU(a);
    T det = 1.0;

    for (MultiArrayIndex j = 0; j < n; ++j) 
    {
        // Apply previous transformations.
        for (MultiArrayIndex i = 0; i < m; ++i) 
        {
            MultiArrayIndex end = std::min(i, j);
            T s = dot(rowVector(LU, Shape(i,0), end), columnVector(LU, Shape(0,j), end));
            LU(i,j) = LU(i,j) -= s;
        }

        // Find pivot and exchange if necessary.
        MultiArrayIndex p = j + argMax(abs(columnVector(LU, Shape(j,j), m)));
        if (p != j) 
        {
            rowVector(LU, p).swapData(rowVector(LU, j));
            det = -det;
        }
        
        det *= LU(j,j);

        // Compute multipliers.
        if (LU(j,j) != 0.0)
            columnVector(LU, Shape(j+1,j), m) /= LU(j,j);
        else
            break; // det is zero
    }
    return det;
}

// returns the new value of 'a' (when this Givens rotation is applied to 'a' and 'b')
// the new value of 'b' is zero, of course
template <class T>
T givensCoefficients(T a, T b, T & c, T & s)
{
    if(abs(a) < abs(b))
    {
        T t = a/b, 
          r = std::sqrt(1.0 + t*t);
        s = 1.0 / r;
        c = t*s;
        return r*b;
    }
    else if(a != 0.0)
    {
        T t = b/a, 
          r = std::sqrt(1.0 + t*t);
        c = 1.0 / r;
        s = t*c;
        return r*a;
    }
    else // a == b == 0.0
    {
        c = 1.0;
        s = 0.0;
        return 0.0;
    }
}

// see Golub, van Loan: Algorithm 5.1.3 (p. 216)
template <class T>
bool givensRotationMatrix(T a, T b, Matrix<T> & gTranspose)
{
    if(b == 0.0)
        return false; // no rotation needed
    givensCoefficients(a, b, gTranspose(0,0), gTranspose(0,1));
    gTranspose(1,1) = gTranspose(0,0);
    gTranspose(1,0) = -gTranspose(0,1);
    return true;
}

// reflections are symmetric matrices and can thus be applied to rows
// and columns in the same way => code simplification relative to rotations
template <class T>
inline bool 
givensReflectionMatrix(T a, T b, Matrix<T> & g)
{
    if(b == 0.0)
        return false; // no reflection needed
    givensCoefficients(a, b, g(0,0), g(0,1));
    g(1,1) = -g(0,0);
    g(1,0) = g(0,1);
    return true;
}

// see Golub, van Loan: Algorithm 5.2.2 (p. 227) and Section 12.5.2 (p. 608)
template <class T, class C1, class C2>
bool 
qrGivensStepImpl(MultiArrayIndex i, MultiArrayView<2, T, C1> r, MultiArrayView<2, T, C2> rhs)
{
    typedef typename Matrix<T>::difference_type Shape;
    
    const MultiArrayIndex m = rowCount(r);
    const MultiArrayIndex n = columnCount(r);
    const MultiArrayIndex rhsCount = columnCount(rhs);
    vigra_precondition(m == rowCount(rhs),
                       "qrGivensStepImpl(): Matrix shape mismatch.");

    Matrix<T> givens(2,2);
    for(int k=m-1; k>(int)i; --k)
    {
        if(!givensReflectionMatrix(r(k-1,i), r(k,i), givens))
            continue; // r(k,i) was already zero

        r(k-1,i) = givens(0,0)*r(k-1,i) + givens(0,1)*r(k,i);
        r(k,i) = 0.0;
        
        r.subarray(Shape(k-1,i+1), Shape(k+1,n)) = givens*r.subarray(Shape(k-1,i+1), Shape(k+1,n));
        rhs.subarray(Shape(k-1,0), Shape(k+1,rhsCount)) = givens*rhs.subarray(Shape(k-1,0), Shape(k+1,rhsCount));
    }
    return r(i,i) != 0.0;
}

// see Golub, van Loan: Section 12.5.2 (p. 608)
template <class T, class C1, class C2, class Permutation>
void 
upperTriangularCyclicShiftColumns(MultiArrayIndex i, MultiArrayIndex j, 
                                  MultiArrayView<2, T, C1> &r, MultiArrayView<2, T, C2> &rhs, Permutation & permutation)
{
    typedef typename Matrix<T>::difference_type Shape;
    
    const MultiArrayIndex m = rowCount(r);
    const MultiArrayIndex n = columnCount(r);
    const MultiArrayIndex rhsCount = columnCount(rhs);
    vigra_precondition(i < n && j < n,
                       "upperTriangularCyclicShiftColumns(): Shift indices out of range.");
    vigra_precondition(m == rowCount(rhs),
                       "upperTriangularCyclicShiftColumns(): Matrix shape mismatch.");

    if(j == i)
        return;
    if(j < i)
        std::swap(j,i);
        
    Matrix<T> t = columnVector(r, i);
    MultiArrayIndex ti = permutation[i];
    for(MultiArrayIndex k=i; k<j;++k)
    {
        columnVector(r, k) = columnVector(r, k+1);
        permutation[k] = permutation[k+1];
    }
    columnVector(r, j) = t;
    permutation[j] = ti;
    
    Matrix<T> givens(2,2);
    for(MultiArrayIndex k=i; k<j; ++k)
    {
        if(!givensReflectionMatrix(r(k,k), r(k+1,k), givens))
            continue;  // r(k+1,k) was already zero
        
        r(k,k) = givens(0,0)*r(k,k) + givens(0,1)*r(k+1,k);
        r(k+1,k) = 0.0;
        
        r.subarray(Shape(k,k+1), Shape(k+2,n)) = givens*r.subarray(Shape(k,k+1), Shape(k+2,n));
        rhs.subarray(Shape(k,0), Shape(k+2,rhsCount)) = givens*rhs.subarray(Shape(k,0), Shape(k+2,rhsCount));
    }
}

// see Golub, van Loan: Section 12.5.2 (p. 608)
template <class T, class C1, class C2, class Permutation>
void 
upperTriangularSwapColumns(MultiArrayIndex i, MultiArrayIndex j, 
                           MultiArrayView<2, T, C1> &r, MultiArrayView<2, T, C2> &rhs, Permutation & permutation)
{    
    typedef typename Matrix<T>::difference_type Shape;
    
    const MultiArrayIndex m = rowCount(r);
    const MultiArrayIndex n = columnCount(r);
    const MultiArrayIndex rhsCount = columnCount(rhs);
    vigra_precondition(i < n && j < n,
                       "upperTriangularSwapColumns(): Swap indices out of range.");
    vigra_precondition(m == rowCount(rhs),
                       "upperTriangularSwapColumns(): Matrix shape mismatch.");

    if(j == i)
        return;
    if(j < i)
        std::swap(j,i);

    columnVector(r, i).swapData(columnVector(r, j));
    std::swap(permutation[i], permutation[j]);
    
    Matrix<T> givens(2,2);
    for(int k=m-1; k>(int)i; --k)
    {
        if(!givensReflectionMatrix(r(k-1,i), r(k,i), givens))
            continue; // r(k,i) was already zero

        r(k-1,i) = givens(0,0)*r(k-1,i) + givens(0,1)*r(k,i);
        r(k,i) = 0.0;
        
        r.subarray(Shape(k-1,i+1), Shape(k+1,n)) = givens*r.subarray(Shape(k-1,i+1), Shape(k+1,n));
        rhs.subarray(Shape(k-1,0), Shape(k+1,rhsCount)) = givens*rhs.subarray(Shape(k-1,0), Shape(k+1,rhsCount));
    }
    MultiArrayIndex end = std::min(j, m-1);
    for(MultiArrayIndex k=i+1; k<end; ++k)
    {
        if(!givensReflectionMatrix(r(k,k), r(k+1,k), givens))
            continue;  // r(k+1,k) was already zero
        
        r(k,k) = givens(0,0)*r(k,k) + givens(0,1)*r(k+1,k);
        r(k+1,k) = 0.0;
        
        r.subarray(Shape(k,k+1), Shape(k+2,n)) = givens*r.subarray(Shape(k,k+1), Shape(k+2,n));
        rhs.subarray(Shape(k,0), Shape(k+2,rhsCount)) = givens*rhs.subarray(Shape(k,0), Shape(k+2,rhsCount));
    }
}

// see Lawson & Hanson: Algorithm H1 (p. 57)
template <class T, class C1, class C2, class U>
bool householderVector(MultiArrayView<2, T, C1> const & v, MultiArrayView<2, T, C2> & u, U & vnorm)
{
    vnorm = (v(0,0) > 0.0)
                 ? -norm(v)
                 :  norm(v);
    U f = std::sqrt(vnorm*(vnorm - v(0,0)));
    
    if(f == NumericTraits<U>::zero())
    {
        u.init(NumericTraits<T>::zero());
        return false;
    }
    else
    {
        u(0,0) = (v(0,0) - vnorm) / f;
        for(MultiArrayIndex k=1; k<rowCount(u); ++k)
            u(k,0) = v(k,0) / f;
        return true;
    }
}

// see Lawson & Hanson: Algorithm H1 (p. 57)
template <class T, class C1, class C2, class C3>
bool 
qrHouseholderStepImpl(MultiArrayIndex i, MultiArrayView<2, T, C1> & r, 
                      MultiArrayView<2, T, C2> & rhs, MultiArrayView<2, T, C3> & householderMatrix)
{
    typedef typename Matrix<T>::difference_type Shape;
    
    const MultiArrayIndex m = rowCount(r);
    const MultiArrayIndex n = columnCount(r);
    const MultiArrayIndex rhsCount = columnCount(rhs);

    vigra_precondition(i < n && i < m,
        "qrHouseholderStepImpl(): Index i out of range.");

    Matrix<T> u(m-i,1);
    T vnorm;
    bool nontrivial = householderVector(columnVector(r, Shape(i,i), m), u, vnorm);
    
    r(i,i) = vnorm;
    columnVector(r, Shape(i+1,i), m).init(NumericTraits<T>::zero());

    if(columnCount(householderMatrix) == n)
        columnVector(householderMatrix, Shape(i,i), m) = u;

    if(nontrivial)
    {
        for(MultiArrayIndex k=i+1; k<n; ++k)
            columnVector(r, Shape(i,k), m) -= dot(columnVector(r, Shape(i,k), m), u) * u;
        for(MultiArrayIndex k=0; k<rhsCount; ++k)
            columnVector(rhs, Shape(i,k), m) -= dot(columnVector(rhs, Shape(i,k), m), u) * u;
    }
    return r(i,i) != 0.0;
}

template <class T, class C1, class C2>
bool 
qrColumnHouseholderStep(MultiArrayIndex i, MultiArrayView<2, T, C1> &r, MultiArrayView<2, T, C2> &rhs)
{
    Matrix<T> dontStoreHouseholderVectors; // intentionally empty
    return qrHouseholderStepImpl(i, r, rhs, dontStoreHouseholderVectors);
}

template <class T, class C1, class C2>
bool 
qrRowHouseholderStep(MultiArrayIndex i, MultiArrayView<2, T, C1> &r, MultiArrayView<2, T, C2> & householderMatrix)
{
    Matrix<T> dontTransformRHS; // intentionally empty
    MultiArrayView<2, T, StridedArrayTag> rt = transpose(r),
                                          ht = transpose(householderMatrix);
    return qrHouseholderStepImpl(i, rt, dontTransformRHS, ht);
}

// O(n) algorithm due to Bischof: Incremental Condition Estimation, 1990
template <class T, class C1, class C2, class SNType>
void
incrementalMaxSingularValueApproximation(MultiArrayView<2, T, C1> const & newColumn, 
                                         MultiArrayView<2, T, C2> & z, SNType & v) 
{
    typedef typename Matrix<T>::difference_type Shape;
    MultiArrayIndex n = rowCount(newColumn) - 1;
    
    SNType vneu = squaredNorm(newColumn);
    T yv = dot(columnVector(newColumn, Shape(0,0),n), columnVector(z, Shape(0,0),n));
    // use atan2 as it is robust against overflow/underflow
    T t = 0.5*std::atan2(T(2.0*yv), T(sq(v)-vneu)),
      s = std::sin(t),
      c = std::cos(t);
    v = std::sqrt(sq(c*v) + sq(s)*vneu + 2.0*s*c*yv);
    columnVector(z, Shape(0,0),n) = c*columnVector(z, Shape(0,0),n) + s*columnVector(newColumn, Shape(0,0),n);
    z(n,0) = s*newColumn(n,0);
}

// O(n) algorithm due to Bischof: Incremental Condition Estimation, 1990
template <class T, class C1, class C2, class SNType>
void
incrementalMinSingularValueApproximation(MultiArrayView<2, T, C1> const & newColumn, 
                                         MultiArrayView<2, T, C2> & z, SNType & v, double tolerance) 
{
    typedef typename Matrix<T>::difference_type Shape;

    if(v <= tolerance)
    {
        v = 0.0;
        return;
    }

    MultiArrayIndex n = rowCount(newColumn) - 1;
    
    T gamma = newColumn(n,0);
    if(gamma == 0.0) 
    {
        v = 0.0;
        return;
    }
    
    T yv = dot(columnVector(newColumn, Shape(0,0),n), columnVector(z, Shape(0,0),n));
    // use atan2 as it is robust against overflow/underflow
    T t = 0.5*std::atan2(T(-2.0*yv), T(squaredNorm(gamma / v) + squaredNorm(yv) - 1.0)),
      s = std::sin(t),
      c = std::cos(t);
    columnVector(z, Shape(0,0),n) *= c;
    z(n,0) = (s - c*yv) / gamma;
    v *= norm(gamma) / hypot(c*gamma, v*(s - c*yv));
}

// QR algorithm with optional column pivoting
template <class T, class C1, class C2, class C3>
unsigned int 
qrTransformToTriangularImpl(MultiArrayView<2, T, C1> & r, MultiArrayView<2, T, C2> & rhs, MultiArrayView<2, T, C3> & householder,
                            ArrayVector<MultiArrayIndex> & permutation, double epsilon)
{
    typedef typename Matrix<T>::difference_type Shape;
    typedef typename NormTraits<MultiArrayView<2, T, C1> >::NormType NormType;
    typedef typename NormTraits<MultiArrayView<2, T, C1> >::SquaredNormType SNType;
    
    const MultiArrayIndex m = rowCount(r);
    const MultiArrayIndex n = columnCount(r);
    const MultiArrayIndex maxRank = std::min(m, n);
    
    vigra_precondition(m >= n,
        "qrTransformToTriangularImpl(): Coefficient matrix with at least as many rows as columns required.");

    const MultiArrayIndex rhsCount = columnCount(rhs);
    bool transformRHS = rhsCount > 0;
    vigra_precondition(!transformRHS || m == rowCount(rhs),
                       "qrTransformToTriangularImpl(): RHS matrix shape mismatch.");
                       
    bool storeHouseholderSteps = columnCount(householder) > 0;
    vigra_precondition(!storeHouseholderSteps || r.shape() == householder.shape(),
                       "qrTransformToTriangularImpl(): Householder matrix shape mismatch.");
                       
    bool pivoting = permutation.size() > 0;
    vigra_precondition(!pivoting || n == (MultiArrayIndex)permutation.size(),
                       "qrTransformToTriangularImpl(): Permutation array size mismatch.");

    if(n == 0)
        return 0; // trivial solution
        
    Matrix<SNType> columnSquaredNorms;
    if(pivoting)
    {
        columnSquaredNorms.reshape(Shape(1,n));
        for(MultiArrayIndex k=0; k<n; ++k)
            columnSquaredNorms[k] = squaredNorm(columnVector(r, k));
            
        int pivot = argMax(columnSquaredNorms);
        if(pivot != 0)
        {
            columnVector(r, 0).swapData(columnVector(r, pivot));
            std::swap(columnSquaredNorms[0], columnSquaredNorms[pivot]);
            std::swap(permutation[0], permutation[pivot]);
        }
    }
    
    qrHouseholderStepImpl(0, r, rhs, householder);
    
    MultiArrayIndex rank = 1;
    NormType maxApproxSingularValue = norm(r(0,0)),
             minApproxSingularValue = maxApproxSingularValue;
    
    double tolerance = (epsilon == 0.0)
                          ? m*maxApproxSingularValue*NumericTraits<T>::epsilon()
                          : epsilon;
    
    bool simpleSingularValueApproximation = (n < 4);
    Matrix<T> zmax, zmin;
    if(minApproxSingularValue <= tolerance)
    {
        rank = 0;
        pivoting = false;
        simpleSingularValueApproximation = true;
    }
    if(!simpleSingularValueApproximation)
    {
        zmax.reshape(Shape(m,1));
        zmin.reshape(Shape(m,1));
        zmax(0,0) = r(0,0);
        zmin(0,0) = 1.0 / r(0,0);
    }

    for(MultiArrayIndex k=1; k<maxRank; ++k)
    {
        if(pivoting)
        {
            for(MultiArrayIndex l=k; l<n; ++l)
                columnSquaredNorms[l] -= squaredNorm(r(k, l));
            int pivot = k + argMax(rowVector(columnSquaredNorms, Shape(0,k), n));
            if(pivot != (int)k)
            {
                columnVector(r, k).swapData(columnVector(r, pivot));
                std::swap(columnSquaredNorms[k], columnSquaredNorms[pivot]);
                std::swap(permutation[k], permutation[pivot]);
            }
        }
        
        qrHouseholderStepImpl(k, r, rhs, householder);

        if(simpleSingularValueApproximation)
        {
            NormType nv = norm(r(k,k));        
            maxApproxSingularValue = std::max(nv, maxApproxSingularValue);
            minApproxSingularValue = std::min(nv, minApproxSingularValue);
        }
        else
        {
            incrementalMaxSingularValueApproximation(columnVector(r, Shape(0,k),k+1), zmax, maxApproxSingularValue);
            incrementalMinSingularValueApproximation(columnVector(r, Shape(0,k),k+1), zmin, minApproxSingularValue, tolerance);
        }
        
#if 0
        Matrix<T> u(k+1,k+1), s(k+1, 1), v(k+1,k+1);
        singularValueDecomposition(r.subarray(Shape(0,0), Shape(k+1,k+1)), u, s, v);
        std::cerr << "estimate, svd " << k << ": " << minApproxSingularValue << " " << s(k,0) << "\n";
#endif
        
        if(epsilon == 0.0)
            tolerance = m*maxApproxSingularValue*NumericTraits<T>::epsilon();

        if(minApproxSingularValue > tolerance)
            ++rank;
        else
            pivoting = false; // matrix doesn't have full rank, triangulize the rest without pivoting
    }
    return (unsigned int)rank;
}

template <class T, class C1, class C2>
unsigned int 
qrTransformToUpperTriangular(MultiArrayView<2, T, C1> & r, MultiArrayView<2, T, C2> & rhs, 
                             ArrayVector<MultiArrayIndex> & permutation, double epsilon = 0.0)
{
    Matrix<T> dontStoreHouseholderVectors; // intentionally empty
    return qrTransformToTriangularImpl(r, rhs, dontStoreHouseholderVectors, permutation, epsilon);
}

// QR algorithm with optional row pivoting
template <class T, class C1, class C2, class C3>
unsigned int 
qrTransformToLowerTriangular(MultiArrayView<2, T, C1> & r, MultiArrayView<2, T, C2> & rhs, MultiArrayView<2, T, C3> & householderMatrix, 
                      double epsilon = 0.0)
{
    ArrayVector<MultiArrayIndex> permutation((unsigned int)rowCount(rhs));
    for(MultiArrayIndex k=0; k<(MultiArrayIndex)permutation.size(); ++k)
        permutation[k] = k;
    Matrix<T> dontTransformRHS; // intentionally empty
    MultiArrayView<2, T, StridedArrayTag> rt = transpose(r),
                                          ht = transpose(householderMatrix);
    unsigned int rank = qrTransformToTriangularImpl(rt, dontTransformRHS, ht, permutation, epsilon);
    
    // apply row permutation to RHS
    Matrix<T> tempRHS(rhs);
    for(MultiArrayIndex k=0; k<(MultiArrayIndex)permutation.size(); ++k)
        rowVector(rhs, k) = rowVector(tempRHS, permutation[k]);
    return rank;
}

// QR algorithm without column pivoting
template <class T, class C1, class C2>
inline bool
qrTransformToUpperTriangular(MultiArrayView<2, T, C1> & r, MultiArrayView<2, T, C2> & rhs, 
                      double epsilon = 0.0)
{
    ArrayVector<MultiArrayIndex> noPivoting; // intentionally empty
    
    return (qrTransformToUpperTriangular(r, rhs, noPivoting, epsilon) == 
            (unsigned int)columnCount(r));
}

// QR algorithm without row pivoting
template <class T, class C1, class C2>
inline bool
qrTransformToLowerTriangular(MultiArrayView<2, T, C1> & r, MultiArrayView<2, T, C2> & householder, 
                      double epsilon = 0.0)
{
    Matrix<T> noPivoting; // intentionally empty
    
    return (qrTransformToLowerTriangular(r, noPivoting, householder, epsilon) == 
           (unsigned int)rowCount(r));
}

// restore ordering of result vector elements after QR solution with column pivoting
template <class T, class C1, class C2, class Permutation>
void inverseRowPermutation(MultiArrayView<2, T, C1> &permuted, MultiArrayView<2, T, C2> &res,
                           Permutation const & permutation)
{
    for(MultiArrayIndex k=0; k<columnCount(permuted); ++k)
        for(MultiArrayIndex l=0; l<rowCount(permuted); ++l)
            res(permutation[l], k) = permuted(l,k);
}

template <class T, class C1, class C2>
void applyHouseholderColumnReflections(MultiArrayView<2, T, C1> const &householder, MultiArrayView<2, T, C2> &res)
{
    typedef typename Matrix<T>::difference_type Shape;
    MultiArrayIndex n = rowCount(householder);
    MultiArrayIndex m = columnCount(householder);
    MultiArrayIndex rhsCount = columnCount(res);
    
    for(int k = m-1; k >= 0; --k)
    {
        MultiArrayView<2, T, C1> u = columnVector(householder, Shape(k,k), n);
        for(MultiArrayIndex l=0; l<rhsCount; ++l)
            columnVector(res, Shape(k,l), n) -= dot(columnVector(res, Shape(k,l), n), u) * u;
    }
}

} // namespace detail

template <class T, class C1, class C2, class C3>
unsigned int 
linearSolveQRReplace(MultiArrayView<2, T, C1> &A, MultiArrayView<2, T, C2> &b,
                     MultiArrayView<2, T, C3> & res, 
                     double epsilon = 0.0)
{
    typedef typename Matrix<T>::difference_type Shape;

    MultiArrayIndex n = columnCount(A);
    MultiArrayIndex m = rowCount(A);
    MultiArrayIndex rhsCount = columnCount(res);
    MultiArrayIndex rank = std::min(m,n);
    Shape ul(MultiArrayIndex(0), MultiArrayIndex(0));

    
    vigra_precondition(rhsCount == columnCount(b),
           "linearSolveQR(): RHS and solution must have the same number of columns.");
    vigra_precondition(m == rowCount(b),
           "linearSolveQR(): Coefficient matrix and RHS must have the same number of rows.");
    vigra_precondition(n == rowCount(res),
           "linearSolveQR(): Mismatch between column count of coefficient matrix and row count of solution.");
    vigra_precondition(epsilon >= 0.0,
           "linearSolveQR(): 'epsilon' must be non-negative.");
    
    if(m < n)
    {
        // minimum norm solution of underdetermined system
        Matrix<T> householderMatrix(n, m);
        MultiArrayView<2, T, StridedArrayTag> ht = transpose(householderMatrix);
        rank = (MultiArrayIndex)detail::qrTransformToLowerTriangular(A, b, ht, epsilon);
        res.subarray(Shape(rank,0), Shape(n, rhsCount)).init(NumericTraits<T>::zero());
        if(rank < m)
        {
            // system is also rank-deficient => compute minimum norm least squares solution
            MultiArrayView<2, T, C1> Asub = A.subarray(ul, Shape(m,rank));
            detail::qrTransformToUpperTriangular(Asub, b, epsilon);
            linearSolveUpperTriangular(A.subarray(ul, Shape(rank,rank)), 
                                       b.subarray(ul, Shape(rank,rhsCount)), 
                                       res.subarray(ul, Shape(rank, rhsCount)));
        }
        else
        {
            // system has full rank => compute minimum norm solution
            linearSolveLowerTriangular(A.subarray(ul, Shape(rank,rank)), 
                                       b.subarray(ul, Shape(rank, rhsCount)), 
                                       res.subarray(ul, Shape(rank, rhsCount)));
        }
        detail::applyHouseholderColumnReflections(householderMatrix.subarray(ul, Shape(n, rank)), res);
    }
    else
    {
        // solution of well-determined or overdetermined system
        ArrayVector<MultiArrayIndex> permutation((unsigned int)n);
        for(MultiArrayIndex k=0; k<n; ++k)
            permutation[k] = k;

        rank = (MultiArrayIndex)detail::qrTransformToUpperTriangular(A, b, permutation, epsilon);
        
        Matrix<T> permutedSolution(n, rhsCount);
        if(rank < n)
        {
            // system is rank-deficient => compute minimum norm solution
            Matrix<T> householderMatrix(n, rank);
            MultiArrayView<2, T, StridedArrayTag> ht = transpose(householderMatrix);
            MultiArrayView<2, T, C1> Asub = A.subarray(ul, Shape(rank,n));
            detail::qrTransformToLowerTriangular(Asub, ht, epsilon);
            linearSolveLowerTriangular(A.subarray(ul, Shape(rank,rank)), 
                                       b.subarray(ul, Shape(rank, rhsCount)), 
                                       permutedSolution.subarray(ul, Shape(rank, rhsCount)));
            detail::applyHouseholderColumnReflections(householderMatrix, permutedSolution);
        }
        else
        {
            // system has full rank => compute exact or least squares solution
            linearSolveUpperTriangular(A.subarray(ul, Shape(rank,rank)), 
                                       b.subarray(ul, Shape(rank,rhsCount)), 
                                       permutedSolution);
        }
        detail::inverseRowPermutation(permutedSolution, res, permutation);
    }
    return (unsigned int)rank;
}

template <class T, class C1, class C2, class C3>
unsigned int linearSolveQR(MultiArrayView<2, T, C1> const & A, MultiArrayView<2, T, C2> const & b,
                                  MultiArrayView<2, T, C3> & res)
{
    Matrix<T> r(A), rhs(b);
    return linearSolveQRReplace(r, rhs, res);
}

/** \defgroup MatrixAlgebra Advanced Matrix Algebra
    
    \brief Solution of linear systems, eigen systems, linear least squares etc.
    
    \ingroup LinearAlgebraModule
 */
//@{
    /** Create the inverse or pseudo-inverse of matrix \a v.

        If the matrix \a v is square, \a res must have the same shape and will contain the
        inverse of \a v. If \a v is rectangular, \a res must have the transposed shape 
        of \a v. The inverse is then computed in the least-squares 
        sense, i.e. \a res will be the pseudo-inverse (Moore-Penrose inverse).
        The function returns <tt>true</tt> upon success, and <tt>false</tt> if \a v 
        is not invertible (has not full rank). The inverse is computed by means of QR 
        decomposition. This function can be applied in-place.
        
    <b>\#include</b> \<vigra/linear_solve.hxx\> or<br>
    <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2>
bool inverse(const MultiArrayView<2, T, C1> &v, MultiArrayView<2, T, C2> &res)
{
    typedef typename MultiArrayShape<2>::type Shape;
    
    const MultiArrayIndex n = columnCount(v);
    const MultiArrayIndex m = rowCount(v);
    vigra_precondition(n == rowCount(res) && m == columnCount(res),
       "inverse(): shape of output matrix must be the transpose of the input matrix' shape.");
    
    if(m < n)
    {
        MultiArrayView<2, T, StridedArrayTag> vt = transpose(v);
        Matrix<T> r(vt.shape()), q(n, n);
        if(!qrDecomposition(vt, q, r))
            return false; // a didn't have full rank
        linearSolveUpperTriangular(r.subarray(Shape(0,0), Shape(m,m)), 
                                   transpose(q).subarray(Shape(0,0), Shape(m,n)), 
                                   transpose(res)); 
    }
    else
    {
        Matrix<T> r(v.shape()), q(m, m);
        if(!qrDecomposition(v, q, r))
            return false; // a didn't have full rank
        linearSolveUpperTriangular(r.subarray(Shape(0,0), Shape(n,n)), 
                                   transpose(q).subarray(Shape(0,0), Shape(n,m)), 
                                   res); 
    }
    return true;
}

    /** Create the inverse or pseudo-inverse of matrix \a v.

        The result is returned as a temporary matrix. If the matrix \a v is square, 
        the result will have the same shape and contains the inverse of \a v. 
        If \a v is rectangular, the result will have the transposed shape of \a v. 
        The inverse is then computed in the least-squares 
        sense, i.e. \a res will be the pseudo-inverse (Moore-Penrose inverse).
        The inverse is computed by means of QR decomposition. If \a v
        is not invertible, <tt>vigra::PreconditionViolation</tt> exception is thrown.
        Usage:

        \code
        vigra::Matrix<double> v(n, n);
        v = ...;

        vigra::Matrix<double> m = inverse(v);
        \endcode

    <b>\#include</b> \<vigra/linear_solve.hxx\> or<br>
    <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
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

    /** Compute the determinant of a square matrix.

        \a method must be one of the following:
        <DL>
        <DT>"Cholesky"<DD> Compute the solution by means of Cholesky decomposition. This
                           method is faster than "LU", but requires the matrix \a a 
                           to be symmetric positive definite. If this is 
                           not the case, a <tt>ContractViolation</tt> exception is thrown.
                           
        <DT>"LU"<DD> (default) Compute the solution by means of LU decomposition.
        </DL>

    <b>\#include</b> \<vigra/linear_solve.hxx\> or<br>
    <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1>
T determinant(MultiArrayView<2, T, C1> const & a, std::string method = "LU")
{
    MultiArrayIndex n = columnCount(a);
    vigra_precondition(rowCount(a) == n,
               "determinant(): Square matrix required.");    

    method = tolower(method);
    
    if(n == 1)
        return a(0,0);
    if(n == 2)
        return a(0,0)*a(1,1) - a(0,1)*a(1,0);
    if(method == "lu")
    {
        return detail::determinantByLUDecomposition(a);
    }
    else if(method == "cholesky")
    {
        Matrix<T> L(a.shape());
        vigra_precondition(choleskyDecomposition(a, L),
           "determinant(): Cholesky method requires symmetric positive definite matrix.");
        T det = L(0,0);
        for(MultiArrayIndex k=1; k<n; ++k)
            det *= L(k,k);
        return sq(det);
    }
    else
    {
        vigra_precondition(false, "determinant(): Unknown solution method.");
    }
    return T();
}

    /** Compute the logarithm of the determinant of a symmetric positive definite matrix.

        This is useful to avoid multiplication of very large numbers in big matrices.
        It is implemented by means of Cholesky decomposition.

    <b>\#include</b> \<vigra/linear_solve.hxx\> or<br>
    <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1>
T logDeterminant(MultiArrayView<2, T, C1> const & a)
{
    MultiArrayIndex n = columnCount(a);
    vigra_precondition(rowCount(a) == n,
               "logDeterminant(): Square matrix required.");    
    if(n == 1)
    {
        vigra_precondition(a(0,0) > 0.0,
                   "logDeterminant(): Matrix not positive definite.");    
        return std::log(a(0,0));
    }
    if(n == 2)
    {
        T det = a(0,0)*a(1,1) - a(0,1)*a(1,0);
        vigra_precondition(det > 0.0,
                   "logDeterminant(): Matrix not positive definite.");    
        return std::log(det);
    }
    else
    {
        Matrix<T> L(a.shape());
        vigra_precondition(choleskyDecomposition(a, L),
                "logDeterminant(): Matrix not positive definite.");  
        T logdet = std::log(L(0,0)); 
        for(MultiArrayIndex k=1; k<n; ++k)
            logdet += std::log(L(k,k));  // L(k,k) is guaranteed to be positive
        return 2.0*logdet;
    }
}

    /** Cholesky decomposition.

        \a A must be a symmetric positive definite matrix, and \a L will be a lower
        triangular matrix, such that (up to round-off errors):

        \code
        A == L * transpose(L);
        \endcode

        This implementation cannot be applied in-place, i.e. <tt>&L == &A</tt> is an error.
        If \a A is not symmetric, a <tt>ContractViolation</tt> exception is thrown. If it
        is not positive definite, the function returns <tt>false</tt>.

    <b>\#include</b> \<vigra/linear_solve.hxx\> or<br>
    <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2>
bool choleskyDecomposition(MultiArrayView<2, T, C1> const & A,
                           MultiArrayView<2, T, C2> &L)
{
    MultiArrayIndex n = columnCount(A); 
    vigra_precondition(rowCount(A) == n,
                       "choleskyDecomposition(): Input matrix must be square.");
    vigra_precondition(n == columnCount(L) && n == rowCount(L),
                       "choleskyDecomposition(): Output matrix must have same shape as input matrix.");
    vigra_precondition(isSymmetric(A),
                       "choleskyDecomposition(): Input matrix must be symmetric.");

    for (MultiArrayIndex j = 0; j < n; ++j) 
    {
        T d(0.0);
        for (MultiArrayIndex k = 0; k < j; ++k) 
        {
            T s(0.0);
            for (MultiArrayIndex i = 0; i < k; ++i) 
            {
               s += L(k, i)*L(j, i);
            }
            L(j, k) = s = (A(j, k) - s)/L(k, k);
            d = d + s*s;
        }
        d = A(j, j) - d;
        if(d <= 0.0)
            return false;  // A is not positive definite
        L(j, j) = std::sqrt(d);
        for (MultiArrayIndex k = j+1; k < n; ++k) 
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

        If \a a doesn't have full rank, the function returns <tt>false</tt>. 
        The decomposition is computed by householder transformations. It can be applied in-place,
        i.e. <tt>&a == &q</tt> or <tt>&a == &r</tt> are allowed.

    <b>\#include</b> \<vigra/linear_solve.hxx\> or<br>
    <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
bool qrDecomposition(MultiArrayView<2, T, C1> const & a,
                     MultiArrayView<2, T, C2> &q, MultiArrayView<2, T, C3> &r,
                     double epsilon = 0.0)
{
    const MultiArrayIndex m = rowCount(a);
    const MultiArrayIndex n = columnCount(a);
    vigra_precondition(n == columnCount(r) && m == rowCount(r) &&
                       m == columnCount(q) && m == rowCount(q),
                       "qrDecomposition(): Matrix shape mismatch.");

    q = identityMatrix<T>(m);
    MultiArrayView<2,T, StridedArrayTag> tq = transpose(q);
    r = a;
    ArrayVector<MultiArrayIndex> noPivoting; // intentionally empty
    return ((MultiArrayIndex)detail::qrTransformToUpperTriangular(r, tq, noPivoting, epsilon) == std::min(m,n));
}

    /** Deprecated, use \ref linearSolveUpperTriangular().
     */
template <class T, class C1, class C2, class C3>
inline 
bool reverseElimination(const MultiArrayView<2, T, C1> &r, const MultiArrayView<2, T, C2> &b,
                        MultiArrayView<2, T, C3> x)
{
    return linearSolveUpperTriangular(r, b, x);
}

    /** Solve a linear system with upper-triangular coefficient matrix.

        The square matrix \a r must be an upper-triangular coefficient matrix as can,
        for example, be obtained by means of QR decomposition. If \a r doesn't have full rank
        the function fails and returns <tt>false</tt>, otherwise it returns <tt>true</tt>. The 
        lower triangular part of matrix \a r will not be touched, so it doesn't need to contain zeros.
        
        The column vectors of matrix \a b are the right-hand sides of the equation (several equations
        with the same coefficients can thus be solved in one go). The result is returned
        int \a x, whose columns contain the solutions for the corresponding
        columns of \a b. This implementation can be applied in-place, i.e. <tt>&b == &x</tt> is allowed.
        The following size requirements apply:
        
        \code
        rowCount(r) == columnCount(r);
        rowCount(r) == rowCount(b);
        columnCount(r) == rowCount(x);
        columnCount(b) == columnCount(x);
        \endcode

    <b>\#include</b> \<vigra/linear_solve.hxx\> or<br>
    <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
bool linearSolveUpperTriangular(const MultiArrayView<2, T, C1> &r, const MultiArrayView<2, T, C2> &b,
                                MultiArrayView<2, T, C3> x)
{
    typedef MultiArrayShape<2>::type Shape;
    MultiArrayIndex m = rowCount(r);
    MultiArrayIndex rhsCount = columnCount(b);
    vigra_precondition(m == columnCount(r),
        "linearSolveUpperTriangular(): square coefficient matrix required.");
    vigra_precondition(m == rowCount(b) && m == rowCount(x) && rhsCount == columnCount(x),
        "linearSolveUpperTriangular(): matrix shape mismatch.");

    for(MultiArrayIndex k = 0; k < rhsCount; ++k)
    {
        for(int i=m-1; i>=0; --i)
        {
            if(r(i,i) == NumericTraits<T>::zero())
                return false;  // r doesn't have full rank
            T sum = b(i, k);
            for(MultiArrayIndex j=i+1; j<m; ++j)
                 sum -= r(i, j) * x(j, k);
            x(i, k) = sum / r(i, i);
        }
    }
    return true;
}

    /** Solve a linear system with lower-triangular coefficient matrix.

        The square matrix \a l must be a lower-triangular coefficient matrix. If \a l 
        doesn't have full rank the function fails and returns <tt>false</tt>, 
        otherwise it returns <tt>true</tt>. The upper triangular part of matrix \a l will not be touched, 
        so it doesn't need to contain zeros.
        
        The column vectors of matrix \a b are the right-hand sides of the equation (several equations
        with the same coefficients can thus be solved in one go). The result is returned
        in \a x, whose columns contain the solutions for the corresponding
        columns of \a b. This implementation can be applied in-place, i.e. <tt>&b == &x</tt> is allowed.
        The following size requirements apply:
        
        \code
        rowCount(l) == columnCount(l);
        rowCount(l) == rowCount(b);
        columnCount(l) == rowCount(x);
        columnCount(b) == columnCount(x);
        \endcode

    <b>\#include</b> \<vigra/linear_solve.hxx\> or<br>
    <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
bool linearSolveLowerTriangular(const MultiArrayView<2, T, C1> &l, const MultiArrayView<2, T, C2> &b,
                            MultiArrayView<2, T, C3> x)
{
    MultiArrayIndex m = columnCount(l);
    MultiArrayIndex n = columnCount(b);
    vigra_precondition(m == rowCount(l),
        "linearSolveLowerTriangular(): square coefficient matrix required.");
    vigra_precondition(m == rowCount(b) && m == rowCount(x) && n == columnCount(x),
        "linearSolveLowerTriangular(): matrix shape mismatch.");

    for(MultiArrayIndex k = 0; k < n; ++k)
    {
        for(MultiArrayIndex i=0; i<m; ++i)
        {
            if(l(i,i) == NumericTraits<T>::zero())
                return false;  // l doesn't have full rank
            T sum = b(i, k);
            for(MultiArrayIndex j=0; j<i; ++j)
                 sum -= l(i, j) * x(j, k);
            x(i, k) = sum / l(i, i);
        }
    }
    return true;
}


    /** Solve a linear system when the Cholesky decomposition of the left hand side is given.

        The square matrix \a L must be a lower-triangular matrix resulting from Cholesky
        decomposition of some positive definite coefficient matrix.
        
        The column vectors of matrix \a b are the right-hand sides of the equation (several equations
        with the same matrix \a L can thus be solved in one go). The result is returned
        in \a x, whose columns contain the solutions for the corresponding
        columns of \a b. This implementation can be applied in-place, i.e. <tt>&b == &x</tt> is allowed.
        The following size requirements apply:
        
        \code
        rowCount(L) == columnCount(L);
        rowCount(L) == rowCount(b);
        columnCount(L) == rowCount(x);
        columnCount(b) == columnCount(x);
        \endcode

    <b>\#include</b> \<vigra/linear_solve.hxx\> or<br>
    <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
inline 
void choleskySolve(MultiArrayView<2, T, C1> & L, MultiArrayView<2, T, C2> const & b, MultiArrayView<2, T, C3> & x)
{
    /* Solve L * y = b */
    linearSolveLowerTriangular(L, b, x);
    /* Solve L^T * x = y */
    linearSolveUpperTriangular(transpose(L), x, x);
}

    /** Solve a linear system.

        \a A is the coefficient matrix, and the column vectors
        in \a b are the right-hand sides of the equation (so, several equations
        with the same coefficients can be solved in one go). The result is returned 
        in \a res, whose columns contain the solutions for the corresponding
        columns of \a b. The number of columns of \a A must equal the number of rows of
        both \a b and \a res, and the number of columns of \a b and \a res must match. 
        
        \a method must be one of the following:
        <DL>
        <DT>"Cholesky"<DD> Compute the solution by means of Cholesky decomposition. The 
                           coefficient matrix \a A must by symmetric positive definite. If
                           this is not the case, the function returns <tt>false</tt>.
                           
        <DT>"QR"<DD> (default) Compute the solution by means of QR decomposition.  The 
                           coefficient matrix \a A can be square or rectangular. In the latter case,
                           it must have more rows than columns, and the solution will be computed in the 
                           least squares sense. If \a A doesn't have full rank, the function 
                           returns <tt>false</tt>.

        <DT>"SVD"<DD> Compute the solution by means of singular value decomposition.  The 
                           coefficient matrix \a A can be square or rectangular. In the latter case,
                           it must have more rows than columns, and the solution will be computed in the 
                           least squares sense. If \a A doesn't have full rank, the function 
                           returns <tt>false</tt>.

        <DT>"NE"<DD> Compute the solution by means of the normal equations, i.e. by applying Cholesky
                           decomposition to the equivalent problem <tt>A'*A*x = A'*b</tt>. This only makes sense
                           when the equation is to be solved in the least squares sense, i.e. when \a A is a 
                           rectangular matrix with more rows than columns. If \a A doesn't have full column rank, 
                           the function returns <tt>false</tt>.
        </DL>
        
        This function can be applied in-place, i.e. <tt>&b == &res</tt> or <tt>&A == &res</tt> are allowed
        (provided they have the required shapes).

        The following size requirements apply:
        
        \code
        rowCount(r) == rowCount(b);
        columnCount(r) == rowCount(x);
        columnCount(b) == columnCount(x);
        \endcode

    <b>\#include</b> \<vigra/linear_solve.hxx\> or<br>
    <b>\#include</b> \<vigra/linear_algebra.hxx\><br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
bool linearSolve(const MultiArrayView<2, T, C1> &A, const MultiArrayView<2, T, C2> &b,
                 MultiArrayView<2, T, C3> & res, std::string method = "QR")
{
    typedef typename Matrix<T>::difference_type Shape;
    typedef typename Matrix<T>::view_type SubMatrix;
    
    const MultiArrayIndex n = columnCount(A);
    const MultiArrayIndex m = rowCount(A);

    vigra_precondition(n <= m,
        "linearSolve(): Coefficient matrix A must have at least as many rows as columns.");
    vigra_precondition(n == rowCount(res) && 
                       m == rowCount(b) && columnCount(b) == columnCount(res),
        "linearSolve(): matrix shape mismatch.");

    method = tolower(method);
    if(method == "cholesky")
    {
        vigra_precondition(columnCount(A) == rowCount(A),
            "linearSolve(): Cholesky method requires square coefficient matrix.");
        Matrix<T> L(A.shape());
        if(!choleskyDecomposition(A, L))
            return false; // false if A wasn't symmetric positive definite
        choleskySolve(L, b, res);
    }
    else if(method == "qr")
    {
        return (MultiArrayIndex)linearSolveQR(A, b, res) == n;
    }
    else if(method == "ne")
    {
        return linearSolve(transpose(A)*A, transpose(A)*b, res, "Cholesky");
    }
    else if(method == "svd")
    {
        MultiArrayIndex rhsCount = columnCount(b);
        Matrix<T> u(A.shape()), s(n, 1), v(n, n);

        MultiArrayIndex rank = (MultiArrayIndex)singularValueDecomposition(A, u, s, v);

        Matrix<T> t = transpose(u)*b;
        for(MultiArrayIndex l=0; l<rhsCount; ++l)
        {
            for(MultiArrayIndex k=0; k<rank; ++k)
                t(k,l) /= s(k,0);
            for(MultiArrayIndex k=rank; k<n; ++k)
                t(k,l) = NumericTraits<T>::zero();
        }
        res = v*t;

        return (rank == n);
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
using linalg::logDeterminant;
using linalg::linearSolve;
using linalg::choleskySolve;
using linalg::choleskyDecomposition;
using linalg::qrDecomposition;
using linalg::linearSolveUpperTriangular;
using linalg::linearSolveLowerTriangular;

} // namespace vigra


#endif // VIGRA_LINEAR_SOLVE_HXX
