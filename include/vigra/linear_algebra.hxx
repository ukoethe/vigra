/************************************************************************/
/*                                                                      */
/*        Copyright 2004 by Gunnar Kedenburg and Ullrich Koethe         */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/


#ifndef VIGRA_LINEAR_ALGREBRA_HXX
#define VIGRA_LINEAR_ALGREBRA_HXX

#include <cmath>
#include <algorithm>
#include <iosfwd>
#include <iomanip>
#include "vigra/multi_array.hxx"
#include "vigra/mathutil.hxx"


namespace vigra
{

namespace linalg
{

template <class T, class C>
inline std::size_t rowCount(const MultiArrayView<2, T, C> &x);

template <class T, class C>
inline std::size_t columnCount(const MultiArrayView<2, T, C> &x);

template <class T, class C>
MultiArrayView <2, T, C>
rowVector(MultiArrayView <2, T, C> const & m, int d);

template <class T, class C>
MultiArrayView <2, T, C>
columnVector(MultiArrayView<2, T, C> const & m, int d);

template <class T, class C>
T squaredNorm(const MultiArrayView<2, T, C> &a);

template <class T, class ALLOC>
class TemporaryMatrix;

/********************************************************/
/*                                                      */
/*                        Matrix                        */
/*                                                      */
/********************************************************/

/** Matrix class.

    This is the basic class for all linear algebra computations. Matrices are
    strored in a <i>row-major</i> format, i.e. the row index is varying fastest.
    This is the same format as in the lapack and gmm++ libraries, so it will
    be easy to interface these libraries. In fact, if you need optimized
    high performance code, you should use them. The VIGRA linear algebra
    functionality is provided for smaller problems and rapid prototyping
    (no one wants to spend half the day installing a new library just to 
    discover that the new algorithm idea didn't work anyway).

    <b>See also:</b>
    <ul>
    <li> \ref LinearAlgebraFunctions
    </ul>

    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
*/
template <class T, class ALLOC = std::allocator<T> >
class Matrix
: public MultiArray<2, T, ALLOC>
{
    typedef MultiArray<2, T, ALLOC> BaseType;
    
  public:
    typedef Matrix<T, ALLOC>                        matrix_type;
    typedef TemporaryMatrix<T, ALLOC>               temp_type;
    typedef MultiArrayView<2, T, UnstridedArrayTag> view_type;
    typedef typename BaseType::value_type           value_type;
    typedef typename BaseType::pointer              pointer;
    typedef typename BaseType::const_pointer        const_pointer;
    typedef typename BaseType::reference            reference;
    typedef typename BaseType::const_reference      const_reference;
    typedef typename BaseType::difference_type      difference_type;
    typedef ALLOC                                   allocator_type;
    
        /** default constructor
         */
    Matrix() 
    {}

        /** construct with given allocator
         */
    explicit Matrix(ALLOC const & alloc)
    : BaseType(alloc)
    {}

        /** construct with given shape and init with all 
            elements with zero. Note that the order of the axes is
            <tt>difference_type(rows, columns)</tt> which
            is the opposite of the usual VIGRA convention.
         */
    explicit Matrix(const difference_type &shape,
                    ALLOC const & alloc = allocator_type())
    : BaseType(shape, alloc)
    {}

        /** construct with given shape and init with all 
            elements with zero. Note that the order of the axes is
            <tt>(rows, columns)</tt> which
            is the opposite of the usual VIGRA convention.
         */
    Matrix(std::size_t rows, std::size_t columns,
                    ALLOC const & alloc = allocator_type())
    : BaseType(difference_type(rows, columns), alloc)
    {}

        /** construct with given shape and init with all 
            elements with the constant \a init. Note that the order of the axes is
            <tt>difference_type(rows, columns)</tt> which
            is the opposite of the usual VIGRA convention.
         */
    Matrix(const difference_type &shape, const_reference init,
           allocator_type const & alloc = allocator_type())
    : BaseType(shape, init, alloc)
    {}

        /** construct with given shape and init with all 
            elements with the constant \a init. Note that the order of the axes is
            <tt>(rows, columns)</tt> which
            is the opposite of the usual VIGRA convention.
         */
    Matrix(std::size_t rows, std::size_t columns, const_reference init,
           allocator_type const & alloc = allocator_type())
    : BaseType(difference_type(rows, columns), init, alloc)
    {}

        /** construct with given shape and copy data from C-style array \a init.
            Data in this array are expected to be given in column-major
            order (the C standard order) and will automatically be
            converted to the required row-major format. Note that the order of the axes is
            <tt>difference_type(rows, columns)</tt> which
            is the opposite of the usual VIGRA convention.
         */
    Matrix(const difference_type &shape, const_pointer init,
           allocator_type const & alloc = allocator_type())
    : BaseType(shape, alloc)
    {
        difference_type trans(shape[1], shape[0]);
        linalg::transpose(MultiArrayView<2, T>(trans, const_cast<pointer>(init)), *this);
    }

        /** construct with given shape and copy data from C-style array \a init.
            Data in this array are expected to be given in column-major
            order (the C standard order) and will automatically be
            converted to the required row-major format. Note that the order of 
            the axes is <tt>(rows, columns)</tt> which
            is the opposite of the usual VIGRA convention.
         */
    Matrix(std::size_t rows, std::size_t columns, const_pointer init,
           allocator_type const & alloc = allocator_type())
    : BaseType(difference_type(rows, columns), alloc)
    {
        difference_type trans(columns, rows);
        linalg::transpose(MultiArrayView<2, T>(trans, const_cast<pointer>(init)), *this);
    }

        /** copy constructor. Allocates new memory and 
            copies tha data.
         */
    Matrix(const Matrix &rhs)
    : BaseType(rhs)
    {}

        /** construct from temporary matrix, which looses its data.
            
            This operation is equivalent to
            \code
                TemporaryMatrix<T> temp = ...;
                
                Matrix<T> m;
                m.swap(temp);
            \endcode
         */
    Matrix(const TemporaryMatrix<T, ALLOC> &rhs)
    : BaseType(rhs.allocator())
    {
        this->swap(const_cast<TemporaryMatrix<T, ALLOC> &>(rhs));
    }
    
        /** construct from a MultiArrayView. Allocates new memory and 
            copies tha data. \a rhs is assumed to be in row-major order already.
         */
    template<class U, class C>
    Matrix(const MultiArrayView<2, U, C> &rhs)
    : BaseType(rhs)
    {}

        /** assignment.
            If the size of \a rhs is the same as the matrix's old size, only the data
            are copied. Otherwise, new storage is allocated, which invalidates 
            all objects (array views, iterators) depending on the matrix.
         */
    Matrix & operator=(const Matrix &rhs)
    {
        BaseType::operator=(rhs); // has the correct semantics already
        return *this;
    }

        /** assign a temporary matrix. This is implemented by swapping the data
            between the two matrices, so that all depending objects 
            (array views, iterators) ar invalidated.
         */
    Matrix & operator=(const TemporaryMatrix<T, ALLOC> &rhs)
    {
        this->swap(const_cast<TemporaryMatrix<T, ALLOC> &>(rhs));
        return *this;
    }

        /** assignment from arbitrary 2-dimensional MultiArrayView.<br>
            If the size of \a rhs is the same as the matrix's old size, only the data
            are copied. Otherwise, new storage is allocated, which invalidates 
            all objects (array views, iterators) depending on the matrix. 
            \a rhs is assumed to be in row-major order already.
         */
    template <class U, class C>
    Matrix & operator=(const MultiArrayView<2, U, C> &rhs)
    {
        BaseType::operator=(rhs); // has the correct semantics already
        return *this;
    }
    
        /** Create a matrix view that represents the row vector of row \a d.
         */
    view_type rowVector(std::size_t d) const
    {
        return vigra::linalg::rowVector(*this, d);
    }
    
        /** Create a matrix view that represents the column vector of column \a d.
         */
    view_type columnVector(std::size_t d) const
    {
        return vigra::linalg::columnVector(*this, d);
    }
    
        /** number of rows (height) of the matrix.
        */
    std::size_t rowCount() const
    {
        return this->m_shape[0];
    }
    
        /** number of columns (width) of the matrix.
        */
    std::size_t columnCount() const
    {
        return this->m_shape[1];
    }
    
        /** number of elements (width*height) of the matrix.
        */
    std::size_t elementCount() const
    {
        return rowCount()*columnCount();
    }
    
#ifdef DOXYGEN 
// repeat the index functions for documentation. In real code, they are inherited.

        /** read/write access to matrix element <tt>(row, column)</tt>.
            Note that the order of the argument is the opposite of the usual 
            VIGRA convention due to row-major matrix order.
        */
    value_type & operator()(std::size_t row, std::size_t column);

        /** read access to matrix element <tt>(row, column)</tt>.
            Note that the order of the argument is the opposite of the usual 
            VIGRA convention due to row-major matrix order.
        */
    value_type operator()(std::size_t row, std::size_t column) const;
#endif

        /** squared Frobenius norm. Sum of squares of the matrix elements.
        */
    value_type squaredNorm() const
    {
        return vigra::linalg::squaredNorm(*this);
    }
    
        /** Frobenius norm. Root of sum of squares of the matrix elements.
        */
    value_type norm() const
    {
        return VIGRA_CSTD::sqrt(squaredNorm());
    }
    
        /** transpose matrix in-place (precondition: matrix must be square)
         */
    Matrix & transpose();
    
        /** add \a other to this (sizes must match).
         */
    template <class U, class C>
    Matrix & operator+=(MultiArrayView<2, U, C> const & other);
    
        /** subtract \a other from this (sizes must match).
         */
    template <class U, class C>
    Matrix & operator-=(MultiArrayView<2, U, C> const & other);
    
        /** scalar multiply this with \a other
         */
    Matrix & operator*=(T other);
    
        /** scalar devide this by \a other
         */
    Matrix & operator/=(T other);
};

template <class T, class ALLOC>
Matrix<T, ALLOC> & Matrix<T, ALLOC>::transpose()
{
    const unsigned int cols = columnCount();
    vigra_precondition(cols == rowCount(),
        "Matrix::transpose(): in-place transposition requires square matrix.");
    for(unsigned int i = 0; i < cols; ++i)
        for(unsigned int j = i+1; j < cols; ++j)
            std::swap((*this)(j, i), (*this)(i, j));
    return *this;
}

template <class T, class ALLOC>
template <class U, class C>
Matrix<T, ALLOC> & Matrix<T, ALLOC>::operator+=(MultiArrayView<2, U, C> const & other)
{
    const unsigned int rows = rowCount();
    const unsigned int cols = columnCount();
    vigra_precondition(rows == vigra::linalg::rowCount(other) && cols == vigra::linalg::columnCount(other),
       "Matrix::operator+=(): Shape mismatch.");
    
    for(unsigned int i = 0; i < cols; ++i)
        for(unsigned int j = 0; j < rows; ++j)
            (*this)(j, i) += other(j, i);
    return *this;
}

template <class T, class ALLOC>
template <class U, class C>
Matrix<T, ALLOC> & Matrix<T, ALLOC>::operator-=(MultiArrayView<2, U, C> const & other)
{
    const unsigned int rows = rowCount();
    const unsigned int cols = columnCount();
    vigra_precondition(rows == vigra::linalg::rowCount(other) && cols == vigra::linalg::columnCount(other),
       "Matrix::operator-=(): Shape mismatch.");
    
    for(unsigned int i = 0; i < cols; ++i)
        for(unsigned int j = 0; j < rows; ++j)
            (*this)(j, i) -= other(j, i);
    return *this;
}

template <class T, class ALLOC>
Matrix<T, ALLOC> & Matrix<T, ALLOC>::operator*=(T other)
{
    const unsigned int rows = rowCount();
    const unsigned int cols = columnCount();
    
    for(unsigned int i = 0; i < cols; ++i)
        for(unsigned int j = 0; j < rows; ++j)
            (*this)(j, i) *= other;
    return *this;
}

template <class T, class ALLOC>
Matrix<T, ALLOC> & Matrix<T, ALLOC>::operator/=(T other)
{
    const unsigned int rows = rowCount();
    const unsigned int cols = columnCount();
    
    for(unsigned int i = 0; i < cols; ++i)
        for(unsigned int j = 0; j < rows; ++j)
            (*this)(j, i) /= other;
    return *this;
}

template <class T, class ALLOC = std::allocator<T> >
class TemporaryMatrix
: public Matrix<T, ALLOC>
{
    typedef Matrix<T, ALLOC> BaseType;
  public:
    typedef Matrix<T, ALLOC>                        matrix_type;
    typedef TemporaryMatrix<T, ALLOC>               temp_type;
    typedef MultiArrayView<2, T, UnstridedArrayTag> view_type;
    typedef typename BaseType::value_type           value_type;
    typedef typename BaseType::pointer              pointer;
    typedef typename BaseType::const_pointer        const_pointer;
    typedef typename BaseType::reference            reference;
    typedef typename BaseType::const_reference      const_reference;
    typedef typename BaseType::difference_type      difference_type;
    typedef ALLOC                                   allocator_type;

    TemporaryMatrix(std::size_t rows, std::size_t columns)
    : BaseType(rows, columns, ALLOC())
    {}

    TemporaryMatrix(std::size_t rows, std::size_t columns, const_reference init)
    : BaseType(rows, columns, init, ALLOC())
    {}

    template<class U, class C>
    TemporaryMatrix(const MultiArrayView<2, U, C> &rhs)
    : BaseType(rhs)
    {}
    
    TemporaryMatrix(const TemporaryMatrix &rhs)
    : BaseType()
    {
        this->swap(const_cast<TemporaryMatrix &>(rhs));
    }
    
    TemporaryMatrix & transpose()
    {
        BaseType::transpose();
        return *this;
    }
    
    template <class U, class C>
    TemporaryMatrix & operator+=(MultiArrayView<2, U, C> const & other)
    {
        BaseType::operator+=(other);
        return *this;
    }
    
    template <class U, class C>
    TemporaryMatrix & operator-=(MultiArrayView<2, U, C> const & other)
    {
        BaseType::operator-=(other);
        return *this;
    }

    TemporaryMatrix & operator*=(T other)
    {
        BaseType::operator*=(other);
        return *this;
    }
    
    TemporaryMatrix & operator/=(T other)
    {
        BaseType::operator/=(other);
        return *this;
    }
  private:

    TemporaryMatrix &operator=(const TemporaryMatrix &rhs); // not implemented
};


/** \addtogroup LinearAlgebraFunctions Matrix functions
 */
//@{
 
    /** Number of rows of a matrix represented as a <tt>MultiArrayView&lt;2,...&gt;</tt>
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C>
inline std::size_t rowCount(const MultiArrayView<2, T, C> &x)
{
    return x.shape(0);
}

    /** Number of columns of a matrix represented as a <tt>MultiArrayView&lt;2,...&gt;</tt>
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C>
inline std::size_t columnCount(const MultiArrayView<2, T, C> &x)
{
    return x.shape(1);
}

    /** Create a row vector view for row \a d of the matrix \a m
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C>
MultiArrayView <2, T, C>
rowVector(MultiArrayView <2, T, C> const & m, int d)
{
    typedef typename MultiArrayView <2, T, C>::difference_type Shape;
    return m.subarray(Shape(d, 0), Shape(d+1, columnCount(m)));
}

    /** Create a column vector view for column \a d of the matrix \a m
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C>
MultiArrayView <2, T, C>
columnVector(MultiArrayView<2, T, C> const & m, int d)
{
    typedef typename MultiArrayView <2, T, C>::difference_type Shape;
    return m.subarray(Shape(0, d), Shape(rowCount(m), d+1));
}

    /** calculate the inner product of two matrices representing vectors. 
        That is, matrix \a x must have a single row, and matrix \a y must 
        have a single column, and the other dimensions must match.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C1, class C2>
T dot(const MultiArrayView<2, T, C1> &x, const MultiArrayView<2, T, C2> &y)
{
    const unsigned int n = columnCount(x);
    vigra_precondition(n == rowCount(y) && 1 == rowCount(x) && 1 == columnCount(y),
       "dot(): shape mismatch.");
    T ret = NumericTraits<T>::zero();
    for(unsigned int i = 0; i < n; ++i)
        ret += x(0, i) * y(i, 0);
    return ret;
}

    /** calculate the inner product of two vectors. The vector
        lenths must match.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C1, class C2>
T dot(const MultiArrayView<1, T, C1> &x, const MultiArrayView<1, T, C2> &y)
{
    const unsigned int n = x.elementCount();
    vigra_precondition(n == y.elementCount(),
       "dot(): shape mismatch.");
    T ret = NumericTraits<T>::zero();
    for(unsigned int i = 0; i < n; ++i)
        ret += x(i) * y(i);
    return ret;
}

    /** calculate the outer product of two matrices representing vectors. 
        That is, matrix \a x must have a single column, and matrix \a y must 
        have a single row, and the other dimensions must match. The result
        is written into \a r.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C1, class C2, class C3>
void outer(const MultiArrayView<2, T, C1> &x, const MultiArrayView<2, T, C2> &y,
      MultiArrayView<2, T, C3> &r)
{
    const unsigned int rows = rowCount(r);
    const unsigned int cols = columnCount(r);
    vigra_precondition(rows == rowCount(x) && cols == columnCount(y) && 
                       1 == columnCount(x) && 1 == rowCount(y),
       "outer(): shape mismatch.");
    for(unsigned int i = 0; i < cols; ++i)
        for(unsigned int j = 0; j < rows; ++j)
            r(j, i) = x(j, 0) * y(0, i);
}

    /** calculate the outer product of two matrices representing vectors. 
        That is, matrix \a x must have a single column, and matrix \a y must 
        have a single row, and the other dimensions must match. The result
        is returned as a temporary matrix.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C1, class C2>
TemporaryMatrix<T> 
outer(const MultiArrayView<2, T, C1> &x, const MultiArrayView<2, T, C2> &y)
{
    const unsigned int rows = rowCount(x);
    const unsigned int cols = columnCount(y);
    vigra_precondition(1 == columnCount(x) && 1 == rowCount(y),
       "outer(): shape mismatch.");
    TemporaryMatrix<T> ret(rows, cols);
    outer(x, y, ret);
    return ret;
}

    /** calculate the squared Frobenius norm of a matrix. 
        Equal to the sum of squares of the matrix elements.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C>
T squaredNorm(const MultiArrayView<2, T, C> &a)
{
    const unsigned int rows = rowCount(a);
    const unsigned int cols = columnCount(a);
    T ret = NumericTraits<T>::zero();
    for(unsigned int j = 0; j < cols; ++j)
        for(unsigned int i = 0; i < rows; ++i)
            ret += sq(a(i, j));
    return ret;
}

    /** calculate the squared norm of a vector. 
        Equal to the sum of squares of the vector elements.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C>
T squaredNorm(const MultiArrayView<1, T, C> &a)
{
    const unsigned int size = a.elementCount();
    T ret = NumericTraits<T>::zero();
    for(unsigned int i = 0; i < size; ++i)
        ret += sq(a(i));
    return ret;
}

    /** calculate the Frobenius norm of a matrix or vector. 
        Equal to the square root of sum of squares of the matrix elements.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <unsigned int N, class T, class C>
T norm(const MultiArrayView<N, T, C> &a)
{
    return VIGRA_CSTD::sqrt(squaredNorm(a));
}

    /** initialize the given square matrix as an identity matrix.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C>
void identityMatrix(MultiArrayView<2, T, C> &r)
{
    const unsigned int rows = rowCount(r);
    vigra_precondition(rows == columnCount(r),
       "identityMatrix(): Matrix must be square.");
    for(unsigned int i = 0; i < rows; ++i) {
        for(unsigned int j = 0; j < rows; ++j)
            r(j, i) = NumericTraits<T>::zero();
        r(i, i) = NumericTraits<T>::one();
    }
}

    /** create n identity matrix of the given size.
        Usage:
        
        \code
        vigra::Matrix<double> m = vigra::identityMatrix<double>(size);
        \endcode
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T>
TemporaryMatrix<T> identityMatrix(unsigned int size)
{
    TemporaryMatrix<T> ret(size, size, NumericTraits<T>::zero());
    for(unsigned int i = 0; i < size; ++i)
        ret(i, i) = NumericTraits<T>::one();
    return ret;
}

template <class T, class C1, class C2>
void diagonalMatrixImpl(MultiArrayView<1, T, C1> const & v, MultiArrayView<2, T, C2> &r)
{
    const unsigned int size = v.elementCount();
    vigra_precondition(rowCount(r) == size && columnCount(r) == size,
        "diagonalMatrix(): result must be a square matrix.");
    for(unsigned int i = 0; i < size; ++i)
        r(i, i) = v(i);
}

    /** make a diagonal matrix from a vector.
        The vector is given as matrix \a v, which must either have a single
        row or column. The result is witten into the square matrix \a r.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C1, class C2>
void diagonalMatrix(MultiArrayView<2, T, C1> const & v, MultiArrayView<2, T, C2> &r)
{
    vigra_precondition(rowCount(v) == 1 || columnCount(v) == 1,
        "diagonalMatrix(): input must be a vector.");
    r.init(NumericTraits<T>::zero());
    if(rowCount(v) == 1)
        diagonalMatrixImpl(v.bindInner(0), r);
    else
        diagonalMatrixImpl(v.bindOuter(0), r);
}

    /** create a diagonal matrix from a vector.
        The vector is given as matrix \a v, which must either have a single
        row or column. The result is returned as a temporary matrix.
        Usage:
        
        \code
        vigra::Matrix<double> v(1, len);
        v = ...;
        
        vigra::Matrix<double> m = diagonalMatrix(v);
        \endcode
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C>
TemporaryMatrix<T> diagonalMatrix(MultiArrayView<2, T, C> const & v)
{
    vigra_precondition(rowCount(v) == 1 || columnCount(v) == 1,
        "diagonalMatrix(): input must be a vector.");
    unsigned int size = v.elementCount();
    TemporaryMatrix<T> ret(size, size, NumericTraits<T>::zero());
    if(rowCount(v) == 1)
        diagonalMatrixImpl(v.bindInner(0), ret);
    else
        diagonalMatrixImpl(v.bindOuter(0), ret);
    return ret;
}

    /** transpose matrix \a v.
        The result is written into \a r which must have the correct (i.e.
        transposed) shape.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C1, class C2>
void transpose(const MultiArrayView<2, T, C1> &v, MultiArrayView<2, T, C2> &r)
{
    const unsigned int rows = rowCount(r);
    const unsigned int cols = columnCount(r);
    vigra_precondition(rows == columnCount(v) && cols == rowCount(v),
       "transpose(): arrays must have transposed shapes.");
    for(unsigned int i = 0; i < cols; ++i)
        for(unsigned int j = 0; j < rows; ++j)
            r(j, i) = v(i, j);
}

    /** create the transpose of a matrix \a v.
        The result is returned as a temporary matrix.
        Usage:
        
        \code
        vigra::Matrix<double> v(rows, cols);
        v = ...;
        
        vigra::Matrix<double> m = transpose(v);
        \endcode
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C>
TemporaryMatrix<T> transpose(MultiArrayView<2, T, C> const & v)
{
    TemporaryMatrix<T> ret(columnCount(v), rowCount(v));
    transpose(v, ret);
    return ret;
}

template <class T>
TemporaryMatrix<T> transpose(TemporaryMatrix<T> const & v)
{
    const unsigned int rows = v.rowCount();
    const unsigned int cols = v.columnCount();
    if(rows == cols)
    {
        return const_cast<TemporaryMatrix<T> &>(v).transpose();
    }
    else
    {
        TemporaryMatrix<T> ret(cols, rows);
        transpose(v, ret);
        return ret;
    }
}

    /** invert square matrix \a v.
        The result is written into \a r which must have the same shape.
        The inverse is calculated by means of QR decomposition. If \a v
        is not invertible, <tt>vigra::PreconditionViolation</tt> exception is thrown.
    
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
 
    /** add matrices \a a and \a b.
        The result is written into \a r. All three matrices must have the same shape.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2, class C3>
void add(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b,
              MultiArrayView<2, T, C3> &r)
{
    const unsigned int rrows = rowCount(r);
    const unsigned int rcols = columnCount(r);
    vigra_precondition(rrows == rowCount(a) && rcols == columnCount(a) && 
                       rrows == rowCount(b) && rcols == columnCount(b),
                       "add(): Matrix shapes must agree.");

    for(unsigned int i = 0; i < rcols; ++i) {
        for(unsigned int j = 0; j < rrows; ++j) {
            r(j, i) = a(j, i) + b(j, i);
        }
    }
}

 
    /** add matrices \a a and \a b.
        The two matrices must have the same shape.
        The result is returned as a temporary matrix. 
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2>
inline TemporaryMatrix<T>
operator+(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b)
{
    return TemporaryMatrix<T>(a) += b;
}

template <class T, class C>
inline TemporaryMatrix<T>
operator+(const TemporaryMatrix<T> &a, const MultiArrayView<2, T, C> &b)
{
    return const_cast<TemporaryMatrix<T> &>(a) += b;
}

template <class T, class C>
inline TemporaryMatrix<T>
operator+(const MultiArrayView<2, T, C> &a, const TemporaryMatrix<T> &b)
{
    return const_cast<TemporaryMatrix<T> &>(b) += a;
}

template <class T>
inline TemporaryMatrix<T>
operator+(const TemporaryMatrix<T> &a, const TemporaryMatrix<T> &b)
{
    return const_cast<TemporaryMatrix<T> &>(a) += b;
}

    /** subtract matrix \a b from \a a.
        The result is written into \a r. All three matrices must have the same shape.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2, class C3>
void sub(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b,
              MultiArrayView<2, T, C3> &r)
{
    const unsigned int rrows = rowCount(r);
    const unsigned int rcols = columnCount(r);
    vigra_precondition(rrows == rowCount(a) && rcols == columnCount(a) && 
                       rrows == rowCount(b) && rcols == columnCount(b),
                       "subtract(): Matrix shapes must agree.");

    for(unsigned int i = 0; i < rcols; ++i) {
        for(unsigned int j = 0; j < rrows; ++j) {
            r(j, i) = a(j, i) - b(j, i);
        }
    }
}
  
    /** subtract matrix \a b from \a a.
        The two matrices must have the same shape.
        The result is returned as a temporary matrix. 
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2>
inline TemporaryMatrix<T>
operator-(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b)
{
    return TemporaryMatrix<T>(a) -= b;
}

template <class T, class C>
inline TemporaryMatrix<T>
operator-(const TemporaryMatrix<T> &a, const MultiArrayView<2, T, C> &b)
{
    return const_cast<TemporaryMatrix<T> &>(a) -= b;
}

template <class T, class C>
TemporaryMatrix<T>
operator-(const MultiArrayView<2, T, C> &a, const TemporaryMatrix<T> &b)
{
    const unsigned int rows = rowCount(a);
    const unsigned int cols = columnCount(a);
    vigra_precondition(rows == b.rowCount() && cols == b.columnCount(),
       "Matrix::operator-(): Shape mismatch.");
    
    for(unsigned int i = 0; i < cols; ++i)
        for(unsigned int j = 0; j < rows; ++j)
            const_cast<TemporaryMatrix<T> &>(b)(j, i) = a(j, i) - b(j, i);
    return b;
}

template <class T>
inline TemporaryMatrix<T>
operator-(const TemporaryMatrix<T> &a, const TemporaryMatrix<T> &b)
{
    return const_cast<TemporaryMatrix<T> &>(a) -= b;
}

    /** negate matrix \a a.
        The result is returned as a temporary matrix. 
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C>
inline TemporaryMatrix<T>
operator-(const MultiArrayView<2, T, C> &a)
{
    return TemporaryMatrix<T>(a) *= -NumericTraits<T>::one();
}

template <class T>
inline TemporaryMatrix<T>
operator-(const TemporaryMatrix<T> &a)
{
    return const_cast<TemporaryMatrix<T> &>(a) *= -NumericTraits<T>::one();
}

    /** multiply matrix \a a with scalar \a b.
        The result is written into \a r. \a a and \a r must have the same shape.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2>
void mul(const MultiArrayView<2, T, C1> &a, T b, MultiArrayView<2, T, C2> &r)
{
    const unsigned int rows = rowCount(a);
    const unsigned int cols = columnCount(a);
    vigra_precondition(rows == rowCount(r) && cols == columnCount(r),
                       "mul(): Matrix sizes must agree.");
    
    for(unsigned int i = 0; i < cols; ++i)
        for(unsigned int j = 0; j < rows; ++j)
            r(j, i) = a(j, i) * b;
}

    /** multiply scalar \a a with matrix \a b.
        The result is written into \a r. \a b and \a r must have the same shape.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C2, class C3>
void mul(T a, const MultiArrayView<2, T, C2> &b, MultiArrayView<N, T, C3> &r)
{
    mul(b, a, r);
}

    /** multiply matrices \a a and \a b.
        The result is written into \a r. The three matrices must have matching shapes.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2, class C3>
void mul(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b,
         MultiArrayView<2, T, C3> &r)
{
    const unsigned int rrows = rowCount(r);
    const unsigned int rcols = columnCount(r);
    const unsigned int acols = columnCount(a);
    vigra_precondition(rrows == rowCount(a) && rcols == columnCount(b) && acols == rowCount(b),
                       "mul(): Matrix shapes must agree.");

    for(unsigned int i = 0; i < rcols; ++i) {
        for(unsigned int j = 0; j < rrows; ++j) {
            r(j, i) = 0.0;
            for(unsigned int k = 0; k < acols; ++k) {
                r(j, i) += a(j, k) * b(k, i);
            }
        }
    }
}

    /** multiply matrices \a a and \a b.
        \a a and \a b must have matching shapes.
        The result is returned as a temporary matrix. 
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2>
inline TemporaryMatrix<T>
mul(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b)
{
    TemporaryMatrix<T> ret(rowCount(a), columnCount(b));
    mul(a, b, ret);
    return ret;
}

    /** multiply matrix \a a with scalar \a b.
        The result is returned as a temporary matrix. 
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C>
inline TemporaryMatrix<T>
operator*(const MultiArrayView<2, T, C> &a, T b)
{
    return TemporaryMatrix<T>(a) *= b;
}

template <class T>
inline TemporaryMatrix<T>
operator*(const TemporaryMatrix<T> &a, T b)
{
    return const_cast<TemporaryMatrix<T> &>(a) *= b;
}

    /** multiply scalar \a a with matrix \a b.
        The result is returned as a temporary matrix. 
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C>
inline TemporaryMatrix<T>
operator*(T a, const MultiArrayView<2, T, C> &b)
{
    return TemporaryMatrix<T>(b) *= a;
}

template <class T>
inline TemporaryMatrix<T>
operator*(T a, const TemporaryMatrix<T> &b)
{
    return const_cast<TemporaryMatrix<T> &>(b) *= b;
}

    /** multiply matrices \a a and \a b.
        \a a and \a b must have matching shapes.
        The result is returned as a temporary matrix. 
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2>
inline TemporaryMatrix<T>
operator*(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b)
{
    TemporaryMatrix<T> ret(rowCount(a), columnCount(b));
    mul(a, b, ret);
    return ret;
}

    /** divide matrix \a a by scalar \a b.
        The result is written into \a r. \a a and \a r must have the same shape.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2>
void div(const MultiArrayView<2, T, C1> &a, T b, MultiArrayView<2, T, C2> &r)
{
    const unsigned int rows = rowCount(a);
    const unsigned int cols = columnCount(a);
    vigra_precondition(rows == rowCount(r) && cols == columnCount(r),
                       "div(): Matrix sizes must agree.");
    
    for(unsigned int i = 0; i < cols; ++i)
        for(unsigned int j = 0; j < rows; ++j)
            r(j, i) = a(j, i) / b;
}


    /** divide matrix \a a by scalar \a b.
        The result is returned as a temporary matrix. 
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C>
inline TemporaryMatrix<T>
operator/(const MultiArrayView<2, T, C> &a, T b)
{
    return TemporaryMatrix<T>(a) /= b;
}

template <class T>
inline TemporaryMatrix<T>
operator/(const TemporaryMatrix<T> &a, T b)
{
    return const_cast<TemporaryMatrix<T> &>(a) /= b;
}

    /** QR decomposition.
     
        \a a contains the original matrix, results are returned in \a q and \a r, where
        \a q is a orthogonal matrix, and \a r is an uppr triangular matrix, and
        the following relation holds:
        
        \code
        assert(a == q * r);
        \endcode
        
        This implementation uses householder transformations. It can be applied in-place,
        i.e. <tt>&a == &r</tt> is allowed.
    
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

//@}

namespace detail {

// solves a linear equation system with the
//  defining matrix being right triangular, as obtained by qrDecomposition().
template <class T, class C1, class C2, class C3>
void reverseEliminationQR(const MultiArrayView<2, T, C1> &r, const MultiArrayView<2, T, C2> &b,
                          MultiArrayView<2, T, C3> & x)
{
    const unsigned int m = columnCount(r);
    x(m-1, 0) = b(m-1, 0) / r(m-1, m-1);
    if(m >= 2) {
        for(int i = m-2; i >= 0; --i) {
            // compute the i'th inner product, excluding the diagonal entry.
            T sum = 0.0;
            for(unsigned int j = i+1; j < m; ++j)
                sum += r(i, j) * x(j, 0);
            if(r(i, i) != NumericTraits<T>::zero())
                x(i, 0) = (b(i, 0) - sum) / r(i, i);
            else
                x(i, 0) = 0;
        }
    }
}

// code adapted from JAMA
// a and b will be overwritten
template <class T, class C1, class C2>
void 
housholderTridiagonalMatrixization(MultiArrayView<2, T, C1> &a, MultiArrayView<2, T, C2> &b)
{
    const unsigned int n = rowCount(a);
    vigra_precondition(n == columnCount(a),
        "housholderTridiagonalMatrixization(): matrix must be square.");
    vigra_precondition(n == rowCount(b) && 2 <= columnCount(b),
        "housholderTridiagonalMatrixization(): matrix size mismatch.");

    MultiArrayView<1, T, C2> d = b.bindOuter(0);
    MultiArrayView<1, T, C2> e = b.bindOuter(1);
    
    for(int j = 0; j < n; ++j) 
    {
        d(j) = a(n-1, j);
    }

    // Householder reduction to tridiagonalMatrix form.
 
    for(int i = n-1; i > 0; --i) 
    {
        // Scale to avoid under/overflow.
 
        T scale = 0.0;
        T h = 0.0;
        for(int k = 0; k < i; ++k)
        {
            scale = scale + abs(d(k));
        }
        if(scale == 0.0)
        {
            e(i) = d(i-1);
            for(int j = 0; j < i; ++j)
            {
                d(j) = a(i-1, j);
                a(i, j) = 0.0;
                a(j, i) = 0.0;
            }
        }
        else
        {
            // Generate Householder vector.
 
            for(int k = 0; k < i; ++k)
            {
                d(k) /= scale;
                h += sq(d(k));
            }
            T f = d(i-1);
            T g = VIGRA_CSTD::sqrt(h);
            if(f > 0) {
                g = -g;
            }
            e(i) = scale * g;
            h -= f * g;
            d(i-1) = f - g;
            for(int j = 0; j < i; ++j)
            {
               e(j) = 0.0;
            }
 
            // Apply similarity transformation to remaining columns.
 
            for(int j = 0; j < i; ++j)
            {
               f = d(j);
               a(j, i) = f;
               g = e(j) + a(j, j) * f;
               for(int k = j+1; k <= i-1; ++k)
               {
                   g += a(k, j) * d(k);
                   e(k) += a(k, j) * f;
               }
               e(j) = g;
            }
            f = 0.0;
            for(int j = 0; j < i; ++j)
            {
                e(j) /= h;
                f += e(j) * d(j);
            }
            T hh = f / (h + h);
            for(int j = 0; j < i; ++j)
            {
                e(j) -= hh * d(j);
            }
            for(int j = 0; j < i; ++j)
            {
                f = d(j);
                g = e(j);
                for(int k = j; k <= i-1; ++k)
                {
                    a(k, j) -= (f * e(k) + g * d(k));
                }
                d(j) = a(i-1, j);
                a(i, j) = 0.0;
            }
        }
        d(i) = h;
    }
 
    // Accumulate transformations.
 
    for(int i = 0; i < n-1; ++i) 
    {
        a(n-1, i) = a(i, i);
        a(i, i) = 1.0;
        T h = d(i+1);
        if(h != 0.0) 
        {
            for(int k = 0; k <= i; ++k) 
            {
                d(k) = a(k, i+1) / h;
            }
            for(int j = 0; j <= i; ++j) 
            {
                T g = 0.0;
                for(int k = 0; k <= i; ++k) 
                {
                    g += a(k, i+1) * a(k, j);
                }
                for(int k = 0; k <= i; ++k) 
                {
                    a(k, j) -= g * d(k);
                }
            }
        }
        for(int k = 0; k <= i; ++k) 
        {
            a(k, i+1) = 0.0;
        }
    }
    for(int j = 0; j < n; ++j) 
    {
        d(j) = a(n-1, j);
        a(n-1, j) = 0.0;
    }
    a(n-1, n-1) = 1.0;
    e(0) = 0.0;
}

// code adapted from JAMA
// de and z will be overwritten
template <class T, class C1, class C2>
bool 
tridiagonalMatrixEigensystem(MultiArrayView<2, T, C1> &de, MultiArrayView<2, T, C2> &z) 
{ 
    const unsigned int n = rowCount(z);
    vigra_precondition(n == columnCount(z),
        "tridiagonalMatrixEigensystem(): matrix must be square.");
    vigra_precondition(n == rowCount(de) && 2 <= columnCount(de),
        "tridiagonalMatrixEigensystem(): matrix size mismatch.");

    MultiArrayView<1, T, C2> d = de.bindOuter(0);
    MultiArrayView<1, T, C2> e = de.bindOuter(1);
    
    for(int i = 1; i < n; i++) {
       e(i-1) = e(i);
    }
    e(n-1) = 0.0;
 
    T f = 0.0;
    T tst1 = 0.0;
    T eps = VIGRA_CSTD::pow(2.0,-52.0);
    for(int l = 0; l < n; ++l) 
    {
        // Find small subdiagonalMatrix element
 
        tst1 = std::max(tst1, abs(d(l)) + abs(e(l)));
        int m = l;

        // Original while-loop from Java code
        while(m < n) 
        {
            if(abs(e(m)) <= eps*tst1)
                break;
            ++m;
        }
 
        // If m == l, d(l) is an eigenvalue,
        // otherwise, iterate.
 
        if(m > l) 
        {
            int iter = 0;
            do
            {
                if(++iter > 50)
                   return false; // too many iterations
 
                // Compute implicit shift
 
                T g = d(l);
                T p = (d(l+1) - g) / (2.0 * e(l));
                T r = hypot(p,1.0);
                if(p < 0)
                {
                    r = -r;
                }
                d(l) = e(l) / (p + r);
                d(l+1) = e(l) * (p + r);
                T dl1 = d(l+1);
                T h = g - d(l);
                for(int i = l+2; i < n; ++i)
                {
                   d(i) -= h;
                }
                f = f + h;
 
                // Implicit QL transformation.
 
                p = d(m);
                T c = 1.0;
                T c2 = c;
                T c3 = c;
                T el1 = e(l+1);
                T s = 0.0;
                T s2 = 0.0;
                for(int i = m-1; i >= l; --i)
                {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e(i);
                    h = c * p;
                    r = hypot(p,e(i));
                    e(i+1) = s * r;
                    s = e(i) / r;
                    c = p / r;
                    p = c * d(i) - s * g;
                    d(i+1) = h + s * (c * g + s * d(i));
  
                    // Accumulate transformation.
  
                    for(int k = 0; k < n; ++k)
                    {
                         h = z(k, i+1);
                         z(k, i+1) = s * z(k, i) + c * h;
                        z(k, i) = c * z(k, i) - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e(l) / dl1;
                e(l) = s * p;
                d(l) = c * p;
 
                // Check for convergence.
 
            } while(abs(e(l)) > eps*tst1);
        }
        d(l) = d(l) + f;
        e(l) = 0.0;
    }
   
    // Sort eigenvalues and corresponding vectors.
 
    for(int i = 0; i < n-1; ++i) 
    {
        int k = i;
        T p = abs(d(i));
        for(int j = i+1; j < n; ++j)
        {
            T p1 = abs(d(j));
            if(p < p1)
            {
                k = j;
                p = p1;
            }
        }
        if(k != i)
        {
            std::swap(d(k), d(i));
            for(int j = 0; j < n; ++j)
            {
                std::swap(z(j, i), z(j, k));
            }
        }
    }
    return true;
}

} // namespace detail

/** \addtogroup LinearAlgebraFunctions Matrix functions
 */
//@{
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
    if(r(acols-1, acols-1) == NumericTraits<T>::zero())
         return false; // a didn't have full rank.
    q.transpose();
    Matrix<T> qb = q * b;
    for(unsigned int i = 0; i < bcols; ++i)
    {
        MultiArrayView<2, T, C3> resi = columnVector(res, i);
        detail::reverseEliminationQR(r, columnVector(qb, i), resi);
    }
    return true;
}

    /** Compute the eigensystem of a symmetric matrix.
     
        \a a is a real symmetric matrix, \a ew is a single-column matrix
        holding the eigenvalues, and \a ev is a matrix of the same size as 
        \a a whose columns are the corresponding eigenvectors. Eigenvalues 
        will be sorted from largest to smallest magnitude. 
        The algorithm returns <tt>false</tt> when it doesn't 
        converge. It can be applied in-place, i.e. <tt>&a == &ev</tt> is allowed.
        The code of this function was adapted from JAMA.
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
bool 
eigenSystem(MultiArrayView<2, T, C1> const & a, 
            MultiArrayView<2, T, C2> & ew, MultiArrayView<2, T, C3> & ev) 
{
    unsigned int acols = columnCount(a);
    vigra_precondition(acols == rowCount(a),
        "eigenSystem(): square input matrix required.");
    vigra_precondition(1 == columnCount(ew) && acols == rowCount(ew) && 
                       acols == columnCount(ev) && acols == rowCount(ev),
        "eigenSystem(): matrix shape mismatch.");
    
    ev.copy(a); // does nothing if &e == &a
    Matrix<T> de(acols, 2);
    detail::housholderTridiagonalMatrixization(ev, de);
    if(!detail::tridiagonalMatrixEigensystem(de, ev))
        return false;
    
    ew.copy(columnVector(de, 0));
    return true;
}

//@}

} // namespace linalg

using linalg::Matrix;
using linalg::identityMatrix;
using linalg::diagonalMatrix;
using linalg::squaredNorm;
using linalg::norm;
using linalg::transpose;
using linalg::inverse;
using linalg::dot;
using linalg::outer;
using linalg::rowCount;
using linalg::columnCount;
using linalg::rowVector;
using linalg::columnVector;
using linalg::linearSolve;
using linalg::qrDecomposition;
using linalg::eigenSystem;

} // namespace vigra

namespace std {

/** \addtogroup LinearAlgebraFunctions Matrix functions
 */
//@{
 
    /** print a matrix \a m to the stream \a s. 
    
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: std
     */ 
template <class T, class C>
ostream &
operator<<(ostream & s, const vigra::MultiArrayView<2, T, C> &m)
{
    const unsigned int rows = vigra::linalg::rowCount(m);
    const unsigned int cols = vigra::linalg::columnCount(m);
    ios::fmtflags flags = 
        s.setf(ios::right | ios::fixed, ios::adjustfield | ios::floatfield);
    for(unsigned int j = 0; j < rows; ++j) 
    {
        for(unsigned int i = 0; i < cols; ++i)
        {
            s << setw(7) << setprecision(4) << m(j, i) << " ";
        }
        s << endl;
    }
    s.setf(flags);
    return s;
}

//@}

} // namespace std



#endif // VIGRA_LINEAR_ALGREBRA_HXX
