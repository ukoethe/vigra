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


#ifndef VIGRA_MATRIX_HXX
#define VIGRA_MATRIX_HXX

#include <cmath>
#include <iosfwd>
#include <iomanip>
#include "vigra/multi_array.hxx"
#include "vigra/mathutil.hxx"
#include "vigra/numerictraits.hxx"


namespace vigra
{

namespace linalg
{

template <class T, class C>
inline unsigned int rowCount(const MultiArrayView<2, T, C> &x);

template <class T, class C>
inline unsigned int columnCount(const MultiArrayView<2, T, C> &x);

template <class T, class C>
MultiArrayView <2, T, C>
rowVector(MultiArrayView <2, T, C> const & m, int d);

template <class T, class C>
MultiArrayView <2, T, C>
columnVector(MultiArrayView<2, T, C> const & m, int d);

template <class T, class ALLOC>
class TemporaryMatrix;

template <class T, class C1, class C2>
void transpose(const MultiArrayView<2, T, C1> &v, MultiArrayView<2, T, C2> &r);

template <class T, class C>
bool isSymmetric(const MultiArrayView<2, T, C> &v);

enum RawArrayMemoryLayout { RowMajor, ColumnMajor };

/********************************************************/
/*                                                      */
/*                        Matrix                        */
/*                                                      */
/********************************************************/

/** Matrix class.

    This is the basic class for all linear algebra computations. Matrices are
    strored in a <i>column-major</i> format, i.e. the row index is varying fastest.
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

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    typedef typename BaseType::SquaredNormType      SquaredNormType;
    typedef typename BaseType::NormType             NormType;
    
        /** default constructor
         */
    Matrix() 
    {}

        /** construct with given allocator
         */
    explicit Matrix(ALLOC const & alloc)
    : BaseType(alloc)
    {}

        /** construct with given shape and init all 
            elements with zero. Note that the order of the axes is
            <tt>difference_type(rows, columns)</tt> which
            is the opposite of the usual VIGRA convention.
         */
    explicit Matrix(const difference_type &shape,
                    ALLOC const & alloc = allocator_type())
    : BaseType(shape, alloc)
    {}

        /** construct with given shape and init all 
            elements with zero. Note that the order of the axes is
            <tt>(rows, columns)</tt> which
            is the opposite of the usual VIGRA convention.
         */
    Matrix(unsigned int rows, unsigned int columns,
                    ALLOC const & alloc = allocator_type())
    : BaseType(difference_type(rows, columns), alloc)
    {}

        /** construct with given shape and init all 
            elements with the constant \a init. Note that the order of the axes is
            <tt>difference_type(rows, columns)</tt> which
            is the opposite of the usual VIGRA convention.
         */
    Matrix(const difference_type &shape, const_reference init,
           allocator_type const & alloc = allocator_type())
    : BaseType(shape, init, alloc)
    {}

        /** construct with given shape and init all 
            elements with the constant \a init. Note that the order of the axes is
            <tt>(rows, columns)</tt> which
            is the opposite of the usual VIGRA convention.
         */
    Matrix(unsigned int rows, unsigned int columns, const_reference init,
           allocator_type const & alloc = allocator_type())
    : BaseType(difference_type(rows, columns), init, alloc)
    {}

        /** construct with given shape and copy data from C-style array \a init.
            Unless \a layout is <tt>ColumnMajor</tt>, the elements in this array 
            are assumed to be given in row-major order (the C standard order) and 
            will automatically be converted to the required column-major format. 
            Note that the order of the axes is <tt>difference_type(rows, columns)</tt> which
            is the opposite of the usual VIGRA convention.
         */
    Matrix(const difference_type &shape, const_pointer init, RawArrayMemoryLayout layout = RowMajor,
           allocator_type const & alloc = allocator_type())
    : BaseType(shape, alloc) // FIXME: this function initializes the memory twice
    {
        if(layout == RowMajor)
        {
            difference_type trans(shape[1], shape[0]);
            linalg::transpose(MultiArrayView<2, T>(trans, const_cast<pointer>(init)), *this);
        }
        else
        {
            std::copy(init, init + elementCount(), this->data());
        }
    }

        /** construct with given shape and copy data from C-style array \a init.
            Unless \a layout is <tt>ColumnMajor</tt>, the elements in this array 
            are assumed to be given in row-major order (the C standard order) and 
            will automatically be converted to the required column-major format. 
            Note that the order of the axes is <tt>(rows, columns)</tt> which
            is the opposite of the usual VIGRA convention.
         */
    Matrix(unsigned int rows, unsigned int columns, const_pointer init, RawArrayMemoryLayout layout = RowMajor,
           allocator_type const & alloc = allocator_type())
    : BaseType(difference_type(rows, columns), alloc) // FIXME: this function initializes the memory twice
    {
        if(layout == RowMajor)
        {
            difference_type trans(columns, rows);
            linalg::transpose(MultiArrayView<2, T>(trans, const_cast<pointer>(init)), *this);
        }
        else
        {
            std::copy(init, init + elementCount(), this->data());
        }
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
            copies tha data. \a rhs is assumed to be in column-major order already.
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
            \a rhs is assumed to be in column-major order already.
         */
    template <class U, class C>
    Matrix & operator=(const MultiArrayView<2, U, C> &rhs)
    {
        BaseType::operator=(rhs); // has the correct semantics already
        return *this;
    }
    
        /** Create a matrix view that represents the row vector of row \a d.
         */
    view_type rowVector(unsigned int d) const
    {
        return vigra::linalg::rowVector(*this, d);
    }
    
        /** Create a matrix view that represents the column vector of column \a d.
         */
    view_type columnVector(unsigned int d) const
    {
        return vigra::linalg::columnVector(*this, d);
    }
    
        /** number of rows (height) of the matrix.
        */
    unsigned int rowCount() const
    {
        return this->m_shape[0];
    }
    
        /** number of columns (width) of the matrix.
        */
    unsigned int columnCount() const
    {
        return this->m_shape[1];
    }
    
        /** number of elements (width*height) of the matrix.
        */
    unsigned int elementCount() const
    {
        return rowCount()*columnCount();
    }
       
        /** check whether the matrix is symmetric.
        */
    bool isSymmetric() const
    {
        return vigra::linalg::isSymmetric(*this);
    }
    
#ifdef DOXYGEN 
// repeat the index functions for documentation. In real code, they are inherited.

        /** read/write access to matrix element <tt>(row, column)</tt>.
            Note that the order of the argument is the opposite of the usual 
            VIGRA convention due to column-major matrix order.
        */
    value_type & operator()(unsigned int row, unsigned int column);

        /** read access to matrix element <tt>(row, column)</tt>.
            Note that the order of the argument is the opposite of the usual 
            VIGRA convention due to column-major matrix order.
        */
    value_type operator()(unsigned int row, unsigned int column) const;
#endif

        /** squared Frobenius norm. Sum of squares of the matrix elements.
        */
    SquaredNormType squaredNorm() const
    {
        return BaseType::squaredNorm();
    }
    
        /** Frobenius norm. Root of sum of squares of the matrix elements.
        */
    NormType norm() const
    {
        return BaseType::norm();
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

// TemporaryMatrix is provided as an optimization: Functions returning a matrix can 
// use TemporaryMatrix to make explicit that it was allocated as a temporary data structure.
// Functions receiving a TemporaryMatrix can thus often avoid to allocate new temporary 
// memory.
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

    TemporaryMatrix(unsigned int rows, unsigned int columns)
    : BaseType(rows, columns, ALLOC())
    {}

    TemporaryMatrix(unsigned int rows, unsigned int columns, const_reference init)
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C>
inline unsigned int rowCount(const MultiArrayView<2, T, C> &x)
{
    return x.shape(0);
}

    /** Number of columns of a matrix represented as a <tt>MultiArrayView&lt;2,...&gt;</tt>
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C>
inline unsigned int columnCount(const MultiArrayView<2, T, C> &x)
{
    return x.shape(1);
}

    /** Create a row vector view for row \a d of the matrix \a m
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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

    /** Check whether matrix \a m is symmetric.
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */ 
template <class T, class C>
bool
isSymmetric(MultiArrayView<2, T, C> const & m)
{
    const unsigned int size = rowCount(m);
    if(size != columnCount(m))
        return false;
        
    for(unsigned int i = 0; i < size; ++i)
        for(unsigned int j = i+1; j < size; ++j)
            if(m(j, i) != m(i, j))
                return false;
    return true;
}

#ifdef DOXYGEN // documentation only -- function is already defined in vigra/multi_array.hxx

    /** calculate the squared Frobenius norm of a matrix. 
        Equal to the sum of squares of the matrix elements.
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>"
        Namespace: vigra
     */ 
template <class T, class ALLOC>
typename Matrix<T, ALLLOC>::SquaredNormType
squaredNorm(const Matrix<T, ALLLOC> &a);

    /** calculate the squared Frobenius norm of a matrix. 
        Equal to the sum of squares of the matrix elements.
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>"
        Namespace: vigra
     */ 
template <class T, class ALLOC>
typename Matrix<T, ALLLOC>::NormType
norm(const Matrix<T, ALLLOC> &a);

#endif // DOXYGEN

    /** initialize the given square matrix as an identity matrix.
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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

    /** add matrices \a a and \a b.
        The result is written into \a r. All three matrices must have the same shape.
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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

    /** calculate the inner product of two matrices representing vectors. 
        That is, matrix \a x must have a single row, and matrix \a y must 
        have a single column, and the other dimensions must match.
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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

	/** multiply matrix \a a with scalar \a b.
        The result is written into \a r. \a a and \a r must have the same shape.
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2>
void smul(const MultiArrayView<2, T, C1> &a, T b, MultiArrayView<2, T, C2> &r)
{
    const unsigned int rows = rowCount(a);
    const unsigned int cols = columnCount(a);
    vigra_precondition(rows == rowCount(r) && cols == columnCount(r),
                       "smul(): Matrix sizes must agree.");
    
    for(unsigned int i = 0; i < cols; ++i)
        for(unsigned int j = 0; j < rows; ++j)
            r(j, i) = a(j, i) * b;
}

    /** multiply scalar \a a with matrix \a b.
        The result is written into \a r. \a b and \a r must have the same shape.
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C2, class C3>
void smul(T a, const MultiArrayView<2, T, C2> &b, MultiArrayView<2, T, C3> &r)
{
    smul(b, a, r);
}

    /** perform matrix multiplication of matrices \a a and \a b.
        The result is written into \a r. The three matrices must have matching shapes.
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2, class C3>
void mmul(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b,
         MultiArrayView<2, T, C3> &r)
{
    const unsigned int rrows = rowCount(r);
    const unsigned int rcols = columnCount(r);
    const unsigned int acols = columnCount(a);
    vigra_precondition(rrows == rowCount(a) && rcols == columnCount(b) && acols == rowCount(b),
                       "mmul(): Matrix shapes must agree.");

    for(unsigned int i = 0; i < rcols; ++i) {
        for(unsigned int j = 0; j < rrows; ++j) {
            r(j, i) = 0.0;
            for(unsigned int k = 0; k < acols; ++k) {
                r(j, i) += a(j, k) * b(k, i);
            }
        }
    }
}

    /** perform matrix multiplication of matrices \a a and \a b.
        \a a and \a b must have matching shapes.
        The result is returned as a temporary matrix. 
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2>
inline TemporaryMatrix<T>
mmul(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b)
{
    TemporaryMatrix<T> ret(rowCount(a), columnCount(b));
    mmul(a, b, ret);
    return ret;
}

    /** multiply two matrices \a a and \a b pointwise.
        The result is written into \a r. All three matrices must have the same shape.
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2, class C3>
void pmul(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b,
              MultiArrayView<2, T, C3> &r)
{
    const unsigned int rrows = rowCount(r);
    const unsigned int rcols = columnCount(r);
    vigra_precondition(rrows == rowCount(a) && rcols == columnCount(a) && 
                       rrows == rowCount(b) && rcols == columnCount(b),
                       "pmul(): Matrix shapes must agree.");

    for(unsigned int i = 0; i < rcols; ++i) {
        for(unsigned int j = 0; j < rrows; ++j) {
            r(j, i) = a(j, i) * b(j, i);
        }
    }
}

    /** multiply matrices \a a and \a b pointwise.
        \a a and \a b must have matching shapes.
        The result is returned as a temporary matrix. 
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2>
inline TemporaryMatrix<T>
pmul(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b)
{
    TemporaryMatrix<T> ret(rowCount(a), columnCount(b));
    pmul(a, b, ret);
    return ret;
}

    /** multiply matrix \a a with scalar \a b.
        The result is returned as a temporary matrix. 
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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

    /** multiply matrix \a a with TinyVector \a b.
        \a a must be of size <tt>N x N</tt>. Vector \a b and the result 
        vector are interpreted as column vectors.
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class A, int N, class DATA, class DERIVED>
TinyVector<T, N> 
operator*(const Matrix<T, A> &a, const TinyVectorBase<T, N, DATA, DERIVED> &b)
{
    vigra_precondition(N == rowCount(a) && N == columnCount(a),
         "operator*(Matrix, TinyVector): Shape mismatch.");

    TinyVector<T, N> res = TinyVectorView<T, N>(&a(0,0)) * b[0];
    for(unsigned int i = 1; i < N; ++i)
        res += TinyVectorView<T, N>(&a(0,i)) * b[i];
    return res;
}

    /** multiply TinyVector \a a with matrix \a b.
        \a b must be of size <tt>N x N</tt>. Vector \a a and the result 
        vector are interpreted as row vectors.
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, int N, class DATA, class DERIVED, class A>
TinyVector<T, N> 
operator*(const TinyVectorBase<T, N, DATA, DERIVED> &a, const Matrix<T, A> &b)
{
    vigra_precondition(N == rowCount(b) && N == columnCount(b),
         "operator*(TinyVector, Matrix): Shape mismatch.");

    TinyVector<T, N> res;
    for(unsigned int i = 0; i < N; ++i)
        res[i] = dot(a, TinyVectorView<T, N>(&b(0,i)));
    return res;
}

    /** perform matrix multiplication of matrices \a a and \a b.
        \a a and \a b must have matching shapes.
        The result is returned as a temporary matrix. 
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2>
inline TemporaryMatrix<T>
operator*(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b)
{
    TemporaryMatrix<T> ret(rowCount(a), columnCount(b));
    mmul(a, b, ret);
    return ret;
}

    /** divide matrix \a a by scalar \a b.
        The result is written into \a r. \a a and \a r must have the same shape.
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2>
void sdiv(const MultiArrayView<2, T, C1> &a, T b, MultiArrayView<2, T, C2> &r)
{
    const unsigned int rows = rowCount(a);
    const unsigned int cols = columnCount(a);
    vigra_precondition(rows == rowCount(r) && cols == columnCount(r),
                       "sdiv(): Matrix sizes must agree.");
    
    for(unsigned int i = 0; i < cols; ++i)
        for(unsigned int j = 0; j < rows; ++j)
            r(j, i) = a(j, i) / b;
}

    /** divide two matrices \a a and \a b pointwise.
        The result is written into \a r. All three matrices must have the same shape.
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2, class C3>
void pdiv(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b,
              MultiArrayView<2, T, C3> &r)
{
    const unsigned int rrows = rowCount(r);
    const unsigned int rcols = columnCount(r);
    vigra_precondition(rrows == rowCount(a) && rcols == columnCount(a) && 
                       rrows == rowCount(b) && rcols == columnCount(b),
                       "pdiv(): Matrix shapes must agree.");

    for(unsigned int i = 0; i < rcols; ++i) {
        for(unsigned int j = 0; j < rrows; ++j) {
            r(j, i) = a(j, i) * b(j, i);
        }
    }
}

    /** divide matrices \a a and \a b pointwise.
        \a a and \a b must have matching shapes.
        The result is returned as a temporary matrix. 
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */ 
template <class T, class C1, class C2>
inline TemporaryMatrix<T>
pdiv(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b)
{
    TemporaryMatrix<T> ret(rowCount(a), columnCount(b));
    pdiv(a, b, ret);
    return ret;
}

    /** divide matrix \a a by scalar \a b.
        The result is returned as a temporary matrix. 
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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

//@}

} // namespace linalg

using linalg::RowMajor;
using linalg::ColumnMajor;
using linalg::Matrix;
using linalg::identityMatrix;
using linalg::diagonalMatrix;
using linalg::transpose;
using linalg::dot;
using linalg::outer;
using linalg::rowCount;
using linalg::columnCount;
using linalg::rowVector;
using linalg::columnVector;
using linalg::isSymmetric;

/********************************************************/
/*                                                      */
/*                       NormTraits                     */
/*                                                      */
/********************************************************/

template <class T, class ALLOC>
struct NormTraits<linalg::Matrix<T, ALLOC> >
{
    typedef linalg::Matrix<T, ALLOC> Type;
    typedef typename Type::SquaredNormType SquaredNormType;
    typedef typename Type::NormType NormType;
};

template <class T, class ALLOC>
struct NormTraits<linalg::TemporaryMatrix<T, ALLOC> >
{
    typedef linalg::TemporaryMatrix<T, ALLOC> Type;
    typedef typename Type::SquaredNormType SquaredNormType;
    typedef typename Type::NormType NormType;
};

} // namespace vigra

namespace std {

/** \addtogroup LinearAlgebraFunctions Matrix functions
 */
//@{
 
    /** print a matrix \a m to the stream \a s. 
    
    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
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



#endif // VIGRA_MATRIX_HXX
