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


#ifndef VIGRA_MATRIX_HXX
#define VIGRA_MATRIX_HXX

#include <cmath>
#include <iosfwd>
#include <iomanip>
#include "multi_array.hxx"
#include "mathutil.hxx"
#include "numerictraits.hxx"


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

template <class T, class ALLOC>
class TemporaryMatrix;

template <class T, class C1, class C2>
void transpose(const MultiArrayView<2, T, C1> &v, MultiArrayView<2, T, C2> &r);

template <class T, class C>
bool isSymmetric(const MultiArrayView<2, T, C> &v);


// implement operator+= etc. for const lhs so that it can be applied to temporary matrices
template <class T, class C1, class U, class C2>
MultiArrayView<2, T, C1> & 
operator+=(MultiArrayView<2, T, C1> const & l, MultiArrayView<2, U, C2> const & other)
{
    MultiArrayView<2, T, C1> & left = const_cast<MultiArrayView<2, T, C1> &>(l);
    const std::size_t rows = rowCount(left);
    const std::size_t cols = columnCount(left);
    vigra_precondition(rows == rowCount(other) && cols == columnCount(other),
       "operator+=(MultiArrayView, MultiArrayView): Shape mismatch.");

    for(std::size_t i = 0; i < cols; ++i)
        for(std::size_t j = 0; j < rows; ++j)
            left(j, i) += other(j, i);
    return left;
}

template <class T, class C1, class U, class C2>
MultiArrayView<2, T, C1> & operator-=(MultiArrayView<2, T, C1> const & l, MultiArrayView<2, U, C2> const & other)
{
    MultiArrayView<2, T, C1> & left = const_cast<MultiArrayView<2, T, C1> &>(l);
    const std::size_t rows = rowCount(left);
    const std::size_t cols = columnCount(left);
    vigra_precondition(rows == rowCount(other) && cols == columnCount(other),
       "operator-=(MultiArrayView, MultiArrayView): Shape mismatch.");

    for(std::size_t i = 0; i < cols; ++i)
        for(std::size_t j = 0; j < rows; ++j)
            left(j, i) -= other(j, i);
    return left;
}

template <class T, class C>
MultiArrayView<2, T, C> & operator+=(MultiArrayView<2, T, C> const & l, T other)
{
    MultiArrayView<2, T, C> & left = const_cast<MultiArrayView<2, T, C> &>(l);
    const std::size_t rows = rowCount(left);
    const std::size_t cols = columnCount(left);

    for(std::size_t i = 0; i < cols; ++i)
        for(std::size_t j = 0; j < rows; ++j)
            left(j, i) += other;
    return left;
}

template <class T, class C>
MultiArrayView<2, T, C> & operator-=(MultiArrayView<2, T, C> const & l, T other)
{
    MultiArrayView<2, T, C> & left = const_cast<MultiArrayView<2, T, C> &>(l);
    const std::size_t rows = rowCount(left);
    const std::size_t cols = columnCount(left);

    for(std::size_t i = 0; i < cols; ++i)
        for(std::size_t j = 0; j < rows; ++j)
            left(j, i) -= other;
    return left;
}

template <class T, class C>
MultiArrayView<2, T, C> & operator*=(MultiArrayView<2, T, C> const & l, T other)
{
    MultiArrayView<2, T, C> & left = const_cast<MultiArrayView<2, T, C> &>(l);
    const std::size_t rows = rowCount(left);
    const std::size_t cols = columnCount(left);

    for(std::size_t i = 0; i < cols; ++i)
        for(std::size_t j = 0; j < rows; ++j)
            left(j, i) *= other;
    return left;
}

template <class T, class C>
MultiArrayView<2, T, C> & operator/=(MultiArrayView<2, T, C> const & l, T other)
{
    MultiArrayView<2, T, C> & left = const_cast<MultiArrayView<2, T, C> &>(l);
    const std::size_t rows = rowCount(left);
    const std::size_t cols = columnCount(left);

    for(std::size_t i = 0; i < cols; ++i)
        for(std::size_t j = 0; j < rows; ++j)
            left(j, i) /= other;
    return left;
}

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
    Matrix(std::size_t rows, std::size_t columns,
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
    Matrix(std::size_t rows, std::size_t columns, const_reference init,
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
    Matrix(std::size_t rows, std::size_t columns, const_pointer init, RawArrayMemoryLayout layout = RowMajor,
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

        /** assign a temporary matrix. If the shapes of the two matrices match,
            only the data are copied (in order to not invalidate views and iterators
            depending on this matrix). Otherwise, the memory is swapped
            between the two matrices, so that all depending objects
            (array views, iterators) ar invalidated.
         */
    Matrix & operator=(const TemporaryMatrix<T, ALLOC> &rhs)
    {
        if(this->shape() == rhs.shape())
            this->copy(rhs);
        else
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

        /** reshape to the given shape and initialize with zero.
         */
    void reshape(std::size_t rows, std::size_t columns)
    {
        BaseType::reshape(difference_type(rows, columns));
    }

        /** reshape to the given shape and initialize with \a init.
         */
    void reshape(std::size_t rows, std::size_t columns, const_reference init)
    {
        BaseType::reshape(difference_type(rows, columns), init);
    }

        /** reshape to the given shape and initialize with zero.
         */
    void reshape(difference_type const & shape)
    {
        BaseType::reshape(shape);
    }

        /** reshape to the given shape and initialize with \a init.
         */
    void reshape(difference_type const & shape, const_reference init)
    {
        BaseType::reshape(shape, init);
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

        /** check whether the matrix is symmetric.
        */
    bool isSymmetric() const
    {
        return vigra::linalg::isSymmetric(*this);
    }

#ifdef DOXYGEN
// repeat the following functions for documentation. In real code, they are inherited.

        /** read/write access to matrix element <tt>(row, column)</tt>.
            Note that the order of the argument is the opposite of the usual
            VIGRA convention due to column-major matrix order.
        */
    value_type & operator()(std::size_t row, std::size_t column);

        /** read access to matrix element <tt>(row, column)</tt>.
            Note that the order of the argument is the opposite of the usual
            VIGRA convention due to column-major matrix order.
        */
    value_type operator()(std::size_t row, std::size_t column) const;

        /** squared Frobenius norm. Sum of squares of the matrix elements.
        */
    typename NormTraits<Matrix>::SquaredNormType squaredNorm() const;

        /** Frobenius norm. Root of sum of squares of the matrix elements.
        */
    typename NormTraits<Matrix>::NormType norm() const;

        /** create a transposed view of this matrix.
            No data are copied. If you want to transpose this matrix permanently, 
            you have to assign the transposed view:
            
            \code
            a = a.transpose();
            \endcode
         */
    MultiArrayView<2, vluae_type, StridedArrayTag> transpose() const;
#endif

        /** add \a other to this (sizes must match).
         */
    template <class U, class C>
    Matrix & operator+=(MultiArrayView<2, U, C> const & other)
    {
        vigra::linalg::operator+=(*this, other);
        return *this;
    }

        /** subtract \a other from this (sizes must match).
         */
    template <class U, class C>
    Matrix & operator-=(MultiArrayView<2, U, C> const & other)
    {
        vigra::linalg::operator-=(*this, other);
        return *this;
    }

        /** add \a other to each element of this matrix
         */
    Matrix & operator+=(T other)
    {
        vigra::linalg::operator+=(*this, other);
        return *this;
    }

        /** subtraxt \a other from each element of this matrix
         */
    Matrix & operator-=(T other)
    {
        vigra::linalg::operator-=(*this, other);
        return *this;
    }

        /** scalar multiply this with \a other
         */
    Matrix & operator*=(T other)
    {
        vigra::linalg::operator*=(*this, other);
        return *this;
    }

        /** scalar devide this by \a other
         */
    Matrix & operator/=(T other)
    {
        vigra::linalg::operator/=(*this, other);
        return *this;
    }
};

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

    TemporaryMatrix(difference_type const & shape)
    : BaseType(shape, ALLOC())
    {}

    TemporaryMatrix(difference_type const & shape, const_reference init)
    : BaseType(shape, init, ALLOC())
    {}

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
    
    template <class U, class C>
    TemporaryMatrix & operator+=(MultiArrayView<2, U, C> const & other)
    {
        vigra::linalg::operator+=(*this, other);
        return *this;
    }

    template <class U, class C>
    TemporaryMatrix & operator-=(MultiArrayView<2, U, C> const & other)
    {
        vigra::linalg::operator-=(*this, other);
        return *this;
    }

    TemporaryMatrix & operator+=(T other)
    {
        vigra::linalg::operator+=(*this, other);
        return *this;
    }

    TemporaryMatrix & operator-=(T other)
    {
        vigra::linalg::operator-=(*this, other);
        return *this;
    }

    TemporaryMatrix & operator*=(T other)
    {
        vigra::linalg::operator*=(*this, other);
        return *this;
    }

    TemporaryMatrix & operator/=(T other)
    {
        vigra::linalg::operator/=(*this, other);
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
inline std::size_t rowCount(const MultiArrayView<2, T, C> &x)
{
    return x.shape(0);
}

    /** Number of columns of a matrix represented as a <tt>MultiArrayView&lt;2,...&gt;</tt>

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C>
inline std::size_t columnCount(const MultiArrayView<2, T, C> &x)
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


    /** Create a row vector view of the matrix \a m starting at element \a first and ranging 
        to column \a end (non-inclusive).

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C>
MultiArrayView <2, T, C>
rowVector(MultiArrayView <2, T, C> const & m, MultiArrayShape<2>::type first, int end)
{
    typedef typename MultiArrayView <2, T, C>::difference_type Shape;
    return m.subarray(first, Shape(first[0]+1, end));
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

    /** Create a column vector view of the matrix \a m starting at element \a first and 
        ranging to row \a end (non-inclusive).

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     **/
template <class T, class C>
MultiArrayView <2, T, C>
columnVector(MultiArrayView<2, T, C> const & m, MultiArrayShape<2>::type first, int end)
{
    typedef typename MultiArrayView <2, T, C>::difference_type Shape;
    return m.subarray(first, Shape(end, first[1]+1));
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
    const std::size_t size = rowCount(m);
    if(size != columnCount(m))
        return false;

    for(std::size_t i = 0; i < size; ++i)
        for(std::size_t j = i+1; j < size; ++j)
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

    /** calculate the Frobenius norm of a matrix.
        Equal to the root of the sum of squares of the matrix elements.

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
    const std::size_t rows = rowCount(r);
    vigra_precondition(rows == columnCount(r),
       "identityMatrix(): Matrix must be square.");
    for(std::size_t i = 0; i < rows; ++i) {
        for(std::size_t j = 0; j < rows; ++j)
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
TemporaryMatrix<T> identityMatrix(std::size_t size)
{
    TemporaryMatrix<T> ret(size, size, NumericTraits<T>::zero());
    for(std::size_t i = 0; i < size; ++i)
        ret(i, i) = NumericTraits<T>::one();
    return ret;
}

template <class T, class C1, class C2>
void diagonalMatrixImpl(MultiArrayView<1, T, C1> const & v, MultiArrayView<2, T, C2> &r)
{
    const std::size_t size = v.elementCount();
    vigra_precondition(rowCount(r) == size && columnCount(r) == size,
        "diagonalMatrix(): result must be a square matrix.");
    for(std::size_t i = 0; i < size; ++i)
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
    std::size_t size = v.elementCount();
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
    const std::size_t rows = rowCount(r);
    const std::size_t cols = columnCount(r);
    vigra_precondition(rows == columnCount(v) && cols == rowCount(v),
       "transpose(): arrays must have transposed shapes.");
    for(std::size_t i = 0; i < cols; ++i)
        for(std::size_t j = 0; j < rows; ++j)
            r(j, i) = v(i, j);
}

    /** create the transpose of matrix \a v.
        This does not copy any data, but only creates a transposed view 
        to the original matrix. A copy is only made when the transposed view
        is assigned to another matrix.
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
MultiArrayView<2, T, StridedArrayTag> transpose(MultiArrayView<2, T, C> const & v)
{
    return v.transpose();
}

    /** Create new matrix by concatenating two matrices \a a and \a b vertically, i.e. on top of each other.
        The two matrices must have the same number of columns.
        The result is returned as a temporary matrix.

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */
template <class T, class C1, class C2>
inline TemporaryMatrix<T>
joinVertically(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b)
{
    typedef typename TemporaryMatrix<T>::difference_type Shape;
    
    int n = columnCount(a);
    vigra_precondition(n == columnCount(b),
       "joinVertically(): shape mismatch.");
    
    int ma = rowCount(a);
    int mb = rowCount(b);
    TemporaryMatrix<T> t(ma + mb, n, T());
    t.subarray(Shape(0,0), Shape(ma, n)) = a;
    t.subarray(Shape(ma,0), Shape(ma+mb, n)) = b;
    return t;
}

    /** Create new matrix by concatenating two matrices \a a and \a b horizontally, i.e. side by side.
        The two matrices must have the same number of rows.
        The result is returned as a temporary matrix.

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */
template <class T, class C1, class C2>
inline TemporaryMatrix<T>
joinHorizontally(const MultiArrayView<2, T, C1> &a, const MultiArrayView<2, T, C2> &b)
{
    typedef typename TemporaryMatrix<T>::difference_type Shape;
    
    int m = rowCount(a);
    vigra_precondition(m == rowCount(b),
       "joinHorizontally(): shape mismatch.");
    
    int na = columnCount(a);
    int nb = columnCount(b);
    TemporaryMatrix<T> t(m, na + nb, T());
    t.subarray(Shape(0,0), Shape(m, na)) = a;
    t.subarray(Shape(0, na), Shape(m, na + nb)) = b;
    return t;
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
    const std::size_t rrows = rowCount(r);
    const std::size_t rcols = columnCount(r);
    vigra_precondition(rrows == rowCount(a) && rcols == columnCount(a) &&
                       rrows == rowCount(b) && rcols == columnCount(b),
                       "add(): Matrix shapes must agree.");

    for(std::size_t i = 0; i < rcols; ++i) {
        for(std::size_t j = 0; j < rrows; ++j) {
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

    /** add scalar \a b to every element of the matrix \a a.
        The result is returned as a temporary matrix.

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */
template <class T, class C>
inline TemporaryMatrix<T>
operator+(const MultiArrayView<2, T, C> &a, T b)
{
    return TemporaryMatrix<T>(a) += b;
}

template <class T>
inline TemporaryMatrix<T>
operator+(const TemporaryMatrix<T> &a, T b)
{
    return const_cast<TemporaryMatrix<T> &>(a) += b;
}

    /** add scalar \a a to every element of the matrix \a b.
        The result is returned as a temporary matrix.

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */
template <class T, class C>
inline TemporaryMatrix<T>
operator+(T a, const MultiArrayView<2, T, C> &b)
{
    return TemporaryMatrix<T>(b) += a;
}

template <class T>
inline TemporaryMatrix<T>
operator+(T a, const TemporaryMatrix<T> &b)
{
    return const_cast<TemporaryMatrix<T> &>(b) += a;
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
    const std::size_t rrows = rowCount(r);
    const std::size_t rcols = columnCount(r);
    vigra_precondition(rrows == rowCount(a) && rcols == columnCount(a) &&
                       rrows == rowCount(b) && rcols == columnCount(b),
                       "subtract(): Matrix shapes must agree.");

    for(std::size_t i = 0; i < rcols; ++i) {
        for(std::size_t j = 0; j < rrows; ++j) {
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
    const std::size_t rows = rowCount(a);
    const std::size_t cols = columnCount(a);
    vigra_precondition(rows == b.rowCount() && cols == b.columnCount(),
       "Matrix::operator-(): Shape mismatch.");

    for(std::size_t i = 0; i < cols; ++i)
        for(std::size_t j = 0; j < rows; ++j)
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

    /** subtract scalar \a b from every element of the matrix \a a.
        The result is returned as a temporary matrix.

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */
template <class T, class C>
inline TemporaryMatrix<T>
operator-(const MultiArrayView<2, T, C> &a, T b)
{
    return TemporaryMatrix<T>(a) -= b;
}

template <class T>
inline TemporaryMatrix<T>
operator-(const TemporaryMatrix<T> &a, T b)
{
    return const_cast<TemporaryMatrix<T> &>(a) -= b;
}

    /** subtract every element of the matrix \a b from scalar \a a.
        The result is returned as a temporary matrix.

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */
template <class T, class C>
inline TemporaryMatrix<T>
operator-(T a, const MultiArrayView<2, T, C> &b)
{
    return TemporaryMatrix<T>(b.shape(), a) -= b;
}

    /** calculate the inner product of two matrices representing vectors.
        Typically, matrix \a x has a single row, and matrix \a y has
        a single column, and the other dimensions match. In addition, this
        function handles the cases when either or both of the two inputs are 
        transposed (e.g. it can compute the dot product of two column vectors). 
        A <tt>PreconditionViolation</tt> exception is thrown when
        the shape conditions are violated. 

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2>
typename NormTraits<T>::SquaredNormType 
dot(const MultiArrayView<2, T, C1> &x, const MultiArrayView<2, T, C2> &y)
{
    typename NormTraits<T>::SquaredNormType ret = NumericTraits<typename NormTraits<T>::SquaredNormType>::zero();
    if(y.shape(1) == 1)
    {
        unsigned int size = y.shape(0);
        if(x.shape(0) == 1 && x.shape(1) == size) // proper scalar product
            for(std::size_t i = 0; i < size; ++i)
                ret += x(0, i) * y(i, 0);
        else if(x.shape(1) == 1 && x.shape(0) == size) // two column vectors
            for(std::size_t i = 0; i < size; ++i)
                ret += x(i, 0) * y(i, 0);
        else    
            vigra_precondition(false, "dot(): wrong matrix shapes.");
    }
    else if(y.shape(0) == 1)
    {
        unsigned int size = y.shape(1);
        if(x.shape(0) == 1 && x.shape(1) == size) // two row vectors
            for(std::size_t i = 0; i < size; ++i)
                ret += x(0, i) * y(0, i);
        else if(x.shape(1) == 1 && x.shape(0) == size) // column dot row
            for(std::size_t i = 0; i < size; ++i)
                ret += x(i, 0) * y(0, i);
        else    
            vigra_precondition(false, "dot(): wrong matrix shapes.");
    }
    else
            vigra_precondition(false, "dot(): wrong matrix shapes.");
    return ret;
}

    /** calculate the inner product of two vectors. The vector
        lengths must match.

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2>
typename NormTraits<T>::SquaredNormType 
dot(const MultiArrayView<1, T, C1> &x, const MultiArrayView<1, T, C2> &y)
{
    const std::size_t n = x.elementCount();
    vigra_precondition(n == y.elementCount(),
       "dot(): shape mismatch.");
    typename NormTraits<T>::SquaredNormType ret = NumericTraits<typename NormTraits<T>::SquaredNormType>::zero();
    for(std::size_t i = 0; i < n; ++i)
        ret += x(i) * y(i);
    return ret;
}

    /** calculate the cross product of two vectors of length 3.
        The result is written into \a r.

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
void cross(const MultiArrayView<1, T, C1> &x, const MultiArrayView<1, T, C2> &y,
           MultiArrayView<1, T, C3> &r)
{
    vigra_precondition(3 == x.elementCount() && 3 == y.elementCount() && 3 == r.elementCount(),
       "cross(): vectors must have length 3.");
    r(0) = x(1)*y(2) - x(2)*y(1);
    r(1) = x(2)*y(0) - x(0)*y(2);
    r(2) = x(0)*y(1) - x(1)*y(0);
}

    /** calculate the cross product of two matrices representing vectors.
        That is, \a x, \a y, and \a r must have a single column of length 3. The result
        is written into \a r.

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2, class C3>
void cross(const MultiArrayView<2, T, C1> &x, const MultiArrayView<2, T, C2> &y,
           MultiArrayView<2, T, C3> &r)
{
    vigra_precondition(3 == rowCount(x) && 3 == rowCount(y) && 3 == rowCount(r) ,
       "cross(): vectors must have length 3.");
    r(0,0) = x(1,0)*y(2,0) - x(2,0)*y(1,0);
    r(1,0) = x(2,0)*y(0,0) - x(0,0)*y(2,0);
    r(2,0) = x(0,0)*y(1,0) - x(1,0)*y(0,0);
}

    /** calculate the cross product of two matrices representing vectors.
        That is, \a x, and \a y must have a single column of length 3. The result
        is returned as a temporary matrix.

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1, class C2>
TemporaryMatrix<T>
cross(const MultiArrayView<2, T, C1> &x, const MultiArrayView<2, T, C2> &y)
{
    TemporaryMatrix<T> ret(3, 1);
    cross(x, y, ret);
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
    const std::size_t rows = rowCount(r);
    const std::size_t cols = columnCount(r);
    vigra_precondition(rows == rowCount(x) && cols == columnCount(y) &&
                       1 == columnCount(x) && 1 == rowCount(y),
       "outer(): shape mismatch.");
    for(std::size_t i = 0; i < cols; ++i)
        for(std::size_t j = 0; j < rows; ++j)
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
    const std::size_t rows = rowCount(x);
    const std::size_t cols = columnCount(y);
    vigra_precondition(1 == columnCount(x) && 1 == rowCount(y),
       "outer(): shape mismatch.");
    TemporaryMatrix<T> ret(rows, cols);
    outer(x, y, ret);
    return ret;
}

    /** calculate the outer product of a matrix (representing a vector) with itself.
        The result is returned as a temporary matrix.

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespaces: vigra and vigra::linalg
     */
template <class T, class C1>
TemporaryMatrix<T>
outer(const MultiArrayView<2, T, C1> &x)
{
    const std::size_t rows = rowCount(x);
    const std::size_t cols = columnCount(x);
    vigra_precondition(rows == 1 || cols == 1,
       "outer(): matrix does not represent a vector.");
    const std::size_t size = std::max(rows, cols);
    TemporaryMatrix<T> ret(size, size);

    if(rows == 1)
    {
        for(std::size_t i = 0; i < size; ++i)
            for(std::size_t j = 0; j < size; ++j)
                ret(j, i) = x(0, j) * x(0, i);
    }
    else
    {
        for(std::size_t i = 0; i < size; ++i)
            for(std::size_t j = 0; j < size; ++j)
                ret(j, i) = x(j, 0) * x(i, 0);
    }
    return ret;
}

template <class T>
class PointWise
{
  public:
    T const & t;
    
    PointWise(T const & it)
    : t(it)
    {}
};

template <class T>
PointWise<T> pointWise(T const & t)
{
    return PointWise<T>(t);
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
    const std::size_t rows = rowCount(a);
    const std::size_t cols = columnCount(a);
    vigra_precondition(rows == rowCount(r) && cols == columnCount(r),
                       "smul(): Matrix sizes must agree.");

    for(std::size_t i = 0; i < cols; ++i)
        for(std::size_t j = 0; j < rows; ++j)
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
    const std::size_t rrows = rowCount(r);
    const std::size_t rcols = columnCount(r);
    const std::size_t acols = columnCount(a);
    vigra_precondition(rrows == rowCount(a) && rcols == columnCount(b) && acols == rowCount(b),
                       "mmul(): Matrix shapes must agree.");

    // order of loops ensures that inner loop goes down columns
    for(std::size_t i = 0; i < rcols; ++i) 
    {
        for(std::size_t j = 0; j < rrows; ++j) 
            r(j, i) = a(j, 0) * b(0, i);
        for(std::size_t k = 1; k < acols; ++k) 
            for(std::size_t j = 0; j < rrows; ++j) 
                r(j, i) += a(j, k) * b(k, i);
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
    const std::size_t rrows = rowCount(r);
    const std::size_t rcols = columnCount(r);
    vigra_precondition(rrows == rowCount(a) && rcols == columnCount(a) &&
                       rrows == rowCount(b) && rcols == columnCount(b),
                       "pmul(): Matrix shapes must agree.");

    for(std::size_t i = 0; i < rcols; ++i) {
        for(std::size_t j = 0; j < rrows; ++j) {
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
    TemporaryMatrix<T> ret(a.shape());
    pmul(a, b, ret);
    return ret;
}

    /** multiply matrices \a a and \a b pointwise.
        \a a and \a b must have matching shapes.
        The result is returned as a temporary matrix.
        
        Usage:
        
        \code
        Matrix<double> a(m,n), b(m,n);
        
        Matrix<double> c = a * pointWise(b);
        // is equivalent to
        // Matrix<double> c = pmul(a, b);
        \endcode

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */
template <class T, class C, class U>
inline TemporaryMatrix<T>
operator*(const MultiArrayView<2, T, C> &a, PointWise<U> b)
{
    return pmul(a, b.t);
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
    return const_cast<TemporaryMatrix<T> &>(b) *= a;
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
    for(std::size_t i = 1; i < N; ++i)
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
    for(std::size_t i = 0; i < N; ++i)
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
    const std::size_t rows = rowCount(a);
    const std::size_t cols = columnCount(a);
    vigra_precondition(rows == rowCount(r) && cols == columnCount(r),
                       "sdiv(): Matrix sizes must agree.");

    for(std::size_t i = 0; i < cols; ++i)
        for(std::size_t j = 0; j < rows; ++j)
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
    const std::size_t rrows = rowCount(r);
    const std::size_t rcols = columnCount(r);
    vigra_precondition(rrows == rowCount(a) && rcols == columnCount(a) &&
                       rrows == rowCount(b) && rcols == columnCount(b),
                       "pdiv(): Matrix shapes must agree.");

    for(std::size_t i = 0; i < rcols; ++i) {
        for(std::size_t j = 0; j < rrows; ++j) {
            r(j, i) = a(j, i) / b(j, i);
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
    TemporaryMatrix<T> ret(a.shape());
    pdiv(a, b, ret);
    return ret;
}

    /** divide matrices \a a and \a b pointwise.
        \a a and \a b must have matching shapes.
        The result is returned as a temporary matrix.
        
        Usage:
        
        \code
        Matrix<double> a(m,n), b(m,n);
        
        Matrix<double> c = a / pointWise(b);
        // is equivalent to
        // Matrix<double> c = pdiv(a, b);
        \endcode

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */
template <class T, class C, class U>
inline TemporaryMatrix<T>
operator/(const MultiArrayView<2, T, C> &a, PointWise<U> b)
{
    return pdiv(a, b.t);
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

    /** Create a matrix whose elements are the quotients between scalar \a a and
        matrix \a b. The result is returned as a temporary matrix.

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra::linalg
     */
template <class T, class C>
inline TemporaryMatrix<T>
operator/(T a, const MultiArrayView<2, T, C> &b)
{
    return TemporaryMatrix<T>(b.shape(), a) / pointWise(b);
}

    /*! Find the index of the minimum element in a matrix.
    
        The function returns the index in column-major scan-order sense,
        i.e. according to the order used by <tt>MultiArrayView::operator[]</tt>.
        If the matrix is actually a vector, this is just the row or columns index.
        In case of a truely 2-dimensional matrix, the index can be converted to an 
        index pair by calling <tt>MultiArrayView::scanOrderIndexToCoordinate()</tt>
        
        <b>Required Interface:</b>
        
        \code
        bool f = a[0] < NumericTraits<T>::max();
        \endcode

        <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>"<br>
        Namespace: vigra
    */
template <class T, class C>
int argMin(MultiArrayView<2, T, C> const & a)
{
    T vopt = NumericTraits<T>::max();
    int best = -1;
    for(unsigned int k=0; k < a.size(); ++k)
    {
        if(a[k] < vopt)
        {
            vopt = a[k];
            best = k;
        }
    }
    return best;
}

    /*! Find the index of the maximum element in a matrix.
    
        The function returns the index in column-major scan-order sense,
        i.e. according to the order used by <tt>MultiArrayView::operator[]</tt>.
        If the matrix is actually a vector, this is just the row or columns index.
        In case of a truely 2-dimensional matrix, the index can be converted to an 
        index pair by calling <tt>MultiArrayView::scanOrderIndexToCoordinate()</tt>
        
        <b>Required Interface:</b>
        
        \code
        bool f = NumericTraits<T>::min() < a[0];
        \endcode

        <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>"<br>
        Namespace: vigra
    */
template <class T, class C>
int argMax(MultiArrayView<2, T, C> const & a)
{
    T vopt = NumericTraits<T>::min();
    int best = -1;
    for(unsigned int k=0; k < a.size(); ++k)
    {
        if(vopt < a[k])
        {
            vopt = a[k];
            best = k;
        }
    }
    return best;
}

    /*! Find the index of the minimum element in a matrix subject to a condition.
    
        The function returns <tt>-1</tt> if no element conforms to \a condition.
        Otherwise, the index of the maximum element is returned in column-major scan-order sense,
        i.e. according to the order used by <tt>MultiArrayView::operator[]</tt>.
        If the matrix is actually a vector, this is just the row or columns index.
        In case of a truely 2-dimensional matrix, the index can be converted to an 
        index pair by calling <tt>MultiArrayView::scanOrderIndexToCoordinate()</tt>
        
        <b>Required Interface:</b>
        
        \code
        bool c = condition(a[0]);
        bool f = a[0] < NumericTraits<T>::max();
        \endcode

        <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>"<br>
        Namespace: vigra
    */
template <class T, class C, class UnaryFunctor>
int argMinIf(MultiArrayView<2, T, C> const & a, UnaryFunctor condition)
{
    T vopt = NumericTraits<T>::max();
    int best = -1;
    for(unsigned int k=0; k < a.size(); ++k)
    {
        if(condition(a[k]) && a[k] < vopt)
        {
            vopt = a[k];
            best = k;
        }
    }
    return best;
}

    /*! Find the index of the maximum element in a matrix subject to a condition.
    
        The function returns <tt>-1</tt> if no element conforms to \a condition.
        Otherwise, the index of the maximum element is returned in column-major scan-order sense,
        i.e. according to the order used by <tt>MultiArrayView::operator[]</tt>.
        If the matrix is actually a vector, this is just the row or columns index.
        In case of a truely 2-dimensional matrix, the index can be converted to an 
        index pair by calling <tt>MultiArrayView::scanOrderIndexToCoordinate()</tt>
        
        <b>Required Interface:</b>
        
        \code
        bool c = condition(a[0]);
        bool f = NumericTraits<T>::min() < a[0];
        \endcode

        <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>"<br>
        Namespace: vigra
    */
template <class T, class C, class UnaryFunctor>
int argMaxIf(MultiArrayView<2, T, C> const & a, UnaryFunctor condition)
{
    T vopt = NumericTraits<T>::min();
    int best = -1;
    for(unsigned int k=0; k < a.size(); ++k)
    {
        if(condition(a[k]) && vopt < a[k])
        {
            vopt = a[k];
            best = k;
        }
    }
    return best;
}

/** Matrix point-wise power.
*/
template <class T, class C>
linalg::TemporaryMatrix<T> pow(MultiArrayView<2, T, C> const & v, T exponent)
{
    linalg::TemporaryMatrix<T> t(v.shape());
    unsigned int m = rowCount(v), n = columnCount(v);

    for(unsigned int i = 0; i < n; ++i)
        for(unsigned int j = 0; j < m; ++j)
            t(j, i) = vigra::pow(v(j, i), exponent);
    return t;
}

template <class T>
linalg::TemporaryMatrix<T> pow(linalg::TemporaryMatrix<T> const & v, T exponent)
{
    linalg::TemporaryMatrix<T> & t = const_cast<linalg::TemporaryMatrix<T> &>(v);
    unsigned int m = rowCount(t), n = columnCount(t);

    for(unsigned int i = 0; i < n; ++i)
        for(unsigned int j = 0; j < m; ++j)
            t(j, i) = vigra::pow(t(j, i), exponent);
    return t;
}

template <class T, class C>
linalg::TemporaryMatrix<T> pow(MultiArrayView<2, T, C> const & v, int exponent)
{
    linalg::TemporaryMatrix<T> t(v.shape());
    unsigned int m = rowCount(v), n = columnCount(v);

    for(unsigned int i = 0; i < n; ++i)
        for(unsigned int j = 0; j < m; ++j)
            t(j, i) = vigra::pow(v(j, i), exponent);
    return t;
}

template <class T>
linalg::TemporaryMatrix<T> pow(linalg::TemporaryMatrix<T> const & v, int exponent)
{
    linalg::TemporaryMatrix<T> & t = const_cast<linalg::TemporaryMatrix<T> &>(v);
    unsigned int m = rowCount(t), n = columnCount(t);

    for(unsigned int i = 0; i < n; ++i)
        for(unsigned int j = 0; j < m; ++j)
            t(j, i) = vigra::pow(t(j, i), exponent);
    return t;
}

template <class C>
linalg::TemporaryMatrix<int> pow(MultiArrayView<2, int, C> const & v, int exponent)
{
    linalg::TemporaryMatrix<int> t(v.shape());
    unsigned int m = rowCount(v), n = columnCount(v);

    for(unsigned int i = 0; i < n; ++i)
        for(unsigned int j = 0; j < m; ++j)
            t(j, i) = (int)vigra::pow((double)v(j, i), exponent);
    return t;
}

inline
linalg::TemporaryMatrix<int> pow(linalg::TemporaryMatrix<int> const & v, int exponent)
{
    linalg::TemporaryMatrix<int> & t = const_cast<linalg::TemporaryMatrix<int> &>(v);
    unsigned int m = rowCount(t), n = columnCount(t);

    for(unsigned int i = 0; i < n; ++i)
        for(unsigned int j = 0; j < m; ++j)
            t(j, i) = (int)vigra::pow((double)t(j, i), exponent);
    return t;
}

    /** Matrix point-wise sqrt. */
template <class T, class C>
linalg::TemporaryMatrix<T> sqrt(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise exp. */
template <class T, class C>
linalg::TemporaryMatrix<T> exp(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise log. */
template <class T, class C>
linalg::TemporaryMatrix<T> log(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise log10. */
template <class T, class C>
linalg::TemporaryMatrix<T> log10(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise sin. */
template <class T, class C>
linalg::TemporaryMatrix<T> sin(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise asin. */
template <class T, class C>
linalg::TemporaryMatrix<T> asin(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise cos. */
template <class T, class C>
linalg::TemporaryMatrix<T> cos(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise acos. */
template <class T, class C>
linalg::TemporaryMatrix<T> acos(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise tan. */
template <class T, class C>
linalg::TemporaryMatrix<T> tan(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise atan. */
template <class T, class C>
linalg::TemporaryMatrix<T> atan(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise round. */
template <class T, class C>
linalg::TemporaryMatrix<T> round(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise floor. */
template <class T, class C>
linalg::TemporaryMatrix<T> floor(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise ceil. */
template <class T, class C>
linalg::TemporaryMatrix<T> ceil(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise abs. */
template <class T, class C>
linalg::TemporaryMatrix<T> abs(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise square. */
template <class T, class C>
linalg::TemporaryMatrix<T> sq(MultiArrayView<2, T, C> const & v);
    /** Matrix point-wise sign. */
template <class T, class C>
linalg::TemporaryMatrix<T> sign(MultiArrayView<2, T, C> const & v);

#define VIGRA_MATRIX_UNARY_FUNCTION(FUNCTION, NAMESPACE) \
using NAMESPACE::FUNCTION; \
template <class T, class C> \
linalg::TemporaryMatrix<T> FUNCTION(MultiArrayView<2, T, C> const & v) \
{ \
    linalg::TemporaryMatrix<T> t(v.shape()); \
    unsigned int m = rowCount(v), n = columnCount(v); \
 \
    for(unsigned int i = 0; i < n; ++i) \
        for(unsigned int j = 0; j < m; ++j) \
            t(j, i) = NAMESPACE::FUNCTION(v(j, i)); \
    return t; \
} \
 \
template <class T> \
linalg::TemporaryMatrix<T> FUNCTION(linalg::Matrix<T> const & v) \
{ \
    linalg::TemporaryMatrix<T> t(v.shape()); \
    unsigned int m = rowCount(v), n = columnCount(v); \
 \
    for(unsigned int i = 0; i < n; ++i) \
        for(unsigned int j = 0; j < m; ++j) \
            t(j, i) = NAMESPACE::FUNCTION(v(j, i)); \
    return t; \
} \
 \
template <class T> \
linalg::TemporaryMatrix<T> FUNCTION(linalg::TemporaryMatrix<T> const & v) \
{ \
    linalg::TemporaryMatrix<T> & t = const_cast<linalg::TemporaryMatrix<T> &>(v); \
    unsigned int m = rowCount(t), n = columnCount(t); \
 \
    for(unsigned int i = 0; i < n; ++i) \
        for(unsigned int j = 0; j < m; ++j) \
            t(j, i) = NAMESPACE::FUNCTION(t(j, i)); \
    return v; \
}\
}\
using linalg::FUNCTION;\
namespace linalg {

VIGRA_MATRIX_UNARY_FUNCTION(sqrt, std)
VIGRA_MATRIX_UNARY_FUNCTION(exp, std)
VIGRA_MATRIX_UNARY_FUNCTION(log, std)
VIGRA_MATRIX_UNARY_FUNCTION(log10, std)
VIGRA_MATRIX_UNARY_FUNCTION(sin, std)
VIGRA_MATRIX_UNARY_FUNCTION(asin, std)
VIGRA_MATRIX_UNARY_FUNCTION(cos, std)
VIGRA_MATRIX_UNARY_FUNCTION(acos, std)
VIGRA_MATRIX_UNARY_FUNCTION(tan, std)
VIGRA_MATRIX_UNARY_FUNCTION(atan, std)
VIGRA_MATRIX_UNARY_FUNCTION(round, vigra)
VIGRA_MATRIX_UNARY_FUNCTION(floor, vigra)
VIGRA_MATRIX_UNARY_FUNCTION(ceil, vigra)
VIGRA_MATRIX_UNARY_FUNCTION(abs, vigra)
VIGRA_MATRIX_UNARY_FUNCTION(sq, vigra)
VIGRA_MATRIX_UNARY_FUNCTION(sign, vigra)

#undef VIGRA_MATRIX_UNARY_FUNCTION

//@}

} // namespace linalg

using linalg::RowMajor;
using linalg::ColumnMajor;
using linalg::Matrix;
using linalg::identityMatrix;
using linalg::diagonalMatrix;
using linalg::transpose;
using linalg::pointWise;
using linalg::dot;
using linalg::cross;
using linalg::outer;
using linalg::rowCount;
using linalg::columnCount;
using linalg::rowVector;
using linalg::columnVector;
using linalg::isSymmetric;
using linalg::joinHorizontally;
using linalg::joinVertically;
using linalg::argMin;
using linalg::argMinIf;
using linalg::argMax;
using linalg::argMaxIf;

/********************************************************/
/*                                                      */
/*                       NormTraits                     */
/*                                                      */
/********************************************************/

template <class T, class ALLOC>
struct NormTraits<Matrix<T, ALLOC> >
: public NormTraits<MultiArray<2, T, ALLOC> >
{
    typedef NormTraits<MultiArray<2, T, ALLOC> > BaseType;
    typedef Matrix<T, ALLOC>                     Type;
    typedef typename BaseType::SquaredNormType   SquaredNormType;
    typedef typename BaseType::NormType          NormType;
};

template <class T, class ALLOC>
struct NormTraits<linalg::TemporaryMatrix<T, ALLOC> >
: public NormTraits<Matrix<T, ALLOC> >
{
    typedef NormTraits<Matrix<T, ALLOC> >        BaseType;
    typedef linalg::TemporaryMatrix<T, ALLOC>    Type;
    typedef typename BaseType::SquaredNormType   SquaredNormType;
    typedef typename BaseType::NormType          NormType;
};

/** \addtogroup LinearAlgebraFunctions Matrix functions
 */
//@{

    /** print a matrix \a m to the stream \a s.

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: std
     */
template <class T, class C>
std::ostream &
operator<<(std::ostream & s, const vigra::MultiArrayView<2, T, C> &m)
{
    const std::size_t rows = vigra::linalg::rowCount(m);
    const std::size_t cols = vigra::linalg::columnCount(m);
    std::ios::fmtflags flags =
        s.setf(std::ios::right | std::ios::fixed, std::ios::adjustfield | std::ios::floatfield);
    for(std::size_t j = 0; j < rows; ++j)
    {
        for(std::size_t i = 0; i < cols; ++i)
        {
            s << std::setw(7) << std::setprecision(4) << m(j, i) << " ";
        }
        s << std::endl;
    }
    s.setf(flags);
    return s;
}

//@}

namespace linalg {

namespace detail {

template <class T1, class C1, class T2, class C2, class T3, class C3>
void
columnStatisticsImpl(const MultiArrayView<2, T1, C1> & A, 
                 MultiArrayView<2, T2, C2> & mean, MultiArrayView<2, T3, C3> & sumOfSquaredDifferences)
{
    unsigned int m = rowCount(A);
    unsigned int n = columnCount(A);
    vigra_precondition(1 == rowCount(mean) && n == columnCount(mean) &&
                       1 == rowCount(sumOfSquaredDifferences) && n == columnCount(sumOfSquaredDifferences),
                       "columnStatistics(): Shape mismatch between input and output.");

    // West's algorithm for incremental variance computation
    mean.init(NumericTraits<T2>::zero());
    sumOfSquaredDifferences.init(NumericTraits<T3>::zero());
    
    for(unsigned int k=0; k<m; ++k)
    {
        Matrix<T2> t = rowVector(A, k) - mean;
        double f  = 1.0 / (k + 1.0),
               f1 = 1.0 - f;
        mean += f*t;
        sumOfSquaredDifferences += f1*sq(t);
    }
}

} // namespace detail

/** \addtogroup LinearAlgebraFunctions Matrix functions
 */
//@{
    /** compute mean of every column of matrix \a A.
    
        \a mean must be a row vector with as many columns as \a A.

    <b>\#include</b> "<a href="matrix_8hxx-source.html">vigra/matrix.hxx</a>" or<br>
    <b>\#include</b> "<a href="linear__algebra_8hxx-source.html">vigra/linear_algebra.hxx</a>"<br>
        Namespace: vigra
     */
template <class T1, class C1, class T2, class C2>
void
columnStatistics(const MultiArrayView<2, T1, C1> & A, 
                 MultiArrayView<2, T2, C2> & mean)
{
    unsigned int m = rowCount(A);
    unsigned int n = columnCount(A);
    vigra_precondition(1 == rowCount(mean) && n == columnCount(mean),
                       "columnStatistics(): Shape mismatch between input and output.");

    mean.init(NumericTraits<T2>::zero());
    
    for(unsigned int k=0; k<m; ++k)
    {
        mean += rowVector(A, k);
    }
    mean /= T2(m);
}

template <class T1, class C1, class T2, class C2, class T3, class C3>
void
columnStatistics(const MultiArrayView<2, T1, C1> & A, 
                 MultiArrayView<2, T2, C2> & mean, MultiArrayView<2, T3, C3> & stdDev)
{
    detail::columnStatisticsImpl(A, mean, stdDev);
    
#if 0
    if(rowCount(A) > 1)
        stdDev = sqrt(stdDev / T3(rowCount(A) - 1.0));
#else
    stdDev = sqrt(stdDev / T3(rowCount(A)));
#endif
}

template <class T1, class C1, class T2, class C2, class T3, class C3, class T4, class C4>
void
columnStatistics(const MultiArrayView<2, T1, C1> & A, 
                 MultiArrayView<2, T2, C2> & mean, MultiArrayView<2, T3, C3> & stdDev, MultiArrayView<2, T4, C4> & norm)
{
    unsigned int m = rowCount(A);
    unsigned int n = columnCount(A);
    vigra_precondition(1 == rowCount(mean) && n == columnCount(mean) &&
                       1 == rowCount(stdDev) && n == columnCount(stdDev) &&
                       1 == rowCount(norm) && n == columnCount(norm),
                       "columnStatistics(): Shape mismatch between input and output.");

    detail::columnStatisticsImpl(A, mean, stdDev);
    norm = sqrt(stdDev + T2(m) * sq(mean));
#if 0
    stdDev = sqrt(stdDev / T3(m - 1.0));
#else
    stdDev = sqrt(stdDev / T3(m));
#endif
}

template <class T1, class C1, class T2, class C2>
void
rowStatistics(const MultiArrayView<2, T1, C1> & A, 
                 MultiArrayView<2, T2, C2> & mean)
{
    vigra_precondition(1 == columnCount(mean) && rowCount(A) == rowCount(mean),
                       "rowStatistics(): Shape mismatch between input and output.");
    MultiArrayView<2, T2, StridedArrayTag> tm = transpose(mean);
    columnStatistics(transpose(A), tm);
}

template <class T1, class C1, class T2, class C2, class T3, class C3>
void
rowStatistics(const MultiArrayView<2, T1, C1> & A, 
                 MultiArrayView<2, T2, C2> & mean, MultiArrayView<2, T3, C3> & stdDev)
{
    vigra_precondition(1 == columnCount(mean) && rowCount(A) == rowCount(mean) &&
                       1 == columnCount(stdDev) && rowCount(A) == rowCount(stdDev),
                       "rowStatistics(): Shape mismatch between input and output.");
    MultiArrayView<2, T2, StridedArrayTag> tm = transpose(mean);
    MultiArrayView<2, T3, StridedArrayTag> ts = transpose(stdDev);
    columnStatistics(transpose(A), tm, ts);
}

template <class T1, class C1, class T2, class C2, class T3, class C3, class T4, class C4>
void
rowStatistics(const MultiArrayView<2, T1, C1> & A, 
                 MultiArrayView<2, T2, C2> & mean, MultiArrayView<2, T3, C3> & stdDev, MultiArrayView<2, T4, C4> & norm)
{
    vigra_precondition(1 == columnCount(mean) && rowCount(A) == rowCount(mean) &&
                       1 == columnCount(stdDev) && rowCount(A) == rowCount(stdDev) &&
                       1 == columnCount(norm) && rowCount(A) == rowCount(norm),
                       "rowStatistics(): Shape mismatch between input and output.");
    MultiArrayView<2, T2, StridedArrayTag> tm = transpose(mean);
    MultiArrayView<2, T3, StridedArrayTag> ts = transpose(stdDev); 
    MultiArrayView<2, T4, StridedArrayTag> tn = transpose(norm);
    columnStatistics(transpose(A), tm, ts, tn);
}

template <class T1, class C1, class U, class T2, class C2, class T3, class C3>
void updateCovarianceMatrix(MultiArrayView<2, T1, C1> const & features,
                       U & count, MultiArrayView<2, T2, C2> & mean, MultiArrayView<2, T3, C3> & covariance)
{
    unsigned int n = std::max(rowCount(features), columnCount(features));
    vigra_precondition(std::min(rowCount(features), columnCount(features)) == 1,
          "updateCovarianceMatrix(): Features must be a row or column vector.");
    vigra_precondition(mean.shape() == features.shape(),
          "updateCovarianceMatrix(): Shape mismatch between feature vector and mean vector.");
    vigra_precondition(n == rowCount(covariance) && n == columnCount(covariance),
          "updateCovarianceMatrix(): Shape mismatch between feature vector and covariance matrix.");
    
    // West's algorithm for incremental covariance matrix computation
    Matrix<T2> t = features - mean;
    ++count;
    double f  = 1.0 / count,
           f1 = 1.0 - f;
    mean += f*t;
    
    if(rowCount(features) == 1)
    {
        for(unsigned int k=0; k<n; ++k)
        {
            covariance(k, k) += f*sq(t(0, k));
            for(unsigned int l=k+1; l<n; ++l)
            {
                covariance(k, l) += f*t(0, k)*t(0, l);
                covariance(l, k) = covariance(k, l);
            }
        }
    }
    else
    {
        for(unsigned int k=0; k<n; ++k)
        {
            covariance(k, k) += f*sq(t(k, 0));
            for(unsigned int l=k+1; l<n; ++l)
            {
                covariance(k, l) += f*t(k, 0)*t(l, 0);
                covariance(l, k) = covariance(k, l);
            }
        }
    }
}

template <class T1, class C1, class T2, class C2>
void covarianceMatrixOfColumns(MultiArrayView<2, T1, C1> const & features,
                               MultiArrayView<2, T2, C2> & covariance)
{
    unsigned int m = rowCount(features), n = columnCount(features);
    vigra_precondition(n == rowCount(covariance) && n == columnCount(covariance),
          "covarianceMatrixOfColumns(): Shape mismatch between feature matrix and covariance matrix.");
    unsigned int count = 0;
    Matrix<T2> means(1, n);
    for(unsigned int k=0; k<m; ++k)
        updateCovarianceMatrix(rowVector(features, k), count, means, covariance);
}

template <class T, class C>
TemporaryMatrix<T> 
covarianceMatrixOfColumns(MultiArrayView<2, T, C> const & features)
{
    TemporaryMatrix<T> res(columnCount(features), columnCount(features));
    covarianceMatrixOfColumns(features, res);
    return res;
}

template <class T1, class C1, class T2, class C2>
void covarianceMatrixOfRows(MultiArrayView<2, T1, C1> const & features,
                            MultiArrayView<2, T2, C2> & covariance)
{
    unsigned int m = rowCount(features), n = columnCount(features);
    vigra_precondition(m == rowCount(covariance) && m == columnCount(covariance),
          "covarianceMatrixOfRows(): Shape mismatch between feature matrix and covariance matrix.");
    unsigned int count = 0;
    Matrix<T2> means(m, 1);
    for(unsigned int k=0; k<n; ++k)
        updateCovarianceMatrix(columnVector(features, k), count, means, covariance);
}

template <class T, class C>
TemporaryMatrix<T> 
covarianceMatrixOfRows(MultiArrayView<2, T, C> const & features)
{
    TemporaryMatrix<T> res(rowCount(features), rowCount(features));
    covarianceMatrixOfRows(features, res);
    return res;
}

enum DataPreparationGoals { ZeroMean = 1, UnitVariance = 2, UnitNorm = 4 };

inline DataPreparationGoals operator|(DataPreparationGoals l, DataPreparationGoals r)
{
    return DataPreparationGoals(int(l) | int(r));
}

namespace detail {

template <class T, class C1, class C2, class C3, class C4>
void
prepareDataImpl(const MultiArrayView<2, T, C1> & A, 
               MultiArrayView<2, T, C2> & res, MultiArrayView<2, T, C3> & offset, MultiArrayView<2, T, C4> & scaling, 
               DataPreparationGoals goals)
{
    unsigned int m = rowCount(A);
    unsigned int n = columnCount(A);
    vigra_precondition(A.shape() == res.shape() && 
                       n == columnCount(offset) && 1 == rowCount(offset) &&
                       n == columnCount(scaling) && 1 == rowCount(scaling),
                       "prepareDataImpl(): Shape mismatch between input and output.");

    if(!goals)
    {
        res = A;
        offset.init(0.0);
        scaling.init(1.0);
        return;
    }
    
    bool zeroMean = (goals & ZeroMean) != 0;
    bool unitVariance = (goals & UnitVariance) != 0;
    bool unitNorm = (goals & UnitNorm) != 0;

    Matrix<T> mean(1, n), sumOfSquaredDifferences(1, n);
    detail::columnStatisticsImpl(A, mean, sumOfSquaredDifferences);
    
    for(unsigned int k=0; k<n; ++k)
    {
        T stdDev = std::sqrt(sumOfSquaredDifferences(0, k) / T(m));
        T norm = std::sqrt(sumOfSquaredDifferences(0, k) + T(m) * sq(mean(0,k)));
        if(zeroMean && stdDev > 0.0)
        {
            columnVector(res, k) = columnVector(A, k) - mean(0,k);
            offset(0, k) = mean(0, k);
            mean(0, k) = 0.0;
        }
        else
        {
            columnVector(res, k) = columnVector(A, k);
            if(closeAtTolerance(mean(0,k), T(0.0)))
                mean(0,k) = 0.0;
            offset(0, k) = mean(0, k);
        }
        vigra_precondition(!(unitVariance && unitNorm && mean(0,k) != 0.0),
            "prepareDataImpl(): If mean is non-zero, unit variance and unit norm cannot be achieved at the same time.");
        if(unitNorm)
        {
            if(mean(0,k) != 0.0)
            {
                columnVector(res, k) /= norm;
                scaling(0, k) = 1.0 / norm;
            }
            else if(stdDev != 0.0)
            {
                columnVector(res, k) /= stdDev;
                scaling(0, k) = 1.0 / stdDev;
            }
            else
            {
                scaling(0, k) = 1.0;
            }
        }
        else if(unitVariance && stdDev != 0.0)
        {
            columnVector(res, k) /= stdDev;
            scaling(0, k) = 1.0 / stdDev;
        }
        else
        {
            scaling(0, k) = 1.0;
        }
    }
}

} // namespace detail

template <class T, class C1, class C2, class C3, class C4>
inline void
prepareColumns(const MultiArrayView<2, T, C1> & A, 
               MultiArrayView<2, T, C2> & res, MultiArrayView<2, T, C3> & offset, MultiArrayView<2, T, C4> & scaling, 
               DataPreparationGoals goals = ZeroMean | UnitVariance)
{
    detail::prepareDataImpl(A, res, offset, scaling, goals);
}

template <class T, class C1, class C2>
inline void
prepareColumns(const MultiArrayView<2, T, C1> & A, MultiArrayView<2, T, C2> & res, 
               DataPreparationGoals goals = ZeroMean | UnitVariance)
{
    Matrix<T> offset(1, columnCount(A)), scaling(1, columnCount(A));
    detail::prepareDataImpl(A, res, offset, scaling, goals);
}

template <class T, class C1, class C2, class C3, class C4>
inline void
prepareRows(const MultiArrayView<2, T, C1> & A, 
               MultiArrayView<2, T, C2> & res, MultiArrayView<2, T, C3> & offset, MultiArrayView<2, T, C4> & scaling, 
               DataPreparationGoals goals = ZeroMean | UnitVariance)
{
    MultiArrayView<2, T, StridedArrayTag> tr = transpose(res), to = transpose(offset), ts = transpose(scaling);
    detail::prepareDataImpl(transpose(A), tr, to, ts, goals);
}

template <class T, class C1, class C2>
inline void
prepareRows(const MultiArrayView<2, T, C1> & A, MultiArrayView<2, T, C2> & res, 
               DataPreparationGoals goals = ZeroMean | UnitVariance)
{
    MultiArrayView<2, T, StridedArrayTag> tr = transpose(res);
    Matrix<T> offset(rowCount(A), 1), scaling(rowCount(A), 1);
    detail::prepareDataImpl(transpose(A), tr, offset, scaling, goals);
}

//@}

} // namespace linalg

using linalg::columnStatistics;
using linalg::prepareColumns;
using linalg::rowStatistics;
using linalg::prepareRows;

}  // namespace vigra



#endif // VIGRA_MATRIX_HXX
