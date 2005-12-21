/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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


#ifndef VIGRA_TINYVECTOR_HXX
#define VIGRA_TINYVECTOR_HXX

#include <cmath>    // abs(double)
#include <cstdlib>  // abs(int)
#include <iosfwd>   // ostream
#include "vigra/config.hxx"
#include "vigra/error.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/mathutil.hxx"

namespace vigra {

using VIGRA_CSTD::abs;
using VIGRA_CSTD::ceil;
using VIGRA_CSTD::floor;


template <class V1, int SIZE, class D1, class D2>
class TinyVectorBase;

template <class V1, int SIZE, class D1, class D2>
inline
typename TinyVectorBase<V1, SIZE, D1, D2>::SquaredNormType
squaredNorm(TinyVectorBase<V1, SIZE, D1, D2> const & t);


namespace detail {

#define VIGRA_EXEC_LOOP(NAME, OPER) \
    template <class T1, class T2>  \
    static void NAME(T1 * left, T2 const * right)  \
    {  \
        for(int i=0; i<LEVEL; ++i)  \
            (left[i]) OPER (right[i]);  \
    }

#define VIGRA_EXEC_LOOP_SCALAR(NAME, OPER) \
    template <class T1, class T2>  \
    static void NAME(T1 * left, T2 right)  \
    {  \
        for(int i=0; i<LEVEL; ++i)  \
            (left[i]) OPER (right);  \
    }

template <int LEVEL>
struct ExecLoop
{
    template <class T1, class T2>
    static void assignCast(T1 * left, T2 const * right)
    {
        for(int i=0; i<LEVEL; ++i)
            left[i] = detail::RequiresExplicitCast<T1>::cast(right[i]);
    }

    VIGRA_EXEC_LOOP(assign, =)
    VIGRA_EXEC_LOOP(add, +=)
    VIGRA_EXEC_LOOP(sub, -=)
    VIGRA_EXEC_LOOP(mul, *=)
    VIGRA_EXEC_LOOP(neg, = -)
    VIGRA_EXEC_LOOP(abs, = vigra::abs)
    VIGRA_EXEC_LOOP(floor, = vigra::floor)
    VIGRA_EXEC_LOOP(ceil, = vigra::ceil)
    VIGRA_EXEC_LOOP(fromPromote, = NumericTraits<T1>::fromPromote)
    VIGRA_EXEC_LOOP(fromRealPromote, = NumericTraits<T1>::fromRealPromote)
    VIGRA_EXEC_LOOP_SCALAR(assignScalar, =)
    VIGRA_EXEC_LOOP_SCALAR(mulScalar, *=)
    VIGRA_EXEC_LOOP_SCALAR(divScalar, /=)

    template <class T1, class T2>
    static bool notEqual(T1 const * left, T2 const * right)
    {
        for(int i=0; i<LEVEL; ++i)
            if(left[i] != right[i])
                return true;
        return false;
    }

    template <class T>
    static typename NumericTraits<T>::Promote
    dot(T const * d)
    {
        typename NumericTraits<T>::Promote  res(*d * *d);
        for(int i=1; i<LEVEL; ++i)
            res += d[i] * d[i];
        return res;
    }

    template <class T1, class T2>
    static typename PromoteTraits<T1, T2>::Promote
    dot(T1 const * left, T2 const * right)
    {
        typename PromoteTraits<T1, T2>::Promote res(*left * *right);
        for(int i=1; i<LEVEL; ++i)
            res += left[i] * right[i];
        return res;
    }

    template <class T>
    static typename NormTraits<T>::SquaredNormType
    squaredNorm(T const * d)
    {
        typename NormTraits<T>::SquaredNormType  res = vigra::squaredNorm(*d);
        for(int i=1; i<LEVEL; ++i)
            res += vigra::squaredNorm(d[i]);
        return res;
    }
};

template <int LEVEL>
struct UnrollDot
{
    template <class T>
    static typename NumericTraits<T>::Promote
    dot(T const * d)
    {
        return *d * *d + UnrollDot<LEVEL-1>::dot(d+1);
    }

    template <class T1, class T2>
    static typename PromoteTraits<T1, T2>::Promote
    dot(T1 const * left, T2 const * right)
    {
        return *left * *right + UnrollDot<LEVEL-1>::dot(left+1, right+1);
    }
};

template <>
struct UnrollDot<1>
{
    template <class T>
    static typename NumericTraits<T>::Promote
    dot(T const * d)
    {
        return *d * *d ;
    }

    template <class T1, class T2>
    static typename PromoteTraits<T1, T2>::Promote
    dot(T1 const * left, T2 const * right)
    {
        return *left * *right;
    }
};

template <int LEVEL>
struct UnrollSquaredNorm
{
    template <class T>
    static typename NormTraits<T>::SquaredNormType
    squaredNorm(T const * d)
    {
        return vigra::squaredNorm(*d) + UnrollSquaredNorm<LEVEL-1>::squaredNorm(d+1);
    }
};

template <>
struct UnrollSquaredNorm<1>
{
    template <class T>
    static typename NormTraits<T>::SquaredNormType
    squaredNorm(T const * d)
    {
        return vigra::squaredNorm(*d);
    }
};

#undef VIGRA_EXEC_LOOP
#undef VIGRA_EXEC_LOOP_SCALAR

#define VIGRA_UNROLL_LOOP(NAME, OPER) \
    template <class T1, class T2>  \
    static void NAME(T1 * left, T2 const * right)  \
    {  \
        (*left) OPER (*right);  \
        UnrollLoop<LEVEL-1>::NAME(left+1, right+1); \
    }

#define VIGRA_UNROLL_LOOP_SCALAR(NAME, OPER) \
    template <class T1, class T2>  \
    static void NAME(T1 * left, T2 right)  \
    {  \
        (*left) OPER (right);  \
        UnrollLoop<LEVEL-1>::NAME(left+1, right); \
    }


template <int LEVEL>
struct UnrollLoop
{
    template <class T1, class T2>
    static void assignCast(T1 * left, T2 const * right)
    {
        *left = detail::RequiresExplicitCast<T1>::cast(*right);
        UnrollLoop<LEVEL-1>::assignCast(left+1, right+1);
    }

    VIGRA_UNROLL_LOOP(assign, =)
    VIGRA_UNROLL_LOOP(add, +=)
    VIGRA_UNROLL_LOOP(sub, -=)
    VIGRA_UNROLL_LOOP(mul, *=)
    VIGRA_UNROLL_LOOP(neg, = -)
    VIGRA_UNROLL_LOOP(abs, = vigra::abs)
    VIGRA_UNROLL_LOOP(floor, = vigra::floor)
    VIGRA_UNROLL_LOOP(ceil, = vigra::ceil)
    VIGRA_UNROLL_LOOP(fromPromote, = NumericTraits<T1>::fromPromote)
    VIGRA_UNROLL_LOOP(fromRealPromote, = NumericTraits<T1>::fromRealPromote)
    VIGRA_UNROLL_LOOP_SCALAR(assignScalar, =)
    VIGRA_UNROLL_LOOP_SCALAR(mulScalar, *=)
    VIGRA_UNROLL_LOOP_SCALAR(divScalar, /=)

    template <class T1, class T2>
    static bool notEqual(T1 const * left, T2 const * right)
    {
        return (*left != *right) || UnrollLoop<LEVEL - 1>::notEqual(left+1, right+1);
    }

    template <class T>
    static typename NumericTraits<T>::Promote
    dot(T const * d)
    {
        return UnrollDot<LEVEL>::dot(d);
    }

    template <class T1, class T2>
    static typename PromoteTraits<T1, T2>::Promote
    dot(T1 const * left, T2 const * right)
    {
        return UnrollDot<LEVEL>::dot(left, right);
    }

    template <class T>
    static typename NormTraits<T>::SquaredNormType
    squaredNorm(T const * d)
    {
        return UnrollSquaredNorm<LEVEL>::squaredNorm(d);
    }
};

#undef VIGRA_UNROLL_LOOP
#undef VIGRA_UNROLL_LOOP_SCALAR

template <>
struct UnrollLoop<0>
{
    template <class T1, class T2>
    static void assignCast(T1, T2) {}
    template <class T1, class T2>
    static void assign(T1, T2) {}
    template <class T1, class T2>
    static void assignScalar(T1, T2) {}
    template <class T1, class T2>
    static void add(T1, T2) {}
    template <class T1, class T2>
    static void sub(T1, T2) {}
    template <class T1, class T2>
    static void mul(T1, T2) {}
    template <class T1, class T2>
    static void mulScalar(T1, T2) {}
    template <class T1, class T2>
    static void div(T1, T2) {}
    template <class T1, class T2>
    static void divScalar(T1, T2) {}
    template <class T1, class T2>
    static void fromPromote(T1, T2) {}
    template <class T1, class T2>
    static void fromRealPromote(T1, T2) {}
    template <class T1, class T2>
    static void neg(T1, T2) {}
    template <class T1, class T2>
    static void abs(T1, T2) {}
    template <class T1, class T2>
    static void floor(T1, T2) {}
    template <class T1, class T2>
    static void ceil(T1, T2) {}
    template <class T1, class T2>
    static bool notEqual(T1, T2) { return false; }
};

template <bool PREDICATE>
struct TinyVectorIf
{
    template <class T, class F>
    struct res
    {
        typedef T type;
    };
};

template <>
struct TinyVectorIf<false>
{
    template <class T, class F>
    struct res
    {
        typedef F type;
    };
};

template <int SIZE>
struct LoopType
{
    typedef typename TinyVectorIf<SIZE < 5>::
            template res<UnrollLoop<SIZE>, ExecLoop<SIZE> >::type type;
};

struct DontInit {};

inline DontInit dontInit() {return DontInit(); }

} // namespace detail

template <class T, int SIZE>
class TinyVector;

template <class T, int SIZE>
class TinyVectorView;

/********************************************************/
/*                                                      */
/*                    TinyVectorBase                    */
/*                                                      */
/********************************************************/

/** \brief Base class for fixed size vectors.

    This class contains functionality shared by
    \ref TinyVector and \ref TinyVectorView, and enables these classes
    to be freely mixed within expressions. It is typically not used directly.

    <b>\#include</b> "<a href="tinyvector_8hxx-source.html">vigra/tinyvector.hxx</a>"<br>
    Namespace: vigra
**/
template <class VALUETYPE, int SIZE, class DATA, class DERIVED>
class TinyVectorBase
{
    TinyVectorBase(TinyVectorBase const &); // do not use

    TinyVectorBase & operator=(TinyVectorBase const & other); // do not use

  protected:

    typedef typename detail::LoopType<SIZE>::type Loop;

    TinyVectorBase()
    {}

  public:
        /** STL-compatible definition of valuetype
        */
    typedef VALUETYPE value_type;

        /** reference (return of operator[]).
        */
    typedef VALUETYPE & reference;

        /** const reference (return of operator[] const).
        */
    typedef VALUETYPE const & const_reference;

        /** pointer (return of operator->).
        */
    typedef VALUETYPE * pointer;

        /** const pointer (return of operator-> const).
        */
    typedef VALUETYPE const * const_pointer;

        /** STL-compatible definition of iterator
        */
    typedef value_type * iterator;

        /** STL-compatible definition of const iterator
        */
    typedef value_type const * const_iterator;

        /** STL-compatible definition of size_type
        */
    typedef unsigned int size_type;

        /** STL-compatible definition of difference_type
        */
    typedef int difference_type;

        /** the scalar type for the outer product
        */
    typedef double scalar_multiplier;

        /** the vector's squared norm type
        */
    typedef typename NormTraits<VALUETYPE>::SquaredNormType SquaredNormType;

        /** the vector's norm type
        */
    typedef typename SquareRootTraits<SquaredNormType>::SquareRootResult NormType;

        /** the vector's size
        */
    enum { static_size = SIZE };

        /** Initialize from another sequence (must have length SIZE!)
        */
    template <class Iterator>
    void init(Iterator i, Iterator end)
    {
		vigra_precondition(end-i == SIZE,
            "TinyVector::init(): Sequence has wrong size.");
        Loop::assignCast(data_, i);
    }

        /** Component-wise add-assignment
        */
    template <class T1, class D1, class D2>
    DERIVED & operator+=(TinyVectorBase<T1, SIZE, D1, D2> const & r)
    {
        Loop::add(data_, r.begin());
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise subtract-assignment
        */
    template <class T1, class D1, class D2>
    DERIVED & operator-=(TinyVectorBase<T1, SIZE, D1, D2> const & r)
    {
        Loop::sub(data_, r.begin());
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise multiply-assignment
        */
    template <class T1, class D1, class D2>
    DERIVED & operator*=(TinyVectorBase<T1, SIZE, D1, D2> const & r)
    {
        Loop::mul(data_, r.begin());
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise scalar multiply-assignment
        */
    DERIVED & operator*=(double r)
    {
        Loop::mulScalar(data_, r);
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise scalar divide-assignment
        */
    DERIVED & operator/=(double r)
    {
        Loop::divScalar(data_, r);
        return static_cast<DERIVED &>(*this);
    }

        /** Calculate magnitude.
        */
    NormType magnitude() const
    {
         return sqrt(static_cast<typename 
              SquareRootTraits<SquaredNormType>::SquareRootArgument>(squaredMagnitude()));
    }

        /** Calculate squared magnitude.
        */
    SquaredNormType squaredMagnitude() const
    {
        return Loop::squaredNorm(data_);
    }

        /** Access component by index.
        */
    reference operator[](difference_type i) { return data_[i]; }

        /** Get component by index.
        */
    const_reference operator[](difference_type i) const { return data_[i]; }

        /** Get random access iterator to begin of vector.
        */
    iterator begin() { return data_; }
        /** Get random access iterator past-the-end of vector.
        */
    iterator end() { return data_ + SIZE; }

        /** Get const random access iterator to begin of vector.
        */
    const_iterator begin() const { return data_; }

        /** Get const random access iterator past-the-end of vector.
        */
    const_iterator end() const { return data_ + SIZE; }

        /** Size of TinyVector vector always equals the template parameter SIZE.
        */
    size_type size() const { return SIZE; }

    pointer data() { return data_; }

    const_pointer data() const { return data_; }


  protected:
    DATA data_;
};

/** \brief Class for fixed size vectors.

    This class contains an array of size SIZE of the specified VALUETYPE.
    The interface conforms to STL vector, except that there are no functions
    that change the size of a TinyVector.

    \ref TinyVectorOperators "Arithmetic operations"
    on TinyVectors are defined as component-wise applications of these
    operations. Addition and subtraction of two TinyVectors
    (+=, -=, +, -, unary -), multiplication and division of an
    TinyVector with a double, and NumericTraits/PromoteTraits are defined,
    so that TinyVector fulfills the requirements of \ref LinearAlgebra.

    VIGRA algorithms typically use \ref vigra::VectorAccessor to access
    TinyVectors as a whole, or specific components of them.

    See also:<br>
    <DL>
        <DT>
            <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
            \ref vigra::TinyVectorBase
            <DD>
        <DT>
            <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
            \ref vigra::TinyVectorView
            <DD>
        <DT>
            <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
            \ref TinyVectorTraits
            <DD>
        <DT>
            <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
            \ref TinyVectorOperators
            <DD>
    </DL>

    <b>\#include</b> "<a href="tinyvector_8hxx-source.html">vigra/tinyvector.hxx</a>"<br>
    Namespace: vigra
**/
template <class T, int SIZE>
class TinyVector
: public TinyVectorBase<T, SIZE, T[SIZE], TinyVector<T, SIZE> >
{
    typedef TinyVectorBase<T, SIZE, T[SIZE], TinyVector<T, SIZE> > BaseType;
    typedef typename BaseType::Loop Loop;

  public:

    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::reference reference;
    typedef typename BaseType::const_reference const_reference;
    typedef typename BaseType::pointer pointer;
    typedef typename BaseType::const_pointer const_pointer;
    typedef typename BaseType::iterator iterator;
    typedef typename BaseType::const_iterator const_iterator;
    typedef typename BaseType::size_type size_type;
    typedef typename BaseType::difference_type difference_type;
    typedef typename BaseType::scalar_multiplier scalar_multiplier;
    typedef typename BaseType::SquaredNormType SquaredNormType;
    typedef typename BaseType::NormType NormType;

        /** Construction with constant value
        */
    explicit TinyVector(value_type const & initial)
    : BaseType()
    {
        Loop::assignScalar(BaseType::begin(), initial);
    }

        /** Construction with explicit values.
            Call only if SIZE == 2
        */
    TinyVector(value_type const & i1, value_type const & i2)
    : BaseType()
    {
        BaseType::data_[0] = i1;
        BaseType::data_[1] = i2;
    }

        /** Construction with explicit values.
            Call only if SIZE == 3
        */
    TinyVector(value_type const & i1, value_type const & i2, value_type const & i3)
    : BaseType()
    {
        BaseType::data_[0] = i1;
        BaseType::data_[1] = i2;
        BaseType::data_[2] = i3;
    }

        /** Construction with explicit values.
            Call only if SIZE == 4
        */
    TinyVector(value_type const & i1, value_type const & i2,
               value_type const & i3, value_type const & i4)
    : BaseType()
    {
        BaseType::data_[0] = i1;
        BaseType::data_[1] = i2;
        BaseType::data_[2] = i3;
        BaseType::data_[3] = i4;
    }

       /** Default constructor (initializes all components with zero)
        */
    TinyVector()
    : BaseType()
    {
        Loop::assignScalar(BaseType::data_, NumericTraits<value_type>::zero());
    }

        /** Copy constructor.
        */
    TinyVector(TinyVector const & r)
    : BaseType()
    {
        Loop::assign(BaseType::data_, r.data_);
    }

        /** Constructor from C array.
        */
    explicit TinyVector(const_pointer data)
    : BaseType()
    {
        Loop::assign(BaseType::data_, data);
    }

        /** Copy assignment.
        */
    TinyVector & operator=(TinyVector const & r)
    {
        Loop::assign(BaseType::data_, r.data_);
        return *this;
    }

        /** Copy with type conversion.
        */
    template <class U, class DATA, class DERIVED>
    TinyVector(TinyVectorBase<U, SIZE, DATA, DERIVED> const & r)
    : BaseType()
    {
        Loop::assignCast(BaseType::data_, r.begin());
    }

        /** Copy assignment with type conversion.
        */
    template <class U, class DATA, class DERIVED>
    TinyVector & operator=(TinyVectorBase<U, SIZE, DATA, DERIVED> const & r)
    {
        Loop::assignCast(BaseType::data_, r.begin());
        return *this;
    }

    explicit TinyVector(detail::DontInit)
    : BaseType()
    {}
};

/** \brief Wrapper for fixed size vectors.

    This class wraps an array of size SIZE of the specified VALUETYPE.
    Thus, the array can be accessed with an interface similar to
    that of std::vector (except that there are no functions
    that change the size of a TinyVectorView). The TinyVectorView
    does <em>not</em> assume ownership of the given memory.

    \ref TinyVectorOperators "Arithmetic operations"
    on TinyVectorViews are defined as component-wise applications of these
    operations. Addition and subtraction of two TinyVectorViews
    (+=, -=, +, -, unary -), multiplication and division of an
    TinyVectorViews with a double, and NumericTraits/PromoteTraits are defined,
    so that TinyVectorView fulfills the requirements of \ref LinearAlgebra.

    VIGRA algorithms typically use \ref vigra::VectorAccessor to access
    TinyVectorViews as a whole, or specific components of them.

    <b>See also:</b>
    <ul>
        <li> \ref vigra::TinyVectorBase
        <li> \ref vigra::TinyVector
        <li> \ref TinyVectorTraits
        <li> \ref TinyVectorOperators
    </ul>

    <b>\#include</b> "<a href="tinyvector_8hxx-source.html">vigra/tinyvector.hxx</a>"<br>
    Namespace: vigra
**/
template <class T, int SIZE>
class TinyVectorView
: public TinyVectorBase<T, SIZE, T *, TinyVectorView<T, SIZE> >
{
    typedef TinyVectorBase<T, SIZE, T *, TinyVectorView<T, SIZE> > BaseType;
    typedef typename BaseType::Loop Loop;

  public:

    typedef typename BaseType::value_type value_type;
    typedef typename BaseType::reference reference;
    typedef typename BaseType::const_reference const_reference;
    typedef typename BaseType::pointer pointer;
    typedef typename BaseType::const_pointer const_pointer;
    typedef typename BaseType::iterator iterator;
    typedef typename BaseType::const_iterator const_iterator;
    typedef typename BaseType::size_type size_type;
    typedef typename BaseType::difference_type difference_type;
    typedef typename BaseType::scalar_multiplier scalar_multiplier;
    typedef typename BaseType::SquaredNormType SquaredNormType;
    typedef typename BaseType::NormType NormType;

        /** Default constructor
            (pointer to wrapped data is NULL).
        */
    TinyVectorView()
    : BaseType()
    {
        BaseType::data_ = 0;
    }

        /** Construct view for given data array
        */
    TinyVectorView(const_pointer data)
    : BaseType()
    {
        BaseType::data_ = const_cast<pointer>(data);
    }

        /** Copy constructor (shallow copy).
        */
    TinyVectorView(TinyVectorView const & other)
    : BaseType()
    {
        BaseType::data_ = const_cast<pointer>(other.data_);
    }

        /** Construct view from other TinyVector.
        */
    template <class DATA, class DERIVED>
    TinyVectorView(TinyVectorBase<T, SIZE, DATA, DERIVED> const & other)
    : BaseType()
    {
        BaseType::data_ = const_cast<pointer>(other.data());
    }

        /** Copy the data (not the pointer) of the rhs.
        */
   TinyVectorView & operator=(TinyVectorView const & r)
    {
        Loop::assign(BaseType::data_, r.begin());
        return *this;
    }

        /** Copy the data of the rhs with cast.
        */
    template <class U, class DATA, class DERIVED>
    TinyVectorView & operator=(TinyVectorBase<U, SIZE, DATA, DERIVED> const & r)
    {
        Loop::assignCast(BaseType::data_, r.begin());
        return *this;
    }
};

/********************************************************/
/*                                                      */
/*                     TinyVector Comparison            */
/*                                                      */
/********************************************************/

/** \addtogroup TinyVectorOperators Functions for TinyVector

    \brief <b>\#include</b> "<a href="tinyvector_8hxx-source.html">vigra/tinyvector.hxx</a>

    These functions fulfill the requirements of a Linear Space (vector space).
    Return types are determined according to \ref TinyVectorTraits.

    Namespace: vigra
    <p>

 */
//@{
    /// component-wise equal
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline bool
operator==(TinyVectorBase<V1, SIZE, D1, D2> const & l,
           TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    return !(l != r);
}

    /// component-wise not equal
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline bool
operator!=(TinyVectorBase<V1, SIZE, D1, D2> const & l,
           TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    return ltype::notEqual(l.begin(), r.begin());
}

/********************************************************/
/*                                                      */
/*                     TinyVector Output                */
/*                                                      */
/********************************************************/

    /// stream output
template <class V1, int SIZE, class DATA, class DERIVED>
std::ostream &
operator<<(std::ostream & out, TinyVectorBase<V1, SIZE, DATA, DERIVED> const & l)
{
    out << "(";
    int i;
    for(i=0; i<SIZE-1; ++i)
        out << l[i] << ", ";
    out << l[i] << ")";
    return out;
}
//@}

/********************************************************/
/*                                                      */
/*                      TinyVector-Traits               */
/*                                                      */
/********************************************************/

/** \page TinyVectorTraits Numeric and Promote Traits of TinyVector
    The numeric and promote traits for TinyVectors follow
    the general specifications for \ref NumericPromotionTraits.
    They are implemented in terms of the traits of the basic types by
    partial template specialization:

    \code

    template <class T, int SIZE>
    struct NumericTraits<TinyVector<T, SIZE> >
    {
        typedef TinyVector<typename NumericTraits<T>::Promote, SIZE> Promote;
        typedef TinyVector<typename NumericTraits<T>::RealPromote, SIZE> RealPromote;

        typedef typename NumericTraits<T>::isIntegral isIntegral;
        typedef VigraFalseType isScalar;

        // etc.
    };

    template <class T, int SIZE>
    struct NormTraits<TinyVector<T, SIZE> >
    {
        typedef TinyVector<T, SIZE> Type;
        typedef typename Type::SquaredNormType    SquaredNormType;
        typedef typename Type::NormType           NormType;
    };

    template <class T1, class T2, SIZE>
    struct PromoteTraits<TinyVector<T1, SIZE>, TinyVector<T2, SIZE> >
    {
        typedef TinyVector<typename PromoteTraits<T1, T2>::Promote, SIZE> Promote;
    };
    \endcode

    <b>\#include</b> "<a href="tinyvector_8hxx-source.html">vigra/tinyvector.hxx</a>"<br>
    Namespace: vigra

    On compilers that don't support pertial template specialization (e.g.
    MS VisualC++), the traits classes are explicitly specialized for
    <TT>TinyVector<VALUETYPE, SIZE></TT> with
    <TT>VALUETYPE = unsigned char | int | float | double</TT> and <TT>SIZE = 2 | 3 | 4</TT>.

*/

#if !defined(NO_PARTIAL_TEMPLATE_SPECIALIZATION)

template <class T, int SIZE>
struct NumericTraits<TinyVector<T, SIZE> >
{
    typedef TinyVector<T, SIZE> Type;
    typedef TinyVector<typename NumericTraits<T>::Promote, SIZE> Promote;
    typedef TinyVector<typename NumericTraits<T>::RealPromote, SIZE> RealPromote;
    typedef TinyVector<typename NumericTraits<T>::ComplexPromote, SIZE> ComplexPromote;
    typedef T ValueType;

    typedef typename NumericTraits<T>::isIntegral isIntegral;
    typedef VigraFalseType isScalar;
    typedef VigraFalseType isOrdered;
    typedef VigraFalseType isComplex;

    static TinyVector<T, SIZE> zero() {
        return TinyVector<T, SIZE>(NumericTraits<T>::zero());
    }
    static TinyVector<T, SIZE> one() {
        return TinyVector<T, SIZE>(NumericTraits<T>::one());
    }
    static TinyVector<T, SIZE> nonZero() {
        return TinyVector<T, SIZE>(NumericTraits<T>::nonZero());
    }

    template <class D1, class D2>
    static Promote toPromote(TinyVectorBase<T, SIZE, D1, D2> const & v)
    {
        return Promote(v);
    }

    template <class D1, class D2>
    static RealPromote toRealPromote(TinyVectorBase<T, SIZE, D1, D2> const & v)
    {
        return RealPromote(v);
    }

    template <class D1, class D2>
    static TinyVector<T, SIZE>
    fromPromote(TinyVectorBase<typename NumericTraits<T>::Promote, SIZE, D1, D2> const & v)
    {
        TinyVector<T, SIZE> res(detail::dontInit());
        typedef typename detail::LoopType<SIZE>::type ltype;
        ltype::fromPromote(res.begin(), v.begin());
        return res;
    }

    template <class D1, class D2>
    static TinyVector<T, SIZE>
    fromRealPromote(TinyVectorBase<typename NumericTraits<T>::RealPromote, SIZE, D1, D2> const & v)
    {
        TinyVector<T, SIZE> res(detail::dontInit());
        typedef typename detail::LoopType<SIZE>::type ltype;
        ltype::fromRealPromote(res.begin(), v.begin());
        return res;
    }
};

template <class T, int SIZE>
struct NumericTraits<TinyVectorView<T, SIZE> >
: public NumericTraits<TinyVector<T, SIZE> >
{
    typedef TinyVector<T, SIZE> Type;
    typedef TinyVector<typename NumericTraits<T>::Promote, SIZE> Promote;
    typedef TinyVector<typename NumericTraits<T>::RealPromote, SIZE> RealPromote;
    typedef TinyVector<typename NumericTraits<T>::ComplexPromote, SIZE> ComplexPromote;
    typedef T ValueType;

    typedef typename NumericTraits<T>::isIntegral isIntegral;
    typedef VigraFalseType isScalar;
    typedef VigraFalseType isOrdered;
    typedef VigraFalseType isComplex;
};

template <class T, int SIZE>
struct NormTraits<TinyVector<T, SIZE> >
{
    typedef TinyVector<T, SIZE> Type;
    typedef typename Type::SquaredNormType    SquaredNormType;
    typedef typename Type::NormType           NormType;
};

template <class T, int SIZE>
struct NormTraits<TinyVectorView<T, SIZE> >
{
    typedef TinyVector<T, SIZE> Type;
    typedef typename Type::SquaredNormType    SquaredNormType;
    typedef typename Type::NormType           NormType;
};

template <class T1, class T2, int SIZE>
struct PromoteTraits<TinyVector<T1, SIZE>, TinyVector<T2, SIZE> >
{
    typedef TinyVector<typename PromoteTraits<T1, T2>::Promote, SIZE> Promote;
};

template <class T1, class T2, int SIZE>
struct PromoteTraits<TinyVectorView<T1, SIZE>, TinyVectorView<T2, SIZE> >
{
    typedef TinyVector<typename PromoteTraits<T1, T2>::Promote, SIZE> Promote;
};

template <class T1, class T2, int SIZE>
struct PromoteTraits<TinyVectorView<T1, SIZE>, TinyVector<T2, SIZE> >
{
    typedef TinyVector<typename PromoteTraits<T1, T2>::Promote, SIZE> Promote;
};

template <class T1, class T2, int SIZE>
struct PromoteTraits<TinyVector<T1, SIZE>, TinyVectorView<T2, SIZE> >
{
    typedef TinyVector<typename PromoteTraits<T1, T2>::Promote, SIZE> Promote;
};

template <class T, int SIZE>
struct PromoteTraits<TinyVector<T, SIZE>, double >
{
    typedef TinyVector<typename NumericTraits<T>::RealPromote, SIZE> Promote;
};

template <class T, int SIZE>
struct PromoteTraits<double, TinyVector<T, SIZE> >
{
    typedef TinyVector<typename NumericTraits<T>::RealPromote, SIZE> Promote;
};

template <class T, int SIZE>
struct PromoteTraits<TinyVectorView<T, SIZE>, double >
{
    typedef TinyVector<typename NumericTraits<T>::RealPromote, SIZE> Promote;
};

template <class T, int SIZE>
struct PromoteTraits<double, TinyVectorView<T, SIZE> >
{
    typedef TinyVector<typename NumericTraits<T>::RealPromote, SIZE> Promote;
};

#else // NO_PARTIAL_TEMPLATE_SPECIALIZATION


#define TINYVECTOR_NUMTRAITS(T, SIZE) \
template<>\
struct NumericTraits<TinyVector<T, SIZE> >\
{\
    typedef TinyVector<T, SIZE> Type;\
    typedef TinyVector<NumericTraits<T>::Promote, SIZE> Promote;\
    typedef TinyVector<NumericTraits<T>::RealPromote, SIZE> RealPromote;\
    typedef TinyVector<NumericTraits<T>::ComplexPromote, SIZE> ComplexPromote;\
    typedef T ValueType; \
    typedef NumericTraits<T>::isIntegral isIntegral;\
    typedef VigraFalseType isScalar;\
    typedef VigraFalseType isOrdered;\
    typedef VigraFalseType isComplex;\
    \
    static TinyVector<T, SIZE> zero() { \
        return TinyVector<T, SIZE>(NumericTraits<T>::zero()); \
    }\
    static TinyVector<T, SIZE> one() { \
        return TinyVector<T, SIZE>(NumericTraits<T>::one()); \
    }\
    static TinyVector<T, SIZE> nonZero() { \
        return TinyVector<T, SIZE>(NumericTraits<T>::nonZero()); \
    }\
    \
    static Promote toPromote(TinyVector<T, SIZE> const & v) { \
        return Promote(v); \
    }\
    static RealPromote toRealPromote(TinyVector<T, SIZE> const & v) { \
        return RealPromote(v); \
    }\
    static TinyVector<T, SIZE> fromPromote(Promote const & v) { \
        TinyVector<T, SIZE> res;\
        TinyVector<T, SIZE>::iterator d = res.begin(), dend = res.end();\
        Promote::const_iterator s = v.begin();\
        for(; d != dend; ++d, ++s)\
            *d = NumericTraits<T>::fromPromote(*s);\
        return res;\
    }\
    static TinyVector<T, SIZE> fromRealPromote(RealPromote const & v) {\
        TinyVector<T, SIZE> res;\
        TinyVector<T, SIZE>::iterator d = res.begin(), dend = res.end();\
        RealPromote::const_iterator s = v.begin();\
        for(; d != dend; ++d, ++s)\
            *d = NumericTraits<T>::fromRealPromote(*s);\
        return res;\
    }\
}; \
template<>\
struct NormTraits<TinyVector<T, SIZE> >\
{\
    typedef TinyVector<T, SIZE> Type;\
    typedef Type::SquaredNormType           SquaredNormType; \
    typedef Type::NormType NormType; \
};

#define TINYVECTOR_PROMTRAITS1(type1, SIZE) \
template<> \
struct PromoteTraits<TinyVector<type1, SIZE>, TinyVector<type1, SIZE> > \
{ \
    typedef TinyVector<PromoteTraits<type1, type1>::Promote, SIZE> Promote; \
    static Promote toPromote(TinyVector<type1, SIZE> const & v) { \
        return static_cast<Promote>(v); } \
};

#define TINYVECTOR_PROMTRAITS2(type1, type2, SIZE) \
template<> \
struct PromoteTraits<TinyVector<type1, SIZE>, TinyVector<type2, SIZE> > \
{ \
    typedef TinyVector<PromoteTraits<type1, type2>::Promote, SIZE> Promote; \
    static Promote toPromote(TinyVector<type1, SIZE> const & v) { \
        return static_cast<Promote>(v); } \
    static Promote toPromote(TinyVector<type2, SIZE> const & v) { \
       return static_cast<Promote>(v); } \
};

#define TINYVECTOR_TRAITS(SIZE) \
TINYVECTOR_NUMTRAITS(unsigned char, SIZE)\
TINYVECTOR_NUMTRAITS(int, SIZE)\
TINYVECTOR_NUMTRAITS(float, SIZE)\
TINYVECTOR_NUMTRAITS(double, SIZE)\
TINYVECTOR_PROMTRAITS1(unsigned char, SIZE)\
TINYVECTOR_PROMTRAITS1(int, SIZE)\
TINYVECTOR_PROMTRAITS1(float, SIZE)\
TINYVECTOR_PROMTRAITS1(double, SIZE)\
TINYVECTOR_PROMTRAITS2(float, unsigned char, SIZE)\
TINYVECTOR_PROMTRAITS2(unsigned char, float, SIZE)\
TINYVECTOR_PROMTRAITS2(int, unsigned char, SIZE)\
TINYVECTOR_PROMTRAITS2(unsigned char, int, SIZE)\
TINYVECTOR_PROMTRAITS2(int, float, SIZE)\
TINYVECTOR_PROMTRAITS2(float, int, SIZE)\
TINYVECTOR_PROMTRAITS2(double, unsigned char, SIZE)\
TINYVECTOR_PROMTRAITS2(unsigned char, double, SIZE)\
TINYVECTOR_PROMTRAITS2(int, double, SIZE)\
TINYVECTOR_PROMTRAITS2(double, int, SIZE)\
TINYVECTOR_PROMTRAITS2(double, float, SIZE)\
TINYVECTOR_PROMTRAITS2(float, double, SIZE)

TINYVECTOR_TRAITS(2)
TINYVECTOR_TRAITS(3)
TINYVECTOR_TRAITS(4)

#undef TINYVECTOR_NUMTRAITS
#undef TINYVECTOR_PROMTRAITS1
#undef TINYVECTOR_PROMTRAITS2
#undef TINYVECTOR_TRAITS

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION


/********************************************************/
/*                                                      */
/*                      TinyVector-Arithmetic           */
/*                                                      */
/********************************************************/

/** \addtogroup TinyVectorOperators
 */
//@{

    /// component-wise addition
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2, SIZE> >::Promote
operator+(TinyVectorBase<V1, SIZE, D1, D2> const & l,
          TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    return typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2 , SIZE> >::Promote(l) += r;
}

    /// component-wise subtraction
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2, SIZE> >::Promote
operator-(TinyVectorBase<V1, SIZE, D1, D2> const & l,
          TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    return typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2 , SIZE> >::Promote(l) -= r;
}

    /// component-wise multiplication
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2, SIZE> >::Promote
operator*(TinyVectorBase<V1, SIZE, D1, D2> const & l,
          TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    return typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2 , SIZE> >::Promote(l) *= r;
}

    /// component-wise left scalar multiplication
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator*(double v, TinyVectorBase<V, SIZE, D1, D2> const & r)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(r) *= v;
}

    /// component-wise right scalar multiplication
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator*(TinyVectorBase<V, SIZE, D1, D2> const & l, double v)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(l) *= v;
}

    /// component-wise scalar division
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator/(TinyVectorBase<V, SIZE, D1, D2> const & l, double v)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(l) /= v;
}


    /** Unary negation (construct TinyVector with negative values)
    */
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<V, SIZE>
operator-(TinyVectorBase<V, SIZE, D1, D2> const & v)
{
    TinyVector<V, SIZE> res(detail::dontInit());
    typedef typename detail::LoopType<SIZE>::type ltype;
    ltype::neg(res.begin(), v.begin());
    return res;
}

    /// component-wise absolute value
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<V, SIZE>
abs(TinyVectorBase<V, SIZE, D1, D2> const & v)
{
    TinyVector<V, SIZE> res(detail::dontInit());
    typedef typename detail::LoopType<SIZE>::type ltype;
    ltype::abs(res.begin(), v.begin());
    return res;
}

    /** Apply ceil() function to each vector component.
    */
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<V, SIZE>
ceil(TinyVectorBase<V, SIZE, D1, D2> const & v)
{
    TinyVector<V, SIZE> res(detail::dontInit());
    typedef typename detail::LoopType<SIZE>::type ltype;
    ltype::ceil(res.begin(), v.begin());
    return res;
}

    /** Apply floor() function to each vector component.
    */
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<V, SIZE>
floor(TinyVectorBase<V, SIZE, D1, D2> const & v)
{
    TinyVector<V, SIZE> res(detail::dontInit());
    typedef typename detail::LoopType<SIZE>::type ltype;
    ltype::floor(res.begin(), v.begin());
    return res;
}

    /// cross product
template <class V1, class D1, class D2, class V2, class D3, class D4>
inline
TinyVector<typename PromoteTraits<V1, V2>::Promote, 3>
cross(TinyVectorBase<V1, 3, D1, D2> const & r1,
      TinyVectorBase<V2, 3, D3, D4> const & r2)
{
    typedef TinyVector<typename PromoteTraits<V1, V2>::Promote, 3>
            Res;
    return  Res(r1[1]*r2[2] - r1[2]*r2[1],
                r1[2]*r2[0] - r1[0]*r2[2],
                r1[0]*r2[1] - r1[1]*r2[0]);
}

    /// dot product
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
typename PromoteTraits<V1, V2>::Promote
dot(TinyVectorBase<V1, SIZE, D1, D2> const & l,
    TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    return ltype::dot(l.begin(), r.begin());
}


    /// squared norm
template <class V1, int SIZE, class D1, class D2>
inline
typename TinyVectorBase<V1, SIZE, D1, D2>::SquaredNormType
squaredNorm(TinyVectorBase<V1, SIZE, D1, D2> const & t)
{
    return t.squaredMagnitude();
}
//@}


} // namespace vigra

#endif // VIGRA_TINYVECTOR_HXX
