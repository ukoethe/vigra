/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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


#ifndef VIGRA_TINYVECTOR_HXX
#define VIGRA_TINYVECTOR_HXX

namespace lemon {

struct Invalid;

} // namespace lemon

#include <cmath>    // abs(double)
#include <cstdlib>  // abs(int)
#include <iosfwd>   // ostream
#include <algorithm>
#include "config.hxx"
#include "error.hxx"
#include "metaprogramming.hxx"
#include "numerictraits.hxx"
#include "memory.hxx"
#include "mathutil.hxx"
#include "diff2d.hxx"

#ifdef VIGRA_CHECK_BOUNDS
#define VIGRA_ASSERT_INSIDE(diff) \
  vigra_precondition(diff >= 0, "Index out of bounds");\
  vigra_precondition(diff < SIZE, "Index out of bounds");
#else
#define VIGRA_ASSERT_INSIDE(diff)
#endif

namespace vigra {

// mask cl.exe shortcomings [begin]
#if defined(_MSC_VER)
#pragma warning( push )
#pragma warning( disable : 4503 )
#endif

using VIGRA_CSTD::abs;
using VIGRA_CSTD::ceil;
using VIGRA_CSTD::floor;
using VIGRA_CSTD::sqrt;


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

#define VIGRA_EXEC_LOOP_MINMAX(NAME, OPER) \
    template <class T1, class T2>  \
    static void NAME(T1 * left, T2 const * right)  \
    {  \
        for(int i=0; i<LEVEL; ++i)  \
            if(left[i] OPER right[i]) \
                left[i] = right[i];  \
    }

#define VIGRA_EXEC_LOOP_SCALAR(NAME, OPER) \
    template <class T1, class T2>  \
    static void NAME(T1 * left, T2 right)  \
    {  \
        for(int i=0; i<LEVEL; ++i)  \
            (left[i]) = detail::RequiresExplicitCast<T1>::cast((left[i]) OPER (right));  \
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

    template <class T1, class T2>
    static void reverseAssign(T1 * left, T2 const * right)
    {
        for(int i=0; i<LEVEL; ++i)
            left[i] = right[-i];
    }

    template <class T1, class T2>
    static void assignScalar(T1 * left, T2 right)
    {
        for(int i=0; i<LEVEL; ++i)
            left[i] = detail::RequiresExplicitCast<T1>::cast(right);
    }

    template <class T1, class T2>
    static void power(T1 * left, T2 right)
    {
        for(int i=0; i<LEVEL; ++i)
            left[i] = detail::RequiresExplicitCast<T1>::cast(pow(left, right));
    }

    VIGRA_EXEC_LOOP(assign, =)
    VIGRA_EXEC_LOOP(add, +=)
    VIGRA_EXEC_LOOP(sub, -=)
    VIGRA_EXEC_LOOP(mul, *=)
    VIGRA_EXEC_LOOP(div, /=)
    VIGRA_EXEC_LOOP(mod, %=)
    VIGRA_EXEC_LOOP(neg, = -)
    VIGRA_EXEC_LOOP(abs, = vigra::abs)
    VIGRA_EXEC_LOOP(floor, = vigra::floor)
    VIGRA_EXEC_LOOP(ceil, = vigra::ceil)
    VIGRA_EXEC_LOOP(round, = vigra::round)
    VIGRA_EXEC_LOOP(sqrt, = vigra::sqrt)
    VIGRA_EXEC_LOOP(fromPromote, = NumericTraits<T1>::fromPromote)
    VIGRA_EXEC_LOOP(fromRealPromote, = NumericTraits<T1>::fromRealPromote)
    VIGRA_EXEC_LOOP_SCALAR(addScalar, +)
    VIGRA_EXEC_LOOP_SCALAR(subScalar, -)
    VIGRA_EXEC_LOOP_SCALAR(mulScalar, *)
    VIGRA_EXEC_LOOP_SCALAR(divScalar, /)
    
    VIGRA_EXEC_LOOP_MINMAX(min, >)
    VIGRA_EXEC_LOOP_MINMAX(max, <)

    template <class T>
    static T const & minimum(T const * p)
    {
        return *std::min_element(p, p+LEVEL);
    }

    template <class T>
    static T const & maximum(T const * p)
    {
        return *std::max_element(p, p+LEVEL);
    }

    template <class T>
    static bool all(T const * p, T const & zero)
    {
        for(int i=0; i<LEVEL; ++i)
            if(p[i] == zero)
                return false;
        return true;
    }

    template <class T>
    static bool any(T const * p, T const & zero)
    {
        for(int i=0; i<LEVEL; ++i)
            if(p[i] != zero)
                return true;
        return false;
    }

    template <class T1, class T2>
    static bool notEqual(T1 const * left, T2 const * right)
    {
        for(int i=0; i<LEVEL; ++i)
            if(left[i] != right[i])
                return true;
        return false;
    }

    template <class T1, class T2>
    static bool lexicographicLessThan(T1 const * left, T2 const * right)
    {
        for(int i=0; i<LEVEL; ++i)
        {
            if(left[i] < right[i])
                return true;
            if(right[i] < left[i])
                return false;
        }
        return false;
    }
    
    template <class T>
    static bool closeAtTolerance(T const * left, T const * right, T epsilon)
    {
        bool res = true;
        for(int i=0; i<LEVEL; ++i)
        {
            res = res && vigra::closeAtTolerance(left[i], right[i], epsilon);
        }
        return res;
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
struct UnrollScalarResult
{
    template <class T>
    static typename NumericTraits<T>::Promote
    dot(T const * d)
    {
        return *d * *d + UnrollScalarResult<LEVEL-1>::dot(d+1);
    }

    template <class T1, class T2>
    static typename PromoteTraits<T1, T2>::Promote
    dot(T1 const * left, T2 const * right)
    {
        return *left * *right + UnrollScalarResult<LEVEL-1>::dot(left+1, right+1);
    }
    
    template <class T>
    static typename NormTraits<T>::SquaredNormType
    squaredNorm(T const * d)
    {
        return vigra::squaredNorm(*d) + UnrollScalarResult<LEVEL-1>::squaredNorm(d+1);
    }

    static std::ptrdiff_t
    squaredNorm(std::ptrdiff_t const * d)
    {
        return (*d)*(*d) + UnrollScalarResult<LEVEL-1>::squaredNorm(d+1);
    }
    
    template <class T>
    static T const & minimum(T const * p)
    {
        T const & m = UnrollScalarResult<LEVEL - 1>::minimum(p+1);
        return *p < m
                    ? *p
                    : m;
    }

    template <class T>
    static T const & maximum(T const * p)
    {
        T const & m = UnrollScalarResult<LEVEL - 1>::maximum(p+1);
        return *p > m
                    ? *p
                    : m;
    }

    template <class T>
    static bool all(T const * p, T const & zero)
    {
        return *p != zero && UnrollScalarResult<LEVEL - 1>::all(p+1, zero);
    }

    template <class T>
    static bool any(T const * p, T const & zero)
    {
        return *p != zero || UnrollScalarResult<LEVEL - 1>::any(p+1, zero);
    }
};

template <>
struct UnrollScalarResult<1>
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

    template <class T>
    static typename NormTraits<T>::SquaredNormType
    squaredNorm(T const * d)
    {
        return vigra::squaredNorm(*d);
    }

    static std::ptrdiff_t
    squaredNorm(std::ptrdiff_t const * d)
    {
        return (*d)*(*d);
    }

    template <class T>
    static T const & minimum(T const * p)
    {
        return *p;
    }

    template <class T>
    static T const & maximum(T const * p)
    {
        return *p;
    }

    template <class T>
    static bool all(T const * p, T const & zero)
    {
        return *p != zero;
    }

    template <class T>
    static bool any(T const * p, T const & zero)
    {
        return *p != zero;
    }
};

#undef VIGRA_EXEC_LOOP
#undef VIGRA_EXEC_LOOP_MINMAX
#undef VIGRA_EXEC_LOOP_SCALAR

#define VIGRA_UNROLL_LOOP(NAME, OPER) \
    template <class T1, class T2>  \
    static void NAME(T1 * left, T2 const * right)  \
    {  \
        (*left) OPER (*right);  \
        UnrollLoop<LEVEL-1>::NAME(left+1, right+1); \
    }

#define VIGRA_UNROLL_LOOP_MINMAX(NAME, OPER) \
    template <class T1, class T2>  \
    static void NAME(T1 * left, T2 const * right)  \
    {  \
        if(*left OPER *right) \
            *left = *right;  \
        UnrollLoop<LEVEL-1>::NAME(left+1, right+1); \
    }

#define VIGRA_UNROLL_LOOP_SCALAR(NAME, OPER) \
    template <class T1, class T2>  \
    static void NAME(T1 * left, T2 right)  \
    {  \
        (*left) = detail::RequiresExplicitCast<T1>::cast((*left) OPER (right));  \
        UnrollLoop<LEVEL-1>::NAME(left+1, right); \
    }


template <int LEVEL>
struct UnrollLoop
{
    template <class T1, class T2>
    static void reverseAssign(T1 * left, T2 const * right)
    {
        *left = *right;
        UnrollLoop<LEVEL-1>::reverseAssign(left+1, right-1);
    }

    template <class T1, class T2>
    static void assignCast(T1 * left, T2 const * right)
    {
        *left = detail::RequiresExplicitCast<T1>::cast(*right);
        UnrollLoop<LEVEL-1>::assignCast(left+1, right+1);
    }

    template <class T1, class T2>
    static void assignScalar(T1 * left, T2 right)
    {
        *left = detail::RequiresExplicitCast<T1>::cast(right);
        UnrollLoop<LEVEL-1>::assignScalar(left+1, right);
    }

    template <class T1, class T2>
    static void power(T1 * left, T2 right)
    {
        *left = detail::RequiresExplicitCast<T1>::cast(pow(*left, right));
        UnrollLoop<LEVEL-1>::power(left+1, right);
    }

    VIGRA_UNROLL_LOOP(assign, =)
    VIGRA_UNROLL_LOOP(add, +=)
    VIGRA_UNROLL_LOOP(sub, -=)
    VIGRA_UNROLL_LOOP(mul, *=)
    VIGRA_UNROLL_LOOP(div, /=)
    VIGRA_UNROLL_LOOP(mod, %=)
    VIGRA_UNROLL_LOOP(neg, = -)
    VIGRA_UNROLL_LOOP(abs, = vigra::abs)
    VIGRA_UNROLL_LOOP(floor, = vigra::floor)
    VIGRA_UNROLL_LOOP(ceil, = vigra::ceil)
    VIGRA_UNROLL_LOOP(round, = vigra::round)
    VIGRA_UNROLL_LOOP(sqrt, = vigra::sqrt)
    VIGRA_UNROLL_LOOP(fromPromote, = NumericTraits<T1>::fromPromote)
    VIGRA_UNROLL_LOOP(fromRealPromote, = NumericTraits<T1>::fromRealPromote)
    VIGRA_UNROLL_LOOP_SCALAR(addScalar, +)
    VIGRA_UNROLL_LOOP_SCALAR(subScalar, -)
    VIGRA_UNROLL_LOOP_SCALAR(mulScalar, *)
    VIGRA_UNROLL_LOOP_SCALAR(divScalar, /)
    
    VIGRA_UNROLL_LOOP_MINMAX(min, >)
    VIGRA_UNROLL_LOOP_MINMAX(max, <)

    template <class T>
    static T const & minimum(T const * p)
    {
        return UnrollScalarResult<LEVEL>::minimum(p);
    }

    template <class T>
    static T const & maximum(T const * p)
    {
        return UnrollScalarResult<LEVEL>::maximum(p);
    }

    template <class T>
    static bool all(T const * p, T const & zero)
    {
        return UnrollScalarResult<LEVEL>::all(p, zero);
    }

    template <class T>
    static bool any(T const * p, T const & zero)
    {
        return UnrollScalarResult<LEVEL>::any(p, zero);
    }

    template <class T1, class T2>
    static bool notEqual(T1 const * left, T2 const * right)
    {
        return (*left != *right) || UnrollLoop<LEVEL - 1>::notEqual(left+1, right+1);
    }

    template <class T1, class T2>
    static bool lexicographicLessThan(T1 const * left, T2 const * right)
    {
        if(*left < *right)
            return true;
        if(*right < *left)
            return false;
        return UnrollLoop<LEVEL - 1>::lexicographicLessThan(left+1, right+1);
    }

    template <class T>
    static bool closeAtTolerance(T const * left, T const * right, T epsilon)
    {
        return vigra::closeAtTolerance(*left, *right, epsilon) && 
                  UnrollLoop<LEVEL - 1>::closeAtTolerance(left+1, right+1, epsilon);
    }
    
    template <class T>
    static typename NumericTraits<T>::Promote
    dot(T const * d)
    {
        return UnrollScalarResult<LEVEL>::dot(d);
    }

    template <class T1, class T2>
    static typename PromoteTraits<T1, T2>::Promote
    dot(T1 const * left, T2 const * right)
    {
        return UnrollScalarResult<LEVEL>::dot(left, right);
    }

    template <class T>
    static typename NormTraits<T>::SquaredNormType
    squaredNorm(T const * d)
    {
        return UnrollScalarResult<LEVEL>::squaredNorm(d);
    }
};

#undef VIGRA_UNROLL_LOOP
#undef VIGRA_UNROLL_LOOP_MINMAX
#undef VIGRA_UNROLL_LOOP_SCALAR

template <>
struct UnrollLoop<0>
{
    template <class T1, class T2>
    static void reverseAssign(T1, T2) {}
    template <class T1, class T2>
    static void assignCast(T1, T2) {}
    template <class T1, class T2>
    static void assign(T1, T2) {}
    template <class T1, class T2>
    static void assignScalar(T1, T2) {}
    template <class T1, class T2>
    static void power(T1, T2) {}
    template <class T1, class T2>
    static void add(T1, T2) {}
    template <class T1, class T2>
    static void addScalar(T1, T2) {}
    template <class T1, class T2>
    static void sub(T1, T2) {}
    template <class T1, class T2>
    static void subScalar(T1, T2) {}
    template <class T1, class T2>
    static void mul(T1, T2) {}
    template <class T1, class T2>
    static void mulScalar(T1, T2) {}
    template <class T1, class T2>
    static void div(T1, T2) {}
    template <class T1, class T2>
    static void mod(T1, T2) {}
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
    static void round(T1, T2) {}
    template <class T1, class T2>
    static void sqrt(T1, T2) {}
    template <class T1, class T2>
    static bool notEqual(T1, T2) { return false; }
    template <class T1, class T2>
    static bool lexicographicLessThan(T1, T2) { return false; }
    template <class T1, class T2>
    static void min(T1, T2) {}
    template <class T1, class T2>
    static void max(T1, T2) {}
    template <class T>
    static T minimum(T const *) { return NumericTraits<T>::max(); }
    template <class T>
    static T maximum(T const *) { return NumericTraits<T>::min(); }
    template <class T>
    static bool all(T const *, T const &) { return true; }
    template <class T>
    static bool any(T const *, T const &) { return false; }
    template <class T>
    static bool closeAtTolerance(T const *, T const *, T) { return true; }
};

template <int SIZE>
struct LoopType
{
    static const int MaxUnrollSize = 5;
    typedef typename IfBool<(SIZE <= MaxUnrollSize), UnrollLoop<SIZE>, ExecLoop<SIZE> >::type type;

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

    <b>\#include</b> \<vigra/tinyvector.hxx\><br>
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
    typedef std::ptrdiff_t difference_type;

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

        /** Initialize with constant value
        */
    void init(value_type initial)
    {
        Loop::assignScalar(data_, initial);
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

        /** Component-wise divide-assignment
        */
    template <class T1, class D1, class D2>
    DERIVED & operator/=(TinyVectorBase<T1, SIZE, D1, D2> const & r)
    {
        Loop::div(data_, r.begin());
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise modulo-assignment
        */
    template <class T1, class D1, class D2>
    DERIVED & operator%=(TinyVectorBase<T1, SIZE, D1, D2> const & r)
    {
        Loop::mod(data_, r.begin());
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise scalar multiply-assignment
        */
    DERIVED & operator+=(double r)
    {
        Loop::addScalar(data_, r);
        return static_cast<DERIVED &>(*this);
    }

        /** Component-wise scalar divide-assignment
        */
    DERIVED & operator-=(double r)
    {
        Loop::subScalar(data_, r);
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

        /** Return the minimal element.
        */
    VALUETYPE const & minimum() const
    {
        return Loop::minimum(data_);
    }

        /** Return the maximal element.
        */
    VALUETYPE const & maximum() const
    {
        return Loop::maximum(data_);
    }
    
        /** Check that all elements of this vector are non-zero (or 'true' if T is bool).
        */
    bool all() const
    {
        return Loop::all(data_, VALUETYPE());
    }
    
        /** Check that at least one element of this vector is non-zero (or 'true' if T is bool).
        */
    bool any() const
    {
        return Loop::any(data_, VALUETYPE());
    }

        /** Access component by index.
        */
    reference operator[](difference_type i)
    {
        VIGRA_ASSERT_INSIDE(i);
        return data_[i];
    }

        /** Get component by index.
        */
    const_reference operator[](difference_type i) const
    {
        VIGRA_ASSERT_INSIDE(i);
        return data_[i];
    }

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
    
        /** Get a view to the subarray with length <tt>(TO-FROM)</tt> starting at <tt>FROM</tt>.
            The bounds must fullfill <tt>0 <= FROM < TO <= SIZE</tt>, but this is only
            checked when <tt>VIGRA_CHECK_BOUNDS</tt> is \#define'd.
        */
    template <int FROM, int TO>
    TinyVectorView<VALUETYPE, TO-FROM> subarray() const
    {
#ifdef VIGRA_CHECK_BOUNDS
        vigra_precondition(FROM >= 0, "Index out of bounds");
        vigra_precondition(FROM < TO, "Index out of bounds");
        vigra_precondition(TO <=SIZE, "Index out of bounds");
#endif
        return TinyVectorView<VALUETYPE, TO-FROM>(data_+FROM);
    }
    
    TinyVector<VALUETYPE, SIZE-1>
    dropIndex(int m) const
    {
#ifdef VIGRA_CHECK_BOUNDS
        vigra_precondition(0 <= m && m < SIZE, "Dimension out of bounds");
#endif
        TinyVector<VALUETYPE, SIZE-1> res(SkipInitialization);
        for(int k=0; k<m; ++k)
            res[k] = data_[k];
        for(int k=m; k<SIZE-1; ++k)
            res[k] = data_[k+1];
        return res;
    }

        /** Size of TinyVector vector always equals the template parameter SIZE.
        */
    size_type size() const { return SIZE; }

    pointer data() { return data_; }

    const_pointer data() const { return data_; }
    
        /** \brief Factory function for a unit vector for dimension \a k.
        */
    static TinyVector<VALUETYPE, SIZE> unitVector(int k)
    {
        VIGRA_ASSERT_INSIDE(k);
        TinyVector<VALUETYPE, SIZE> ret;
        ret[k] = 1;
        return ret;
    }
    
        /** \brief Factory function for a linear sequence.
        
            The result will be initialized as <tt>res[k] = start + k*step</tt>.
        */
    static TinyVector<VALUETYPE, SIZE> linearSequence(VALUETYPE start=VALUETYPE(), VALUETYPE step=VALUETYPE(1))
    {
        TinyVector<VALUETYPE, SIZE> ret(SkipInitialization);
        for(int k=0; k<SIZE; ++k, start+=step)
            ret[k] = start;
        return ret;
    }

  protected:
  
    DATA data_;
};

/** \brief Class for fixed size vectors.
    \ingroup RangesAndPoints

    This class contains an array of size SIZE of the specified VALUETYPE.
    The interface conforms to STL vector, except that there are no functions
    that change the size of a TinyVector.

    \ref TinyVectorOperators "Arithmetic operations"
    on TinyVectors are defined as component-wise applications of these
    operations. Addition and subtraction of two TinyVectors
    (+=, -=, +, -, unary -), multiplication and division of an
    TinyVector with a double, and NumericTraits/PromoteTraits are defined,
    so that TinyVector fulfills the requirements of \ref LinearAlgebraConcept "Linear Algebra".

    VIGRA algorithms typically use \ref vigra::VectorAccessor to access
    TinyVectors as a whole, or specific components of them.

    See also:<br>
    <UL style="list-style-image:url(documents/bullet.gif)">
        <LI> \ref vigra::TinyVectorBase
        <LI> \ref vigra::TinyVectorView
        <LI> \ref TinyVectorTraits
        <LI> \ref TinyVectorOperators
    </UL>

    <b>\#include</b> \<vigra/tinyvector.hxx\><br>
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
    
    enum ReverseCopyTag { ReverseCopy };

        /** Construction with constant value.
        
            Initializes all vector elements with the given value.
        */
    explicit TinyVector(value_type const & initial)
    : BaseType()
    {
        Loop::assignScalar(BaseType::begin(), initial);
    }

        /** Construction from lemon::Invalid.
        
            Initializes all vector elements with -1.
        */
    TinyVector(lemon::Invalid const &)
    : BaseType()
    {
        Loop::assignScalar(BaseType::begin(), -1);
    }

        /** Construction with Diff2D.
        
            Use only when <tt>SIZE == 2</tt>.
        */
    explicit TinyVector(Diff2D const & initial)
    : BaseType()
    {
        BaseType::data_[0] = detail::RequiresExplicitCast<T>::cast(initial.x);
        BaseType::data_[1] = detail::RequiresExplicitCast<T>::cast(initial.y);
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

        /** Construction with explicit values.
            Call only if SIZE == 5
        */
    TinyVector(value_type const & i1, value_type const & i2,
               value_type const & i3, value_type const & i4,
               value_type const & i5)
    : BaseType()
    {
        BaseType::data_[0] = i1;
        BaseType::data_[1] = i2;
        BaseType::data_[2] = i3;
        BaseType::data_[3] = i4;
        BaseType::data_[4] = i5;
    }
    
       /** Default constructor (initializes all elements with zero).
        */
    TinyVector()
    : BaseType()
    {
        Loop::assignScalar(BaseType::data_, value_type());
    }

        /** Construct without initializing the vector elements.
        */
    explicit TinyVector(SkipInitializationTag)
    : BaseType()
    {}

    explicit TinyVector(detail::DontInit)
    : BaseType()
    {}

        /** Copy constructor.
        */
    TinyVector(TinyVector const & r)
    : BaseType()
    {
        Loop::assign(BaseType::data_, r.data_);
    }

        /** Constructor from C array.
        */
    template <class U>
    explicit TinyVector(U const * data)
    : BaseType()
    {
        Loop::assign(BaseType::data_, data);
    }

        /** Constructor by reverse copy from C array.
            
            Usage:
            \code
            TinyVector<int, 3> v(1,2,3);
            TinyVector<int, 3> reversed(v.begin(), TinyVector<int, 3>::ReverseCopy);
            \endcode
        */
    explicit TinyVector(const_pointer data, ReverseCopyTag)
    : BaseType()
    {
        Loop::reverseAssign(BaseType::data_, data+SIZE-1);
    }

        /** Copy with type conversion.
        */
    template <class U, class DATA, class DERIVED>
    TinyVector(TinyVectorBase<U, SIZE, DATA, DERIVED> const & r)
    : BaseType()
    {
        Loop::assignCast(BaseType::data_, r.begin());
    }

        /** Copy assignment.
        */
    TinyVector & operator=(TinyVector const & r)
    {
        Loop::assign(BaseType::data_, r.data_);
        return *this;
    }

        /** Copy assignment with type conversion.
        */
    template <class U, class DATA, class DERIVED>
    TinyVector & operator=(TinyVectorBase<U, SIZE, DATA, DERIVED> const & r)
    {
        Loop::assignCast(BaseType::data_, r.begin());
        return *this;
    }

        /** Assignment from Diff2D.
        
            Use only when <tt>SIZE == 2</tt>.
        */
    TinyVector & operator=(Diff2D const & r)
    {
        BaseType::data_[0] = detail::RequiresExplicitCast<T>::cast(r.x);
        BaseType::data_[1] = detail::RequiresExplicitCast<T>::cast(r.y);
        return *this;
    }

        /** Assignment from scalar. Will set all entries to the given value.
        */
    TinyVector & operator=(value_type const & v)
    {
        Loop::assignScalar(BaseType::begin(), v);
        return *this;
    }

        /** Copy from a TinyVector with a different number of elements.
        
            Only the first <tt>min(SIZE, USIZE)</tt> elements are copied.
        */
    template <class U, int USIZE, class DATA, class DERIVED>
    TinyVector & copy(TinyVectorBase<U, USIZE, DATA, DERIVED> const & r)
    {
        static const int minSize = USIZE < SIZE
                                        ? USIZE
                                        : SIZE;
        
        typedef typename detail::LoopType<minSize>::type MinLoop;
        MinLoop::assignCast(BaseType::data_, r.begin());
        return *this;
    }
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
    so that TinyVectorView fulfills the requirements of \ref LinearAlgebraConcept "Linear Algebra".

    VIGRA algorithms typically use \ref vigra::VectorAccessor to access
    TinyVectorViews as a whole, or specific components of them.

    <b>See also:</b>
    <ul>
        <li> \ref vigra::TinyVectorBase
        <li> \ref vigra::TinyVector
        <li> \ref TinyVectorTraits
        <li> \ref TinyVectorOperators
    </ul>

    <b>\#include</b> \<vigra/tinyvector.hxx\><br>
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

    \brief Implement basic arithmetic and equality for TinyVector.

    These functions fulfill the requirements of a Linear Space (vector space).
    Return types are determined according to \ref TinyVectorTraits.

    <b>\#include</b> \<vigra/tinyvector.hxx\><br>
    Namespace: vigra
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

    /// lexicographical comparison
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline bool
operator<(TinyVectorBase<V1, SIZE, D1, D2> const & l,
                      TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    return ltype::lexicographicLessThan(l.begin(), r.begin());
}


    /// pointwise less-than
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline bool
allLess(TinyVectorBase<V1, SIZE, D1, D2> const & l,
        TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    for(int k=0; k < SIZE; ++k)
        if (l[k] >= r[k])
            return false;
    return true;
}

    /// pointwise greater-than
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline bool
allGreater(TinyVectorBase<V1, SIZE, D1, D2> const & l,
           TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    for(int k=0; k < SIZE; ++k)
        if(l[k] <= r[k])
            return false;
    return true;
}

    /// pointwise less-equal
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline bool
allLessEqual(TinyVectorBase<V1, SIZE, D1, D2> const & l,
             TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    for(int k=0; k < SIZE; ++k)
        if (l[k] > r[k])
            return false;
    return true;
}

    /// pointwise greater-equal
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline bool
allGreaterEqual(TinyVectorBase<V1, SIZE, D1, D2> const & l,
                TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    for(int k=0; k < SIZE; ++k)
        if (l[k] < r[k])
            return false;
    return true;
}

template <class V, int SIZE, class D1, class D2, class D3, class D4>
bool 
closeAtTolerance(TinyVectorBase<V, SIZE, D1, D2> const & l,
                 TinyVectorBase<V, SIZE, D3, D4> const & r, 
                 V epsilon = NumericTraits<V>::epsilon())
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    return ltype::closeAtTolerance(l.begin(), r.begin(), epsilon);
}

template <class V, int SIZE>
bool 
closeAtTolerance(TinyVector<V, SIZE> const & l,
                 TinyVector<V, SIZE> const & r, 
                 V epsilon = NumericTraits<V>::epsilon())
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    return ltype::closeAtTolerance(l.begin(), r.begin(), epsilon);
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
        typedef typename NumericTraits<T>::isSigned isSigned;

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

    <b>\#include</b> \<vigra/tinyvector.hxx\><br>
    Namespace: vigra

    On compilers that don't support partial template specialization (e.g.
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
    typedef typename NumericTraits<T>::isSigned isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;

    static TinyVector<T, SIZE> zero()
    {
        return TinyVector<T, SIZE>(NumericTraits<T>::zero());
    }
    static TinyVector<T, SIZE> one()
    {
        return TinyVector<T, SIZE>(NumericTraits<T>::one());
    }
    static TinyVector<T, SIZE> nonZero()
    {
        return TinyVector<T, SIZE>(NumericTraits<T>::nonZero());
    }

    static TinyVector<T, SIZE> min()
    {
        return TinyVector<T, SIZE>(NumericTraits<T>::min());
    }
    static TinyVector<T, SIZE> max()
    {
        return TinyVector<T, SIZE>(NumericTraits<T>::max());
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
    typedef typename NumericTraits<T>::isSigned isSigned;
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

template<class T, int SIZE>
struct CanSkipInitialization<TinyVectorView<T, SIZE> >
{
    typedef typename CanSkipInitialization<T>::type type;
    static const bool value = type::asBool;
};

template<class T, int SIZE>
struct CanSkipInitialization<TinyVector<T, SIZE> >
{
    typedef typename CanSkipInitialization<T>::type type;
    static const bool value = type::asBool;
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
    typedef NumericTraits<T>::isSigned isSigned; \
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

    /// component-wise division
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2, SIZE> >::Promote
operator/(TinyVectorBase<V1, SIZE, D1, D2> const & l,
          TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    return typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2 , SIZE> >::Promote(l) /= r;
}

    /// component-wise modulo
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2, SIZE> >::Promote
operator%(TinyVectorBase<V1, SIZE, D1, D2> const & l,
          TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    return typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2 , SIZE> >::Promote(l) %= r;
}

    /// component-wise left scalar addition
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator+(double v, TinyVectorBase<V, SIZE, D1, D2> const & r)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(r) += v;
}

    /// component-wise right scalar addition
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator+(TinyVectorBase<V, SIZE, D1, D2> const & l, double v)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(l) += v;
}

    /// component-wise left scalar subtraction
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator-(double v, TinyVectorBase<V, SIZE, D1, D2> const & r)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(v) -= r;
}

    /// component-wise right scalar subtraction
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator-(TinyVectorBase<V, SIZE, D1, D2> const & l, double v)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(l) -= v;
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

    /// component-wise left scalar division
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator/(double v, TinyVectorBase<V, SIZE, D1, D2> const & r)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(v) /= r;
}

    /// component-wise right scalar division
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator/(TinyVectorBase<V, SIZE, D1, D2> const & l, double v)
{
    return typename NumericTraits<TinyVector<V, SIZE> >::RealPromote(l) /= v;
}

    /// component-wise scalar division without type promotion
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<V, SIZE>
div(TinyVectorBase<V, SIZE, D1, D2> const & l, V v)
{
    TinyVector<V, SIZE> result(l);
    typedef typename detail::LoopType<SIZE>::type Loop;
    Loop::divScalar(result.data(), v);
    return result;
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

    /** Apply round() function to each vector component.
    */
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<V, SIZE>
round(TinyVectorBase<V, SIZE, D1, D2> const & v)
{
    TinyVector<V, SIZE> res(detail::dontInit());
    typedef typename detail::LoopType<SIZE>::type ltype;
    ltype::round(res.begin(), v.begin());
    return res;
}

    /** Apply sqrt() function to each vector component.
    */
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<V, SIZE>
sqrt(TinyVectorBase<V, SIZE, D1, D2> const & v)
{
    TinyVector<V, SIZE> res(detail::dontInit());
    typedef typename detail::LoopType<SIZE>::type ltype;
    ltype::sqrt(res.begin(), v.begin());
    return res;
}

using std::pow;

    /** Apply pow() function to each vector component.
    */
template <class V, int SIZE, class D1, class D2, class E>
inline
TinyVector<V, SIZE>
pow(TinyVectorBase<V, SIZE, D1, D2> const & v, E exponent)
{
    TinyVector<V, SIZE> res(v);
    typedef typename detail::LoopType<SIZE>::type ltype;
    ltype::power(res.begin(), exponent);
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

    /// sum of the vector's elements
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<V>::Promote
sum(TinyVectorBase<V, SIZE, D1, D2> const & l)
{
    typename NumericTraits<V>::Promote res = l[0];
    for(int k=1; k<SIZE; ++k)
        res += l[k];
    return res;
}

    /// cumulative sum of the vector's elements
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<typename NumericTraits<V>::Promote, SIZE>
cumsum(TinyVectorBase<V, SIZE, D1, D2> const & l)
{
    TinyVector<typename NumericTraits<V>::Promote, SIZE> res(l);
    for(int k=1; k<SIZE; ++k)
        res[k] += res[k-1];
    return res;
}

    /// product of the vector's elements
template <class V, int SIZE, class D1, class D2>
inline
typename NumericTraits<V>::Promote
prod(TinyVectorBase<V, SIZE, D1, D2> const & l)
{
    typename NumericTraits<V>::Promote res = l[0];
    for(int k=1; k<SIZE; ++k)
        res *= l[k];
    return res;
}

    /// cumulative product of the vector's elements
template <class V, int SIZE, class D1, class D2>
inline
TinyVector<typename NumericTraits<V>::Promote, SIZE>
cumprod(TinyVectorBase<V, SIZE, D1, D2> const & l)
{
    TinyVector<typename NumericTraits<V>::Promote, SIZE> res(l);
    for(int k=1; k<SIZE; ++k)
        res[k] *= res[k-1];
    return res;
}

using std::min;

    /// element-wise minimum
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
TinyVector<typename PromoteTraits<V1, V2>::Promote, SIZE>
min(TinyVectorBase<V1, SIZE, D1, D2> const & l,
    TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    TinyVector<typename PromoteTraits<V1, V2>::Promote, SIZE> res(l);
    ltype::min(res.begin(), r.begin());
    return res;
}

// we also have to overload min for like-typed argument to prevent match of std::min()
template <class V1, int SIZE, class D1, class D2>
inline
TinyVector<V1, SIZE>
min(TinyVectorBase<V1, SIZE, D1, D2> const & l,
    TinyVectorBase<V1, SIZE, D1, D2> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    TinyVector<V1, SIZE> res(l);
    ltype::min(res.begin(), r.begin());
    return res;
}

template <class V1, int SIZE>
inline
TinyVector<V1, SIZE>
min(TinyVector<V1, SIZE> const & l,
    TinyVector<V1, SIZE> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    TinyVector<V1, SIZE> res(l);
    ltype::min(res.begin(), r.begin());
    return res;
}

    /// minimum element
template <class V, int SIZE, class D1, class D2>
inline
V const &
min(TinyVectorBase<V, SIZE, D1, D2> const & l)
{
    return l.minimum();
}

using std::max;

    /// element-wise maximum
template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
inline
TinyVector<typename PromoteTraits<V1, V2>::Promote, SIZE>
max(TinyVectorBase<V1, SIZE, D1, D2> const & l,
    TinyVectorBase<V2, SIZE, D3, D4> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    TinyVector<typename PromoteTraits<V1, V2>::Promote, SIZE> res(l);
    ltype::max(res.begin(), r.begin());
    return res;
}

// we also have to overload max for like-typed argument to prevent match of std::max()
template <class V1, int SIZE, class D1, class D2>
inline
TinyVector<V1, SIZE>
max(TinyVectorBase<V1, SIZE, D1, D2> const & l,
    TinyVectorBase<V1, SIZE, D1, D2> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    TinyVector<V1, SIZE> res(l);
    ltype::max(res.begin(), r.begin());
    return res;
}

template <class V1, int SIZE>
inline
TinyVector<V1, SIZE>
max(TinyVector<V1, SIZE> const & l,
    TinyVector<V1, SIZE> const & r)
{
    typedef typename detail::LoopType<SIZE>::type ltype;
    TinyVector<V1, SIZE> res(l);
    ltype::max(res.begin(), r.begin());
    return res;
}

    /// maximum element
template <class V, int SIZE, class D1, class D2>
inline
V const &
max(TinyVectorBase<V, SIZE, D1, D2> const & l)
{
    return l.maximum();
}

    /// squared norm
template <class V1, int SIZE, class D1, class D2>
inline
typename TinyVectorBase<V1, SIZE, D1, D2>::SquaredNormType
squaredNorm(TinyVectorBase<V1, SIZE, D1, D2> const & t)
{
    return t.squaredMagnitude();
}

    /// squared norm
template <class V, int SIZE>
inline
typename TinyVector<V, SIZE>::SquaredNormType
squaredNorm(TinyVector<V, SIZE> const & t)
{
    return t.squaredMagnitude();
}

using std::reverse;

    /// reversed copy
template <class V, int SIZE>
inline
TinyVector<V, SIZE>
reverse(TinyVector<V, SIZE> const & t)
{
    return TinyVector<V, SIZE>(t.begin(), TinyVector<V, SIZE>::ReverseCopy);
}

    /** \brief transposed copy
    
        Elements are arranged such that <tt>res[k] = t[permutation[k]]</tt>.
    */
template <class V, int SIZE, class T>
inline
TinyVector<V, SIZE>
transpose(TinyVector<V, SIZE> const & t, TinyVector<T, SIZE> const & permutation)
{
    TinyVector<V, SIZE> res(SkipInitialization);
    for(int k=0; k<SIZE; ++k)
    {
        VIGRA_ASSERT_INSIDE(permutation[k]);
        res[k] = t[permutation[k]];
    }
    return res;
}

    /** \brief transposed copy
    
        All elements smaller 0 are clipped to zero.
    */
template<class V,int SIZE>
inline
TinyVector<V, SIZE> clipLower(TinyVector<V, SIZE> const & t){
    TinyVector<V, SIZE> res(SkipInitialization);
    for(int k=0; k<SIZE; ++k){
        res[k]=t[k]< static_cast<V>(0) ?  static_cast<V>(0) :  t[k];
    }
    return res;
}

    /** \brief transposed copy
    
        All elements smaller val are clipped to val.
    */
template<class V,int SIZE>
inline
TinyVector<V, SIZE> clipLower(TinyVector<V, SIZE> const & t,const V val){
    TinyVector<V, SIZE> res(SkipInitialization);
    for(int k=0; k<SIZE; ++k){
        res[k]=t[k]< val ? val :  t[k];
    }
    return res;
}
    /** \brief transposed copy
    
        All elements bigger val are clipped to val.
    */
template<class V,int SIZE>
inline
TinyVector<V, SIZE> clipUpper(TinyVector<V, SIZE> const & t,const V val){
    TinyVector<V, SIZE> res(SkipInitialization);
    for(int k=0; k<SIZE; ++k){
        res[k]=t[k]> val ? val :  t[k];
    }
    return res;
}


template<class V,int SIZE>
inline
TinyVector<V, SIZE> clip(TinyVector<V, SIZE> const & t,const V valLow,const V valUpper){
    TinyVector<V, SIZE> res(SkipInitialization);
    for(int k=0; k<SIZE; ++k){
        res[k]=t[k]< valLow   ? valLow :  t[k];
        res[k]=t[k]> valUpper ? valUpper :  t[k];
    }
    return res;
}

template<class V,int SIZE>
inline
bool isZero(TinyVector<V, SIZE> const & t){
    for(int k=0; k<SIZE; ++k){
        if(t[k]!=static_cast<V>(0))
            return false;
    }
    return true;
}

template<class V,int SIZE>
inline typename NumericTraits<V>::RealPromote
mean(TinyVector<V, SIZE> const & t){
    const V sumVal = sum(t);
    return static_cast< typename NumericTraits<V>::RealPromote>(sumVal)/SIZE;
}


template<class V,int SIZE>
inline typename NumericTraits<V>::RealPromote
sizeDividedSquaredNorm(TinyVector<V, SIZE> const & t){
    return squaredNorm(t)/SIZE;
}

template<class V,int SIZE>
inline typename NumericTraits<V>::RealPromote
sizeDividedNorm(TinyVector<V, SIZE> const & t){
    return norm(t)/SIZE;
}



//@}

// mask cl.exe shortcomings [end]
#if defined(_MSC_VER)
#pragma warning( pop )
#endif

} // namespace vigra
#undef VIGRA_ASSERT_INSIDE
#endif // VIGRA_TINYVECTOR_HXX
