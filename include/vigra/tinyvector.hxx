/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
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
 
 
#ifndef VIGRA_TINYVECTOR_HXX
#define VIGRA_TINYVECTOR_HXX

#include <cmath>    // abs(double)
#include <cstdlib>  // abs(int)
#include "vigra/config.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/accessor.hxx"

namespace vigra {

/********************************************************/
/*                                                      */
/*                      TinyVector                      */
/*                                                      */
/********************************************************/

/** \brief Class for fixed size vectors.

    This class contains an array of size SIZE of the specified VALUETYPE. 
    The interface conforms with STL vector, except that there are no functions
    that change the size of a TinyVector. 
    
    \ref TinyVectorOperators "Arithmetic operations"
    on TinyVectors are defined as component-wise applications of these 
    operations. Addition and subtraction of two TinyVectors 
    (+=, -=, +, -, unary -), multiplication and division of an
    TinyVector with a double, and NumericTraits/PromoteTraits are defined, 
    so that TinyVector fulfills the requirements of a Linear Space. 
    
    VIGRA algorithms typically use \ref vigra::VectorAccessor to access 
    TinyVectors as a whole, or specific components of them.
    
    <b>\#include</b> "<a href="tinyvector_8hxx-source.html">vigra/tinyvector.hxx</a>"<br>
    Namespace: vigra
**/
template <class VALUETYPE, int SIZE>
class TinyVector
{
  public:
        /** STL-compatible definition of valuetype
        */
    typedef VALUETYPE value_type;
        /** STL-compatible definition of iterator
        */
    typedef value_type * iterator;
        /** STL-compatible definition of const iterator
        */
    typedef value_type const * const_iterator;

        /** Construct from another sequence 
        */
    template <class Iterator>   
    TinyVector(Iterator i, Iterator end)
    {
        for(iterator p = begin(); i != end; ++i, ++p)
            *p = static_cast<value_type>(*i);
    }
    
        /** Construction with constant value 
        */    
    explicit TinyVector(value_type initial)
    {
        for(iterator p = begin(); p != end(); ++p)
            *p = initial;
    }
    
        /** Default constructor (sets all components to 0)  
        */    
    TinyVector()
    {
        for(iterator p = begin(); p != end(); ++p)
            *p = NumericTraits<value_type>::zero();
    }

#if !defined(TEMPLATE_COPY_CONSTRUCTOR_BUG)
        
    TinyVector(TinyVector const & r)
    {
        iterator p = begin();
        const_iterator q = r.begin();
        for(; p != end(); ++p, ++q)
            *p = *q;
    }

    TinyVector & operator=(TinyVector const & r)
    {
        iterator p = begin();
        const_iterator q = r.begin();
        for(; p != end(); ++p, ++q)
            *p = *q;
        return *this;
    }

#endif // TEMPLATE_COPY_CONSTRUCTOR_BUG

    
        /** Copy constructor.
            */
    template <class U>   
    TinyVector(TinyVector<U, SIZE> const & r)
    {
        iterator p = begin();
        typename TinyVector<U, SIZE>::const_iterator q = r.begin();
        for(; p != end(); ++p, ++q)
            *p = static_cast<value_type>(*q);
    }

        /** Copy assignment.
        */    
    template <class U>   
    TinyVector & operator=(TinyVector<U, SIZE> const & r)
    {
        iterator p = begin();
        typename TinyVector<U, SIZE>::const_iterator q = r.begin();
        for(; p != end(); ++p, ++q)
            *p = static_cast<value_type>(*q);
        return *this;
    }

        /** Unary negation (construct TinyVector with negative values)
        */
    TinyVector operator-() const
    {
        TinyVector r;
        const_iterator s = begin();
        iterator d = r.begin();
        for(; s != end(); ++s, ++d)
            *d = -(*s);
        return r;
    }
    
        /** Calculate magnitude.
        */
    typename NumericTraits<VALUETYPE>::RealPromote
    magnitude() const 
    {
#ifndef CMATH_NOT_IN_STD
         return std::sqrt(squaredMagnitude());
#else
         return sqrt(squaredMagnitude());
#endif
    }
    
        /** Calculate squared magnitude.
        */
    typename NumericTraits<VALUETYPE>::Promote
    squaredMagnitude() const 
    { 
         const_iterator i = begin();
         typename NumericTraits<VALUETYPE>::Promote sum = *i * *i;
         for(++i; i < end(); ++i)
            sum += *i * *i;
         return sum;
   }
    
        /** Access component by index.
        */
    value_type & operator[](int const & i) { return data_[i]; }
    
        /** Get component by index.
        */
    value_type operator[](int const & i) const { return data_[i]; }
    
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
    int size() const { return SIZE; }
    
  private:
    value_type data_[SIZE];
};

/********************************************************/
/*                                                      */
/*                     TinyVector Comparison            */
/*                                                      */
/********************************************************/

/** \addtogroup TinyVectorOperators Functions for TinyVector

    \brief <b>\#include</b> "<a href="tinyvector_8hxx-source.html">vigra/tinyvector.hxx</a>
    
    These functions fulfill the requirements of a Linear Space (vector space).

    Namespace: vigra
    <p>
    
 */
//@{
    /// component-wise equal
template <class V1, class V2, int SIZE>
inline
bool 
operator==(TinyVector<V1, SIZE> const & l, TinyVector<V2, SIZE> const & r)
{
    return !(l != r);
}

    /// component-wise not equal
template <class V1, class V2, int SIZE>
inline
bool 
operator!=(TinyVector<V1, SIZE> const & l, TinyVector<V2, SIZE> const & r)
{
    typename TinyVector<V1, SIZE>::const_iterator i1 = l.begin();
    typename TinyVector<V2, SIZE>::const_iterator i2 = r.begin();
    for(; i1 != l.end(); ++i1, ++i2)
        if(*i1 != *i2)
            return true;
    return false;
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

    typedef typename NumericTraits<T>::isIntegral isIntegral;
    typedef VigraFalseType isScalar;
    typedef VigraFalseType isOrdered;
    
    static TinyVector<T, SIZE> zero() { 
        return TinyVector<T, SIZE>(NumericTraits<T>::zero()); 
    }
    static TinyVector<T, SIZE> one() { 
        return TinyVector<T, SIZE>(NumericTraits<T>::one()); 
    }
    static TinyVector<T, SIZE> nonZero() { 
        return TinyVector<T, SIZE>(NumericTraits<T>::nonZero()); 
    }
    
    static Promote toPromote(TinyVector<T, SIZE> const & v) { 
        return Promote(v); 
    }
    static RealPromote toRealPromote(TinyVector<T, SIZE> const & v) { 
        return RealPromote(v); 
    }
    static TinyVector<T, SIZE> fromPromote(Promote const & v) {
        TinyVector<T, SIZE> res;
        TinyVector<T, SIZE>::iterator d = res.begin();
        typename Promote::const_iterator s = v.begin();
        for(; d != res.end(); ++d, ++s)
            *d = NumericTraits<T>::fromPromote(*s);
        return res;
    }
    static TinyVector<T, SIZE> fromRealPromote(RealPromote const & v) {
        TinyVector<T, SIZE> res;
        TinyVector<T, SIZE>::iterator d = res.begin();
        typename RealPromote::const_iterator s = v.begin();
        for(; d != res.end(); ++d, ++s)
            *d = NumericTraits<T>::fromRealPromote(*s);
        return res;
    }
};

template <class T1, class T2, int SIZE>
struct PromoteTraits<TinyVector<T1, SIZE>, TinyVector<T2, SIZE> >
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

#else // NO_PARTIAL_TEMPLATE_SPECIALIZATION


#define TINYVECTOR_NUMTRAITS(T, SIZE) \
template<>\
struct NumericTraits<TinyVector<T, SIZE> >\
{\
    typedef TinyVector<T, SIZE> Type;\
    typedef TinyVector<NumericTraits<T>::Promote, SIZE> Promote;\
    typedef TinyVector<NumericTraits<T>::RealPromote, SIZE> RealPromote;\
    typedef NumericTraits<T>::isIntegral isIntegral;\
    typedef VigraFalseType isScalar;\
    typedef VigraFalseType isOrdered;\
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
        TinyVector<T, SIZE>::iterator d = res.begin();\
        Promote::const_iterator s = v.begin();\
        for(; d != res.end(); ++d, ++s)\
            *d = NumericTraits<T>::fromPromote(*s);\
        return res;\
    }\
    static TinyVector<T, SIZE> fromRealPromote(RealPromote const & v) {\
        TinyVector<T, SIZE> res;\
        TinyVector<T, SIZE>::iterator d = res.begin();\
        RealPromote::const_iterator s = v.begin();\
        for(; d != res.end(); ++d, ++s)\
            *d = NumericTraits<T>::fromRealPromote(*s);\
        return res;\
    }\
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
    /// componentwise add-assignment
template <class V1, class V2, int SIZE>
inline
TinyVector<V1, SIZE> & 
operator+=(TinyVector<V1, SIZE> & l, TinyVector<V2, SIZE> const & r)
{
    typename TinyVector<V1, SIZE>::iterator i1 = l.begin();
    typename TinyVector<V2, SIZE>::const_iterator i2 = r.begin();
    for(; i1 != l.end(); ++i1, ++i2)
        *i1 += *i2;
    return l;
}

    /// componentwise subtract-assignment
template <class V1, class V2, int SIZE>
inline
TinyVector<V1, SIZE> & 
operator-=(TinyVector<V1, SIZE> & l, TinyVector<V2, SIZE> const & r)
{
    typename TinyVector<V1, SIZE>::iterator i1 = l.begin();
    typename TinyVector<V2, SIZE>::const_iterator i2 = r.begin();
    for(; i1 != l.end(); ++i1, ++i2)
        *i1 -= *i2;
    return l;
}


    /// componentwise scalar multiply-assignment
template <class V, int SIZE>
inline
TinyVector<V, SIZE> & 
operator*=(TinyVector<V, SIZE> & l, double r)
{
    typename TinyVector<V, SIZE>::iterator i = l.begin();
    for(; i != l.end(); ++i)
        *i *= r;
    return l;
}

    /// componentwise scalar divide-assignment 
template <class V, int SIZE>
inline
TinyVector<V, SIZE> & 
operator/=(TinyVector<V, SIZE> & l, double r)
{
    typename TinyVector<V, SIZE>::iterator i = l.begin();
    for(; i != l.end(); ++i)
        *i /= r;
    return l;
}

#ifndef CMATH_NOT_IN_STD
using std::abs;
#else
using ::abs;
#endif

    /// component-wise absolute value
template <class T, int SIZE>
inline
TinyVector<T, SIZE> abs(TinyVector<T, SIZE> const & v) { 
    TinyVector<T, SIZE> res;
    typename TinyVector<T, SIZE>::iterator d = res.begin();
    typename TinyVector<T, SIZE>::const_iterator s = v.begin();
    for(; d != res.end(); ++d, ++s)
        *d = abs(*s);
    return res;
}



    /// component-wise addition
template <class V1, class V2, int SIZE>
inline 
typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2, SIZE> >::Promote
operator+(TinyVector<V1, SIZE> const & r1, TinyVector<V2, SIZE> const & r2)
{
    typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2 , SIZE> >::Promote res(r1);
    
    res += r2;
    
    return res;
}

    /// component-wise subtraction
template <class V1, class V2, int SIZE>
inline 
typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2, SIZE> >::Promote
operator-(TinyVector<V1, SIZE> const & r1, TinyVector<V2, SIZE> const & r2)
{
    typename PromoteTraits<TinyVector<V1, SIZE>, TinyVector<V2 , SIZE> >::Promote res(r1);
    
    res -= r2;
    
    return res;
}


    /// component-wise left scalar multiplication
template <class V, int SIZE>
inline 
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator*(double v, TinyVector<V, SIZE> const & r)
{
    typename NumericTraits<TinyVector<V, SIZE> >::RealPromote res(r);
    
    res *= v;
    
    return res;
}

    /// component-wise right scalar multiplication
template <class V, int SIZE>
inline 
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator*(TinyVector<V, SIZE> const & r, double v)
{
    typename NumericTraits<TinyVector<V, SIZE> >::RealPromote res(r);
    
    res *= v;
    
    return res;
}

    /// component-wise scalar division
template <class V, int SIZE>
inline 
typename NumericTraits<TinyVector<V, SIZE> >::RealPromote
operator/(TinyVector<V, SIZE> const & r, double v)
{
    typename NumericTraits<TinyVector<V, SIZE> >::RealPromote res(r);
    
    res /= v;
    
    return res;
}

    /// dot product
template <class V1, class V2, int SIZE>
inline 
typename PromoteTraits<V1, V2>::Promote
dot(TinyVector<V1, SIZE> const & r1, TinyVector<V2, SIZE> const & r2)
{
    typename TinyVector<V1, SIZE>::const_iterator i1 = r1.begin();
    typename TinyVector<V2, SIZE>::const_iterator i2 = r2.begin();
    typename PromoteTraits<V1, V2>::Promote sum = *i1 * *i2;
    for(++i1, ++i2; i1 < r1.end(); ++i1, ++i2)
        sum += *i1 * *i2;
    return sum;
}

using std::ceil;

    /** Apply ceil() function to each vector component.
    */
template <class T, int SIZE>
inline
TinyVector<T, SIZE>
ceil(TinyVector<T, SIZE> const & v)
{
    TinyVector<T, SIZE> res;
    typename TinyVector<T, SIZE>::iterator d = res.begin();
    typename TinyVector<T, SIZE>::const_iterator s = v.begin();
    for(; d != res.end(); ++d, ++s)
        *d = ceil(*s);
    return res;
};

using std::floor;

    /** Apply floor() function to each vector component.
    */
template <class T, int SIZE>
inline
TinyVector<T, SIZE>
floor(TinyVector<T, SIZE> const & v)
{
    TinyVector<T, SIZE> res;
    typename TinyVector<T, SIZE>::iterator d = res.begin();
    typename TinyVector<T, SIZE>::const_iterator s = v.begin();
    for(; d != res.end(); ++d, ++s)
        *d = floor(*s);
    return res;
};

//@}


} // namespace vigra

#endif // VIGRA_TINYVECTOR_HXX
