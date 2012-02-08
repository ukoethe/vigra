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


#ifndef VIGRA_RGBVALUE_HXX
#define VIGRA_RGBVALUE_HXX

#include <cmath>    // abs(double)
#include <cstdlib>  // abs(int)
#include "config.hxx"
#include "numerictraits.hxx"
#include "accessor.hxx"
#include "tinyvector.hxx"
#include "static_assert.hxx"

namespace vigra {

namespace detail {

template <unsigned int I, unsigned int R, unsigned int G, unsigned int B>
struct SelectColorIndexRHS;

template <unsigned int R, unsigned int G, unsigned int B>
struct SelectColorIndexRHS<0, R, G, B>
{
    enum { res = R };
};

template <unsigned int R, unsigned int G, unsigned int B>
struct SelectColorIndexRHS<1, R, G, B>
{
    enum { res = G };
};

template <unsigned int R, unsigned int G, unsigned int B>
struct SelectColorIndexRHS<2, R, G, B>
{
    enum { res = B };
};

} // namespace detail

#ifndef DOXYGEN

template <unsigned int R, unsigned int G, unsigned int B>
struct RGBValue_bad_color_indices
: staticAssert::AssertBool<(R < 3 && G < 3 && B < 3 &&
                           ((1 << R) + (1 << G) + (1 << B) == 7))>
{};

#endif /* DOXYGEN */


/********************************************************/
/*                                                      */
/*                      RGBValue                        */
/*                                                      */
/********************************************************/

/** \brief Class for a single RGB value.

    This class contains three values (of the specified type) that represent
    red, green, and blue color channels. By means of the template parameters
    <tt>RED_IDX, GREEN_IDX, BLUE_IDX</tt>, the indices 0, 1, 2 can be assigned to
    the three colors arbitrarily, so that, for example, a BGR type can be created
    as

    \code
    typedef RGBValue<unsigned char, 2,1,0> BGRValue;
    \endcode

    The standard order red=0, green=1, blue=2 is the default. There are three possibilities
    to access the color values: accessor functions (\ref red(), \ref green(),
    \ref blue()), index operator (operator[](dx), where the <tt>rgb[RED_IDX]</tt>
    returns red etc.) and iterator (STL-compatible random access
    iterator that references the three colors in turn). The latter two
    methods, together with the necessary embedded typedefs, ensure
    compatibility of a RGBValue with a STL vector.

    \ref RGBValueOperators "Arithmetic operations" are defined as component-wise applications of these
    operations. Addition, subtraction, and multiplication of two RGBValues
    (+=, -=, *=, +, -, *, unary -), multiplication and division of an
    RGBValue with a double, and NumericTraits/PromoteTraits are defined,
    so that RGBValue fulfills the requirements of a \ref LinearAlgebraConcept "Linear Algebra".

    A number of \ref RGBValueAccessors "accessors" are provided
    that support access to RGBValues as a whole, to a selected
    color component, or to the luminance value.

    <b>\#include</b> \<vigra/rgbvalue.hxx\><br>
    Namespace: vigra
*/
template <class VALUETYPE, unsigned int RED_IDX = 0, unsigned int GREEN_IDX = 1, unsigned int BLUE_IDX = 2>
class RGBValue
: public TinyVector<VALUETYPE, 3>
{
    typedef TinyVector<VALUETYPE, 3> Base;

        // inverse mapping from index to color
    enum {
      IDX0 = (RED_IDX == 0) ? 0 : (GREEN_IDX == 0) ? 1 : 2,
      IDX1 = (RED_IDX == 1) ? 0 : (GREEN_IDX == 1) ? 1 : 2,
      IDX2 = (RED_IDX == 2) ? 0 : (GREEN_IDX == 2) ? 1 : 2
    };

  public:
        /** STL-compatible definition of valuetype
        */
    typedef typename Base::value_type value_type;
        /** STL-compatible definition of iterator
        */
    typedef typename Base::iterator iterator;
        /** STL-compatible definition of const iterator
        */
    typedef typename Base::const_iterator const_iterator;
        /** squared norm type (result of squaredManitude())
        */
    typedef typename Base::SquaredNormType SquaredNormType;
        /** norm type (result of magnitude())
        */
    typedef typename Base::NormType NormType;

    typedef typename Base::reference reference;
    typedef typename Base::const_reference const_reference;
    typedef typename Base::pointer pointer;
    typedef typename Base::const_pointer const_pointer;
    typedef typename Base::size_type size_type;
    typedef typename Base::difference_type difference_type;
    typedef typename Base::scalar_multiplier scalar_multiplier;
    typedef typename Base::ReverseCopyTag ReverseCopyTag;

        /** Color index positions
        */
    enum
    {
      RedIdx = RED_IDX,
      GreenIdx = GREEN_IDX,
      BlueIdx = BLUE_IDX
    };

        /** Construct from explicit color values.
            \a first, \a second, \a third are written in this order,
            irrespective of how the color indices are specified.
        */
    RGBValue(value_type first, value_type second, value_type third)
    : Base(first, second, third)
    {
        VIGRA_STATIC_ASSERT((RGBValue_bad_color_indices<RED_IDX, GREEN_IDX, BLUE_IDX>));
    }

        /** Construct gray value.
        */
    RGBValue(value_type gray)
    : Base(gray, gray, gray)
    {
        VIGRA_STATIC_ASSERT((RGBValue_bad_color_indices<RED_IDX, GREEN_IDX, BLUE_IDX>));
    }

        /** Copy from raw memory. The order is preserved,
            irrespective of how the color indices are specified.
        */
    explicit RGBValue(const_pointer i)
    : Base(i)
    {
        VIGRA_STATIC_ASSERT((RGBValue_bad_color_indices<RED_IDX, GREEN_IDX, BLUE_IDX>));
    }

        /** Construct by reverse copying from raw memory.
        */
    RGBValue(const_pointer i, ReverseCopyTag reverse)
    : Base(i, reverse)
    {
        VIGRA_STATIC_ASSERT((RGBValue_bad_color_indices<RED_IDX, GREEN_IDX, BLUE_IDX>));
    }

        /** Default constructor (sets all components to 0)
        */
    RGBValue()
    : Base(0, 0, 0)
    {
        VIGRA_STATIC_ASSERT((RGBValue_bad_color_indices<RED_IDX, GREEN_IDX, BLUE_IDX>));
    }

#if !defined(TEMPLATE_COPY_CONSTRUCTOR_BUG)

    RGBValue(RGBValue const & r)
    : Base((Base const &)r)
    {
        VIGRA_STATIC_ASSERT((RGBValue_bad_color_indices<RED_IDX, GREEN_IDX, BLUE_IDX>));
    }

    RGBValue & operator=(RGBValue const & r)
    {
        Base::operator=(r);
        return *this;
    }

#endif // TEMPLATE_COPY_CONSTRUCTOR_BUG

        /** Copy constructor.
        */
    template <class U, unsigned int R, unsigned int G, unsigned int B>
    RGBValue(RGBValue<U, R, G, B> const & r)
    : Base(detail::RequiresExplicitCast<value_type>::cast(r[detail::SelectColorIndexRHS<IDX0, R, G, B>::res]),
           detail::RequiresExplicitCast<value_type>::cast(r[detail::SelectColorIndexRHS<IDX1, R, G, B>::res]),
           detail::RequiresExplicitCast<value_type>::cast(r[detail::SelectColorIndexRHS<IDX2, R, G, B>::res]))
    {
        VIGRA_STATIC_ASSERT((RGBValue_bad_color_indices<RED_IDX, GREEN_IDX, BLUE_IDX>));
    }

        /** Copy assignment.
        */
    template <class U, unsigned int R, unsigned int G, unsigned int B>
    RGBValue & operator=(RGBValue<U, R, G, B> const & r)
    {
        setRed(detail::RequiresExplicitCast<value_type>::cast(r.red()));
        setGreen(detail::RequiresExplicitCast<value_type>::cast(r.green()));
        setBlue(detail::RequiresExplicitCast<value_type>::cast(r.blue()));
        return *this;
    }

        /** construct from TinyVector
        */
    RGBValue(TinyVector<value_type, 3> const & r)
    : Base(r)
    {
        VIGRA_STATIC_ASSERT((RGBValue_bad_color_indices<RED_IDX, GREEN_IDX, BLUE_IDX>));
    }

        /** assign TinyVector.
        */
    RGBValue & operator=(TinyVector<value_type, 3> const & r)
    {
        Base::operator=(r);
        return *this;
    }

        /** Unary negation (construct RGBValue with negative values)
        */
    RGBValue operator-() const
    {
        return RGBValue(-(*this)[0], -(*this)[1], -(*this)[2]);
    }

        /** Access red component.
        */
    value_type & red() { return (*this)[RED_IDX]; }

        /** Access green component.
        */
    value_type & green() { return (*this)[GREEN_IDX]; }

        /** Access blue component.
        */
    value_type & blue() { return (*this)[BLUE_IDX]; }

        /** Get red component.
        */
    value_type const & red() const { return (*this)[RED_IDX]; }

        /** Get green component.
        */
    value_type const & green() const { return (*this)[GREEN_IDX]; }

        /** Get blue component.
        */
    value_type const & blue() const { return (*this)[BLUE_IDX]; }

        /** Calculate luminance.
        */
    value_type luminance() const {
         return detail::RequiresExplicitCast<value_type>::cast(0.3*red() + 0.59*green() + 0.11*blue()); }

        /** Calculate magnitude.
        */
    NormType magnitude() const {
         return Base::magnitude();
    }

        /** Calculate squared magnitude.
        */
    SquaredNormType squaredMagnitude() const {
         return Base::squaredMagnitude();
    }

        /** Set red component. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>VALUETYPE</TT>.
        */
    template <class V>
    void setRed(V value) { (*this)[RED_IDX] = detail::RequiresExplicitCast<value_type>::cast(value); }

        /** Set green component.The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>VALUETYPE</TT>.
        */
    template <class V>
    void setGreen(V value) { (*this)[GREEN_IDX] = detail::RequiresExplicitCast<value_type>::cast(value); }

        /** Set blue component.The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>VALUETYPE</TT>.
        */
    template <class V>
    void setBlue(V value) { (*this)[BLUE_IDX] = detail::RequiresExplicitCast<value_type>::cast(value); }


    template <class V>
    void setRGB(V r, V g, V b)
    {
        (*this)[RED_IDX] = detail::RequiresExplicitCast<value_type>::cast(r);
        (*this)[GREEN_IDX] = detail::RequiresExplicitCast<value_type>::cast(g);
        (*this)[BLUE_IDX] = detail::RequiresExplicitCast<value_type>::cast(b);
    }
};

/********************************************************/
/*                                                      */
/*                     RGBValue Comparison              */
/*                                                      */
/********************************************************/

/** \addtogroup RGBValueOperators Functions for RGBValue

    \brief Implement basic arithmetic and equality for RGBValue.

    These functions fulfill the requirements of a Linear Algebra.
    Return types are determined according to \ref RGBValueTraits.

    <b>\#include</b> \<vigra/rgbvalue.hxx\><br>
    Namespace: vigra
    <p>

 */
//@{
    /// component-wise equal
template <class V1, unsigned int RIDX1, unsigned int GIDX1, unsigned int BIDX1,
          class V2, unsigned int RIDX2, unsigned int GIDX2, unsigned int BIDX2>
inline
bool
operator==(RGBValue<V1, RIDX1, GIDX1, BIDX1> const & l,
           RGBValue<V2, RIDX2, GIDX2, BIDX2> const & r)
{
    return (l.red() == r.red()) &&
           (l.green() == r.green()) &&
           (l.blue() == r.blue());
}

    /// component-wise not equal
template <class V1, unsigned int RIDX1, unsigned int GIDX1, unsigned int BIDX1,
          class V2, unsigned int RIDX2, unsigned int GIDX2, unsigned int BIDX2>
inline
bool
operator!=(RGBValue<V1, RIDX1, GIDX1, BIDX1> const & l,
           RGBValue<V2, RIDX2, GIDX2, BIDX2> const & r)
{
    return (l.red() != r.red()) ||
           (l.green() != r.green()) ||
           (l.blue() != r.blue());
}


//@}

/********************************************************/
/*                                                      */
/*                      RGBValue-Traits                 */
/*                                                      */
/********************************************************/

/** \page RGBValueTraits Numeric and Promote Traits of RGBValue
    The numeric and promote traits for RGBValues follow
    the general specifications for \ref NumericPromotionTraits.
    They are implemented in terms of the traits of the basic types by
    partial template specialization. Note that PromoteTraits are only defined
    for the case that the color indices are the same in both RGBValues.

    \code

    template <class T, unsigned int R, unsigned int G, unsigned int B>
    struct NumericTraits<RGBValue<T, R, G, B> >
    {
        typedef RGBValue<T, R, G, B> Type;
        typedef RGBValue<typename NumericTraits<T>::Promote, R, G, B> Promote;
        typedef RGBValue<typename NumericTraits<T>::RealPromote, R, G, B> RealPromote;
        typedef RGBValue<typename NumericTraits<T>::ComplexPromote, R, G, B> ComplexPromote;
        typedef T ValueType;

        typedef typename NumericTraits<T>::isIntegral isIntegral;
        typedef VigraFalseType isScalar;
        typedef typename NumericTraits<T>::isSigned isSigned;

        // etc.
    };

    template <class T, unsigned int R, unsigned int G, unsigned int B>
    struct NormTraits<RGBValue<T, R, G, B> >
    {
        typedef RGBValue<T, R, G, B> Type;
        typedef typename Type::SquaredNormType    SquaredNormType;
        typedef typename Type::NormType           NormType;
    };

    template <class T1, unsigned int R, unsigned int G, unsigned int B, class T2>
    struct PromoteTraits<RGBValue<T1, R, G, B>, RGBValue<T2, R, G, B> >
    {
        typedef RGBValue<typename PromoteTraits<T1, T2>::Promote, R, G, B> Promote;
    };

    template <class T, unsigned int R, unsigned int G, unsigned int B>
    struct PromoteTraits<RGBValue<T, R, G, B>, double >
    {
        typedef RGBValue<typename NumericTraits<T>::RealPromote, R, G, B> Promote;
    };

    template <class T, unsigned int R, unsigned int G, unsigned int B>
    struct PromoteTraits<double, RGBValue<T, R, G, B> >
    {
        typedef RGBValue<typename NumericTraits<T>::RealPromote, R, G, B> Promote;
    };
    \endcode

    <b>\#include</b> \<vigra/rgbvalue.hxx\><br>
    Namespace: vigra

*/

#if !defined(NO_PARTIAL_TEMPLATE_SPECIALIZATION)

template <class T, unsigned int R, unsigned int G, unsigned int B>
struct NumericTraits<RGBValue<T, R, G, B> >
{
    typedef RGBValue<T, R, G, B> Type;
    typedef RGBValue<typename NumericTraits<T>::Promote, R, G, B> Promote;
    typedef RGBValue<typename NumericTraits<T>::RealPromote, R, G, B> RealPromote;
    typedef RGBValue<typename NumericTraits<T>::ComplexPromote, R, G, B> ComplexPromote;
    typedef T ValueType;

    typedef typename NumericTraits<T>::isIntegral isIntegral;
    typedef VigraFalseType isScalar;
    typedef typename NumericTraits<T>::isSigned isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;

    static Type zero()
    {
        return Type(NumericTraits<T>::zero());
    }
    static Type one()
    {
        return Type(NumericTraits<T>::one());
    }
    static Type nonZero()
    {
        return Type(NumericTraits<T>::nonZero());
    }

    static Type min()
    {
        return Type(NumericTraits<T>::min());
    }
    static Type max()
    {
        return Type(NumericTraits<T>::max());
    }

    static Promote toPromote(Type const & v)
    {
        return Promote(v);
    }
    static RealPromote toRealPromote(Type const & v)
    {
        return RealPromote(v);
    }
    static Type fromPromote(Promote const & v)
    {
      return Type(NumericTraits<T>::fromPromote(v.red()),
                  NumericTraits<T>::fromPromote(v.green()),
                  NumericTraits<T>::fromPromote(v.blue()));
    }
    static Type fromRealPromote(RealPromote const & v)
    {
        return Type(NumericTraits<T>::fromRealPromote(v.red()),
                    NumericTraits<T>::fromRealPromote(v.green()),
                    NumericTraits<T>::fromRealPromote(v.blue()));
    }
};

template <class T, unsigned int R, unsigned int G, unsigned int B>
struct NormTraits<RGBValue<T, R, G, B> >
{
    typedef RGBValue<T, R, G, B> Type;
    typedef typename Type::SquaredNormType    SquaredNormType;
    typedef typename Type::NormType           NormType;
};

template <class T1, unsigned int R, unsigned int G, unsigned int B, class T2>
struct PromoteTraits<RGBValue<T1, R, G, B>, RGBValue<T2, R, G, B> >
{
    typedef RGBValue<typename PromoteTraits<T1, T2>::Promote, R, G, B> Promote;
};

template <class T, unsigned int R, unsigned int G, unsigned int B>
struct PromoteTraits<RGBValue<T, R, G, B>, double >
{
    typedef RGBValue<typename NumericTraits<T>::RealPromote, R, G, B> Promote;
};

template <class T, unsigned int R, unsigned int G, unsigned int B>
struct PromoteTraits<double, RGBValue<T, R, G, B> >
{
    typedef RGBValue<typename NumericTraits<T>::RealPromote, R, G, B> Promote;
};

template<class T, unsigned int R, unsigned int G, unsigned int B>
struct CanSkipInitialization<RGBValue<T, R, G, B> >
{
    typedef typename CanSkipInitialization<T>::type type;
    static const bool value = type::asBool;
};


#else // NO_PARTIAL_TEMPLATE_SPECIALIZATION

#define RGBVALUE_NUMTRAITS(T) \
template<>\
struct NumericTraits<RGBValue<T, 0, 1, 2> >\
{\
    typedef RGBValue<T> Type; \
    typedef RGBValue<NumericTraits<T>::Promote> Promote; \
    typedef RGBValue<NumericTraits<T>::RealPromote> RealPromote; \
    typedef RGBValue<NumericTraits<T>::ComplexPromote> ComplexPromote; \
    typedef T ValueType; \
    \
    typedef NumericTraits<T>::isIntegral isIntegral; \
    typedef VigraFalseType isScalar; \
    typedef NumericTraits<T>::isSigned isSigned; \
    typedef VigraFalseType isOrdered; \
    typedef VigraFalseType isComplex; \
    \
    static RGBValue<T> zero() { \
        return RGBValue<T>(NumericTraits<T>::zero()); \
    }\
    static RGBValue<T> one() { \
        return RGBValue<T>(NumericTraits<T>::one()); \
    }\
    static RGBValue<T> nonZero() { \
        return RGBValue<T>(NumericTraits<T>::nonZero()); \
    }\
    \
    static Promote toPromote(RGBValue<T> const & v) { \
        return Promote(v); \
    }\
    static RealPromote toRealPromote(RGBValue<T> const & v) { \
        return RealPromote(v); \
    }\
    static RGBValue<T> fromPromote(Promote const & v) { \
        RGBValue<T> res;\
        RGBValue<T>::iterator d = res.begin();\
        Promote::const_iterator s = v.begin();\
        for(; d != res.end(); ++d, ++s)\
            *d = NumericTraits<T>::fromPromote(*s);\
        return res;\
    }\
    static RGBValue<T> fromRealPromote(RealPromote const & v) {\
        RGBValue<T> res;\
        RGBValue<T>::iterator d = res.begin();\
        RealPromote::const_iterator s = v.begin();\
        for(; d != res.end(); ++d, ++s)\
            *d = NumericTraits<T>::fromRealPromote(*s);\
        return res;\
    }\
}; \
template<>\
struct NormTraits<RGBValue<T, 0, 1, 2> >\
{\
    typedef RGBValue<T> Type;\
    typedef Type::SquaredNormType           SquaredNormType; \
    typedef Type::NormType NormType; \
};

#define RGBVALUE_PROMTRAITS1(type1) \
template<> \
struct PromoteTraits<RGBValue<type1, 0, 1, 2>, RGBValue<type1, 0, 1, 2> > \
{ \
    typedef RGBValue<PromoteTraits<type1, type1>::Promote> Promote; \
    static Promote toPromote(RGBValue<type1> const & v) { \
        return static_cast<Promote>(v); } \
}; \
template <> \
struct PromoteTraits<RGBValue<type1, 0, 1, 2>, double > \
{ \
    typedef RGBValue<typename NumericTraits<type1>::RealPromote> Promote; \
}; \
template <> \
struct PromoteTraits<double, RGBValue<type1, 0, 1, 2> > \
{ \
    typedef RGBValue<typename NumericTraits<type1>::RealPromote> Promote; \
};

#define RGBVALUE_PROMTRAITS2(type1, type2) \
template<> \
struct PromoteTraits<RGBValue<type1, 0, 1, 2>, RGBValue<type2, 0, 1, 2> > \
{ \
    typedef RGBValue<PromoteTraits<type1, type2>::Promote> Promote; \
    static Promote toPromote(RGBValue<type1> const & v) { \
        return static_cast<Promote>(v); } \
    static Promote toPromote(RGBValue<type2> const & v) { \
        return static_cast<Promote>(v); } \
};

RGBVALUE_NUMTRAITS(unsigned char)
RGBVALUE_NUMTRAITS(int)
RGBVALUE_NUMTRAITS(float)
RGBVALUE_NUMTRAITS(double)
RGBVALUE_PROMTRAITS1(unsigned char)
RGBVALUE_PROMTRAITS1(int)
RGBVALUE_PROMTRAITS1(float)
RGBVALUE_PROMTRAITS1(double)
RGBVALUE_PROMTRAITS2(float, unsigned char)
RGBVALUE_PROMTRAITS2(unsigned char, float)
RGBVALUE_PROMTRAITS2(int, unsigned char)
RGBVALUE_PROMTRAITS2(unsigned char, int)
RGBVALUE_PROMTRAITS2(int, float)
RGBVALUE_PROMTRAITS2(float, int)
RGBVALUE_PROMTRAITS2(double, unsigned char)
RGBVALUE_PROMTRAITS2(unsigned char, double)
RGBVALUE_PROMTRAITS2(int, double)
RGBVALUE_PROMTRAITS2(double, int)
RGBVALUE_PROMTRAITS2(double, float)
RGBVALUE_PROMTRAITS2(float, double)

#undef RGBVALUE_NUMTRAITS
#undef RGBVALUE_PROMTRAITS1
#undef RGBVALUE_PROMTRAITS2

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION


/********************************************************/
/*                                                      */
/*                      RGBValue-Arithmetic             */
/*                                                      */
/********************************************************/

/** \addtogroup RGBValueOperators
 */
//@{
    /// componentwise add-assignment
template <class V1, unsigned int RIDX1, unsigned int GIDX1, unsigned int BIDX1,
          class V2, unsigned int RIDX2, unsigned int GIDX2, unsigned int BIDX2>
inline
RGBValue<V1, RIDX1, GIDX1, BIDX1> &
operator+=(RGBValue<V1, RIDX1, GIDX1, BIDX1> & l,
           RGBValue<V2, RIDX2, GIDX2, BIDX2> const & r)
{
    l.red() += r.red();
    l.green() += r.green();
    l.blue() += r.blue();
    return l;
}

    /// componentwise subtract-assignment
template <class V1, unsigned int RIDX1, unsigned int GIDX1, unsigned int BIDX1,
          class V2, unsigned int RIDX2, unsigned int GIDX2, unsigned int BIDX2>
inline
RGBValue<V1, RIDX1, GIDX1, BIDX1> &
operator-=(RGBValue<V1, RIDX1, GIDX1, BIDX1> & l,
           RGBValue<V2, RIDX2, GIDX2, BIDX2> const & r)
{
    l.red() -= r.red();
    l.green() -= r.green();
    l.blue() -= r.blue();
    return l;
}

    /// componentwise multiply-assignment
template <class V1, unsigned int RIDX1, unsigned int GIDX1, unsigned int BIDX1,
          class V2, unsigned int RIDX2, unsigned int GIDX2, unsigned int BIDX2>
inline
RGBValue<V1, RIDX1, GIDX1, BIDX1> &
operator*=(RGBValue<V1, RIDX1, GIDX1, BIDX1> & l,
           RGBValue<V2, RIDX2, GIDX2, BIDX2> const & r)
{
    l.red() = V1(l.red() * r.red());
    l.green() = V1(l.green() * r.green());
    l.blue() = V1(l.blue() * r.blue());
    return l;
}

    /// componentwise scalar multiply-assignment
template <class V, unsigned int RIDX, unsigned int GIDX, unsigned int BIDX>
inline
RGBValue<V, RIDX, GIDX, BIDX> &
operator*=(RGBValue<V, RIDX, GIDX, BIDX> & l, double r)
{
    l.red() = V(l.red() * r);
    l.green() = V(l.green() * r);
    l.blue() = V(l.blue() * r);
    return l;
}

    /// componentwise divide-assignment
template <class V1, unsigned int RIDX1, unsigned int GIDX1, unsigned int BIDX1,
          class V2, unsigned int RIDX2, unsigned int GIDX2, unsigned int BIDX2>
inline
RGBValue<V1, RIDX1, GIDX1, BIDX1> &
operator/=(RGBValue<V1, RIDX1, GIDX1, BIDX1> & l,
           RGBValue<V2, RIDX2, GIDX2, BIDX2> const & r)
{
    l.red() = V1(l.red() / r.red());
    l.green() = V1(l.green() / r.green());
    l.blue() = V1(l.blue() / r.blue());
    return l;
}

    /// componentwise scalar divide-assignment
template <class V, unsigned int RIDX, unsigned int GIDX, unsigned int BIDX>
inline
RGBValue<V, RIDX, GIDX, BIDX> &
operator/=(RGBValue<V, RIDX, GIDX, BIDX> & l, double r)
{
    l.red() = V(l.red() / r);
    l.green() = V(l.green() / r);
    l.blue() = V(l.blue() / r);
    return l;
}

using VIGRA_CSTD::abs;

    /// component-wise absolute value
template <class T, unsigned int RIDX, unsigned int GIDX, unsigned int BIDX>
inline
RGBValue<T, RIDX, GIDX, BIDX>
abs(RGBValue<T, RIDX, GIDX, BIDX> const & v)
{
  return RGBValue<T, RIDX, GIDX, BIDX>(abs(v.red()), abs(v.green()), abs(v.blue()));
}

    /// component-wise addition
template <class V1, unsigned int R, unsigned int G, unsigned int B, class V2>
inline
typename PromoteTraits<RGBValue<V1, R, G, B>,
                       RGBValue<V2, R, G, B> >::Promote
operator+(RGBValue<V1, R, G, B> const & r1,
          RGBValue<V2, R, G, B> const & r2)
{
    typename PromoteTraits<RGBValue<V1, R, G, B>,
                           RGBValue<V2, R, G, B> >::Promote res(r1);

    res += r2;

    return res;
}

    /// component-wise subtraction
template <class V1, unsigned int R, unsigned int G, unsigned int B, class V2>
inline
typename PromoteTraits<RGBValue<V1, R, G, B>,
                       RGBValue<V2, R, G, B> >::Promote
operator-(RGBValue<V1, R, G, B> const & r1,
          RGBValue<V2, R, G, B> const & r2)
{
    typename PromoteTraits<RGBValue<V1, R, G, B>,
                           RGBValue<V2, R, G, B> >::Promote res(r1);

    res -= r2;

    return res;
}

    /// component-wise multiplication
template <class V1, unsigned int R, unsigned int G, unsigned int B, class V2>
inline
typename PromoteTraits<RGBValue<V1, R, G, B>,
                       RGBValue<V2, R, G, B> >::Promote
operator*(RGBValue<V1, R, G, B> const & r1,
          RGBValue<V2, R, G, B> const & r2)
{
    typename PromoteTraits<RGBValue<V1, R, G, B>,
                           RGBValue<V2, R, G, B> >::Promote res(r1);

    res *= r2;

    return res;
}

    /// component-wise left scalar multiplication
template <class V, unsigned int R, unsigned int G, unsigned int B>
inline
typename NumericTraits<RGBValue<V, R, G, B> >::RealPromote
operator*(double v, RGBValue<V, R, G, B> const & r)
{
    typename NumericTraits<RGBValue<V, R, G, B> >::RealPromote res(r);

    res *= v;

    return res;
}

    /// component-wise right scalar multiplication
template <class V, unsigned int R, unsigned int G, unsigned int B>
inline
typename NumericTraits<RGBValue<V, R, G, B> >::RealPromote
operator*(RGBValue<V, R, G, B> const & r, double v)
{
    typename NumericTraits<RGBValue<V, R, G, B> >::RealPromote res(r);

    res *= v;

    return res;
}

    /// component-wise division
template <class V1, unsigned int R, unsigned int G, unsigned int B, class V2>
inline
typename PromoteTraits<RGBValue<V1, R, G, B>,
                       RGBValue<V2, R, G, B> >::Promote
operator/(RGBValue<V1, R, G, B> const & r1,
          RGBValue<V2, R, G, B> const & r2)
{
    typename PromoteTraits<RGBValue<V1, R, G, B>,
                           RGBValue<V2, R, G, B> >::Promote res(r1);

    res /= r2;

    return res;
}

    /// component-wise scalar division
template <class V, unsigned int R, unsigned int G, unsigned int B>
inline
typename NumericTraits<RGBValue<V, R, G, B> >::RealPromote
operator/(RGBValue<V, R, G, B> const & r, double v)
{
    typename NumericTraits<RGBValue<V, R, G, B> >::RealPromote res(r);

    res /= v;

    return res;
}

    /// cross product
template <class V1, unsigned int R, unsigned int G, unsigned int B, class V2>
inline
typename PromoteTraits<RGBValue<V1, R, G, B>,
                       RGBValue<V2, R, G, B> >::Promote
cross(RGBValue<V1, R, G, B> const & r1,
      RGBValue<V2, R, G, B> const & r2)
{
    typedef typename PromoteTraits<RGBValue<V1, R, G, B>,
                                   RGBValue<V2, R, G, B> >::Promote
            Res;

    return  Res(r1.green()*r2.blue() - r1.blue()*r2.green(),
                r1.blue()*r2.red() - r1.red()*r2.blue(),
                r1.red()*r2.green() - r1.green()*r2.red());
}

    /// dot product
template <class V1, unsigned int RIDX1, unsigned int GIDX1, unsigned int BIDX1,
          class V2, unsigned int RIDX2, unsigned int GIDX2, unsigned int BIDX2>
inline
typename PromoteTraits<V1, V2>::Promote
dot(RGBValue<V1, RIDX1, GIDX1, BIDX1> const & r1,
    RGBValue<V2, RIDX2, GIDX2, BIDX2> const & r2)
{
    return r1.red()*r2.red() + r1.green()*r2.green() + r1.blue()*r2.blue();
}

using VIGRA_CSTD::ceil;

    /** Apply ceil() function to each RGB component.
    */
template <class V, unsigned int RIDX, unsigned int GIDX, unsigned int BIDX>
inline
RGBValue<V, RIDX, GIDX, BIDX>
ceil(RGBValue<V, RIDX, GIDX, BIDX> const & r)
{
    return RGBValue<V, RIDX, GIDX, BIDX>(ceil(r.red()),
                                         ceil(r.green()),
                                         ceil(r.blue()));
}

using VIGRA_CSTD::floor;

    /** Apply floor() function to each RGB component.
    */
template <class V, unsigned int RIDX, unsigned int GIDX, unsigned int BIDX>
inline
RGBValue<V, RIDX, GIDX, BIDX>
floor(RGBValue<V, RIDX, GIDX, BIDX> const & r)
{
    return RGBValue<V, RIDX, GIDX, BIDX>(floor(r.red()),
                                         floor(r.green()),
                                         floor(r.blue()));
}

//@}

/********************************************************/
/*                                                      */
/*                      RGBValue-Accessors              */
/*                                                      */
/********************************************************/

/** \addtogroup DataAccessors
*/
//@{
/** \defgroup RGBValueAccessors Accessors for RGBValue */
//@{
    /** Encapsulate access to rgb values.

    <b>\#include</b> \<vigra/rgbvalue.hxx\><br>
    Namespace: vigra
    */
template <class RGBVALUE>
class RGBAccessor
: public VectorAccessor<RGBVALUE>
{
  public:

    typedef typename RGBVALUE::value_type component_type;

        /** Get value of the red component
        */
    template <class RGBIterator>
    component_type const & red(RGBIterator const & rgb) const
    {
        return (*rgb).red();
    }

    template <class V, class RGBIterator>
    void setRGB(V r, V g, V b, RGBIterator const & rgb) const
    {
        (*rgb).setRGB( r, g, b );
    }


        /** Set value of the red component. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>component_type</TT>.
        */
    template <class V, class RGBIterator>
    void setRed(V value, RGBIterator const & rgb) const
    {
        (*rgb).setRed(value);
    }

        /** Get value of the red component at an offset
        */
    template <class RGBIterator, class DIFFERENCE>
    component_type const & red(RGBIterator const & rgb, DIFFERENCE diff) const
    {
        return rgb[diff].red();
    }

        /** Set value of the red component at an offset. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>component_type</TT>.
        */
    template <class V, class RGBIterator, class DIFFERENCE>
    void setRed(V value, RGBIterator const & rgb, DIFFERENCE diff) const
    {
        rgb[diff].setRed(value);
    }

        /** Get value of the green component
        */
    template <class RGBIterator>
    component_type const & green(RGBIterator const & rgb) const
    {
        return (*rgb).green();
    }

        /** Set value of the green component. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>component_type</TT>.
        */
    template <class V, class RGBIterator>
    void setGreen(V value, RGBIterator const & rgb) const
    {
        (*rgb).setGreen(value);
    }

        /** Get value of the green component at an offset
        */
    template <class RGBIterator, class DIFFERENCE>
    component_type const & green(RGBIterator const & rgb, DIFFERENCE d) const
    {
        return rgb[d].green();
    }

        /** Set value of the green component at an offset. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>component_type</TT>.
        */
    template <class V, class RGBIterator, class DIFFERENCE>
    void setGreen(V value, RGBIterator const & rgb, DIFFERENCE d) const
    {
        rgb[d].setGreen(value);
    }

        /** Get value of the blue component
        */
    template <class RGBIterator>
    component_type const & blue(RGBIterator const & rgb) const
    {
        return (*rgb).blue();
    }

        /** Set value of the blue component. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>component_type</TT>.
        */
    template <class V, class RGBIterator>
    void setBlue(V value, RGBIterator const & rgb) const
    {
        (*rgb).setBlue(value);
    }

        /** Get value of the blue component at an offset
        */
    template <class RGBIterator, class DIFFERENCE>
    component_type const & blue(RGBIterator const & rgb, DIFFERENCE d) const
    {
        return rgb[d].blue();
    }

        /** Set value of the blue component at an offset. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>component_type</TT>.
        */
    template <class V, class RGBIterator, class DIFFERENCE>
    void setBlue(V value, RGBIterator const & rgb, DIFFERENCE d) const
    {
        rgb[d].setBlue(value);
    }

};


/********************************************************/
/*                                                      */
/*                       RedAccessor                    */
/*                                                      */
/********************************************************/

    /** Encapsulate access to red band of an rgb value.

    <b>\#include</b> \<vigra/rgbvalue.hxx\><br>
    Namespace: vigra
    */
template <class RGBVALUE>
class RedAccessor
{
  public:
    typedef typename RGBVALUE::value_type value_type;

        /** Get value of the red component
        */
    template <class ITERATOR>
    value_type const & operator()(ITERATOR const & i) const {
        return (*i).red();
    }

        /** Get value of the red component at an offset
        */
    template <class ITERATOR, class DIFFERENCE>
    value_type const & operator()(ITERATOR const & i, DIFFERENCE d) const
    {
        return i[d].red();
    }

        /** Set value of the red component. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>value_type</TT>.
        */
    template <class V, class ITERATOR>
    void set(V value, ITERATOR const & i) const {
        (*i).setRed(value);
    }


        /** Set value of the red component at an offset. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>value_type</TT>.
        */
    template <class V, class ITERATOR, class DIFFERENCE>
    void set(V value, ITERATOR const & i, DIFFERENCE d) const
    {
        i[d].setRed(value);
    }
};

/********************************************************/
/*                                                      */
/*                     GreenAccessor                    */
/*                                                      */
/********************************************************/

    /** Encapsulate access to green band of an rgb value.

    <b>\#include</b> \<vigra/rgbvalue.hxx\><br>
    Namespace: vigra
    */
template <class RGBVALUE>
class GreenAccessor
{
  public:
    typedef typename RGBVALUE::value_type value_type;

        /** Get value of the green component
        */
    template <class ITERATOR>
    value_type const & operator()(ITERATOR const & i) const {
        return (*i).green();
    }

        /** Get value of the green component at an offset
        */
    template <class ITERATOR, class DIFFERENCE>
    value_type const & operator()(ITERATOR const & i, DIFFERENCE d) const
    {
        return i[d].green();
    }

        /** Set value of the green component. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>value_type</TT>.
        */
    template <class V, class ITERATOR>
    void set(V value, ITERATOR const & i) const {
        (*i).setGreen(value);
    }


        /** Set value of the green component at an offset. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>value_type</TT>.
        */
    template <class V, class ITERATOR, class DIFFERENCE>
    void set(V value, ITERATOR const & i, DIFFERENCE d) const
    {
        i[d].setGreen(value);
    }
};

/********************************************************/
/*                                                      */
/*                     BlueAccessor                     */
/*                                                      */
/********************************************************/

    /** Encapsulate access to blue band of an rgb value.

    <b>\#include</b> \<vigra/rgbvalue.hxx\><br>
    Namespace: vigra
    */
template <class RGBVALUE>
class BlueAccessor
{
  public:
    typedef typename RGBVALUE::value_type value_type;

        /** Get value of the blue component
        */
    template <class ITERATOR>
    value_type const & operator()(ITERATOR const & i) const {
        return (*i).blue();
    }

        /** Get value of the blue component at an offset
        */
    template <class ITERATOR, class DIFFERENCE>
    value_type const & operator()(ITERATOR const & i, DIFFERENCE d) const
    {
        return i[d].blue();
    }

        /** Set value of the blue component. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>value_type</TT>.
        */
    template <class V, class ITERATOR>
    void set(V value, ITERATOR const & i) const {
        (*i).setBlue(value);
    }


        /** Set value of the blue component at an offset. The type <TT>V</TT> of the passed
            in <TT>value</TT> is automatically converted to <TT>value_type</TT>.
        */
    template <class V, class ITERATOR, class DIFFERENCE>
    void set(V value, ITERATOR const & i, DIFFERENCE d) const
    {
        i[d].setBlue(value);
    }
};

/********************************************************/
/*                                                      */
/*                  RGBToGrayAccessor                   */
/*                                                      */
/********************************************************/

    /** Encapsulate access to luminance of an rgb value.

    <b>\#include</b> \<vigra/rgbvalue.hxx\><br>
    Namespace: vigra
    */
template <class RGBVALUE>
class RGBToGrayAccessor
{
  public:
    typedef typename RGBVALUE::value_type value_type;

        /** Get value of the luminance
        */
    template <class ITERATOR>
    value_type operator()(ITERATOR const & i) const {
                return (*i).luminance(); }

        /** Get value of the luminance at an offset
        */
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR const & i, DIFFERENCE d) const
    {
        return i[d].luminance();
    }
};


/********************************************************/
/*                                                      */
/*                  GrayToRGBAccessor                   */
/*                                                      */
/********************************************************/

    /** Create an RGB view for a grayscale image by making all three channels
        equal.

    <b>\#include</b> \<vigra/rgbvalue.hxx\><br>
    Namespace: vigra
    */
template <class VALUETYPE>
class GrayToRGBAccessor
{
   public:
     typedef typename vigra::RGBValue<VALUETYPE> value_type;

         /** Get RGB value for the given pixel.
         */
     template <class ITERATOR>
     value_type operator()(ITERATOR const & i) const {
                 return value_type(*i,*i,*i); }

         /** Get RGB value at an offset
         */
     template <class ITERATOR, class DIFFERENCE>
     value_type operator()(ITERATOR const & i, DIFFERENCE d) const
     {
         return value_type(i[d],i[d],i[d]);
     }
};


//@}
//@}


} // namespace vigra

#endif // VIGRA_RGBVALUE_HXX
