/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2004 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    It was adapted from the file boost/rational.hpp of the            */
/*    boost library.                                                    */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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

// this file is based on work by Paul Moore:
//
//  (C) Copyright Paul Moore 1999. Permission to copy, use, modify, sell and
//  distribute this software is granted provided this copyright notice appears
//  in all copies. This software is provided "as is" without express or
//  implied warranty, and with no claim as to its suitability for any purpose.
//
//  See http://www.boost.org/libs/rational for documentation.


#ifndef VIGRA_RATIONAL_HPP
#define VIGRA_RATIONAL_HPP

#include <cmath>
#include <stdexcept>
#include <iosfwd>
#include "config.hxx"
#include "mathutil.hxx"
#include "numerictraits.hxx"
#include "metaprogramming.hxx"

namespace vigra {

/** \addtogroup MathFunctions Mathematical Functions
*/
//@{


/********************************************************/
/*                                                      */
/*                         gcd                          */
/*                                                      */
/********************************************************/

/*! Calculate the greatest common divisor.

    This function works for arbitrary integer types, including user-defined
    (e.g. infinite precision) ones.

    <b>\#include</b> \<<a href="rational_8hxx-source.html">vigra/rational.hxx</a>\><br>
    Namespace: vigra
*/
template <typename IntType>
IntType gcd(IntType n, IntType m)
{
    // Avoid repeated construction
    IntType zero(0);

    // This is abs() - given the existence of broken compilers with Koenig
    // lookup issues and other problems, I code this explicitly. (Remember,
    // IntType may be a user-defined type).
    if (n < zero)
        n = -n;
    if (m < zero)
        m = -m;

    // As n and m are now positive, we can be sure that %= returns a
    // positive value (the standard guarantees this for built-in types,
    // and we require it of user-defined types).
    for(;;) {
      if(m == zero)
        return n;
      n %= m;
      if(n == zero)
        return m;
      m %= n;
    }
}

/********************************************************/
/*                                                      */
/*                         lcm                          */
/*                                                      */
/********************************************************/

/*! Calculate the lowest common multiple.

    This function works for arbitrary integer types, including user-defined
    (e.g. infinite precision) ones.

    <b>\#include</b> \<<a href="rational_8hxx-source.html">vigra/rational.hxx</a>\><br>
    Namespace: vigra
*/
template <typename IntType>
IntType lcm(IntType n, IntType m)
{
    // Avoid repeated construction
    IntType zero(0);

    if (n == zero || m == zero)
        return zero;

    n /= gcd(n, m);
    n *= m;

    if (n < zero)
        n = -n;
    return n;
}

//@}

class bad_rational : public std::domain_error
{
public:
    explicit bad_rational() : std::domain_error("bad rational: zero denominator") {}
};

template <typename IntType>
class Rational;

template <typename IntType>
Rational<IntType> abs(const Rational<IntType>& r);
template <typename IntType>
Rational<IntType> pow(const Rational<IntType>& r, int n);
template <typename IntType>
Rational<IntType> floor(const Rational<IntType>& r);
template <typename IntType>
Rational<IntType> ceil(const Rational<IntType>& r);

/********************************************************/
/*                                                      */
/*                       Rational                       */
/*                                                      */
/********************************************************/

/** Template for rational numbers.

    This template can make use of arbitrary integer types, including
    user-defined (e.g. infinite precision) ones. Note, however,
    that overflow in either the numerator or denominator is not
    detected during calculations -- the standard behavior of the integer type
    (e.g. wrap around) applies.

    The class can represent and handle positive and negative infinity
    resulting from division by zero. Indeterminate expressions such as 0/0
    are signaled by a <tt>bad_rational</tt> exception which is derived from
    <tt>std::domain_error</tt>.

    <tt>Rational</tt> implements the required interface of an
    \ref AlgebraicField and the required \ref RationalTraits "numeric and
    promotion traits". All arithmetic and comparison operators, as well
    as the relevant algebraic functions are supported .

    <b>See also:</b>
    <ul>
    <li> \ref RationalTraits
    <li> \ref RationalOperations
    </ul>


    <b>\#include</b> \<<a href="rational_8hxx-source.html">vigra/rational.hxx</a>\><br>
    Namespace: vigra
*/
template <typename IntType>
class Rational
{
public:
        /** The type of numerator and denominator
        */
    typedef IntType value_type;

        /** Determine whether arguments should be passed as
            <tt>IntType</tt> or <tt>IntType const &</tt>.
        */
    typedef typename If<typename TypeTraits<IntType>::isBuiltinType,
                        IntType, IntType const &>::type param_type;

        /** Default constructor: creates zero (<tt>0/1</tt>)
        */
    Rational()
    : num(0),
      den(1)
    {}

        /** Copy constructor
        */
    template <class U>
    Rational(Rational<U> const & r)
    : num(r.numerator()),
      den(r.denominator())
    {}

        /** Integer constructor: creates <tt>n/1</tt>
        */
    Rational(param_type n)
    : num(n),
      den(IntType(1))
    {}

        /** Ratio constructor: creates <tt>n/d</tt>.

            The ratio will be normalized unless <tt>doNormalize = false</tt>.
            Since the internal representation is assumed to be normalized,
            <tt>doNormalize = false</tt> must only be used as an optimization
            if <tt>n</tt> and <tt>d</tt> are known to be already normalized
            (i.e. have 1 as their greatest common divisor).
        */
    Rational(param_type n, param_type d, bool doNormalize = true)
    : num(n),
      den(d)
    {
        if(doNormalize)
            normalize();
    }

        /** Construct as an approximation of a real number.

            The maximal allowed relative error is given by <tt>epsilon</tt>.
        */
    explicit Rational(double v, double epsilon = 1e-4)
    : num(IntType(v < 0.0 ?
                    v/epsilon - 0.5
                  : v/epsilon + 0.5)),
      den(IntType(1.0/epsilon + 0.5))
    {
        normalize();
    }

    // Default copy constructor and assignment are fine

        /** Assignment from <tt>IntType</tt>.
        */
    Rational& operator=(param_type n) { return assign(n, 1); }

        /** Assignment from <tt>IntType</tt> pair.
        */
    Rational& assign(param_type n, param_type d, bool doNormalize = true);

        /** Access numerator.
        */
    param_type numerator() const { return num; }

        /** Access denominator.
        */
    param_type denominator() const { return den; }

        /** Add-assignment from <tt>Rational</tt>

            <tt>throws bad_rational</tt> if indeterminate expression.
        */
    Rational& operator+= (const Rational& r);

        /** Subtract-assignment from <tt>Rational</tt>

            <tt>throws bad_rational</tt> if indeterminate expression.
        */
    Rational& operator-= (const Rational& r);

        /** Multiply-assignment from <tt>Rational</tt>

            <tt>throws bad_rational</tt> if indeterminate expression.
        */
    Rational& operator*= (const Rational& r);

        /** Divide-assignment from <tt>Rational</tt>

            <tt>throws bad_rational</tt> if indeterminate expression.
        */
    Rational& operator/= (const Rational& r);

        /** Add-assignment from <tt>IntType</tt>

            <tt>throws bad_rational</tt> if indeterminate expression.
        */
    Rational& operator+= (param_type i);

        /** Subtract-assignment from <tt>IntType</tt>

            <tt>throws bad_rational</tt> if indeterminate expression.
        */
    Rational& operator-= (param_type i);

        /** Multiply-assignment from <tt>IntType</tt>

            <tt>throws bad_rational</tt> if indeterminate expression.
        */
    Rational& operator*= (param_type i);

        /** Divide-assignment from <tt>IntType</tt>

            <tt>throws bad_rational</tt> if indeterminate expression.
        */
    Rational& operator/= (param_type i);

        /** Pre-increment.
        */
    Rational& operator++();
        /** Pre-decrement.
        */
    Rational& operator--();

        /** Post-increment.
        */
    Rational operator++(int) { Rational res(*this); operator++(); return res; }
        /** Post-decrement.
        */
    Rational operator--(int) { Rational res(*this); operator--(); return res; }

        /** Check for zero by calling <tt>!numerator()</tt>
        */
    bool operator!() const { return !num; }

        /** Check whether we have positive infinity.
        */
    bool is_pinf() const
    {
        IntType zero(0);
        return den == zero && num > zero;
    }

        /** Check whether we have negative infinity.
        */
    bool is_ninf() const
    {
        IntType zero(0);
        return den == zero && num < zero;
    }

        /** Check whether we have positive or negative infinity.
        */
    bool is_inf() const
    {
        IntType zero(0);
        return den == zero && num != zero;
    }

        /** Check the sign.

            Gives 1 if the number is positive, -1 if negative, and 0 otherwise.
        */
    int sign() const
    {
        IntType zero(0);
        return num == zero ? 0 : num < zero ? -1 : 1;
    }

private:
    // Implementation - numerator and denominator (normalized).
    // Other possibilities - separate whole-part, or sign, fields?
    IntType num;
    IntType den;

    // Representation note: Fractions are kept in normalized form at all
    // times. normalized form is defined as gcd(num,den) == 1 and den > 0.
    // In particular, note that the implementation of abs() below relies
    // on den always being positive.
    void normalize();
};

// Assign in place
template <typename IntType>
inline Rational<IntType>&
Rational<IntType>::assign(param_type n, param_type d, bool doNormalize)
{
    num = n;
    den = d;
    if(doNormalize)
        normalize();
    return *this;
}

// Arithmetic assignment operators
template <typename IntType>
Rational<IntType>& Rational<IntType>::operator+= (const Rational<IntType>& r)
{
    IntType zero(0);

    // handle the Inf and NaN cases
    if(den == zero)
    {
        if(r.den == zero && sign()*r.sign() < 0)
            throw bad_rational();
        return *this;
    }
    if(r.den == zero)
    {
        assign(r.num, zero, false); // Inf or -Inf
        return *this;
    }

    // This calculation avoids overflow, and minimises the number of expensive
    // calculations. Thanks to Nickolay Mladenov for this algorithm.
    //
    // Proof:
    // We have to compute a/b + c/d, where gcd(a,b)=1 and gcd(b,c)=1.
    // Let g = gcd(b,d), and b = b1*g, d=d1*g. Then gcd(b1,d1)=1
    //
    // The result is (a*d1 + c*b1) / (b1*d1*g).
    // Now we have to normalize this ratio.
    // Let's assume h | gcd((a*d1 + c*b1), (b1*d1*g)), and h > 1
    // If h | b1 then gcd(h,d1)=1 and hence h|(a*d1+c*b1) => h|a.
    // But since gcd(a,b1)=1 we have h=1.
    // Similarly h|d1 leads to h=1.
    // So we have that h | gcd((a*d1 + c*b1) , (b1*d1*g)) => h|g
    // Finally we have gcd((a*d1 + c*b1), (b1*d1*g)) = gcd((a*d1 + c*b1), g)
    // Which proves that instead of normalizing the result, it is better to
    // divide num and den by gcd((a*d1 + c*b1), g)

    // Protect against self-modification
    IntType r_num = r.num;
    IntType r_den = r.den;

    IntType g = gcd(den, r_den);
    den /= g;  // = b1 from the calculations above
    num = num * (r_den / g) + r_num * den;
    g = gcd(num, g);
    num /= g;
    den *= r_den/g;

    return *this;
}

template <typename IntType>
Rational<IntType>& Rational<IntType>::operator-= (const Rational<IntType>& r)
{
    IntType zero(0);

    // handle the Inf and NaN cases
    if(den == zero)
    {
        if(r.den == zero && sign()*r.sign() > 0)
            throw bad_rational();
        return *this;
    }
    if(r.den == zero)
    {
        assign(-r.num, zero, false); // Inf or -Inf
        return *this;
    }

    // Protect against self-modification
    IntType r_num = r.num;
    IntType r_den = r.den;

    // This calculation avoids overflow, and minimises the number of expensive
    // calculations. It corresponds exactly to the += case above
    IntType g = gcd(den, r_den);
    den /= g;
    num = num * (r_den / g) - r_num * den;
    g = gcd(num, g);
    num /= g;
    den *= r_den/g;

    return *this;
}

template <typename IntType>
Rational<IntType>& Rational<IntType>::operator*= (const Rational<IntType>& r)
{
    IntType zero(0);

    // handle the Inf and NaN cases
    if(den == zero)
    {
        if(r.num == zero)
            throw bad_rational();
        num *= r.sign();
        return *this;
    }
    if(r.den == zero)
    {
        if(num == zero)
            throw bad_rational();
        num = r.num * sign();
        den = zero;
        return *this;
    }

    // Protect against self-modification
    IntType r_num = r.num;
    IntType r_den = r.den;

    // Avoid overflow and preserve normalization
    IntType gcd1 = gcd<IntType>(num, r_den);
    IntType gcd2 = gcd<IntType>(r_num, den);
    num = (num/gcd1) * (r_num/gcd2);
    den = (den/gcd2) * (r_den/gcd1);
    return *this;
}

template <typename IntType>
Rational<IntType>& Rational<IntType>::operator/= (const Rational<IntType>& r)
{
    IntType zero(0);

    // handle the Inf and NaN cases
    if(den == zero)
    {
        if(r.den == zero)
            throw bad_rational();
        if(r.num != zero)
            num *= r.sign();
        return *this;
    }
    if(r.num == zero)
    {
        if(num == zero)
            throw bad_rational();
        num = IntType(sign());  // normalized inf!
        den = zero;
        return *this;
    }

    if (num == zero)
        return *this;

    // Protect against self-modification
    IntType r_num = r.num;
    IntType r_den = r.den;

    // Avoid overflow and preserve normalization
    IntType gcd1 = gcd<IntType>(num, r_num);
    IntType gcd2 = gcd<IntType>(r_den, den);
    num = (num/gcd1) * (r_den/gcd2);
    den = (den/gcd2) * (r_num/gcd1);

    if (den < zero) {
        num = -num;
        den = -den;
    }
    return *this;
}

// Mixed-mode operators -- implement explicitly to save gcd() calculations
template <typename IntType>
inline Rational<IntType>&
Rational<IntType>::operator+= (param_type i)
{
    num += i * den;
    return *this;
}

template <typename IntType>
inline Rational<IntType>&
Rational<IntType>::operator-= (param_type i)
{
    num -= i * den;
    return *this;
}

template <typename IntType>
Rational<IntType>&
Rational<IntType>::operator*= (param_type i)
{
    if(i == IntType(1))
        return *this;
    IntType zero(0);
    if(i == zero)
    {
        if(den == zero)
        {
            throw bad_rational();
        }
        else
        {
            num = zero;
            den = IntType(1);
            return *this;
        }
    }

    IntType g = gcd(i, den);
    den /= g;
    num *= i / g;
    return *this;
}

template <typename IntType>
Rational<IntType>&
Rational<IntType>::operator/= (param_type i)
{
    if(i == IntType(1))
        return *this;

    IntType zero(0);
    if(i == zero)
    {
        if(num == zero)
            throw bad_rational();
        num = IntType(sign()); // normalized inf!
        den = zero;
        return *this;
    }

    IntType g = gcd(i, num);
    if(i < zero)
    {
        num /= -g;
        den *= -i / g;
    }
    else
    {
        num /= g;
        den *= i / g;
    }
    return *this;
}

// Increment and decrement
template <typename IntType>
inline Rational<IntType>& Rational<IntType>::operator++()
{
    // This can never denormalise the fraction
    num += den;
    return *this;
}

template <typename IntType>
inline Rational<IntType>& Rational<IntType>::operator--()
{
    // This can never denormalise the fraction
    num -= den;
    return *this;
}

// Normalisation
template <typename IntType>
void Rational<IntType>::normalize()
{
    // Avoid repeated construction
    IntType zero(0);

    if (den == zero)
    {
        if(num == zero)
            throw bad_rational();
        if(num < zero)
            num = IntType(-1);
        else
            num = IntType(1);
        return;
    }

    // Handle the case of zero separately, to avoid division by zero
    if (num == zero) {
        den = IntType(1);
        return;
    }

    IntType g = gcd<IntType>(num, den);

    num /= g;
    den /= g;

    // Ensure that the denominator is positive
    if (den < zero) {
        num = -num;
        den = -den;
    }
}

/********************************************************/
/*                                                      */
/*                      Rational-Traits                 */
/*                                                      */
/********************************************************/

/** \page RationalTraits Numeric and Promote Traits of Rational

    The numeric and promote traits for Rational follow
    the general specifications for \ref NumericPromotionTraits and
    \ref AlgebraicField. They are implemented in terms of the traits of the basic types by
    partial template specialization:

    \code

    template <class T>
    struct NumericTraits<Rational<T> >
    {
        typedef Rational<typename NumericTraits<T>::Promote> Promote;
        typedef Rational<typename NumericTraits<T>::RealPromote> RealPromote;

        typedef typename NumericTraits<T>::isIntegral isIntegral;
        typedef VigraTrueType isScalar;
        typedef typename NumericTraits<T>::isSigned isSigned;
        typedef VigraTrueType isOrdered;

        // etc.
    };

    template<class T>
    struct NormTraits<Rational<T> >
    {
        typedef Rational<T>                           Type;
        typedef typename NumericTraits<Type>::Promote SquaredNormType;
        typedef Type                                  NormType;
    };

    template <class T1, class T2>
    struct PromoteTraits<Rational<T1>, Rational<T2> >
    {
        typedef Rational<typename PromoteTraits<T1, T2>::Promote> Promote;
    };
    \endcode

    <b>\#include</b> \<<a href="rational_8hxx-source.html">vigra/rational.hxx</a>\><br>
    Namespace: vigra

*/
#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template<class T>
struct NumericTraits<Rational<T> >
{
    typedef Rational<T> Type;
    typedef Rational<typename NumericTraits<T>::Promote> Promote;
    typedef Rational<typename NumericTraits<T>::RealPromote> RealPromote;
    typedef std::complex<Rational<RealPromote> > ComplexPromote;
    typedef T ValueType;

    typedef typename NumericTraits<T>::isIntegral isIntegral;
    typedef VigraTrueType isScalar;
    typedef typename NumericTraits<T>::isSigned isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;

    static Type zero() { return Type(0); }
    static Type one() { return Type(1); }
    static Type nonZero() { return one(); }

    static Promote toPromote(Type const & v)
        { return Promote(v.numerator(), v.denominator(), false); }
    static RealPromote toRealPromote(Type const & v)
        { return RealPromote(v.numerator(), v.denominator(), false); }
    static Type fromPromote(Promote const & v)
        { return Type(NumericTraits<T>::fromPromote(v.numerator()),
                      NumericTraits<T>::fromPromote(v.denominator()), false); }
    static Type fromRealPromote(RealPromote const & v)
        { return Type(NumericTraits<T>::fromRealPromote(v.numerator()),
                      NumericTraits<T>::fromRealPromote(v.denominator()), false); }
};

template<class T>
struct NormTraits<Rational<T> >
{
    typedef Rational<T>                           Type;
    typedef typename NumericTraits<Type>::Promote SquaredNormType;
    typedef Type                                  NormType;
};

template <class T>
struct PromoteTraits<Rational<T>, Rational<T> >
{
    typedef Rational<typename PromoteTraits<T, T>::Promote> Promote;
    static Promote toPromote(Rational<T> const & v) { return v; }
};

template <class T1, class T2>
struct PromoteTraits<Rational<T1>, Rational<T2> >
{
    typedef Rational<typename PromoteTraits<T1, T2>::Promote> Promote;
    static Promote toPromote(Rational<T1> const & v) { return v; }
    static Promote toPromote(Rational<T2> const & v) { return v; }
};

template <class T1, class T2>
struct PromoteTraits<Rational<T1>, T2 >
{
    typedef Rational<typename PromoteTraits<T1, T2>::Promote> Promote;
    static Promote toPromote(Rational<T1> const & v) { return v; }
    static Promote toPromote(T2 const & v) { return Promote(v); }
};

template <class T1, class T2>
struct PromoteTraits<T1, Rational<T2> >
{
    typedef Rational<typename PromoteTraits<T1, T2>::Promote> Promote;
    static Promote toPromote(T1 const & v) { return Promote(v); }
    static Promote toPromote(Rational<T2> const & v) { return v; }
};

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

/********************************************************/
/*                                                      */
/*                   RationalOperations                 */
/*                                                      */
/********************************************************/

/** \addtogroup RationalOperations Functions for Rational

    \brief     <b>\#include</b> \<<a href="rational_8hxx-source.html">vigra/rational.hxx</a>\><br>

    These functions fulfill the requirements of an \ref AlgebraicField.

    Namespace: vigra
    <p>

 */
//@{

/********************************************************/
/*                                                      */
/*                        arithmetic                    */
/*                                                      */
/********************************************************/

    /// unary plus
template <typename IntType>
inline Rational<IntType> operator+ (const Rational<IntType>& r)
{
    return r;
}

    /// unary minus (negation)
template <typename IntType>
inline Rational<IntType> operator- (const Rational<IntType>& r)
{
    return Rational<IntType>(-r.numerator(), r.denominator(), false);
}

    /// addition
template <typename IntType>
inline Rational<IntType> operator+(Rational<IntType> l, Rational<IntType> const & r)
{
    return l += r;
}

    /// subtraction
template <typename IntType>
inline Rational<IntType> operator-(Rational<IntType> l, Rational<IntType> const & r)
{
    return l -= r;
}

    /// multiplication
template <typename IntType>
inline Rational<IntType> operator*(Rational<IntType> l, Rational<IntType> const & r)
{
    return l *= r;
}

    /// division
template <typename IntType>
inline Rational<IntType> operator/(Rational<IntType> l, Rational<IntType> const & r)
{
    return l /= r;
}

    /// addition of right-hand <tt>IntType</tt> argument
template <typename IntType>
inline Rational<IntType>
operator+(Rational<IntType> l, typename Rational<IntType>::param_type r)
{
    return l += r;
}

    /// subtraction of right-hand <tt>IntType</tt> argument
template <typename IntType>
inline Rational<IntType>
operator-(Rational<IntType> l, typename Rational<IntType>::param_type r)
{
    return l -= r;
}

    /// multiplication with right-hand <tt>IntType</tt> argument
template <typename IntType>
inline Rational<IntType>
operator*(Rational<IntType> l, typename Rational<IntType>::param_type r)
{
    return l *= r;
}

    /// division by right-hand <tt>IntType</tt> argument
template <typename IntType>
inline Rational<IntType>
operator/(Rational<IntType> l, typename Rational<IntType>::param_type r)
{
    return l /= r;
}

    /// addition of left-hand <tt>IntType</tt> argument
template <typename IntType>
inline Rational<IntType>
operator+(typename Rational<IntType>::param_type l, Rational<IntType> r)
{
    return r += l;
}

    /// subtraction from left-hand <tt>IntType</tt> argument
template <typename IntType>
inline Rational<IntType>
operator-(typename Rational<IntType>::param_type l, Rational<IntType> const & r)
{
    return (-r) += l;
}

    /// multiplication with left-hand <tt>IntType</tt> argument
template <typename IntType>
inline Rational<IntType>
operator*(typename Rational<IntType>::param_type l, Rational<IntType> r)
{
    return r *= l;
}

    /// division of left-hand <tt>IntType</tt> argument
template <typename IntType>
inline Rational<IntType>
operator/(typename Rational<IntType>::param_type l, Rational<IntType> const & r)
{
    if(r.numerator() < IntType(0))
        return Rational<IntType>(-r.denominator(), -r.numerator(), false) *= l;
    else
        return Rational<IntType>(r.denominator(), r.numerator(), false) *= l;
}

/********************************************************/
/*                                                      */
/*                        comparison                    */
/*                                                      */
/********************************************************/


    /// equality
template <typename IntType1, typename IntType2>
inline bool
operator== (const Rational<IntType1> & l, const Rational<IntType2>& r)
{
    return l.denominator() == r.denominator() &&
           l.numerator() == r.numerator(); // works since numbers are normalized, even
                                           // if they represent +-infinity
}

    /// equality with right-hand <tt>IntType2</tt> argument
template <typename IntType1, typename IntType2>
inline bool
operator== (const Rational<IntType1> & l, IntType2 const & i)
{
    return ((l.denominator() == IntType1(1)) && (l.numerator() == i));
}

    /// equality with left-hand <tt>IntType1</tt> argument
template <typename IntType1, typename IntType2>
inline bool
operator==(IntType1 const & l, Rational<IntType2> const & r)
{
    return r == l;
}

    /// inequality
template <typename IntType1, typename IntType2>
inline bool
operator!=(Rational<IntType1> const & l, Rational<IntType2> const & r)
{
    return l.denominator() != r.denominator() ||
           l.numerator() != r.numerator(); // works since numbers are normalized, even
                                           // if they represent +-infinity
}

    /// inequality with right-hand <tt>IntType2</tt> argument
template <typename IntType1, typename IntType2>
inline bool
operator!= (const Rational<IntType1> & l, IntType2 const & i)
{
    return ((l.denominator() != IntType1(1)) || (l.numerator() != i));
}

    /// inequality with left-hand <tt>IntType1</tt> argument
template <typename IntType1, typename IntType2>
inline bool
operator!=(IntType1 const & l, Rational<IntType2> const & r)
{
    return r != l;
}

    /// less-than
template <typename IntType1, typename IntType2>
bool
operator< (const Rational<IntType1> & l, const Rational<IntType2>& r)
{
    // Avoid repeated construction
    typedef typename PromoteTraits<IntType1, IntType2>::Promote IntType;
    IntType zero(0);

    // Handle the easy cases. Take advantage of the fact
    // that the denominator is never negative.
    if(l.denominator() == zero)
    {
        if(r.denominator() == zero)
            // -inf < inf, !(-inf < -inf), !(inf < -inf), !(inf < inf)
            return l.numerator() < r.numerator();
        else
            // -inf < -1, -inf < 0, -inf < 1
            // !(inf < -1), !(inf < 0), !(inf < 1)
            return l.numerator() < zero;
    }
    if(r.denominator() == zero)
        // -1 < inf, 0 < inf, 1 < inf
        // !(-1 < -inf), !(0 < -inf), !(1 < -inf)
        return r.numerator() > zero;
    // !(1 < -1), !(1 < 0), !(0 < -1), !(0 < 0)
    if(l.numerator() >= zero && r.numerator() <= zero)
        return false;
    // -1 < 0, -1 < 1, 0 < 1 (note: !(0 < 0) was already handled!)
    if(l.numerator() <= zero && r.numerator() >= zero)
        return true;

    // both numbers have the same sign (and are neither zero or +-infinity)
    // => calculate result, avoid overflow
    IntType gcd1 = gcd<IntType>(l.numerator(), r.numerator());
    IntType gcd2 = gcd<IntType>(r.denominator(), l.denominator());
    return (l.numerator()/gcd1) * (r.denominator()/gcd2) <
           (l.denominator()/gcd2) * (r.numerator()/gcd1);
}

    /// less-than with right-hand <tt>IntType2</tt> argument
template <typename IntType1, typename IntType2>
bool
operator< (const Rational<IntType1> & l, IntType2 const & i)
{
    // Avoid repeated construction
    typedef typename PromoteTraits<IntType1, IntType2>::Promote IntType;
    IntType zero(0);

    // Handle the easy cases. Take advantage of the fact
    // that the denominator is never negative.
    if(l.denominator() == zero)
        // -inf < -1, -inf < 0, -inf < 1
        // !(inf < -1), !(inf < 0), !(inf < 1)
        return l.numerator() < zero;
    // !(1 < -1), !(1 < 0), !(0 < -1), !(0 < 0)
    if(l.numerator() >= zero && i <= zero)
        return false;
    // -1 < 0, -1 < 1, 0 < 1 (note: !(0 < 0) was already handled!)
    if(l.numerator() <= zero && i >= zero)
        return true;

    // Now, use the fact that n/d truncates towards zero as long as n and d
    // are both positive.
    // Divide instead of multiplying to avoid overflow issues. Of course,
    // division may be slower, but accuracy is more important than speed...
    if (l.numerator() > zero)
        return (l.numerator()/l.denominator()) < i;
    else
        return -i < (-l.numerator()/l.denominator());
}

    /// less-than with left-hand <tt>IntType1</tt> argument
template <typename IntType1, typename IntType2>
inline bool
operator<(IntType1 const & l, Rational<IntType2> const & r)
{
    return r > l;
}

    /// greater-than
template <typename IntType1, typename IntType2>
inline bool
operator>(Rational<IntType1> const & l, Rational<IntType2> const & r)
{
    return r < l;
}

    /// greater-than with right-hand <tt>IntType2</tt> argument
template <typename IntType1, typename IntType2>
bool
operator> (const Rational<IntType1> & l, IntType2 const & i)
{
    // Trap equality first
    if (l.numerator() == i && l.denominator() == IntType1(1))
        return false;

    // Otherwise, we can use operator<
    return !(l < i);
}

    /// greater-than with left-hand <tt>IntType1</tt> argument
template <typename IntType1, typename IntType2>
inline bool
operator>(IntType1 const & l, Rational<IntType2> const & r)
{
    return r < l;
}

    /// less-equal
template <typename IntType1, typename IntType2>
inline bool
operator<=(Rational<IntType1> const & l, Rational<IntType2> const & r)
{
    return !(r < l);
}

    /// less-equal with right-hand <tt>IntType2</tt> argument
template <typename IntType1, typename IntType2>
inline bool
operator<=(Rational<IntType1> const & l, IntType2 const & r)
{
    return !(l > r);
}

    /// less-equal with left-hand <tt>IntType1</tt> argument
template <typename IntType1, typename IntType2>
inline bool
operator<=(IntType1 const & l, Rational<IntType2> const & r)
{
    return r >= l;
}

    /// greater-equal
template <typename IntType1, typename IntType2>
inline bool
operator>=(Rational<IntType1> const & l, Rational<IntType2> const & r)
{
    return !(l < r);
}

    /// greater-equal with right-hand <tt>IntType2</tt> argument
template <typename IntType1, typename IntType2>
inline bool
operator>=(Rational<IntType1> const & l, IntType2 const & r)
{
    return !(l < r);
}

    /// greater-equal with left-hand <tt>IntType1</tt> argument
template <typename IntType1, typename IntType2>
inline bool
operator>=(IntType1 const & l, Rational<IntType2> const & r)
{
    return r <= l;
}

/********************************************************/
/*                                                      */
/*                 algebraic functions                  */
/*                                                      */
/********************************************************/

    /// absolute value
template <typename IntType>
inline Rational<IntType>
abs(const Rational<IntType>& r)
{
    if (r.numerator() >= IntType(0))
        return r;

    return Rational<IntType>(-r.numerator(), r.denominator(), false);
}

    /// norm (same as <tt>abs(r)</tt>)
template <typename IntType>
inline Rational<IntType>
norm(const Rational<IntType>& r)
{
    return abs(r);
}

    /// squared norm
template <typename IntType>
inline typename NormTraits<Rational<IntType> >::SquaredNormType
squaredNorm(const Rational<IntType>& r)
{
    return typename NormTraits<Rational<IntType> >::SquaredNormType(sq(r.numerator()), sq(r.denominator()), false);
}

    /** integer powers

        <tt>throws bad_rational</tt> if indeterminate expression.
    */
template <typename IntType>
Rational<IntType>
pow(const Rational<IntType>& r, int e)
{
    IntType zero(0);
    int ae;
    if(e == 0)
    {
        if(r.denominator() == zero)
                throw bad_rational();
        return Rational<IntType>(IntType(1));
    }
    else if(e < 0)
    {
        if(r.numerator() == zero)
            return Rational<IntType>(IntType(1), zero, false);
        if(r.denominator() == zero)
            return Rational<IntType>(zero);
        ae = -e;
    }
    else
    {
        if(r.denominator() == zero || r.numerator() == zero)
            return r;
        ae = e;
    }

    IntType nold = r.numerator(), dold = r.denominator(),
            nnew = IntType(1), dnew = IntType(1);
    for(; ae != 0; ae >>= 1, nold *= nold, dold *= dold)
    {
        if(ae % 2 != 0)
        {
            nnew *= nold;
            dnew *= dold;
        }
    }
    if(e < 0)
    {
        if(nnew < zero)
            return Rational<IntType>(-dnew, -nnew, false);
        else
            return Rational<IntType>(dnew, nnew, false);
    }
    else
        return Rational<IntType>(nnew, dnew, false);
}

    /// largest integer not larger than <tt>r</tt>
template <typename IntType>
Rational<IntType>
floor(const Rational<IntType>& r)
{
    IntType zero(0), one(1);
    if(r.denominator() == zero || r.denominator() == one)
        return r;
    return r.numerator() < zero ?
                   Rational<IntType>(r.numerator() / r.denominator() - one)
                 : Rational<IntType>(r.numerator() / r.denominator());
}

    /// smallest integer not smaller than <tt>r</tt>
template <typename IntType>
Rational<IntType>
ceil(const Rational<IntType>& r)
{
    IntType zero(0), one(1);
    if(r.denominator() == zero || r.denominator() == one)
        return r;
    return r.numerator() < IntType(0) ?
                   Rational<IntType>(r.numerator() / r.denominator())
                 : Rational<IntType>(r.numerator() / r.denominator() + one);
}


    /** Type conversion

        Executes <tt>static_cast<T>(numerator()) / denominator()</tt>.

        <b>Usage:</b>

        \code
        Rational<int> r;
        int i;
        double d;
        i = rational_cast<int>(r);     // round r downwards
        d = rational_cast<double>(r);  // represent rational as a double
        r = rational_cast<Rational<int> >(r);   // no change
        \endcode
    */
template <typename T, typename IntType>
inline T rational_cast(const Rational<IntType>& src)
{
    return static_cast<T>(src.numerator())/src.denominator();
}

template <class T>
inline T const & rational_cast(T const & v)
{
    return v;
}

//@}

template <typename IntType>
std::ostream& operator<< (std::ostream& os, const vigra::Rational<IntType>& r)
{
    os << r.numerator() << '/' << r.denominator();
    return os;
}

} // namespace vigra

#endif  // VIGRA_RATIONAL_HPP

