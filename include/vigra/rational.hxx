/************************************************************************/
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    It was adapted from the file boost/rational.hpp of the            */
/*    boost library.                                                    */
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
#include "vigra/config.hxx"
#include "vigra/mathutil.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/metaprogramming.hxx"

namespace vigra {

// Note: We use n and m as temporaries in this function, so there is no value
// in using const IntType& as we would only need to make a copy anyway...
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

template <typename IntType>
class Rational
{
public:
    typedef IntType value_type;
    typedef typename If<typename TypeTraits<IntType>::isBuiltinType,
                        IntType, IntType const &>::type param_type;

    Rational() : num(0), den(1) {}
    
    template <class U>
    Rational(Rational<U> const & r)
    : num(r.numerator()),
      den(r.denominator())
    {}
    
    Rational(param_type n) 
    : num(n), 
      den(IntType(1)) 
    {}
    
    Rational(param_type n, param_type d, bool doNormalize = true) 
    : num(n), 
      den(d) 
    { 
        if(doNormalize)
            normalize();
    }
    
    explicit Rational(double v, double epsilon = 1e-4) 
    : num(IntType(v < 0.0 ?
                      IntType(1.0/epsilon)*v - 0.5
                    : IntType(1.0/epsilon)*v + 0.5)), 
      den(IntType(1.0/epsilon)) 
    {
        normalize();
    }
    
    // Default copy constructor and assignment are fine

    // Add assignment from IntType
    Rational& operator=(param_type n) { return assign(n, 1); }

    // Assign in place
    Rational& assign(param_type n, param_type d);

    // Access to representation
    IntType numerator() const { return num; }
    IntType denominator() const { return den; }

    // Arithmetic assignment operators
    Rational& operator+= (const Rational& r);
    Rational& operator-= (const Rational& r);
    Rational& operator*= (const Rational& r);
    Rational& operator/= (const Rational& r);

    Rational& operator+= (param_type i);
    Rational& operator-= (param_type i);
    Rational& operator*= (param_type i);
    Rational& operator/= (param_type i);
    
    // Increment and decrement
    const Rational& operator++();
    const Rational& operator--();
    
    Rational operator++(int) { Rational res(*this); operator++(); return res; }
    Rational operator--(int) { Rational res(*this); operator--(); return res; }

    // Operator not
    bool operator!() const { return !num; }

    // Comparison operators
    bool operator== (const Rational& r) const;    
    bool operator< (const Rational& r) const;

    bool operator==(param_type i) const;
    bool operator< (param_type i) const;
    bool operator> (param_type i) const;

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

template <typename IntType>
inline Rational<IntType> operator+(Rational<IntType> l, Rational<IntType> const & r)
{
    return l += r;
}

template <typename IntType>
inline Rational<IntType> operator-(Rational<IntType> l, Rational<IntType> const & r)
{
    return l -= r;
}

template <typename IntType>
inline Rational<IntType> operator*(Rational<IntType> l, Rational<IntType> const & r)
{
    return l *= r;
}

template <typename IntType>
inline Rational<IntType> operator/(Rational<IntType> l, Rational<IntType> const & r)
{
    return l /= r;
}

template <typename IntType>
inline Rational<IntType> 
operator+(Rational<IntType> l, typename Rational<IntType>::param_type r)
{
    return l += r;
}

template <typename IntType>
inline Rational<IntType> 
operator-(Rational<IntType> l, typename Rational<IntType>::param_type r)
{
    return l -= r;
}

template <typename IntType>
inline Rational<IntType> 
operator*(Rational<IntType> l, typename Rational<IntType>::param_type r)
{
    return l *= r;
}

template <typename IntType>
inline Rational<IntType> 
operator/(Rational<IntType> l, typename Rational<IntType>::param_type r)
{
    return l /= r;
}

template <typename IntType>
inline Rational<IntType> 
operator+(typename Rational<IntType>::param_type l, Rational<IntType> r)
{
    return r += l;
}

template <typename IntType>
inline Rational<IntType> 
operator-(typename Rational<IntType>::param_type l, Rational<IntType> const & r)
{
    return (-r) += l;
}

template <typename IntType>
inline Rational<IntType> 
operator*(typename Rational<IntType>::param_type l, Rational<IntType> r)
{
    return r *= l;
}

template <typename IntType>
inline Rational<IntType> 
operator/(typename Rational<IntType>::param_type l, Rational<IntType> const & r)
{
    IntType zero(0);
    if(r.numerator() == zero)
        throw bad_rational();
    if(r.numerator() < zero)
        return Rational<IntType>(-r.denominator(), -r.numerator(), false) *= l;
    else
        return Rational<IntType>(r.denominator(), r.numerator(), false) *= l;
}

template <typename IntType1, typename IntType2>
inline bool operator!=(Rational<IntType1> const & l, Rational<IntType2> const & r)
{
    return !(r == l);
}

template <typename IntType1, typename IntType2>
inline bool operator>(Rational<IntType1> const & l, Rational<IntType2> const & r)
{
    return r < l;
}

template <typename IntType1, typename IntType2>
inline bool operator<=(Rational<IntType1> const & l, Rational<IntType2> const & r)
{
    return !(r < l);
}

template <typename IntType1, typename IntType2>
inline bool operator>=(Rational<IntType1> const & l, Rational<IntType2> const & r)
{
    return !(l < r);
}

template <typename IntType1, typename IntType2>
inline bool operator!=(Rational<IntType1> const & l, IntType2 r)
{
    return !(l == r);
}

template <typename IntType1, typename IntType2>
inline bool operator<=(Rational<IntType1> const & l, IntType2 r)
{
    return !(l > r);
}

template <typename IntType1, typename IntType2>
inline bool operator>=(Rational<IntType1> const & l, IntType2 r)
{
    return !(l < r);
}

template <typename IntType1, typename IntType2>
inline bool operator==(IntType1 l, Rational<IntType2> const & r)
{
    return r == l;
}

template <typename IntType1, typename IntType2>
inline bool operator!=(IntType1 l, Rational<IntType2> const & r)
{
    return r != l;
}

template <typename IntType1, typename IntType2>
inline bool operator<(IntType1 l, Rational<IntType2> const & r)
{
    return r > l;
}

template <typename IntType1, typename IntType2>
inline bool operator<=(IntType1 l, Rational<IntType2> const & r)
{
    return r >= l;
}

template <typename IntType1, typename IntType2>
inline bool operator>(IntType1 l, Rational<IntType2> const & r)
{
    return r < l;
}

template <typename IntType1, typename IntType2>
inline bool operator>=(IntType1 l, Rational<IntType2> const & r)
{
    return r <= l;
}


// Assign in place
template <typename IntType>
inline Rational<IntType>& Rational<IntType>::assign(param_type n, param_type d)
{
    num = n;
    den = d;
    normalize();
    return *this;
}

// Unary plus and minus
template <typename IntType>
inline Rational<IntType> operator+ (const Rational<IntType>& r)
{
    return r;
}

template <typename IntType>
inline Rational<IntType> operator- (const Rational<IntType>& r)
{
    return Rational<IntType>(-r.numerator(), r.denominator(), false);
}

// Arithmetic assignment operators
template <typename IntType>
Rational<IntType>& Rational<IntType>::operator+= (const Rational<IntType>& r)
{
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
    // Protect against self-modification
    IntType r_num = r.num;
    IntType r_den = r.den;

    // Avoid repeated construction
    IntType zero(0);

    // Trap division by zero
    if (r_num == zero)
        throw bad_rational();
    if (num == zero)
        return *this;

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
    IntType g = gcd(i, num);
    if(i < IntType(0))
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
inline const Rational<IntType>& Rational<IntType>::operator++()
{
    // This can never denormalise the fraction
    num += den;
    return *this;
}

template <typename IntType>
inline const Rational<IntType>& Rational<IntType>::operator--()
{
    // This can never denormalise the fraction
    num -= den;
    return *this;
}

// Comparison operators
template <typename IntType>
bool Rational<IntType>::operator< (const Rational<IntType>& r) const
{
    // Avoid repeated construction
    IntType zero(0);

    // If the two values have different signs, we don't need to do the
    // expensive calculations below. We take advantage here of the fact
    // that the denominator is always positive.
    if (num < zero && r.num >= zero) // -ve < +ve
        return true;
    if (num >= zero && r.num <= zero) // +ve or zero is not < -ve or zero
        return false;

    // Avoid overflow
    IntType gcd1 = gcd<IntType>(num, r.num);
    IntType gcd2 = gcd<IntType>(r.den, den);
    return (num/gcd1) * (r.den/gcd2) < (den/gcd2) * (r.num/gcd1);
}

template <typename IntType>
bool Rational<IntType>::operator< (param_type i) const
{
    // Avoid repeated construction
    IntType zero(0);

    // If the two values have different signs, we don't need to do the
    // expensive calculations below. We take advantage here of the fact
    // that the denominator is always positive.
    if (num < zero && i >= zero) // -ve < +ve
        return true;
    if (num >= zero && i <= zero) // +ve or zero is not < -ve or zero
        return false;

    // Now, use the fact that n/d truncates towards zero as long as n and d
    // are both positive.
    // Divide instead of multiplying to avoid overflow issues. Of course,
    // division may be slower, but accuracy is more important than speed...
    if (num > zero)
        return (num/den) < i;
    else
        return -i < (-num/den);
}

template <typename IntType>
bool Rational<IntType>::operator> (param_type i) const
{
    // Trap equality first
    if (num == i && den == IntType(1))
        return false;

    // Otherwise, we can use operator<
    return !operator<(i);
}

template <typename IntType>
inline bool Rational<IntType>::operator== (const Rational<IntType>& r) const
{
    return ((num == r.num) && (den == r.den));
}

template <typename IntType>
inline bool Rational<IntType>::operator== (param_type i) const
{
    return ((den == IntType(1)) && (num == i));
}

// Normalisation
template <typename IntType>
void Rational<IntType>::normalize()
{
    // Avoid repeated construction
    IntType zero(0);

    if (den == zero)
        throw bad_rational();

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

template <typename IntType>
std::ostream& operator<< (std::ostream& os, const Rational<IntType>& r)
{
    os << r.numerator() << '/' << r.denominator();
    return os;
}

// Type conversion
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

// Do not use any abs() defined on IntType - it isn't worth it, given the
// difficulties involved (Koenig lookup required, there may not *be* an abs()
// defined, etc etc).
template <typename IntType>
inline Rational<IntType> abs(const Rational<IntType>& r)
{
    if (r.numerator() >= IntType(0))
        return r;

    return Rational<IntType>(-r.numerator(), r.denominator(), false);
}

template <typename IntType>
Rational<IntType> 
pow(const Rational<IntType>& r, int e)
{
    IntType zero(0);
    int ae;
    if(e < 0)
    {
        if(r.numerator() == zero)
            throw bad_rational();
        ae = -e;
    }
    else
        ae = e;
   
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

template <typename IntType>
Rational<IntType> floor(const Rational<IntType>& r)
{
    IntType one(1);
    return r.denominator() == one ?
               r
             : r.numerator() < IntType(0) ?
                   Rational<IntType>(r.numerator() / r.denominator() - one)
                 : Rational<IntType>(r.numerator() / r.denominator());
}

template <typename IntType>
Rational<IntType> ceil(const Rational<IntType>& r)
{
    IntType one(1);
    return r.denominator() == one ?
               r
             : r.numerator() < IntType(0) ?
                   Rational<IntType>(r.numerator() / r.denominator())
                 : Rational<IntType>(r.numerator() / r.denominator() + one);
}


} // namespace vigra

#endif  // VIGRA_RATIONAL_HPP

