/************************************************************************/
/*                                                                      */
/*               Copyright 2004-2005 by Ullrich Koethe                  */
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

#ifndef VIGRA_FIXEDPOINT_HXX
#define VIGRA_FIXEDPOINT_HXX

#include "vigra/mathutil.hxx"
#include "vigra/static_assert.hxx"
#include "vigra/error.hxx"
#include "vigra/numerictraits.hxx"

namespace vigra {

template <unsigned IntBits, unsigned FractionalBits>
class FixedPoint;

struct Error_FixedPointTraits_not_specialized_for_this_case;

template <class T1, class T2>
class FixedPointTraits
{
public:
    typedef Error_FixedPointTraits_not_specialized_for_this_case PlusType;
    typedef Error_FixedPointTraits_not_specialized_for_this_case MinusType;
    typedef Error_FixedPointTraits_not_specialized_for_this_case MultipliesType;
//    typedef Error_FixedPointTraits_not_specialized_for_this_case DividesType;
};

// return type policy: 
//     * try to allocate enough bits to represent the biggest possible result
//     * in case of ass/subtract: if all bits of the internal int are used up, 
//                                keep the representation
template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
class FixedPointTraits<FixedPoint<IntBits1, FracBits1>, FixedPoint<IntBits2, FracBits2> >
{
    enum { MaxIntBits  = (IntBits1 < IntBits2) ? IntBits2 : IntBits1,
           MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1,
           PlusMinusIntBits = (MaxIntBits + 1 + MaxFracBits < 32) ?
                               MaxIntBits + 1 : MaxIntBits};
public:
    typedef FixedPoint<PlusMinusIntBits, MaxFracBits>               PlusType;
    typedef FixedPoint<PlusMinusIntBits, MaxFracBits>               MinusType;
    typedef FixedPoint<IntBits1 + IntBits2, FracBits1 + FracBits2>  MultipliesType;
//    typedef FixedPoint<IntBits1 + FracBits2, FracBits1 + IntBits2>  DividesType;
};

template <int N>
struct FixedPoint_overflow_error__More_than_31_bits_requested
: staticAssert::AssertBool<(N < 32)>
{};


template <bool Predicate>
struct FixedPoint_assignment_error__Target_object_has_too_few_integer_bits
: staticAssert::AssertBool<Predicate>
{};

enum FixedPointNoShift { FPNoShift };

namespace detail {

template <bool MustRound>
struct FPAssignWithRound;

template <>
struct FPAssignWithRound<false>
{
    template <int N>
    static inline int exec(int v) { return v << (-N); }
};

template <>
struct FPAssignWithRound<true>
{
    template <int N>
    static inline int exec(int const v)
    {
        return (v + (1 << (N - 1))) >> (N);
    }
};

template <bool MustRound>
struct FPMulImplementation;

template <>
struct FPMulImplementation<false>
{
    template <int N>
    static inline int exec(int l, int r) { return (l * r) << (-N); }
};

template <>
struct FPMulImplementation<true>
{
    template <int N>
    static inline int exec(int l, int r)
    {
        // there is not enough space in the result
        // => perform calculations that preserve as much accuracy as possible
        enum { diffl = N / 2, diffr = N  - diffl, maskl = (1 << diffl) - 1, maskr = (1 << diffr) - 1 };
        int shiftl = l >> diffl;
        int shiftr = r >> diffr;

        return shiftl * shiftr + (((l & maskl) * shiftr) >> diffl) +
                                 (((r & maskr) * shiftl) >> diffr);
    }
};

} // namespace detail

template <unsigned IntBits, unsigned FractionalBits>
class FixedPoint
{
public:
    enum {
        INT_BITS        = IntBits,
        FRACTIONAL_BITS = FractionalBits,
        TOTAL_BITS      = IntBits + FractionalBits,
        MAX             = (int)(((unsigned)1 << TOTAL_BITS) - 1),
        ONE             = 1 << FractionalBits,
        ONE_HALF        = ONE >> 1,
        FRACTIONAL_MASK = ONE - 1,
        INT_MASK        = MAX ^ FRACTIONAL_MASK
    };

    int value;

    FixedPoint()
    {}

        /**
        */
    explicit FixedPoint(int v)
    : value(v << FractionalBits)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_overflow_error__More_than_31_bits_requested<(IntBits + FractionalBits)>));
    }

        /**
        */
    FixedPoint(int v, FixedPointNoShift)
    : value(v)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_overflow_error__More_than_31_bits_requested<(IntBits + FractionalBits)>));
    }

    explicit FixedPoint(double rhs)
    : value((int)round(rhs * ONE))
    {
        vigra_precondition(abs(rhs * ONE) <= (double)MAX,
            "FixedPoint(double rhs): Too few integer bits to convert rhs.");
    }


    FixedPoint(const FixedPoint &other)
    : value(other.value)
    {}

    template <unsigned Int2, unsigned Frac2>
    FixedPoint(const FixedPoint<Int2, Frac2> &other)
    : value(detail::FPAssignWithRound<(Frac2 > FractionalBits)>::template exec<Frac2 - FractionalBits>(other.value))
    {
        VIGRA_STATIC_ASSERT((FixedPoint_assignment_error__Target_object_has_too_few_integer_bits<(IntBits >= Int2)>));
    }

    FixedPoint &operator=(int rhs)
    {
        vigra_precondition(abs(rhs) < (1 << IntBits),
            "FixedPoint::operator=(int rhs): Too few integer bits to represent rhs.");
        value = rhs << FractionalBits;
        return *this;
    }

    FixedPoint &operator=(double rhs)
    {
        vigra_precondition(abs(rhs) <= ((1 << IntBits) - 1),
            "FixedPoint::operator=(double rhs): Too few integer bits to convert rhs.");
        value = (int)round(rhs * ONE);
        return *this;
    }

    FixedPoint & operator=(const FixedPoint &other)
    {
        value = other.value;
        return *this;
    }

    template <unsigned Int2, unsigned Frac2>
    FixedPoint & operator=(const FixedPoint<Int2, Frac2> &other)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_assignment_error__Target_object_has_too_few_integer_bits<(IntBits >= Int2)>));
        value = detail::FPAssignWithRound<(Frac2 > FractionalBits)>::template exec<Frac2 - FractionalBits>(other.value);
        return *this;
    }

    operator double() const
    {
        return (double)value / ONE;
    }

    FixedPoint operator-() const
    {
        return FixedPoint(-value, FPNoShift);
    }

    FixedPoint & operator++()
    {
        value += ONE;
        return *this;
    }

    FixedPoint operator++(int)
    {
        FixedPoint old(*this);
        value += ONE;
        return old;
    }

    FixedPoint & operator--()
    {
        value -= ONE;
        return *this;
    }

    FixedPoint operator--(int)
    {
        FixedPoint old(*this);
        value -= ONE;
        return old;
    }

    template <unsigned Int2, unsigned Frac2>
    FixedPoint & operator+=(const FixedPoint<Int2, Frac2> &other)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_assignment_error__Target_object_has_too_few_integer_bits<(IntBits >= Int2)>));
        value += detail::FPAssignWithRound<(Frac2 > FractionalBits)>::template exec<Frac2 - FractionalBits>(other.value);
        return *this;
    }

    template <unsigned Int2, unsigned Frac2>
    FixedPoint & operator-=(const FixedPoint<Int2, Frac2> &other)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_assignment_error__Target_object_has_too_few_integer_bits<(IntBits >= Int2)>));
        value -= detail::FPAssignWithRound<(Frac2 > FractionalBits)>::template exec<Frac2 - FractionalBits>(other.value);
        return *this;
    }
    
    template <unsigned Int2, unsigned Frac2>
    FixedPoint & operator*=(const FixedPoint<Int2, Frac2> &other)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_assignment_error__Target_object_has_too_few_integer_bits<(IntBits >= Int2)>));
        value = detail::FPMulImplementation<(Frac2 > 0)>::template exec<Frac2>(value, other.value);
        return *this;
    }
};

#define VIGRA_FIXED_POINT_FACTORY(T, INTBITS) \
    inline FixedPoint<INTBITS, 0> fixedPoint(T t) \
    { \
        return FixedPoint<INTBITS, 0>(t, FPNoShift); \
    }

VIGRA_FIXED_POINT_FACTORY(unsigned char, 8)
VIGRA_FIXED_POINT_FACTORY(signed char, 7)
VIGRA_FIXED_POINT_FACTORY(unsigned short, 16)
VIGRA_FIXED_POINT_FACTORY(signed short, 15)
VIGRA_FIXED_POINT_FACTORY(int, 31)

#undef VIGRA_FIXED_POINT_FACTORY

template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
bool operator==(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    enum { MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1 };
    return (l.value << (MaxFracBits - FracBits1)) == (r.value << (MaxFracBits - FracBits2));
}

template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
bool operator!=(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    enum { MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1 };
    return (l.value << (MaxFracBits - FracBits1)) != (r.value << (MaxFracBits - FracBits2));
}

template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
bool operator<(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    enum { MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1 };
    return (l.value << (MaxFracBits - FracBits1)) < (r.value << (MaxFracBits - FracBits2));
}

template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
bool operator<=(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    enum { MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1 };
    return (l.value << (MaxFracBits - FracBits1)) <= (r.value << (MaxFracBits - FracBits2));
}

template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
bool operator>(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    enum { MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1 };
    return (l.value << (MaxFracBits - FracBits1)) > (r.value << (MaxFracBits - FracBits2));
}

template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
bool operator>=(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    enum { MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1 };
    return (l.value << (MaxFracBits - FracBits1)) >= (r.value << (MaxFracBits - FracBits2));
}

template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
typename FixedPointTraits<FixedPoint<IntBits1, FracBits1>, FixedPoint<IntBits2, FracBits2> >::PlusType
operator+(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    enum { MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1 };
    return typename
        FixedPointTraits<FixedPoint<IntBits1, FracBits1>, FixedPoint<IntBits2, FracBits2> >::
        PlusType((l.value << (MaxFracBits - FracBits1)) + (r.value << (MaxFracBits - FracBits2)), FPNoShift);
}

template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2,
          unsigned IntBits3, unsigned FracBits3>
inline void
add(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r,
    FixedPoint<IntBits3, FracBits3> & result)
{
    result = l + r;
}

template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
typename FixedPointTraits<FixedPoint<IntBits1, FracBits1>, FixedPoint<IntBits2, FracBits2> >::MinusType
operator-(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    enum { MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1 };
    return typename
        FixedPointTraits<FixedPoint<IntBits1, FracBits1>, FixedPoint<IntBits2, FracBits2> >::
        MinusType((l.value << (MaxFracBits - FracBits1)) - (r.value << (MaxFracBits - FracBits2)), FPNoShift);
}

template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2,
          unsigned IntBits3, unsigned FracBits3>
inline void
sub(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r,
    FixedPoint<IntBits3, FracBits3> & result)
{
    result = l - r;
}

template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
typename FixedPointTraits<FixedPoint<IntBits1, FracBits1>, FixedPoint<IntBits2, FracBits2> >::MultipliesType
operator*(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    return typename
        FixedPointTraits<FixedPoint<IntBits1, FracBits1>, FixedPoint<IntBits2, FracBits2> >::
        MultipliesType(l.value * r.value, FPNoShift);
}

template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2,
          unsigned IntBits3, unsigned FracBits3>
inline void
mul(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r,
    FixedPoint<IntBits3, FracBits3> & result)
{
    VIGRA_STATIC_ASSERT((FixedPoint_assignment_error__Target_object_has_too_few_integer_bits<(IntBits1 + IntBits2 <= IntBits3)>));
    enum { diff = FracBits1 + FracBits2 - FracBits3 };
    result.value = detail::FPMulImplementation<(diff > 0)>::template exec<diff>(l.value, r.value);
}

template <unsigned IntBits, unsigned FracBits>
inline FixedPoint<IntBits, FracBits>
abs(FixedPoint<IntBits, FracBits> v)
{
    return FixedPoint<IntBits, FracBits>(abs(v.value), FPNoShift);
}

template <unsigned IntBits, unsigned FracBits>
inline FixedPoint<0, FracBits>
frac(FixedPoint<IntBits, FracBits> v)
{
    return FixedPoint<0, FracBits>(v.value & FixedPoint<IntBits, FracBits>::FRACTIONAL_MASK, FPNoShift);
}

template <unsigned IntBits, unsigned FracBits>
inline FixedPoint<0, FracBits>
dual_frac(FixedPoint<IntBits, FracBits> v)
{
    return FixedPoint<0, FracBits>(FixedPoint<0, FracBits>::ONE - 
                                   (v.value & FixedPoint<IntBits, FracBits>::FRACTIONAL_MASK), FPNoShift);
}

template <unsigned IntBits, unsigned FracBits>
inline int
floor(FixedPoint<IntBits, FracBits> v)
{
    return(v.value >> FracBits);
}

template <unsigned IntBits, unsigned FracBits>
inline int
ceil(FixedPoint<IntBits, FracBits> v)
{
    return((v.value + FixedPoint<IntBits, FracBits>::FRACTIONAL_MASK) >> FracBits);
}

template <unsigned IntBits, unsigned FracBits>
inline int
round(FixedPoint<IntBits, FracBits> v)
{
    return((v.value + FixedPoint<IntBits, FracBits>::ONE_HALF) >> FracBits);
}

template <unsigned IntBits, unsigned FracBits>
struct NumericTraits<FixedPoint<IntBits, FracBits> >
{
    typedef FixedPoint<IntBits, FracBits> Type;
        //typedef FixedPoint<IntBits, FracBits> Promote;
        //typedef FixedPoint<IntBits, FracBits> RealPromote;
        //typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

        //typedef VigraFalseType isIntegral;
    typedef VigraTrueType  isScalar;
    typedef VigraTrueType  isOrdered;
    typedef VigraFalseType isComplex;

    static Type zero() { return Type(0, FPNoShift); }
    static Type one() { return Type(Type::ONE, FPNoShift); }
    static Type nonZero() { return one(); }
    static Type epsilon() { return Type(1, FPNoShift); }
    static Type smallestPositive() { return Type(1, FPNoShift); }
    static Type max() { return Type( Type::MAX, FPNoShift); }
    static Type min() { return -max(); }
};

template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
struct PromoteTraits<FixedPoint<IntBits1, FracBits1>,
                     FixedPoint<IntBits2, FracBits2> >
{
    typedef typename 
        FixedPointTraits<FixedPoint<IntBits1, FracBits1>, FixedPoint<IntBits2, FracBits2> >::PlusType 
        Promote;
};

} // namespace vigra

#endif // VIGRA_FIXEDPOINT_HXX
