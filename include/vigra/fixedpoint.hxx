/************************************************************************/
/*                                                                      */
/*               Copyright 2004-2005 by Ullrich Koethe                  */
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

#ifndef VIGRA_FIXEDPOINT_HXX
#define VIGRA_FIXEDPOINT_HXX

#include "mathutil.hxx"
#include "static_assert.hxx"
#include "error.hxx"
#include "numerictraits.hxx"

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
//     * in case of add/subtract: if all bits of the internal int are used up, 
//                                keep the representation
template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
class FixedPointTraits<FixedPoint<IntBits1, FracBits1>, FixedPoint<IntBits2, FracBits2> >
{
    enum { MaxIntBits  = (IntBits1 < IntBits2) ? IntBits2 : IntBits1,
           MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1,
           PlusMinusIntBits = (MaxIntBits + 1 + MaxFracBits < 32) ?
                               MaxIntBits + 1 : MaxIntBits,
           MultipliesFracBits = (IntBits1 + IntBits2 < 31) 
                                    ? (FracBits1 + FracBits2) > (31 - IntBits1 - IntBits2) 
                                            ? 31 - IntBits1 - IntBits2
                                            : FracBits1 + FracBits2
                                    : 0
         };
public:
    typedef FixedPoint<PlusMinusIntBits, MaxFracBits>               PlusType;
    typedef FixedPoint<PlusMinusIntBits, MaxFracBits>               MinusType;
    typedef FixedPoint<IntBits1 + IntBits2, MultipliesFracBits>  MultipliesType;
//    typedef FixedPoint<IntBits1 + FracBits2, FracBits1 + IntBits2>  DividesType;
};

template <unsigned IntBits, unsigned FracBits>
struct SquareRootTraits<FixedPoint<IntBits, FracBits> >
{
    enum { SRTotalBits = (IntBits + FracBits + 1) / 2,
           SRIntBits   = (IntBits + 1) / 2,
           SRFracBits  = SRTotalBits - SRIntBits
         };
public:
    typedef FixedPoint<IntBits, FracBits>      Type;
    typedef FixedPoint<SRIntBits, SRFracBits>  SquareRootResult;
    typedef Type                               SquareRootArgument;
};


#ifndef DOXYGEN

template <int N>
struct FixedPoint_overflow_error__More_than_31_bits_requested
: staticAssert::AssertBool<(N < 32)>
{};

#endif /* DOXYGEN */



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

/********************************************************/
/*                                                      */
/*                      FixedPoint                      */
/*                                                      */
/********************************************************/

/** Template for fixed point arithmetic.

    Fixed point arithmetic is used when computations with fractional accuracy
    must be made at the highest speed possible (e.g. in the inner loop
    of a volume rendering routine). The speed-up relative to floating
    point arithmetic can be dramatic, especially when one can avoid
    conversions between integer and floating point numbers (these are 
    very expensive because integer and floating point arithmetic
    resides in different pipelines). 
    
    The template wraps an <tt>int</tt> and uses <tt>IntBits</tt> to
    represent the integral part of a number, and <tt>FractionalBits</tt>
    for the fractional part, where <tt>IntBits + FractionalBits < 32</tt>.
    (The 32rd bit is reserved because FixedPoint is a signed type).
    These numbers will be automatically allocated in an intelligent way
    in the result of an arithmetic operation. For example, when two 
    fixed point numbers are multiplied, the required number of integer
    bits in the result is the sum of the number of integer bits of the
    arguments, but only when so many bits are available. This is figured out
    by means of FixedPointTraits, and a compile-time error is raised
    when no suitable representation can be found. The idea is that the right
    thing happens automatically as often as possible.

    <tt>FixedPoint</tt> implements the required interface of an
    \ref AlgebraicRing and the required numeric and
    promotion traits. In addition, it supports functions <tt>add</tt>, 
    <tt>sub</tt>, and <tt>mul</tt>, where a particular layout of the result can
    be enforced. 
    
    <tt>unsigned char, signed char, unsigned short, signed short, int</tt> can be
    transformed into a FixedPoint with appropriate layout by means of the factory
    function <tt>fixedPoint()</tt>.

    <b>See also:</b>
    <ul>
    <li> \ref FixedPointOperations
    <li> \ref FixedPointTraits
    </ul>

    <b>\#include</b> \<vigra/fixedpoint.hxx\><br>
    Namespace: vigra
*/
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

    Int32 value;

    FixedPoint()
    {
        VIGRA_STATIC_ASSERT((FixedPoint_overflow_error__More_than_31_bits_requested<(IntBits + FractionalBits)>));
    }

        /** Construct from an int (fractional part will become zero).
        */
    explicit FixedPoint(int v)
    : value(v << FractionalBits)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_overflow_error__More_than_31_bits_requested<(IntBits + FractionalBits)>));
    }

        /** Construct from an int by a bitwise copy. This is normally only used internally.
        */
    FixedPoint(int v, FixedPointNoShift)
    : value(v)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_overflow_error__More_than_31_bits_requested<(IntBits + FractionalBits)>));
    }

        /** Construct from an double and round the fractional part to 
            <tt>FractionalBits</tt> accuracy. A PreconditionViolation exception is raised when
            the integer part is too small to represent the number.
        */
    explicit FixedPoint(double rhs)
    : value((int)round(rhs * ONE))
    {
        VIGRA_STATIC_ASSERT((FixedPoint_overflow_error__More_than_31_bits_requested<(IntBits + FractionalBits)>));
        vigra_precondition(abs(rhs * ONE) <= (double)MAX,
            "FixedPoint(double rhs): Too few integer bits to convert rhs.");
    }

        /** Copy constructor.
        */
    FixedPoint(const FixedPoint &other)
    : value(other.value)
    {}

        /** Construct from a FixedPoint with different layout. It rounds as appropriate and raises
            a compile-time error when the target type has too few integer bits.
        */
    template <unsigned Int2, unsigned Frac2>
    FixedPoint(const FixedPoint<Int2, Frac2> &other)
    : value(detail::FPAssignWithRound<(Frac2 > FractionalBits)>::template exec<int(Frac2) - int(FractionalBits)>(other.value))
    {
        VIGRA_STATIC_ASSERT((FixedPoint_overflow_error__More_than_31_bits_requested<(IntBits + FractionalBits)>));
        VIGRA_STATIC_ASSERT((FixedPoint_assignment_error__Target_object_has_too_few_integer_bits<(IntBits >= Int2)>));
    }

        /** Assignment from int. The fractional part will become zero.  
            A PreconditionViolation exception is raised when
            the integer part is too small to represent the number.
        */
    FixedPoint &operator=(int rhs)
    {
        vigra_precondition(abs(rhs) < (1 << IntBits),
            "FixedPoint::operator=(int rhs): Too few integer bits to represent rhs.");
        value = rhs << FractionalBits;
        return *this;
    }

        /** Assignment form double. The fractional part is rounded, and a 
            PreconditionViolation exception is raised when
            the integer part is too small to represent the number.
        */
    FixedPoint &operator=(double rhs)
    {
        vigra_precondition(abs(rhs) <= ((1 << IntBits) - 1),
            "FixedPoint::operator=(double rhs): Too few integer bits to convert rhs.");
        value = (int)round(rhs * ONE);
        return *this;
    }

        /** Copy assignment.
        */
    FixedPoint & operator=(const FixedPoint &other)
    {
        value = other.value;
        return *this;
    }

        /** Assignment from a FixedPoint with different layout. It rounds as appropriate and raises
            a compile-time error when the target type has too few integer bits.
        */
    template <unsigned Int2, unsigned Frac2>
    FixedPoint & operator=(const FixedPoint<Int2, Frac2> &other)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_assignment_error__Target_object_has_too_few_integer_bits<(IntBits >= Int2)>));
        value = detail::FPAssignWithRound<(Frac2 > FractionalBits)>::template exec<int(Frac2) - int(FractionalBits)>(other.value);
        return *this;
    }

        /** Negation.
        */
    FixedPoint operator-() const
    {
        return FixedPoint(-value, FPNoShift);
    }

        /** Pre-increment.
        */
    FixedPoint & operator++()
    {
        value += ONE;
        return *this;
    }

        /** Post-increment.
        */
    FixedPoint operator++(int)
    {
        FixedPoint old(*this);
        value += ONE;
        return old;
    }

        /** Pre-decrement.
        */
    FixedPoint & operator--()
    {
        value -= ONE;
        return *this;
    }

        /** Post-decrement.
        */
    FixedPoint operator--(int)
    {
        FixedPoint old(*this);
        value -= ONE;
        return old;
    }

        /** Add-assignment from a FixedPoint with different layout. It rounds as appropriate and raises
            a compile-time error when the target type has too few integer bits.
        */
    template <unsigned Int2, unsigned Frac2>
    FixedPoint & operator+=(const FixedPoint<Int2, Frac2> &other)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_assignment_error__Target_object_has_too_few_integer_bits<(IntBits >= Int2)>));
        value += detail::FPAssignWithRound<(Frac2 > FractionalBits)>::template exec<Frac2 - FractionalBits>(other.value);
        return *this;
    }

        /** Subtract-assignment from a FixedPoint with different layout. It rounds as appropriate and raises
            a compile-time error when the target type has too few integer bits.
        */
    template <unsigned Int2, unsigned Frac2>
    FixedPoint & operator-=(const FixedPoint<Int2, Frac2> &other)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_assignment_error__Target_object_has_too_few_integer_bits<(IntBits >= Int2)>));
        value -= detail::FPAssignWithRound<(Frac2 > FractionalBits)>::template exec<Frac2 - FractionalBits>(other.value);
        return *this;
    }
    
        /** Multiply-assignment from a FixedPoint with different layout. It rounds as appropriate and raises
            a compile-time error when the target type has too few integer bits.
        */
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

template <class T>
struct FixedPointCast;

#define VIGRA_FIXED_POINT_CAST(type) \
template <> \
struct FixedPointCast<type> \
{ \
    template <unsigned IntBits, unsigned FracBits> \
    static type cast(FixedPoint<IntBits, FracBits> v) \
    { \
        return round(v); \
    } \
};

VIGRA_FIXED_POINT_CAST(Int8)
VIGRA_FIXED_POINT_CAST(UInt8)
VIGRA_FIXED_POINT_CAST(Int16)
VIGRA_FIXED_POINT_CAST(UInt16)
VIGRA_FIXED_POINT_CAST(Int32)
VIGRA_FIXED_POINT_CAST(UInt32)

#undef VIGRA_FIXED_POINT_CAST

template <>
struct FixedPointCast<float>
{
    template <unsigned IntBits, unsigned FracBits>
    static float cast(FixedPoint<IntBits, FracBits> v)
    {
        return (float)v.value / FixedPoint<IntBits, FracBits>::ONE;
    }
};

template <>
struct FixedPointCast<double>
{
    template <unsigned IntBits, unsigned FracBits>
    static double cast(FixedPoint<IntBits, FracBits> v)
    {
        return (double)v.value / FixedPoint<IntBits, FracBits>::ONE;
    }
};

/********************************************************/
/*                                                      */
/*                 FixedPointOperations                 */
/*                                                      */
/********************************************************/

/** \addtogroup FixedPointOperations Functions for FixedPoint

    \brief     <b>\#include</b> \<vigra/fixedpoint.hxx\><br>

    These functions fulfill the requirements of an \ref AlgebraicRing.

    Namespace: vigra
    <p>

 */
//@{

    /** Convert a FixedPoint to a built-in type.
        If the target is integral, the value is rounded.<br>
        Usage:
        \code
        FixedPoint<16,15> fp(...);
        
        double d = fixed_point_cast<double>(fp);
        \endcode
    */
template <class TARGET, unsigned IntBits, unsigned FracBits>
TARGET fixed_point_cast(FixedPoint<IntBits, FracBits> v)
{
    return FixedPointCast<TARGET>::cast(v);
}

    /// equal
template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
bool operator==(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    enum { MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1 };
    return (l.value << (MaxFracBits - FracBits1)) == (r.value << (MaxFracBits - FracBits2));
}

    /// not equal
template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
bool operator!=(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    enum { MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1 };
    return (l.value << (MaxFracBits - FracBits1)) != (r.value << (MaxFracBits - FracBits2));
}

    /// less than
template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
bool operator<(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    enum { MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1 };
    return (l.value << (MaxFracBits - FracBits1)) < (r.value << (MaxFracBits - FracBits2));
}

    /// less or equal
template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
bool operator<=(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    enum { MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1 };
    return (l.value << (MaxFracBits - FracBits1)) <= (r.value << (MaxFracBits - FracBits2));
}

    /// greater
template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
bool operator>(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    enum { MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1 };
    return (l.value << (MaxFracBits - FracBits1)) > (r.value << (MaxFracBits - FracBits2));
}

    /// greater or equal
template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
bool operator>=(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    enum { MaxFracBits = (FracBits1 < FracBits2) ? FracBits2 : FracBits1 };
    return (l.value << (MaxFracBits - FracBits1)) >= (r.value << (MaxFracBits - FracBits2));
}

    /// addition with automatic determination of the appropriate result type.
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

    /// addition with enforced result type.
template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2,
          unsigned IntBits3, unsigned FracBits3>
inline void
add(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r,
    FixedPoint<IntBits3, FracBits3> & result)
{
    result = l + r;
}

    /// subtraction with automatic determination of the appropriate result type.
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

    /// subtraction with enforced result type.
template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2,
          unsigned IntBits3, unsigned FracBits3>
inline void
sub(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r,
    FixedPoint<IntBits3, FracBits3> & result)
{
    result = l - r;
}

    /// multiplication with automatic determination of the appropriate result type.
template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
inline
typename FixedPointTraits<FixedPoint<IntBits1, FracBits1>, FixedPoint<IntBits2, FracBits2> >::MultipliesType
operator*(FixedPoint<IntBits1, FracBits1> l, FixedPoint<IntBits2, FracBits2> r)
{
    typename FixedPointTraits<FixedPoint<IntBits1, FracBits1>, FixedPoint<IntBits2, FracBits2> >::
        MultipliesType res;
    mul(l, r, res);
    return res;
}

    /// multiplication with enforced result type.
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

    /// square root.
template <unsigned IntBits, unsigned FracBits>
inline typename SquareRootTraits<FixedPoint<IntBits, FracBits> >::SquareRootResult
sqrt(FixedPoint<IntBits, FracBits> v)
{
    return typename SquareRootTraits<FixedPoint<IntBits, FracBits> >::SquareRootResult(sqrti(v.value), FPNoShift);
}

    /// absolute value.
template <unsigned IntBits, unsigned FracBits>
inline FixedPoint<IntBits, FracBits>
abs(FixedPoint<IntBits, FracBits> v)
{
    return FixedPoint<IntBits, FracBits>(abs(v.value), FPNoShift);
}

    /// squared norm (same as v*v).
template <unsigned IntBits, unsigned FracBits>
inline
typename FixedPointTraits<FixedPoint<IntBits, FracBits>, FixedPoint<IntBits, FracBits> >::MultipliesType
squaredNorm(FixedPoint<IntBits, FracBits> v)
{
    return v*v;
}

    /// norm (same as abs).
template <unsigned IntBits, unsigned FracBits>
inline
FixedPoint<IntBits, FracBits>
norm(FixedPoint<IntBits, FracBits> const & v)
{
    return abs(v);
}

    /// fractional part.
template <unsigned IntBits, unsigned FracBits>
inline FixedPoint<0, FracBits>
frac(FixedPoint<IntBits, FracBits> v)
{
    return FixedPoint<0, FracBits>(v.value & FixedPoint<IntBits, FracBits>::FRACTIONAL_MASK, FPNoShift);
}

    /// dual fractional part: <tt>1 - frac(v)</tt>.
template <unsigned IntBits, unsigned FracBits>
inline FixedPoint<0, FracBits>
dual_frac(FixedPoint<IntBits, FracBits> v)
{
    return FixedPoint<0, FracBits>(FixedPoint<0, FracBits>::ONE - 
                                   (v.value & FixedPoint<IntBits, FracBits>::FRACTIONAL_MASK), FPNoShift);
}

    /// rounding down.
template <unsigned IntBits, unsigned FracBits>
inline int
floor(FixedPoint<IntBits, FracBits> v)
{
    return(v.value >> FracBits);
}

    /// rounding up.
template <unsigned IntBits, unsigned FracBits>
inline int
ceil(FixedPoint<IntBits, FracBits> v)
{
    return((v.value + FixedPoint<IntBits, FracBits>::FRACTIONAL_MASK) >> FracBits);
}

    /// rounding to the nearest integer.
template <unsigned IntBits, unsigned FracBits>
inline int
round(FixedPoint<IntBits, FracBits> v)
{
    return((v.value + FixedPoint<IntBits, FracBits>::ONE_HALF) >> FracBits);
}

//@}

/********************************************************/
/*                                                      */
/*                     FixedPoint-Traits                */
/*                                                      */
/********************************************************/

/** \page FixedPointTraits Numeric and Promote Traits of FixedPoint

    The numeric and promote traits for FixedPoint follow
    the general specifications for \ref NumericPromotionTraits and
    \ref AlgebraicRing. They are implemented in terms of the traits of the basic types by
    partial template specialization:

    \code

    template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
    class FixedPointTraits<FixedPoint<IntBits1, FracBits1>, FixedPoint<IntBits2, FracBits2> >
    {
        typedef FixedPoint<PlusMinusIntBits, MaxFracBits>               PlusType;
        typedef FixedPoint<PlusMinusIntBits, MaxFracBits>               MinusType;
        typedef FixedPoint<IntBits1 + IntBits2, FracBits1 + FracBits2>  MultipliesType;
    };

    template <unsigned IntBits, unsigned FracBits>
    struct NumericTraits<FixedPoint<IntBits, FracBits> >
    {
        typedef FixedPoint<IntBits, FracBits> Type;
            // Promote undefined because it depends on the layout, use FixedPointTraits
            // RealPromote in AlgebraicRing -- multiplication with double is not supported.
            // ComplexPromote in AlgebraicRing -- multiplication with double is not supported.
        typedef Type ValueType;

        typedef VigraFalseType isIntegral;
        typedef VigraTrueType  isScalar;
        typedef VigraTrueType  isSigned;
        typedef VigraTrueType  isOrdered;
        typedef VigraFalseType isComplex;

        ... // etc.
    };

    template <unsigned IntBits, unsigned FracBits>
    struct SquareRootTraits<FixedPoint<IntBits, FracBits> >
    {
        typedef FixedPoint<IntBits, FracBits>      Type;
        typedef FixedPoint<SRIntBits, SRFracBits>  SquareRootResult;
        typedef Type                               SquareRootArgument;
    };
    
    template <unsigned IntBits, unsigned FracBits>
    struct NormTraits<FixedPoint<IntBits, FracBits> >
    {
        typedef FixedPoint<IntBits, FracBits>         Type;
        typedef typename 
            FixedPointTraits<FixedPoint<IntBits, FracBits>, FixedPoint<IntBits, FracBits> >::MultipliesType
                                                      SquaredNormType;
        typedef Type                                  NormType;
    };

    template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
    struct PromoteTraits<FixedPoint<IntBits1, FracBits1>,
                         FixedPoint<IntBits2, FracBits2> >
    {
        typedef typename 
            FixedPointTraits<FixedPoint<IntBits1, FracBits1>, FixedPoint<IntBits2, FracBits2> >::PlusType 
            Promote;
    };
    \endcode

    <b>\#include</b> \<vigra/fixedpoint.hxx\><br>
    Namespace: vigra

*/
template <unsigned IntBits, unsigned FracBits>
struct NumericTraits<FixedPoint<IntBits, FracBits> >
{
    typedef FixedPoint<IntBits, FracBits> Type;
        //typedef FixedPoint<IntBits, FracBits> Promote;
        //typedef FixedPoint<IntBits, FracBits> RealPromote;
        //typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    typedef VigraFalseType isIntegral;
    typedef VigraTrueType  isScalar;
    typedef VigraTrueType  isSigned;
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

template <unsigned IntBits, unsigned FracBits>
struct NormTraits<FixedPoint<IntBits, FracBits> >
{
    typedef FixedPoint<IntBits, FracBits>         Type;
    typedef typename 
        FixedPointTraits<FixedPoint<IntBits, FracBits>, FixedPoint<IntBits, FracBits> >::MultipliesType
                                                  SquaredNormType;
    typedef Type                                  NormType;
};

template <unsigned IntBits1, unsigned FracBits1, unsigned IntBits2, unsigned FracBits2>
struct PromoteTraits<FixedPoint<IntBits1, FracBits1>,
                     FixedPoint<IntBits2, FracBits2> >
{
    typedef typename 
        FixedPointTraits<FixedPoint<IntBits1, FracBits1>, FixedPoint<IntBits2, FracBits2> >::PlusType 
        Promote;
};

/***********************************************************************************/

enum FPOverflowHandling { FPOverflowIgnore, FPOverflowSaturate, FPOverflowError };

template <int IntBits, FPOverflowHandling OverflowHandling = FPOverflowIgnore>
class FixedPoint16;

/********************************************************/
/*                                                      */
/*                     FixedPoint16-Traits              */
/*                                                      */
/********************************************************/

/** \page FixedPoint16Traits Numeric and Promote Traits of FixedPoint16

    The numeric and promote traits for FixedPoint16 follow
    the general specifications for \ref NumericPromotionTraits and
    \ref AlgebraicRing. They are implemented in terms of the traits of the basic types by
    partial template specialization:

    \code
    template <int IntBits, FPOverflowHandling OverflowHandling>
    struct NumericTraits<FixedPoint16<IntBits, OverflowHandling> >
    {
        typedef FixedPoint16<IntBits, OverflowHandling> Type;
        typedef Type                                    Promote;
            // RealPromote undefined -- multiplication with double is not supported.
            // ComplexPromote undefined -- multiplication with double is not supported.
        typedef Type ValueType;

        typedef VigraFalseType isIntegral;
        typedef VigraTrueType  isScalar;
        typedef VigraTrueType  isSigned;
        typedef VigraTrueType  isOrdered;
        typedef VigraFalseType isComplex;

        ... // etc.
    };

    template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2>
    struct PromoteTraits<FixedPoint16<IntBits1, OverflowHandling>,
                         FixedPoint16<IntBits2, OverflowHandling> >
    {
        typedef FixedPoint16<MetaMax<IntBits1, IntBits2>::value, OverflowHandling> Promote;
        ... // etc.
    };

    template <int IntBits, FPOverflowHandling OverflowHandling>
    struct NormTraits<FixedPoint16<IntBits, OverflowHandling> >
    {
        typedef FixedPoint16<IntBits, OverflowHandling>     Type;
        typedef typename PromoteTraits<Type, Type>::Promote SquaredNormType;
        typedef Type                                        NormType;
    };

    template <int IntBits, FPOverflowHandling OverflowHandling>
    struct SquareRootTraits<FixedPoint16<IntBits, OverflowHandling> >
    {
        typedef FixedPoint16<IntBits, OverflowHandling>            Type;
        typedef FixedPoint16<(IntBits + 1) / 2, OverflowHandling>  SquareRootResult;
        typedef Type                                               SquareRootArgument;
    };
    \endcode

    <b>\#include</b> \<vigra/fixedpoint.hxx\><br>
    Namespace: vigra

*/
template <int IntBits, FPOverflowHandling OverflowHandling>
struct NumericTraits<FixedPoint16<IntBits, OverflowHandling> >
{
    typedef FixedPoint16<IntBits, OverflowHandling> Type;
    typedef Type                                    Promote;
        // RealPromote undefined -- multiplication with double is not supported.
        // ComplexPromote undefined -- multiplication with double is not supported.
    typedef Type ValueType;

    typedef VigraFalseType isIntegral;
    typedef VigraTrueType  isScalar;
    typedef VigraTrueType  isSigned;
    typedef VigraTrueType  isOrdered;
    typedef VigraFalseType isComplex;

    static Type zero() { return Type(0, FPNoShift); }
    static Type one() { return Type(Type::ONE, FPNoShift); }
    static Type nonZero() { return one(); }
    static Type epsilon() { return Type(1, FPNoShift); }
    static Type smallestPositive() { return Type(1, FPNoShift); }
    static Type max() { return Type( Type::MAX, FPNoShift); }
    static Type min() { return Type( Type::MIN, FPNoShift); }

    static Promote toPromote(Type v) { return v; }
    static Type fromPromote(Promote v) { return v; }; 
};

template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2>
struct PromoteTraits<FixedPoint16<IntBits1, OverflowHandling>,
                     FixedPoint16<IntBits2, OverflowHandling> >
{
    typedef FixedPoint16<MetaMax<IntBits1, IntBits2>::value, OverflowHandling> Promote;
    static Promote toPromote(FixedPoint16<IntBits1, OverflowHandling> v) { return Promote(v); }
    static Promote toPromote(FixedPoint16<IntBits2, OverflowHandling> v) { return Promote(v); }
};

template <int IntBits, FPOverflowHandling OverflowHandling>
struct PromoteTraits<FixedPoint16<IntBits, OverflowHandling>,
                     FixedPoint16<IntBits, OverflowHandling> >
{
    typedef FixedPoint16<IntBits, OverflowHandling> Promote;
    static Promote toPromote(FixedPoint16<IntBits, OverflowHandling> v) { return v; }
};

template <int IntBits, FPOverflowHandling OverflowHandling>
struct NormTraits<FixedPoint16<IntBits, OverflowHandling> >
{
    typedef FixedPoint16<IntBits, OverflowHandling>     Type;
    typedef typename PromoteTraits<Type, Type>::Promote SquaredNormType;
    typedef Type                                        NormType;
};

template <int IntBits, FPOverflowHandling OverflowHandling>
struct SquareRootTraits<FixedPoint16<IntBits, OverflowHandling> >
{
    typedef FixedPoint16<IntBits, OverflowHandling>            Type;
    typedef FixedPoint16<(IntBits + 1) / 2, OverflowHandling>  SquareRootResult;
    typedef Type                                               SquareRootArgument;
};

#ifndef DOXYGEN

template <bool Compatible>
struct FixedPoint_error__Right_shift_operator_has_unsupported_semantics
: staticAssert::AssertBool<Compatible>
{};

#endif /* DOXYGEN */

template <bool Predicate>
struct FixedPoint16_assignment_error__Target_object_has_too_few_integer_bits
: staticAssert::AssertBool<Predicate>
{};

namespace detail {

template<int BeforeIntBits, int AfterIntBits, 
         bool Round = false,
         bool RightShift = (AfterIntBits >= BeforeIntBits)>
struct FP16Align;

template<int BeforeIntBits>
struct FP16Align<BeforeIntBits, BeforeIntBits, true, true>
{
    static inline Int32 exec(Int32 v)
    {
        return v;
    }
};

template<int BeforeIntBits>
struct FP16Align<BeforeIntBits, BeforeIntBits, false, true>
{
    static inline Int32 exec(Int32 v)
    {
        return v;
    }
};

template<int BeforeIntBits, int AfterIntBits>
struct FP16Align<BeforeIntBits, AfterIntBits, false, true>
{
    static inline Int32 exec(Int32 v)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_error__Right_shift_operator_has_unsupported_semantics<((-1 >> 8) == -1)>));
        return v >> (AfterIntBits - BeforeIntBits);
    }
};

template<int BeforeIntBits, int AfterIntBits>
struct FP16Align<BeforeIntBits, AfterIntBits, true, true>
{
    enum { ONE_HALF = 1 << (AfterIntBits - BeforeIntBits - 1) };
    static inline Int32 exec(Int32 v)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_error__Right_shift_operator_has_unsupported_semantics<((-1 >> 8) == -1)>));
        return (v + ONE_HALF) >> (AfterIntBits - BeforeIntBits);
    }
};

template<int BeforeIntBits, int AfterIntBits, bool Round>
struct FP16Align<BeforeIntBits, AfterIntBits, Round, false>
{
    static inline Int32 exec(Int32 v)
    {
        return v << (BeforeIntBits - AfterIntBits);
    }
};

template <FPOverflowHandling OverflowHandling = FPOverflowIgnore>
struct FP16OverflowHandling
{
    static inline Int32 exec(Int32 v)
    {
        return v;
    }
    
    static inline Int32 exec(UInt32 v)
    {
        return v;
    }
};

template <>
struct FP16OverflowHandling<FPOverflowSaturate>
{
    static inline Int32 exec(Int32 v)
    {
        if(v >= 1 << 15)
            return (1 << 15) - 1;
        if(v < -(1 << 15))
            return -(1 << 15);
        return v;
    }
    static inline Int32 exec(UInt32 v)
    {
        if(v >= 1 << 15)
            return (1 << 15) - 1;
        return v;
    }
};

template <>
struct FP16OverflowHandling<FPOverflowError>
{
    static inline Int32 exec(Int32 v)
    {
        vigra_precondition(v < (1 << 15) && v >= -(1 << 15),
             "FixedPoint16: Operation overflows.");
        return v;
    }
    static inline Int32 exec(UInt32 v)
    {
        vigra_precondition(v < (1 << 15),
             "FixedPoint16: Operation overflows.");
        return v;
    }
};


template <int IntBits1, int IntBits2, int IntBitsOut, 
          FPOverflowHandling OverflowHandling >
struct FP16AddImpl
{
    enum { MinIntBits = MetaMin<IntBits1, IntBits2>::value };
    static inline Int32 exec(Int32 t1, Int32 t2)
    {
       return FP16OverflowHandling<OverflowHandling>::exec(
                  FP16Align<MinIntBits, IntBitsOut, /*Round*/ true>::exec(
                      FP16Align<IntBits1, MinIntBits, /*Round*/ false>::exec(t1) + 
                      FP16Align<IntBits2, MinIntBits, /*Round*/ false>::exec(t2)));
    }
};

template <int IntBits1, int IntBits2, int IntBitsOut, 
          FPOverflowHandling OverflowHandling >
struct FP16SubImpl
{
    enum { MinIntBits = MetaMin<IntBits1, IntBits2>::value };
    static inline Int32 exec(Int32 t1, Int32 t2)
    {
       return FP16OverflowHandling<OverflowHandling>::exec(
                  FP16Align<MinIntBits, IntBitsOut, /*Round*/ true>::exec(
                      FP16Align<IntBits1, MinIntBits, /*Round*/ false>::exec(t1) - 
                      FP16Align<IntBits2, MinIntBits, /*Round*/ false>::exec(t2)));
    }
};

template <int IntBits1, int IntBits2, int IntBitsOut, 
          FPOverflowHandling OverflowHandling >
struct FP16MulImpl
{
    static inline Int32 exec(Int32 t1, Int32 t2)
    {
        return FP16OverflowHandling<OverflowHandling>::exec(
                   FP16Align<IntBits1+IntBits2, IntBitsOut+15, /*Round*/ true>::exec(t1*t2));
    }
};

template <int IntBits1, int IntBits2, int IntBitsOut, 
          FPOverflowHandling OverflowHandling >
struct FP16DivImpl
{
    static inline Int32 exec(Int32 t1, Int32 t2)
    {
        if(t2 == 0)
            return (t1 >= 0)
                       ?  (1 << 15) - 1
                       : -(1 << 15);
        return FP16OverflowHandling<OverflowHandling>::exec(
                   FP16Align<IntBits1-IntBits2, IntBitsOut+1, /*Round*/ true>::exec((t1<<16)/t2));
    }
};

} // namespace detail

/********************************************************/
/*                                                      */
/*                      FixedPoint16                    */
/*                                                      */
/********************************************************/

template <class TARGET, int IntBits, FPOverflowHandling OverflowHandling>
TARGET fixed_point_cast(FixedPoint16<IntBits, OverflowHandling> v);

/** Template for 16-bit signed fixed point arithmetic.

    Fixed point arithmetic is used when computations with fractional accuracy
    must be made at the highest speed possible (e.g. in the inner loop
    of a volume rendering routine). The speed-up relative to floating
    point arithmetic can be dramatic, especially when one can avoid
    conversions between integer and floating point numbers (these are 
    very expensive because integer and floating point arithmetic
    resides in different pipelines). 
    
    The template wraps an <tt>Int16</tt> and uses <tt>IntBits</tt> to
    represent the integral part of a number, and <tt>15 - IntBits</tt>
    for the fractional part. The 16th bit is reserved because FixedPoint16 
    is a signed type. Results of expressions with mixed types will preserve
    larger number of <tt>IntBits</tt> of the results, in order to minimize
    the possibility for overflow. Nonetheless, overflow can occur, and the 
    template parameter <tt>OverflowHandling</tt> determines how this will be
    handled:
    
    <DL>
    <DT>FPOverflowIgnore<DD> (default) Ignore overflow, i.e. use the usual modulo behavior of the
                             built-in integer types.
                       
    <DT>FPOverflowSaturate<DD> Use the largest or smallest representable number (depending on sign)
                               in case of overflow.

    <DT>FPOverflowError<DD> Throw <tt>PreconditionViolation</tt> upon overflow. This is useful for 
                            debugging.
    </DL>
    
    The implementation relies on Int32-arithmetic and requires that the right-shift operator
    preserves signedness. Although not enforced by the C++ standard, this is implemented
    by most of today's processors. This property is checked by a 
    VIGRA_STATIC_ASSERT(FixedPoint_error__Right_shift_operator_has_unsupported_semantics).

    <tt>FixedPoint16</tt> implements the required interface of an
    \ref AlgebraicRing and the required numeric and
    promotion traits. In addition, it supports functions <tt>add</tt>, 
    <tt>sub</tt>, <tt>mul</tt>, and <tt>div</tt>, where a particular layout 
    of the result can be enforced. 
    
    Built-in numeric types can be converted into <tt>FixedPoint16</tt> by the 
    appropriate constructors, and from <tt>FixedPoint16</tt> by means of
    <tt>fixed_point_cast&lt;TargetType&gt;(fixedPoint)</tt>.

    <b>See also:</b>
    <ul>
    <li> \ref FixedPoint16Operations
    <li> \ref FixedPoint16Traits
    </ul>

    <b>\#include</b> \<vigra/fixedpoint.hxx\><br>
    Namespace: vigra
*/
template <int IntBits, FPOverflowHandling OverflowHandling>
class FixedPoint16
{
public:
    static const Int32 TOTAL_BITS      = 15; // bit 16 is sign
    static const Int32 INT_BITS        = IntBits;
    static const Int32 FRACTIONAL_BITS = TOTAL_BITS - INT_BITS;
    static const Int32 MAX             = (Int32)((1u << TOTAL_BITS) - 1);
    static const Int32 MIN             = -(Int32)(1u << TOTAL_BITS);
    static const Int32 ONE             = 1 << FRACTIONAL_BITS;
    static const Int32 ONE_HALF        = ONE >> 1;
    static const Int32 FRACTIONAL_MASK = (1u << FRACTIONAL_BITS) - 1;
    static const Int32 INT_MASK        = 0xffffffffu ^ FRACTIONAL_MASK;
    
    static const FixedPoint16 zero, pi, pi_2, mpi_2; 

    Int16 value;

    FixedPoint16()
    : value(0)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_error__Right_shift_operator_has_unsupported_semantics<((-1 >> 8) == -1)>));
    }

        /** Construct from an int (fractional part will become zero).
            Possible overflow is handled according to the target type's <tt>OverflowHandling</tt>.
        */
    explicit FixedPoint16(Int32 v)
    : value(detail::FP16OverflowHandling<OverflowHandling>::exec(v << FRACTIONAL_BITS))
    {
        VIGRA_STATIC_ASSERT((FixedPoint_error__Right_shift_operator_has_unsupported_semantics<((-1 >> 8) == -1)>));
    }

        /** Construct from an int by a bitwise copy. This is normally only used internally.
        */
    FixedPoint16(Int32 v, FixedPointNoShift)
    : value(detail::FP16OverflowHandling<OverflowHandling>::exec(v))
    {
        VIGRA_STATIC_ASSERT((FixedPoint_error__Right_shift_operator_has_unsupported_semantics<((-1 >> 8) == -1)>));
    }

        /** Construct from a double and round the fractional part to 
            <tt>FRACTIONAL_BITS</tt> accuracy. Possible overflow is handled according 
            to the target type's <tt>OverflowHandling</tt>.
        */
    explicit FixedPoint16(double rhs)
    : value(detail::FP16OverflowHandling<OverflowHandling>::exec((Int32)roundi(rhs * ONE)))
    {
        VIGRA_STATIC_ASSERT((FixedPoint_error__Right_shift_operator_has_unsupported_semantics<((-1 >> 8) == -1)>));
    }

        /** Copy constructor.
        */
    FixedPoint16(const FixedPoint16 &other)
    : value(other.value)
    {
        VIGRA_STATIC_ASSERT((FixedPoint_error__Right_shift_operator_has_unsupported_semantics<((-1 >> 8) == -1)>));
    }

        /** Construct from a FixedPoint16 with different layout. It rounds as appropriate and 
            handles possible overflow according to the target type's <tt>OverflowHandling</tt>.
        */
    template <int IntBits2, FPOverflowHandling OverflowHandling2>
    FixedPoint16(const FixedPoint16<IntBits2, OverflowHandling2> &other)
    : value(detail::FP16OverflowHandling<OverflowHandling>::exec(
            detail::FP16Align<IntBits2, IntBits, /*Round*/true>::exec(other.value)))
    {
        VIGRA_STATIC_ASSERT((FixedPoint_error__Right_shift_operator_has_unsupported_semantics<((-1 >> 8) == -1)>));
    }

        /** Assignment from int. The fractional part will become zero.  
            Possible overflow is handled according to the target type's <tt>OverflowHandling</tt>.
        */
    FixedPoint16 &operator=(Int32 rhs)
    {
        value = detail::FP16OverflowHandling<OverflowHandling>::exec(rhs << FRACTIONAL_BITS);
        return *this;
    }

        /** Assignment form double. The fractional part is rounded, and possible overflow is 
            handled according to the target type's <tt>OverflowHandling</tt>.
        */
    FixedPoint16 &operator=(double rhs)
    {
        value = detail::FP16OverflowHandling<OverflowHandling>::exec(roundi(rhs * ONE));
        return *this;
    }

        /** Copy assignment.
        */
    FixedPoint16 & operator=(const FixedPoint16 &other)
    {
        value = other.value;
        return *this;
    }

        /** Assignment from a FixedPoint16 with different layout. It rounds as appropriate, and possible overflow is 
            handled according to the target type's <tt>OverflowHandling</tt>.
        */
    template <int IntBits2>
    FixedPoint16 & operator=(const FixedPoint16<IntBits2, OverflowHandling> &other)
    {
        value = detail::FP16OverflowHandling<OverflowHandling>::exec(
                detail::FP16Align<IntBits2, IntBits, /*Round*/true>::exec(other.value));
        return *this;
    }

        /** Conversion to float
        */
    operator float() const
    {
        return fixed_point_cast<float>(*this);
    }

        /** Conversion to double
        */
    operator double() const
    {
        return fixed_point_cast<double>(*this);
    }

        /** Unary plus.
        */
    FixedPoint16 operator+() const
    {
        return *this;
    }

        /** Negation.
        */
    FixedPoint16 operator-() const
    {
        return FixedPoint16(-value, FPNoShift);
    }

        /** Pre-increment.
        */
    FixedPoint16 & operator++()
    {
        value = detail::FP16OverflowHandling<OverflowHandling>::exec(value+ONE);
        return *this;
    }

        /** Post-increment.
        */
    FixedPoint16 operator++(int)
    {
        FixedPoint16 old(*this);
        ++(*this);
        return old;
    }

        /** Pre-decrement.
        */
    FixedPoint16 & operator--()
    {
        value = detail::FP16OverflowHandling<OverflowHandling>::exec(value-ONE);
        return *this;
    }

        /** Post-decrement.
        */
    FixedPoint16 operator--(int)
    {
        FixedPoint16 old(*this);
        --(*this);
        return old;
    }

        /** Add-assignment from a FixedPoint16 with different layout. It rounds as appropriate, and possible overflow is 
            handled according to the target type's <tt>OverflowHandling</tt>.
        */
    template <int IntBits2>
    FixedPoint16 & operator+=(const FixedPoint16<IntBits2, OverflowHandling> &other)
    {
        value = detail::FP16AddImpl<IntBits, IntBits2, IntBits, OverflowHandling>::exec(value, other.value);
        return *this;
    }

        /** Subtract-assignment from a FixedPoint16 with different layout. It rounds as appropriate, and possible overflow is 
            handled according to the target type's <tt>OverflowHandling</tt>.
        */
    template <int IntBits2>
    FixedPoint16 & operator-=(const FixedPoint16<IntBits2, OverflowHandling> &other)
    {
        value = detail::FP16SubImpl<IntBits, IntBits2, IntBits, OverflowHandling>::exec(value, other.value);
        return *this;
    }
    
        /** Multiply-assignment from a FixedPoint16 with different layout. It rounds as appropriate, and possible overflow is 
            handled according to the target type's <tt>OverflowHandling</tt>.
        */
    template <int IntBits2>
    FixedPoint16 & operator*=(const FixedPoint16<IntBits2, OverflowHandling> &other)
    {
        value = detail::FP16MulImpl<IntBits, IntBits2, IntBits, OverflowHandling>::exec(value, other.value);
        return *this;
    }
    
        /** Divide-assignment from a FixedPoint16 with different layout. It rounds as appropriate, and possible overflow is 
            handled according to the target type's <tt>OverflowHandling</tt>.
        */
    template <int IntBits2>
    FixedPoint16 & operator/=(const FixedPoint16<IntBits2, OverflowHandling> &other)
    {
        value = detail::FP16DivImpl<IntBits, IntBits2, IntBits, OverflowHandling>::exec(value, other.value);
        return *this;
    }
};

template <int IntBits, FPOverflowHandling OverflowHandling>
const FixedPoint16<IntBits, OverflowHandling> FixedPoint16<IntBits, OverflowHandling>::zero(0);

template <int IntBits, FPOverflowHandling OverflowHandling>
const FixedPoint16<IntBits, OverflowHandling> FixedPoint16<IntBits, OverflowHandling>::pi(M_PI);

template <int IntBits, FPOverflowHandling OverflowHandling>
const FixedPoint16<IntBits, OverflowHandling> FixedPoint16<IntBits, OverflowHandling>::pi_2(0.5 * M_PI);

template <int IntBits, FPOverflowHandling OverflowHandling>
const FixedPoint16<IntBits, OverflowHandling> FixedPoint16<IntBits, OverflowHandling>::mpi_2(-0.5 * M_PI); 

namespace detail {

template <class T>
struct FixedPoint16Cast;

#define VIGRA_FIXED_POINT_CAST(type) \
template <> \
struct FixedPoint16Cast<type> \
{ \
    template <int IntBits, FPOverflowHandling OverflowHandling> \
    static type cast(FixedPoint16<IntBits, OverflowHandling> v) \
    { \
        return round(v); \
    } \
};

VIGRA_FIXED_POINT_CAST(Int8)
VIGRA_FIXED_POINT_CAST(UInt8)
VIGRA_FIXED_POINT_CAST(Int16)
VIGRA_FIXED_POINT_CAST(UInt16)
VIGRA_FIXED_POINT_CAST(Int32)
VIGRA_FIXED_POINT_CAST(UInt32)
VIGRA_FIXED_POINT_CAST(Int64)
VIGRA_FIXED_POINT_CAST(UInt64)

#undef VIGRA_FIXED_POINT_CAST

template <>
struct FixedPoint16Cast<float>
{
    template <int IntBits, FPOverflowHandling OverflowHandling>
    static float cast(FixedPoint16<IntBits, OverflowHandling> v)
    {
        return (float)v.value / FixedPoint16<IntBits, OverflowHandling>::ONE;
    }
};

template <>
struct FixedPoint16Cast<double>
{
    template <int IntBits, FPOverflowHandling OverflowHandling>
    static double cast(FixedPoint16<IntBits, OverflowHandling> v)
    {
        return (double)v.value / FixedPoint16<IntBits, OverflowHandling>::ONE;
    }
};

} // namespace detail

/********************************************************/
/*                                                      */
/*                 FixedPoint16Operations               */
/*                                                      */
/********************************************************/

/** \addtogroup FixedPoint16Operations Functions for FixedPoint16

    \brief     <b>\#include</b> \<vigra/fixedpoint.hxx\><br>

    These functions fulfill the requirements of an \ref AlgebraicRing.

    Namespace: vigra
    <p>

 */
//@{

    /** Convert a FixedPoint16 to a built-in type.
        If the target is integral, the value is rounded.<br>
        Usage:
        \code
        FixedPoint16<16,15> fp(...);
        
        double d = fixed_point_cast<double>(fp);
        \endcode
    */
template <class TARGET, int IntBits, FPOverflowHandling OverflowHandling>
TARGET fixed_point_cast(FixedPoint16<IntBits, OverflowHandling> v)
{
    return detail::FixedPoint16Cast<TARGET>::cast(v);
}


    /// equal
template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2>
inline
bool operator==(FixedPoint16<IntBits1, OverflowHandling> l, FixedPoint16<IntBits2, OverflowHandling> r)
{
    enum { MinIntBits = MetaMin<IntBits1, IntBits2>::value };
    return (l.value << (IntBits1 - MinIntBits)) == (r.value << (IntBits2 - MinIntBits));
}

    /// not equal
template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2>
inline
bool operator!=(FixedPoint16<IntBits1, OverflowHandling> l, FixedPoint16<IntBits2, OverflowHandling> r)
{
    enum { MinIntBits = MetaMin<IntBits1, IntBits2>::value };
    return (l.value << (IntBits1 - MinIntBits)) != (r.value << (IntBits2 - MinIntBits));
}

    /// less than
template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2>
inline
bool operator<(FixedPoint16<IntBits1, OverflowHandling> l, FixedPoint16<IntBits2, OverflowHandling> r)
{
    enum { MinIntBits = MetaMin<IntBits1, IntBits2>::value };
    return (l.value << (IntBits1 - MinIntBits)) < (r.value << (IntBits2 - MinIntBits));
}

    /// less or equal
template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2>
inline
bool operator<=(FixedPoint16<IntBits1, OverflowHandling> l, FixedPoint16<IntBits2, OverflowHandling> r)
{
    enum { MinIntBits = MetaMin<IntBits1, IntBits2>::value };
    return (l.value << (IntBits1 - MinIntBits)) <= (r.value << (IntBits2 - MinIntBits));
}

    /// greater
template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2>
inline
bool operator>(FixedPoint16<IntBits1, OverflowHandling> l, FixedPoint16<IntBits2, OverflowHandling> r)
{
    enum { MinIntBits = MetaMin<IntBits1, IntBits2>::value };
    return (l.value << (IntBits1 - MinIntBits)) > (r.value << (IntBits2 - MinIntBits));
}

    /// greater or equal
template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2>
inline
bool operator>=(FixedPoint16<IntBits1, OverflowHandling> l, FixedPoint16<IntBits2, OverflowHandling> r)
{
    enum { MinIntBits = MetaMin<IntBits1, IntBits2>::value };
    return (l.value << (IntBits1 - MinIntBits)) >= (r.value << (IntBits2 - MinIntBits));
}

    /// addition with automatic determination of the appropriate result type.
template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2>
inline
typename PromoteTraits<FixedPoint16<IntBits1, OverflowHandling>, FixedPoint16<IntBits2, OverflowHandling> >::Promote
operator+(FixedPoint16<IntBits1, OverflowHandling> l, FixedPoint16<IntBits2, OverflowHandling> r)
{
    typedef typename
        PromoteTraits<FixedPoint16<IntBits1, OverflowHandling>, FixedPoint16<IntBits2, OverflowHandling> >::Promote
        Result;
    return Result(detail::FP16AddImpl<IntBits1, IntBits2, Result::INT_BITS, OverflowHandling>::exec(l.value, r.value), FPNoShift);
}

    /// addition with enforced result type.
template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2, int IntBits3>
inline 
FixedPoint16<IntBits3, OverflowHandling> &
add(FixedPoint16<IntBits1, OverflowHandling> l, FixedPoint16<IntBits2, OverflowHandling> r,
    FixedPoint16<IntBits3, OverflowHandling> & result)
{
    result.value = detail::FP16AddImpl<IntBits1, IntBits2, IntBits3, OverflowHandling>::exec(l.value, r.value);
    return result;
}

    /// subtraction with automatic determination of the appropriate result type.
template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2>
inline
typename PromoteTraits<FixedPoint16<IntBits1, OverflowHandling>, FixedPoint16<IntBits2, OverflowHandling> >::Promote
operator-(FixedPoint16<IntBits1, OverflowHandling> l, FixedPoint16<IntBits2, OverflowHandling> r)
{
    typedef typename
        PromoteTraits<FixedPoint16<IntBits1, OverflowHandling>, FixedPoint16<IntBits2, OverflowHandling> >::Promote
        Result;
    return Result(detail::FP16SubImpl<IntBits1, IntBits2, Result::INT_BITS, OverflowHandling>::exec(l.value, r.value), FPNoShift);
}

    /// subtraction with enforced result type.
template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2, int IntBits3>
inline FixedPoint16<IntBits3, OverflowHandling> &
sub(FixedPoint16<IntBits1, OverflowHandling> l, FixedPoint16<IntBits2, OverflowHandling> r,
    FixedPoint16<IntBits3, OverflowHandling> & result)
{
    result.value = detail::FP16SubImpl<IntBits1, IntBits2, IntBits3, OverflowHandling>::exec(l.value, r.value);
    return result;
}

    /// multiplication with automatic determination of the appropriate result type.
template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2>
inline
typename PromoteTraits<FixedPoint16<IntBits1, OverflowHandling>, FixedPoint16<IntBits2, OverflowHandling> >::Promote
operator*(FixedPoint16<IntBits1, OverflowHandling> l, FixedPoint16<IntBits2, OverflowHandling> r)
{
    typedef typename
        PromoteTraits<FixedPoint16<IntBits1, OverflowHandling>, FixedPoint16<IntBits2, OverflowHandling> >::Promote
        Result;
    return Result(detail::FP16MulImpl<IntBits1, IntBits2, Result::INT_BITS, OverflowHandling>::exec(l.value, r.value), FPNoShift);
}

    /// multiplication with enforced result type.
template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2, int IntBits3>
inline 
FixedPoint16<IntBits3, OverflowHandling> &
mul(FixedPoint16<IntBits1, OverflowHandling> l, FixedPoint16<IntBits2, OverflowHandling> r,
    FixedPoint16<IntBits3, OverflowHandling> & result)
{
    result.value = detail::FP16MulImpl<IntBits1, IntBits2, IntBits3, OverflowHandling>::exec(l.value, r.value);
    return result;
}

    /// division with automatic determination of the appropriate result type.
template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2>
inline
typename PromoteTraits<FixedPoint16<IntBits1, OverflowHandling>, FixedPoint16<IntBits2, OverflowHandling> >::Promote
operator/(FixedPoint16<IntBits1, OverflowHandling> l, FixedPoint16<IntBits2, OverflowHandling> r)
{
    typedef typename
        PromoteTraits<FixedPoint16<IntBits1, OverflowHandling>, FixedPoint16<IntBits2, OverflowHandling> >::Promote
        Result;
    return Result(detail::FP16DivImpl<IntBits1, IntBits2, Result::INT_BITS, OverflowHandling>::exec(l.value, r.value), FPNoShift);
}

    /// division with enforced result type.
template <int IntBits1, FPOverflowHandling OverflowHandling, int IntBits2, int IntBits3>
inline 
FixedPoint16<IntBits3, OverflowHandling> &
div(FixedPoint16<IntBits1, OverflowHandling> l, FixedPoint16<IntBits2, OverflowHandling> r,
    FixedPoint16<IntBits3, OverflowHandling> & result)
{
    result.value = detail::FP16DivImpl<IntBits1, IntBits2, IntBits3, OverflowHandling>::exec(l.value, r.value);
    return result;
}

    /// square root.
template <int IntBits, FPOverflowHandling OverflowHandling>
inline typename SquareRootTraits<FixedPoint16<IntBits, OverflowHandling> >::SquareRootResult
sqrt(FixedPoint16<IntBits, OverflowHandling> v)
{
    typedef typename SquareRootTraits<FixedPoint16<IntBits, OverflowHandling> >::SquareRootResult Result;
    enum { Shift = 15 + IntBits - 2*Result::INT_BITS };
    return Result(sqrti(v.value << Shift), FPNoShift);
}

#ifndef VIGRA_NO_HYPOT
    using ::hypot;
#endif

    /// Length of hypotenuse. 
template <int IntBits, FPOverflowHandling OverflowHandling>
inline FixedPoint16<IntBits, OverflowHandling>
hypot(FixedPoint16<IntBits, OverflowHandling> v1, FixedPoint16<IntBits, OverflowHandling> v2)
{
    UInt32 l = abs(v1.value), r = abs(v2.value);   
    // sq(l) + sq(r) <= 2**31, so overflow handling after sqrti is sufficient
    return FixedPoint16<IntBits, OverflowHandling>(
                   detail::FP16OverflowHandling<OverflowHandling>::exec(sqrti(sq(l) + sq(r))), 
                   FPNoShift);
}

using std::atan2;

    /// Arctangent. Accuracy better than 1/3 degree (9 significant bits).
template <int IntBits, FPOverflowHandling OverflowHandling>
FixedPoint16<2, OverflowHandling>
atan2(FixedPoint16<IntBits, OverflowHandling> y, FixedPoint16<IntBits, OverflowHandling> x)
{
    enum { ResIntBits = 2 };
    typedef FixedPoint16<ResIntBits, OverflowHandling> FP;
    static const Int32 Pi_4  = 25736,       // = roundi(0.25 * M_PI * (1 << 15)), at 15 frac bits
                       Pi3_4 = 77208,       // = roundi(0.75 * M_PI * (1 << 15)),
                       c1    = 6497,        // = roundi(0.19826763260224867 * (1 << 15)), 
                       c2    = -1047730238; // = roundi(-0.9757748231899761 * (1 << 30));
                       
    // coefficients c1 and c2 minimize
    //
    // NIntegrate[(c1 r^3 + c2 r + Pi/4 - a)^4 /. r -> (Cos[a] - Sin[a])/(Cos[a] + Sin[a]), {a, 0, Pi/2}]
    //
    // Thanks to Jim Shima, http://www.dspguru.com/comp.dsp/tricks/alg/fxdatan2.htm

    if(x.value == 0)
        return (y.value > 0)
                   ? FP::pi_2
                   : (y.value < 0)
                         ? FP::mpi_2
                         : FP::zero;
                         
    Int32 abs_y = abs(y.value);
    Int32 r, angle;
    if(x.value > 0)
    {
        if(y.value == 0)
            return FP::zero;
        r = ((x.value - abs_y) << 15) / (x.value + abs_y); // 15 frac bits
        angle = Pi_4;
    }
    else
    {
        if(y.value == 0)
            return FP::pi;
        r = ((x.value + abs_y) << 15) / (abs_y - x.value); // 15 frac bits
        angle = Pi3_4;
    }
    
    angle += r*((c2 + c1 * (sq(r) >> 15)) >> 15) >> 15;

    return (y.value > 0)
               ? FP(detail::FP16Align<0, ResIntBits, true>::exec( angle), FPNoShift)
               : FP(detail::FP16Align<0, ResIntBits, true>::exec(-angle), FPNoShift);
}

    /// absolute value.
template <int IntBits, FPOverflowHandling OverflowHandling>
inline FixedPoint16<IntBits, OverflowHandling>
abs(FixedPoint16<IntBits, OverflowHandling> v)
{
    return FixedPoint16<IntBits, OverflowHandling>(abs(v.value), FPNoShift);
}

    /// squared norm (same as v*v).
template <int IntBits, FPOverflowHandling OverflowHandling>
inline
typename NormTraits<FixedPoint16<IntBits, OverflowHandling> >::SquaredNormType
squaredNorm(FixedPoint16<IntBits, OverflowHandling> v)
{
    return v*v;
}

    /// norm (same as abs).
template <int IntBits, FPOverflowHandling OverflowHandling>
inline 
typename NormTraits<FixedPoint16<IntBits, OverflowHandling> >::NormType
norm(FixedPoint16<IntBits, OverflowHandling> const & v)
{
    return abs(v);
}

    /// fractional part. (difference between v and its floor)
template <int IntBits, FPOverflowHandling OverflowHandling>
inline FixedPoint16<IntBits, OverflowHandling>
frac(FixedPoint16<IntBits, OverflowHandling> v)
{
    return FixedPoint16<IntBits, OverflowHandling>(
           v.value - (v.value & FixedPoint16<IntBits, OverflowHandling>::INT_MASK), 
           FPNoShift);
}

    /// dual fractional part. (1 - frac(v))
template <int IntBits, FPOverflowHandling OverflowHandling>
inline FixedPoint16<IntBits, OverflowHandling>
dual_frac(FixedPoint16<IntBits, OverflowHandling> v)
{
    return FixedPoint16<IntBits, OverflowHandling>(
           FixedPoint16<IntBits, OverflowHandling>::ONE - v.value + (v.value & FixedPoint16<IntBits, OverflowHandling>::INT_MASK), 
           FPNoShift);
}

    /// rounding down.
template <int IntBits, FPOverflowHandling OverflowHandling>
inline Int32
floor(FixedPoint16<IntBits, OverflowHandling> v)
{
    return(v.value >> FixedPoint16<IntBits, OverflowHandling>::FRACTIONAL_BITS);
}

    /// rounding up.
template <int IntBits, FPOverflowHandling OverflowHandling>
inline Int32
ceil(FixedPoint16<IntBits, OverflowHandling> v)
{
    return((v.value + FixedPoint16<IntBits, OverflowHandling>::FRACTIONAL_MASK) >> 
                      FixedPoint16<IntBits, OverflowHandling>::FRACTIONAL_BITS);
}

    /// rounding to the nearest integer.
template <int IntBits, FPOverflowHandling OverflowHandling>
inline Int32
round(FixedPoint16<IntBits, OverflowHandling> v)
{
    return((v.value + FixedPoint16<IntBits, OverflowHandling>::ONE_HALF) >> 
                      FixedPoint16<IntBits, OverflowHandling>::FRACTIONAL_BITS);
}

    /// rounding to the nearest integer.
template <int IntBits, FPOverflowHandling OverflowHandling>
inline Int32
roundi(FixedPoint16<IntBits, OverflowHandling> v)
{
    return round(v);
}

//@}

} // namespace vigra

namespace std {

template <int IntBits, vigra::FPOverflowHandling OverflowHandling>
ostream & operator<<(ostream & s, vigra::FixedPoint16<IntBits, OverflowHandling> v)
{
    s << vigra::fixed_point_cast<float>(v);
    return s;
}

} // namespace std

#endif // VIGRA_FIXEDPOINT_HXX
