/************************************************************************/
/*                                                                      */
/*               Copyright 2004-2005 by Ullrich Koethe                  */
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
    conversions between integer anfloating point numbers (these are 
    very expensive because integer and floating point arithmetic
    resides in different pipelines). 
    
    The template wraps an <tt>int</tt> and uses <tt>IntBits</tt> to
    represent the integral part of a number, and <tt>FractionalBits</tt>
    for the fractional part, where <tt>IntBits + FractionalBits &lt; 32</tt>.
    (The 32rd bit is reserved because FixedPoint is a signed type).
    These numbers will be automatically allocated in an intelligent way
    in the result of an arithmetic operation. For example, when two 
    fixed point numbers are multiplied, the required number of integer
    bits in the result is the sum of the number of integer bits of the
    arguments, but only when so many bits are avaiable. This is figured out
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

    <b>\#include</b> "<a href="fixedpoint_8hxx-source.html">vigra/fixedpoint.hxx</a>"<br>
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
    : value(detail::FPAssignWithRound<(Frac2 > FractionalBits)>::template exec<Frac2 - FractionalBits>(other.value))
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
        value = detail::FPAssignWithRound<(Frac2 > FractionalBits)>::template exec<Frac2 - FractionalBits>(other.value);
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

    \brief     <b>\#include</b> "<a href="fixedpoint_8hxx-source.html">vigra/fixedpoint.hxx</a>"<br>

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

    <b>\#include</b> "<a href="fixedpoint_8hxx-source.html">vigra/fixedpoint.hxx</a>"<br>
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

} // namespace vigra

#endif // VIGRA_FIXEDPOINT_HXX
