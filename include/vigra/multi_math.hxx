/************************************************************************/
/*                                                                      */
/*               Copyright 2010-2011 by Ullrich Koethe                  */
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

#ifndef VIGRA_MULTI_MATH_HXX
#define VIGRA_MULTI_MATH_HXX

#include "multi_array.hxx"
#include "tinyvector.hxx"
#include "rgbvalue.hxx"
#include "mathutil.hxx"
#include <complex>

namespace vigra {

/** \defgroup MultiMathModule vigra::multi_math

    Namespace <tt>vigra::multi_math</tt> holds VIGRA's support for efficient arithmetic and algebraic functions on multi-dimensional arrays (that is, \ref MultiArrayView and its subclasses). All <tt>multi_math</tt> functions operate element-wise. If you need matrix multiplication, use \ref LinearAlgebraModule instead.
    
    In order to avoid overload ambiguities, multi-array arithmetic must be explicitly activated by
    \code
    using namespace vigra::multi_math;
    \endcode
    (this should not be done globally, but only in the scope where the functionality is actually used).
    
    You can then use the standard operators in the expected way:
    \code
    MultiArray<2, float> i(Shape2(100, 100)), j(Shape2(100, 100));
    
    MultiArray<2, float> h  = i + 4.0 * j;
                         h += (i.transpose() - j) / 2.0;
    \endcode
    etc. (supported operators are <tt>+ - * / ! ~ % && || == != &lt; &lt;= &gt; &gt;= &lt;&lt; &gt;&gt; & | ^ = += -= *= /=</tt>, with both scalar and array arguments). 
    
    Algebraic functions are available as well:
    \code
    h  = exp(-(sq(i) + sq(j)));
    h *= atan2(-i, j);
    \endcode
    The following functions are implemented: <tt>abs, erf, even, odd, sign, signi, round, roundi, sqrt, sqrti, sq, 
    norm, squaredNorm, gamma, loggamma, exp, log, log10, sin, sin_pi, cos, cos_pi, asin, acos, tan, atan, 
    floor, ceil, conj, real, imag, arg, atan2, pow, fmod, min, max</tt>, 
    provided the array's element type supports the respective function.
    
    Supported element types currently include the built-in numeric types, \ref TinyVector, \ref RGBValue, 
    <tt>std::complex</tt>, and \ref FFTWComplex.

    In addition, <tt>multi_math</tt> supports a number of functions that reduce arrays to scalars:
    \code
    double s = sum<double>(i);  // compute the sum of the elements, using 'double' as accumulator type
    double p = product<double>(abs(i));  // compute the product of the elements' absolute values
    
    bool a = any(i < 0.0);  // check if any element of i is negative
    bool b = all(i > 0.0);  // check if all elements of i are positive
    \endcode
    
    Expressions are expanded so that no temporary arrays have to be created. To optimize cache locality,
    loops are executed in the stride ordering of the left-hand-side array.
    
    <b>\#include</b> \<vigra/multi_math.hxx\>

    Namespace: vigra::multi_math
*/
namespace multi_math {

template <class ARG>
struct MultiMathOperand
{
    typedef typename ARG::result_type result_type;
    
    static const int ndim = ARG::ndim;
    
    MultiMathOperand(ARG const & a)
    : arg_(a)
    {}
        
    // Check if all arrays involved in the expression have compatible shapes
    // (including transparent expansion of singleton axes).
    // 's' is the shape of the LHS array. If 's' is zero (i.e. the LHS is 
    // not yet initialized), it is set to the maximal RHS shape.
    //
    template <class SHAPE>
    bool checkShape(SHAPE & s) const
    {
        return arg_.checkShape(s);
    }
    
    // increment the pointer of all RHS arrays along the given 'axis'
    void inc(unsigned int axis) const
    {
        arg_.inc(axis);
    }
    
    // reset the pointer of all RHS arrays along the given 'axis'
    void reset(unsigned int axis) const
    {
        arg_.reset(axis);
    }
    
    // get the value of the expression at the current pointer location
    result_type operator*() const
    {
        return *arg_;
    }
    
    // get the value of the expression at an offset of the current pointer location
    template <class SHAPE>
    result_type operator[](SHAPE const & s) const
    {
        return arg_[s];
    }
    
    ARG arg_;
};

template <unsigned int N, class T, class C>
struct MultiMathOperand<MultiArrayView<N, T, C> >
{
    typedef MultiMathOperand AllowOverload;    
    typedef typename MultiArrayShape<N>::type Shape;

    typedef T result_type;
    
    static const int ndim = (int)N;
    
    MultiMathOperand(MultiArrayView<N, T, C> const & a)
    : p_(a.data()),
      shape_(a.shape()),
      strides_(a.stride())
    {
        // allow for transparent expansion of singleton dimensions
        for(unsigned int k=0; k<N; ++k)
            if(shape_[k] == 1)
                strides_[k] = 0;
    }
    
    bool checkShape(Shape & s) const
    {
        // support:
        //   * transparent expansion of singleton dimensions
        //   * determining LHS shape in a constructor
        for(unsigned int k=0; k<N; ++k)
        {
            if(shape_[k] == 0)
            {
                return false;
            }
            else if(s[k] <= 1)
            {
                s[k] = shape_[k];
            }
            else if(shape_[k] > 1 && shape_[k] != s[k])
            {
                return false;
            }
        }
        return true;
    }
    
    T const & operator[](Shape const & s) const
    {
        return p_[dot(s, strides_)];
    }
    
    void inc(unsigned int axis) const
    {
        p_ += strides_[axis];
    }
    
    void reset(unsigned int axis) const
    {
        p_ -= shape_[axis]*strides_[axis];
    }
    
    result_type operator*() const
    {
        return *p_;
    }
    
    mutable T const * p_;
    Shape shape_, strides_;
};

template <unsigned int N, class T, class A>
struct MultiMathOperand<MultiArray<N, T, A> >
: public MultiMathOperand<MultiArrayView<N, T, UnstridedArrayTag> >
{
    typedef MultiMathOperand AllowOverload;
    
    MultiMathOperand(MultiArray<N, T, A> const & a)
    : MultiMathOperand<MultiArrayView<N, T, UnstridedArrayTag> >(a)
    {}
};

template <class T>
struct MultiMathScalarOperand
{
    typedef MultiMathOperand<T> AllowOverload;
    typedef T result_type;
    
    static const int ndim = 0;
    
    MultiMathScalarOperand(T const & v)
    : v_(v)
    {}
    
    template <class SHAPE>
    bool checkShape(SHAPE const &) const
    {
        return true;
    }
    
    template <class SHAPE>
    T const & operator[](SHAPE const &) const
    {
        return v_;
    }
    
    void inc(unsigned int /* axis */) const
    {}
    
    void reset(unsigned int /* axis */) const
    {}
    
    T const & operator*() const
    {
        return v_;
    }
    
    T v_;
};

#define VIGRA_CONSTANT_OPERAND(template_dcl, type) \
template template_dcl \
struct MultiMathOperand<type > \
: MultiMathScalarOperand<type > \
{ \
    MultiMathOperand(type const & v) \
    : MultiMathScalarOperand<type >(v) \
    {} \
};

VIGRA_CONSTANT_OPERAND(<>, signed char)
VIGRA_CONSTANT_OPERAND(<>, signed short)
VIGRA_CONSTANT_OPERAND(<>, signed int)
VIGRA_CONSTANT_OPERAND(<>, signed long)
VIGRA_CONSTANT_OPERAND(<>, signed long long)
VIGRA_CONSTANT_OPERAND(<>, unsigned char)
VIGRA_CONSTANT_OPERAND(<>, unsigned short)
VIGRA_CONSTANT_OPERAND(<>, unsigned int)
VIGRA_CONSTANT_OPERAND(<>, unsigned long)
VIGRA_CONSTANT_OPERAND(<>, unsigned long long)
VIGRA_CONSTANT_OPERAND(<>, float)
VIGRA_CONSTANT_OPERAND(<>, double)
VIGRA_CONSTANT_OPERAND(<>, long double)
VIGRA_CONSTANT_OPERAND(<class T>, std::complex<T>)

#define VIGRA_TINYVECTOR_ARGS <class T, int N>
#define VIGRA_TINYVECTOR_DECL TinyVector<T, N>
VIGRA_CONSTANT_OPERAND(VIGRA_TINYVECTOR_ARGS, VIGRA_TINYVECTOR_DECL)
#undef VIGRA_TINYVECTOR_ARGS
#undef VIGRA_TINYVECTOR_DECL

#define VIGRA_RGBVALUE_ARGS <class V, unsigned int R, unsigned int G, unsigned int B>
#define VIGRA_RGBVALUE_DECL RGBValue<V, R, G, B>
VIGRA_CONSTANT_OPERAND(VIGRA_RGBVALUE_ARGS, VIGRA_RGBVALUE_DECL)
#undef VIGRA_RGBVALUE_ARGS
#undef VIGRA_RGBVALUE_DECL

#undef VIGRA_CONSTANT_OPERAND

template <class O, class F>
struct MultiMathUnaryOperator
{
    typedef typename F::template Result<typename O::result_type>::type result_type;
    
    static const int ndim = O::ndim;
                                    
    MultiMathUnaryOperator(O const & o)
    : o_(o)
    {}
    
    template <class SHAPE>
    bool checkShape(SHAPE & s) const
    {
        return o_.checkShape(s);
    }
    
    //
    void inc(unsigned int axis) const
    {
        o_.inc(axis);
    }
    
    void reset(unsigned int axis) const
    {
        o_.reset(axis);
    }
    
    template <class POINT>
    result_type operator[](POINT const & p) const
    {
        return f_(o_[p]);
    }
    
    result_type operator*() const
    {
        return f_(*o_);
    }
    
    O o_;
    F f_;
};

#define VIGRA_MULTIMATH_UNARY_OPERATOR(NAME, FCT, OPNAME, RESTYPE) \
namespace detail { \
struct NAME \
{ \
    template <class T> \
    struct Result \
    { \
        typedef RESTYPE type; \
    }; \
     \
    template <class T> \
    typename Result<T>::type \
    operator()(T const & t) const \
    { \
        return FCT(t); \
    } \
}; \
} \
 \
template <unsigned int N, class T, class C> \
MultiMathOperand<MultiMathUnaryOperator<MultiMathOperand<MultiArrayView<N, T, C> >, \
                                        detail::NAME> > \
OPNAME(MultiArrayView<N, T, C> const & v) \
{ \
    typedef MultiMathOperand<MultiArrayView<N, T, C> > O; \
    typedef MultiMathUnaryOperator<O, detail::NAME> OP; \
    return MultiMathOperand<OP>(OP(v)); \
} \
 \
template <unsigned int N, class T, class A> \
MultiMathOperand<MultiMathUnaryOperator<MultiMathOperand<MultiArray<N, T, A> >, \
                                        detail::NAME> > \
OPNAME(MultiArray<N, T, A> const & v) \
{ \
    typedef MultiMathOperand<MultiArray<N, T, A> > O; \
    typedef MultiMathUnaryOperator<O, detail::NAME> OP; \
    return MultiMathOperand<OP>(OP(v)); \
} \
 \
template <class T> \
MultiMathOperand<MultiMathUnaryOperator<MultiMathOperand<T>, \
                                        detail::NAME> > \
OPNAME(MultiMathOperand<T> const & v) \
{ \
    typedef MultiMathOperand<T> O; \
    typedef MultiMathUnaryOperator<O, detail::NAME> OP; \
    return MultiMathOperand<OP>(OP(v)); \
}

#define VIGRA_REALPROMOTE typename NumericTraits<T>::RealPromote

#ifndef DOXYGEN  // doxygen gets confused by these macros

VIGRA_MULTIMATH_UNARY_OPERATOR(Negate, -, operator-, T)
VIGRA_MULTIMATH_UNARY_OPERATOR(Not, !, operator!, T)
VIGRA_MULTIMATH_UNARY_OPERATOR(BitwiseNot, ~, operator~, T)

using vigra::abs;
VIGRA_MULTIMATH_UNARY_OPERATOR(Abs, vigra::abs, abs, typename NormTraits<T>::NormType)

using vigra::erf;
VIGRA_MULTIMATH_UNARY_OPERATOR(Erf, vigra::erf, erf, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_UNARY_OPERATOR(Even, vigra::even, even, bool)
VIGRA_MULTIMATH_UNARY_OPERATOR(Odd, vigra::odd, odd, bool)
VIGRA_MULTIMATH_UNARY_OPERATOR(Sign, vigra::sign, sign, T)
VIGRA_MULTIMATH_UNARY_OPERATOR(Signi, vigra::signi, signi, int)

using vigra::round;
VIGRA_MULTIMATH_UNARY_OPERATOR(Round, vigra::round, round, T)

VIGRA_MULTIMATH_UNARY_OPERATOR(Roundi, vigra::roundi, roundi, int)
VIGRA_MULTIMATH_UNARY_OPERATOR(Sqrti, vigra::sqrti, sqrti, T)
VIGRA_MULTIMATH_UNARY_OPERATOR(Sq, vigra::sq, sq, typename NumericTraits<T>::Promote)
VIGRA_MULTIMATH_UNARY_OPERATOR(Norm, vigra::norm, norm, typename NormTraits<T>::NormType)
VIGRA_MULTIMATH_UNARY_OPERATOR(SquaredNorm, vigra::squaredNorm, squaredNorm, 
                               typename NormTraits<T>::SquaredNormType)
VIGRA_MULTIMATH_UNARY_OPERATOR(Sin_pi, vigra::sin_pi, sin_pi, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_UNARY_OPERATOR(Cos_pi, vigra::cos_pi, cos_pi, VIGRA_REALPROMOTE)

using vigra::gamma;
VIGRA_MULTIMATH_UNARY_OPERATOR(Gamma, vigra::gamma, gamma, VIGRA_REALPROMOTE)

using vigra::loggamma;
VIGRA_MULTIMATH_UNARY_OPERATOR(Loggamma, vigra::loggamma, loggamma, VIGRA_REALPROMOTE)

VIGRA_MULTIMATH_UNARY_OPERATOR(Sqrt, std::sqrt, sqrt, VIGRA_REALPROMOTE)
using vigra::exp;
VIGRA_MULTIMATH_UNARY_OPERATOR(Exp, vigra::exp, exp, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_UNARY_OPERATOR(Log, std::log, log, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_UNARY_OPERATOR(Log10, std::log10, log10, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_UNARY_OPERATOR(Sin, std::sin, sin, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_UNARY_OPERATOR(Asin, std::asin, asin, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_UNARY_OPERATOR(Cos, std::cos, cos, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_UNARY_OPERATOR(Acos, std::acos, acos, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_UNARY_OPERATOR(Tan, std::tan, tan, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_UNARY_OPERATOR(Atan, std::atan, atan, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_UNARY_OPERATOR(Floor, std::floor, floor, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_UNARY_OPERATOR(Ceil, std::ceil, ceil, VIGRA_REALPROMOTE)

VIGRA_MULTIMATH_UNARY_OPERATOR(Conj, conj, conj, T)
VIGRA_MULTIMATH_UNARY_OPERATOR(Real, real, real, typename T::value_type)
VIGRA_MULTIMATH_UNARY_OPERATOR(Imag, imag, imag, typename T::value_type)
VIGRA_MULTIMATH_UNARY_OPERATOR(Arg, arg, arg, typename T::value_type)

#endif //DOXYGEN

#undef VIGRA_REALPROMOTE
#undef VIGRA_MULTIMATH_UNARY_OPERATOR

template <class O1, class O2, class F>
struct MultiMathBinaryOperator
{
    typedef typename F::template Result<typename O1::result_type,
                                         typename O2::result_type>::type result_type;
                                    
    static const int ndim = O1::ndim > O2::ndim ? O1::ndim : O2::ndim;
    
    MultiMathBinaryOperator(O1 const & o1, O2 const & o2)
    : o1_(o1),
      o2_(o2)
    {}
    
    template <class SHAPE>
    bool checkShape(SHAPE & s) const
    {
        return o1_.checkShape(s) && o2_.checkShape(s);
    }
    
    template <class POINT>
    result_type operator[](POINT const & p) const
    {
        return f_(o1_[p], o2_[p]);
    }
    
    void inc(unsigned int axis) const
    {
        o1_.inc(axis);
        o2_.inc(axis);
    }
    
    void reset(unsigned int axis) const
    {
        o1_.reset(axis);
        o2_.reset(axis);
    }
    
    result_type operator*() const
    {
        return f_(*o1_, *o2_);
    }
    
    O1 o1_;
    O2 o2_;
    F f_;
};


// In the sequel, the nested type 'MultiMathOperand<T>::AllowOverload'
// ensures that template functions only participate in overload
// resolution when this type is defined, i.e. when T is a number 
// or array type. It thus prevents 'ambiguous overload' errors.
//
#define VIGRA_MULTIMATH_BINARY_OPERATOR(NAME, FCT, OPNAME, SEP, RESTYPE) \
\
namespace detail { \
struct NAME \
{ \
    template <class T1, class T2> \
    struct Result \
    { \
        typedef RESTYPE type; \
    }; \
    \
    template <class T1, class T2> \
    typename Result<T1, T2>::type \
    operator()(T1 const & t1, T2 const & t2) const \
    { \
        return FCT(t1 SEP t2); \
    } \
}; \
} \
 \
template <unsigned int N, class T1, class A1, class T2, class A2> \
MultiMathOperand<MultiMathBinaryOperator<MultiMathOperand<MultiArrayView<N, T1> >, \
                                         MultiMathOperand<MultiArrayView<N, T2> >, \
                                         detail::NAME> > \
OPNAME(MultiArray<N, T1, A1> const & v1, MultiArray<N, T2, A2> const & v2) \
{ \
    typedef MultiMathOperand<MultiArrayView<N, T1> > O1; \
    typedef MultiMathOperand<MultiArrayView<N, T2> > O2; \
    typedef MultiMathBinaryOperator<O1, O2, detail::NAME> OP; \
    return MultiMathOperand<OP>(OP((MultiArrayView<N, T1> const &)v1, (MultiArrayView<N, T2> const &)v2)); \
} \
\
template <unsigned int N, class T1, class C1, class T2, class C2> \
MultiMathOperand<MultiMathBinaryOperator<MultiMathOperand<MultiArrayView<N, T1, C1> >, \
                                         MultiMathOperand<MultiArrayView<N, T2, C2> >, \
                                         detail::NAME> > \
OPNAME(MultiArrayView<N, T1, C1> const & v1, MultiArrayView<N, T2, C2> const & v2) \
{ \
    typedef MultiMathOperand<MultiArrayView<N, T1, C1> > O1; \
    typedef MultiMathOperand<MultiArrayView<N, T2, C2> > O2; \
    typedef MultiMathBinaryOperator<O1, O2, detail::NAME> OP; \
    return MultiMathOperand<OP>(OP(v1, v2)); \
} \
\
template <unsigned int N, class T1, class T2, class C2> \
MultiMathOperand<MultiMathBinaryOperator<typename MultiMathOperand<T1>::AllowOverload, \
                                         MultiMathOperand<MultiArrayView<N, T2, C2> >, \
                                         detail::NAME> > \
OPNAME(T1 const & v1, MultiArrayView<N, T2, C2> const & v2) \
{ \
    typedef MultiMathOperand<T1> O1; \
    typedef MultiMathOperand<MultiArrayView<N, T2, C2> > O2; \
    typedef MultiMathBinaryOperator<O1, O2, detail::NAME> OP; \
    return MultiMathOperand<OP>(OP(v1, v2)); \
} \
 \
template <unsigned int N, class T1, class C1, class T2> \
MultiMathOperand<MultiMathBinaryOperator<MultiMathOperand<MultiArrayView<N, T1, C1> >, \
                                         typename MultiMathOperand<T2>::AllowOverload, \
                                         detail::NAME> > \
OPNAME(MultiArrayView<N, T1, C1> const & v1, T2 const & v2) \
{ \
    typedef MultiMathOperand<MultiArrayView<N, T1, C1> > O1; \
    typedef MultiMathOperand<T2> O2; \
    typedef MultiMathBinaryOperator<O1, O2, detail::NAME> OP; \
    return MultiMathOperand<OP>(OP(v1, v2)); \
} \
 \
template <unsigned int N, class T1, class T2, class C2> \
MultiMathOperand<MultiMathBinaryOperator<MultiMathOperand<T1>, \
                                         MultiMathOperand<MultiArrayView<N, T2, C2> >, \
                                         detail::NAME> > \
OPNAME(MultiMathOperand<T1> const & v1, MultiArrayView<N, T2, C2> const & v2) \
{ \
    typedef MultiMathOperand<T1> O1; \
    typedef MultiMathOperand<MultiArrayView<N, T2, C2> > O2; \
    typedef MultiMathBinaryOperator<O1, O2, detail::NAME> OP; \
    return MultiMathOperand<OP>(OP(v1, v2)); \
} \
 \
template <unsigned int N, class T1, class C1, class T2> \
MultiMathOperand<MultiMathBinaryOperator<MultiMathOperand<MultiArrayView<N, T1, C1> >, \
                                         MultiMathOperand<T2>, \
                                         detail::NAME> > \
OPNAME(MultiArrayView<N, T1, C1> const & v1, MultiMathOperand<T2> const & v2) \
{ \
    typedef MultiMathOperand<MultiArrayView<N, T1, C1> > O1; \
    typedef MultiMathOperand<T2> O2; \
    typedef MultiMathBinaryOperator<O1, O2, detail::NAME> OP; \
    return MultiMathOperand<OP>(OP(v1, v2)); \
} \
 \
template <class T1, class T2> \
MultiMathOperand<MultiMathBinaryOperator<MultiMathOperand<T1>, \
                                         MultiMathOperand<T2>, \
                                         detail::NAME> > \
OPNAME(MultiMathOperand<T1> const & v1, MultiMathOperand<T2> const & v2) \
{ \
    typedef MultiMathOperand<T1> O1; \
    typedef MultiMathOperand<T2> O2; \
    typedef MultiMathBinaryOperator<O1, O2, detail::NAME> OP; \
    return MultiMathOperand<OP>(OP(v1, v2)); \
} \
\
template <class T1, class T2> \
MultiMathOperand<MultiMathBinaryOperator<typename MultiMathOperand<T1>::AllowOverload, \
                                         MultiMathOperand<T2>, \
                                         detail::NAME> > \
OPNAME(T1 const & v1, MultiMathOperand<T2> const & v2) \
{ \
    typedef MultiMathOperand<T1> O1; \
    typedef MultiMathOperand<T2> O2; \
    typedef MultiMathBinaryOperator<O1, O2, detail::NAME> OP; \
    return MultiMathOperand<OP>(OP(v1, v2)); \
} \
\
template <class T1, class T2> \
MultiMathOperand<MultiMathBinaryOperator<MultiMathOperand<T1>, \
                                         typename MultiMathOperand<T2>::AllowOverload, \
                                         detail::NAME> > \
OPNAME(MultiMathOperand<T1> const & v1, T2 const & v2) \
{ \
    typedef MultiMathOperand<T1> O1; \
    typedef MultiMathOperand<T2> O2; \
    typedef MultiMathBinaryOperator<O1, O2, detail::NAME> OP; \
    return MultiMathOperand<OP>(OP(v1, v2)); \
}

#define VIGRA_NOTHING
#define VIGRA_COMMA ,
#define VIGRA_PROMOTE typename PromoteTraits<T1, T2>::Promote
#define VIGRA_REALPROMOTE typename PromoteTraits<typename NumericTraits<T1>::RealPromote, \
                                                 typename NumericTraits<T2>::RealPromote>::Promote

VIGRA_MULTIMATH_BINARY_OPERATOR(Plus, VIGRA_NOTHING, operator+, +, VIGRA_PROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(Minus, VIGRA_NOTHING, operator-, -, VIGRA_PROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(Multiplies, VIGRA_NOTHING, operator*, *, VIGRA_PROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(Divides, VIGRA_NOTHING, operator/, /, VIGRA_PROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(Modulo, VIGRA_NOTHING, operator%, %, VIGRA_PROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(And, VIGRA_NOTHING, operator&&, &&, VIGRA_PROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(Or, VIGRA_NOTHING, operator||, ||, VIGRA_PROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(Equal, VIGRA_NOTHING, operator==, ==, bool)
VIGRA_MULTIMATH_BINARY_OPERATOR(NotEqual, VIGRA_NOTHING, operator!=, !=, bool)
VIGRA_MULTIMATH_BINARY_OPERATOR(Less, VIGRA_NOTHING, operator<, <, bool)
VIGRA_MULTIMATH_BINARY_OPERATOR(LessEqual, VIGRA_NOTHING, operator<=, <=, bool)
VIGRA_MULTIMATH_BINARY_OPERATOR(Greater, VIGRA_NOTHING, operator>, >, bool)
VIGRA_MULTIMATH_BINARY_OPERATOR(GreaterEqual, VIGRA_NOTHING, operator>=, >=, bool)
VIGRA_MULTIMATH_BINARY_OPERATOR(Leftshift, VIGRA_NOTHING, operator<<, <<, VIGRA_PROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(Rightshift, VIGRA_NOTHING, operator>>, >>, VIGRA_PROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(BitwiseAnd, VIGRA_NOTHING, operator&, &, VIGRA_PROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(BitwiseOr, VIGRA_NOTHING, operator|, |, VIGRA_PROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(BitwiseXor, VIGRA_NOTHING, operator^, ^, VIGRA_PROMOTE)

VIGRA_MULTIMATH_BINARY_OPERATOR(Atan2, std::atan2, atan2, VIGRA_COMMA, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(Pow, std::pow, pow, VIGRA_COMMA, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(Fmod, std::fmod, fmod, VIGRA_COMMA, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(Min, std::min, min, VIGRA_COMMA, VIGRA_PROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(Max, std::max, max, VIGRA_COMMA, VIGRA_PROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(Minimum, std::min, minimum, VIGRA_COMMA, VIGRA_PROMOTE)
VIGRA_MULTIMATH_BINARY_OPERATOR(Maximum, std::max, maximum, VIGRA_COMMA, VIGRA_PROMOTE)

#undef VIGRA_NOTHING
#undef VIGRA_COMMA
#undef VIGRA_PROMOTE
#undef VIGRA_REALPROMOTE
#undef VIGRA_MULTIMATH_BINARY_OPERATOR

namespace detail {

// We pass 'strideOrder' to the recursion in order to make sure
// that the inner loop iterates over the output's major axis.
// Of course, this does not help when the RHS arrays are ordered 
// differently -- maybe it is better to find the most common order
// among all arguments (both RHS and LHS)?
//
template <unsigned int N, class Assign>
struct MultiMathExec
{
    enum { LEVEL = N-1 };
    
    template <class T, class Shape, class Expression>
    static void exec(T * data, Shape const & shape, Shape const & strides, 
                     Shape const & strideOrder, Expression const & e)
    {
        MultiArrayIndex axis = strideOrder[LEVEL];
        for(MultiArrayIndex k=0; k<shape[axis]; ++k, data += strides[axis], e.inc(axis))
        {
            MultiMathExec<N-1, Assign>::exec(data, shape, strides, strideOrder, e);
        }
        e.reset(axis);
        data -= shape[axis]*strides[axis];
    }
};

template <class Assign>
struct MultiMathExec<1, Assign>
{
    enum { LEVEL = 0 };
    
    template <class T, class Shape, class Expression>
    static void exec(T * data, Shape const & shape, Shape const & strides, 
                     Shape const & strideOrder, Expression const & e)
    {
        MultiArrayIndex axis = strideOrder[LEVEL];
        for(MultiArrayIndex k=0; k<shape[axis]; ++k, data += strides[axis], e.inc(axis))
        {
            Assign::assign(data, e);
        }
        e.reset(axis);
        data -= shape[axis]*strides[axis];
    }
};

#define VIGRA_MULTIMATH_ASSIGN(NAME, OP) \
struct MultiMath##NAME \
{ \
    template <class T, class Expression> \
    static void assign(T * data, Expression const & e) \
    { \
        *data OP vigra::detail::RequiresExplicitCast<T>::cast(*e); \
    } \
}; \
 \
template <unsigned int N, class T, class C, class Expression> \
void NAME(MultiArrayView<N, T, C> a, MultiMathOperand<Expression> const & e) \
{ \
    typename MultiArrayShape<N>::type shape(a.shape()); \
     \
    vigra_precondition(e.checkShape(shape), \
       "multi_math: shape mismatch in expression."); \
        \
    MultiMathExec<N, MultiMath##NAME>::exec(a.data(), a.shape(), a.stride(), \
                                            a.strideOrdering(), e); \
} \
 \
template <unsigned int N, class T, class A, class Expression> \
void NAME##OrResize(MultiArray<N, T, A> & a, MultiMathOperand<Expression> const & e) \
{ \
    typename MultiArrayShape<N>::type shape(a.shape()); \
     \
    vigra_precondition(e.checkShape(shape), \
       "multi_math: shape mismatch in expression."); \
        \
    if(a.size() == 0) \
        a.reshape(shape); \
         \
    MultiMathExec<N, MultiMath##NAME>::exec(a.data(), a.shape(), a.stride(), \
                                            a.strideOrdering(), e); \
}

VIGRA_MULTIMATH_ASSIGN(assign, =)
VIGRA_MULTIMATH_ASSIGN(plusAssign, +=)
VIGRA_MULTIMATH_ASSIGN(minusAssign, -=)
VIGRA_MULTIMATH_ASSIGN(multiplyAssign, *=)
VIGRA_MULTIMATH_ASSIGN(divideAssign, /=)

#undef VIGRA_MULTIMATH_ASSIGN

template <unsigned int N, class Assign>
struct MultiMathReduce
{
    enum { LEVEL = N-1 };
    
    template <class T, class Shape, class Expression>
    static void exec(T & t, Shape const & shape, Expression const & e)
    {
        for(MultiArrayIndex k=0; k<shape[LEVEL]; ++k, e.inc(LEVEL))
        {
            MultiMathReduce<N-1, Assign>::exec(t, shape, e);
        }
        e.reset(LEVEL);
    }
};

template <class Assign>
struct MultiMathReduce<1, Assign>
{
    enum { LEVEL = 0 };
    
    template <class T, class Shape, class Expression>
    static void exec(T & t, Shape const & shape, Expression const & e)
    {
        for(MultiArrayIndex k=0; k<shape[0]; ++k, e.inc(0))
        {
            Assign::assign(&t, e);
        }
        e.reset(0);
    }
};

struct MultiMathReduceAll
{
    template <class T, class Expression>
    static void assign(T * data, Expression const & e)
    {
        *data = *data && (*e != NumericTraits<typename Expression::result_type>::zero());
    }
};

struct MultiMathReduceAny
{
    template <class T, class Expression>
    static void assign(T * data, Expression const & e)
    {
        *data = *data || (*e != NumericTraits<typename Expression::result_type>::zero());
    }
};


} // namespace detail

template <class U, class T>
U
sum(MultiMathOperand<T> const & v, U res = NumericTraits<U>::zero()) 
{ 
    static const int ndim = MultiMathOperand<T>::ndim;
    typename MultiArrayShape<ndim>::type shape;
    v.checkShape(shape);
    detail::MultiMathReduce<ndim, detail::MultiMathplusAssign>::exec(res, shape, v);
    return res;
}

template <class U, unsigned int N, class T, class S>
U
sum(MultiArrayView<N, T, S> const & v, U res = NumericTraits<U>::zero()) 
{ 
    return v.sum<U>() + res;
}

template <class U, class T>
U
product(MultiMathOperand<T> const & v, U res = NumericTraits<U>::one()) 
{ 
    static const int ndim = MultiMathOperand<T>::ndim;
    typename MultiArrayShape<ndim>::type shape;
    v.checkShape(shape);
    detail::MultiMathReduce<ndim, detail::MultiMathmultiplyAssign>::exec(res, shape, v);
    return res;
}

template <class U, unsigned int N, class T, class S>
U
product(MultiArrayView<N, T, S> const & v, U res = NumericTraits<U>::one()) 
{ 
    return v.product<U>() * res;
}

template <class T>
bool
all(MultiMathOperand<T> const & v) 
{ 
    static const int ndim = MultiMathOperand<T>::ndim;
    typename MultiArrayShape<ndim>::type shape;
    v.checkShape(shape);
    bool res = true;
    detail::MultiMathReduce<ndim, detail::MultiMathReduceAll>::exec(res, shape, v);
    return res;
}

template <class T>
bool
any(MultiMathOperand<T> const & v) 
{ 
    static const int ndim = MultiMathOperand<T>::ndim;
    typename MultiArrayShape<ndim>::type shape;
    v.checkShape(shape);
    bool res = false;
    detail::MultiMathReduce<ndim, detail::MultiMathReduceAny>::exec(res, shape, v);
    return res;
}


}} // namespace vigra::multi_math

#endif // VIGRA_MULTI_MATH_HXX
