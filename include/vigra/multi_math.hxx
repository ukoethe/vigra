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

namespace multi_math {

template <class ARG>
struct MultiMathOperand
{
	typedef typename ARG::result_type result_type;
    
    MultiMathOperand(ARG const & a)
    : arg_(a)
    {}
        
    template <class SHAPE>
    bool checkShape(SHAPE & s) const
    {
        return arg_.checkShape(s);
    }
    
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
    
    T const * p_;
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

#define VIGRA_CONSTANT_OPERAND(template_dcl, type) \
template template_dcl \
struct MultiMathOperand<type > \
{ \
    typedef MultiMathOperand<type > AllowOverload; \
	typedef type result_type; \
     \
    MultiMathOperand(type const & v) \
    : v_(v) \
    {} \
     \
    template <class SHAPE> \
    bool checkShape(SHAPE const &) const \
    { \
        return true; \
    } \
     \
    template <class SHAPE> \
    type const & operator[](SHAPE const &) const \
    { \
        return v_; \
    } \
     \
    type v_; \
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
                                    
    MultiMathUnaryOperator(O const & o)
    : o_(o)
    {}
    
    template <class SHAPE>
    bool checkShape(SHAPE & s) const
    {
        return o_.checkShape(s);
    }
    
    template <class POINT>
    result_type operator[](POINT const & p) const
    {
        return f_(o_[p]);
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

VIGRA_MULTIMATH_UNARY_OPERATOR(Negate, -, operator-, T)
VIGRA_MULTIMATH_UNARY_OPERATOR(Not, !, operator!, T)
VIGRA_MULTIMATH_UNARY_OPERATOR(BitwiseNot, ~, operator~, T)

using vigra::abs;
VIGRA_MULTIMATH_UNARY_OPERATOR(Abs, vigra::abs, abs, T)

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

using vigra::gamma;
VIGRA_MULTIMATH_UNARY_OPERATOR(Loggamma, vigra::loggamma, loggamma, VIGRA_REALPROMOTE)

VIGRA_MULTIMATH_UNARY_OPERATOR(Sqrt, std::sqrt, sqrt, VIGRA_REALPROMOTE)
VIGRA_MULTIMATH_UNARY_OPERATOR(Exp, std::exp, exp, VIGRA_REALPROMOTE)
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

#undef VIGRA_REALPROMOTE
#undef VIGRA_MULTIMATH_UNARY_OPERATOR

template <class O1, class O2, class F>
struct MultiMathBinaryOperator
{
    typedef typename F::template Result<typename O1::result_type,
                                         typename O2::result_type>::type result_type;
                                    
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
    
    O1 o1_;
    O2 o2_;
    F f_;
};

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

#define VIGRA_MULTIMATH_ASSIGN(NAME, OP) \
template <unsigned int N, class T, class C, class Expression> \
void NAME(MultiArrayView<N, T, C> a, MultiMathOperand<Expression> const & e) \
{ \
    typename MultiArrayShape<N>::type shape(a.shape()); \
     \
    vigra_precondition(e.checkShape(shape), \
       "multi_math: shape mismatch in expression."); \
        \
    typename MultiArrayView<N, T, C>::iterator i = a.begin(), end = a.end(); \
    for(; i != end; ++i) \
	*i OP (e[i.point()]); \
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
    typename MultiArrayView<N, T>::iterator i = static_cast<MultiArrayView<N, T> &>(a).begin(), \
                                            end = i.getEndIterator(); \
    for(; i != end; ++i) \
        *i OP (e[i.point()]); \
}

VIGRA_MULTIMATH_ASSIGN(assign, = vigra::detail::RequiresExplicitCast<T>::cast)
VIGRA_MULTIMATH_ASSIGN(plusAssign, +=)
VIGRA_MULTIMATH_ASSIGN(minusAssign, -=)
VIGRA_MULTIMATH_ASSIGN(multiplyAssign, *=)
VIGRA_MULTIMATH_ASSIGN(divideAssign, /=)

#undef VIGRA_MULTIMATH_ASSIGN

} // namespace detail

}} // namespace vigra::multi_math

#endif // VIGRA_MULTI_MATH_HXX