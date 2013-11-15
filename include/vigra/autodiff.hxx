/************************************************************************/
/*                                                                      */
/*               Copyright 2012-2013 by Ullrich Koethe                  */
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


#ifndef VIGRA_AUTODIFF_HXX
#define VIGRA_AUTODIFF_HXX

#include "tinyvector.hxx"
#include "mathutil.hxx"
#include <cmath>

namespace vigra {

namespace autodiff {

/** Number type for automatic differentiation.

    <a href="http://en.wikipedia.org/wiki/Automatic_differentiation">Automatic differentiation</a> 
    allows one to compute the value of a numeric expression and its gradient 
    with respect to the expression's arguments automatically and in one go. 
    To support this, one needs a special number type that holds a scalar value
    and the corresponding gradient vector of appropriate length. This is the
    purpose of hte template class <tt>DualVector<T, N></tt>, where <tt>T</tt> is the
    underlying numerical type (usually 'double'), and <tt>N</tt> denotes the 
    length of the gradient vector. 
    
    The standard arithmetic and algebraic functions are overloaded for 
    <tt>DualVector</tt> in order to implement the required arithmetic of 
    dual numbers. When you replace all arguments in a numeric expression
    with the appropriate <tt>DualVector</tt> instances, the result will be
    a <tt>DualVector</tt> that contains the result value and gradient of
    the expression, evaluated at the point defined by the input values.
    
    <b> Usage:</b>
    
    <b>\#include</b> \<vigra/autodiff.hxx\><br>
    Namespace: vigra::autodiff

    \code
    typedef DualVector<double, 2> N;  // for expressions with two arguments
    
    N x(1.0, 0);  // first argument of the expression
    N s(2.0, 1);  // second argument of the expression
    
    N y = exp(-0.5 * sq(x / s));
    
    std::cout << "Evaluated exp(- x^2 / (2 s^2)) at x=1 and s = 2:\n";
    std::cout << "result = " << y.value() <<", gradient = " << y.gradient() << "\n";
    \endcode
    Note that the second argument of the <tt>DualVector</tt> constructors specifies that 
    the derivative w.r.t 'x' shall be the element 0 of the gradient vector, and the
    derivative w.r.t. 's' shall be element 1.
*/
template <class T, int N>
class DualVector
{
  public:
    typedef T                 value_type; ///< type of function values and gradient elements
    typedef TinyVector<T, N>  Gradient;   ///< type of the gradient vector
    
    T v;
    Gradient d;
    
        /** Zero initialization.
        */
    DualVector()
    : v(), d()
    {}
    
        /** Provide a value, but zero-initialize the gradient.
        */
    explicit DualVector(T const & val)
    : v(val), d()
    {}
    
        /** Initialize with given value and gradient.
        */
    DualVector(T const & val, Gradient const & grad)
    : v(val), d(grad)
    {}
    
        /** Shorthand for <tt>DualVector(val, Gradient(g0))</tt> when <tt>N == 1</tt>.
        
            Not to be used when <tt>N != 1</tt>.
        */
    DualVector(T const & val, T const & g0)
    : v(val), d(g0)
    {}
    
        /** Shorthand for <tt>DualVector(val, Gradient(g0, g1))</tt> when <tt>N == 2</tt>.
        
            Not to be used when <tt>N != 2</tt>.
        */
    DualVector(T const & val, T const & g0, T const & g1)
    : v(val), d(g0, g1)
    {}
    
        /** Initialize value to represent the argument number 'targetElement' in an
            expression. 

            The derivative of the expression w.r.t. this variable will be element 'targetElement'
            of the resulting gradient vector.
        */
    DualVector(T const & val, int targetElement)
    : v(val), d()
    {
        d[targetElement] = T(1.0);
    }
    
        /** Get current value.
        */
    T value() const
    {
        return v;
    }
    
        /** Get current gradient.
        */
    Gradient const & gradient() const
    {
        return d;
    }
    
    DualVector operator+() const
    {
        return *this;
    }
    
    DualVector operator-() const
    {
        return DualVector(-v, -d);
    }
    
    DualVector & operator+=(DualVector const & o)
    {
        d += o.d;
        v += o.v;
        return *this;
    }
    
    DualVector & operator+=(T const & o)
    {
        v += o;
        return *this;
    }
    
    DualVector & operator-=(DualVector const & o)
    {
        d -= o.d;
        v -= o.v;
        return *this;
    }
    
    DualVector & operator-=(T const & o)
    {
        v -= o;
        return *this;
    }
    
    DualVector & operator*=(DualVector const & o)
    {
        d = o.v * d + v * o.d;
        v *= o.v;
        return *this;
    }
    
    DualVector & operator*=(T const & o)
    {
        d *= o;
        v *= o;
        return *this;
    }
    
    DualVector & operator/=(DualVector const & o)
    {
        d = (o.v * d - v * o.d) / sq(o.v);
        v /= o.v;
        return *this;
    }
    
    DualVector & operator/=(T const & o)
    {
        d /= o;
        v /= o;
        return *this;
    }
};

    /** Given a vector 'v' of expression arguments, create the corresponding
        vector of dual numbers for automatic differentiation.
    */
template <class T, int N>
TinyVector<DualVector<T, N>, N>
dualMatrix(TinyVector<T, N> const & v)
{
    TinyVector<DualVector<T, N>, N> res;
    for(int k=0; k<N; ++k)
    {
        res[k].v = v[k];
        res[k].d[k] = T(1.0);
    }
    return res;
}

template <class T, int N>
inline DualVector<T, N> operator+(DualVector<T, N> v1, DualVector<T, N> const & v2)
{
    return v1 += v2;
}

template <class T, int N>
inline DualVector<T, N> operator+(DualVector<T, N> v1, T v2)
{
    return v1 += v2;
}

template <class T, int N>
inline DualVector<T, N> operator+(T v1, DualVector<T, N> v2)
{
    return v2 += v1;
}

template <class T, int N>
inline DualVector<T, N> operator-(DualVector<T, N> v1, DualVector<T, N> const & v2)
{
    return v1 -= v2;
}

template <class T, int N>
inline DualVector<T, N> operator-(DualVector<T, N> v1, T v2)
{
    return v1 -= v2;
}

template <class T, int N>
inline DualVector<T, N> operator-(T v1, DualVector<T, N> const & v2)
{
    return DualVector<T, N>(v1 - v2.v, -v2.d);
}

template <class T, int N>
inline DualVector<T, N> operator*(DualVector<T, N> v1, DualVector<T, N> const & v2)
{
    return v1 *= v2;
}

template <class T, int N>
inline DualVector<T, N> operator*(DualVector<T, N> v1, T v2)
{
    return v1 *= v2;
}

template <class T, int N>
inline DualVector<T, N> operator*(T v1, DualVector<T, N> v2)
{
    return v2 *= v1;
}

template <class T, int N>
inline DualVector<T, N> operator/(DualVector<T, N> v1, DualVector<T, N> const & v2)
{
    return v1 /= v2;
}

template <class T, int N>
inline DualVector<T, N> operator/(DualVector<T, N> v1, T v2)
{
    return v1 /= v2;
}

template <class T, int N>
inline DualVector<T, N> operator/(T v1, DualVector<T, N> const & v2)
{
    return DualVector<T, N>(v1 / v2.v, -v1*v2.d / sq(v2.v));
}

using vigra::abs;
// abs(x + h) => x + h or -(x + h)
template <typename T, int N>
inline DualVector<T, N> abs(DualVector<T, N> const & v)
{
    return v.v < T(0.0) ? -v : v;
}

using std::fabs;
// abs(x + h) => x + h or -(x + h)
template <typename T, int N>
inline DualVector<T, N> fabs(DualVector<T, N> const & v)
{
    return v.v < T(0.0) ? -v : v;
}

using std::log;
// log(a + h) => log(a) + h / a
template <typename T, int N>
inline DualVector<T, N> log(DualVector<T, N> v) 
{
    v.d /= v.v;
    v.v = log(v.v);
    return v;
}

using std::exp;
// exp(a + h) => exp(a) + exp(a) h
template <class T, int N>
inline DualVector<T, N> exp(DualVector<T, N> v)
{
    v.v = exp(v.v);
    v.d *= v.v;
    return v;
}

using vigra::sqrt;
// sqrt(a + h) => sqrt(a) + h / (2 sqrt(a))
template <typename T, int N>
inline DualVector<T, N> sqrt(DualVector<T, N> v)
{
    v.v = sqrt(v.v);
    v.d /= T(2.0) * v.v;
    return v;
}

using std::sin;
using std::cos;
// sin(a + h) => sin(a) + cos(a) h
template <typename T, int N>
inline DualVector<T, N> sin(DualVector<T, N> v)
{
    v.d *= cos(v.v);
    v.v = sin(v.v);
    return v;
}

// cos(a + h) => cos(a) - sin(a) h
template <typename T, int N>
inline DualVector<T, N> cos(DualVector<T, N> v)
{
    v.d *= -sin(v.v);
    v.v = cos(v.v);
    return v;
}

using vigra::sin_pi;
using vigra::cos_pi;
// sin_pi(a + h) => sin_pi(a) + pi cos_pi(a) h
template <typename T, int N>
inline DualVector<T, N> sin_pi(DualVector<T, N> v)
{
    v.d *= M_PI*cos_pi(v.v);
    v.v = sin_pi(v.v);
    return v;
}

// cos_pi(a + h) => cos_pi(a) - pi sin_pi(a) h
template <typename T, int N>
inline DualVector<T, N> cos_pi(DualVector<T, N> v)
{
    v.d *= -M_PI*sin_pi(v.v);
    v.v = cos_pi(v.v);
    return v;
}

using std::asin;
// asin(a + h) => asin(a) + 1 / sqrt(1 - a^2) h
template <typename T, int N>
inline DualVector<T, N> asin(DualVector<T, N> v)
{
    v.d /= sqrt(T(1.0) - sq(v.v));
    v.v = asin(v.v);
    return v;
}

using std::acos;
// acos(a + h) => acos(a) - 1 / sqrt(1 - a^2) h
template <typename T, int N>
inline DualVector<T, N> acos(DualVector<T, N> v) 
{
    v.d /= -sqrt(T(1.0) - sq(v.v));
    v.v = acos(v.v);
    return v;
}

using std::tan;
// tan(a + h) => tan(a) + (1 + tan(a)^2) h
template <typename T, int N>
inline DualVector<T, N> tan(DualVector<T, N> v)
{
    v.v = tan(v.v);
    v.d *= T(1.0) + sq(v.v);
    return v;
}

using std::atan;
// atan(a + h) => atan(a) + 1 / (1 + a^2) h
template <typename T, int N>
inline DualVector<T, N> atan(DualVector<T, N> v)
{
    v.d /= T(1.0) + sq(v.v);
    v.v = atan(v.v);
    return v;
}

using std::sinh;
using std::cosh;
// sinh(a + h) => sinh(a) + cosh(a) h
template <typename T, int N>
inline DualVector<T, N> sinh(DualVector<T, N> v)
{
    v.d *= cosh(v.v);
    v.v = sinh(v.v);
    return v;
}

// cosh(a + h) => cosh(a) + sinh(a) h
template <typename T, int N>
inline DualVector<T, N> cosh(DualVector<T, N> v) 
{
    v.d *= sinh(v.v);
    v.v = cosh(v.v);
    return v;
}

using std::tanh;
// tanh(a + h) => tanh(a) + (1 - tanh(a)^2) h
template <typename T, int N>
inline DualVector<T, N> tanh(DualVector<T, N> v)
{
    v.v = tanh(v.v);
    v.d *= T(1.0) - sq(v.v);
    return v;
}

using vigra::sq;
// (a + h)^2 => a^2 + 2 a h
template <class T, int N>
inline DualVector<T, N> sq(DualVector<T, N> v)
{
    v.d *= T(2.0)*v.v;
    v.v *= v.v;
    return v;
}

using std::atan2;
// atan2(b + db, a + da) => atan2(b, a) + (- b da + a db) / (a^2 + b^2)
template <typename T, int N> 
inline DualVector<T, N> atan2(DualVector<T, N> v1, DualVector<T, N> const & v2) 
{
    v1.d = (v2.v * v1.d - v1.v * v2.d) / (sq(v1.v) + sq(v2.v));
    v1.v = atan2(v1.v, v2.v);
    return v1;
}


using vigra::pow;
// (a+da)^p => a^p + p*a^(p-1) da
template <typename T, int N> 
inline DualVector<T, N> pow(DualVector<T, N> v, T p)
{
    T pow_p_1 = pow(v.v, p-T(1.0));
    v.d *= p * pow_p_1;
    v.v *= pow_p_1;
    return v;
}

// (a)^(p+dp) => a^p + a^p log(a) dp
template <typename T, int N> 
inline DualVector<T, N> pow(T v, DualVector<T, N> p)
{
    p.v = pow(v, p.v);
    p.d *= p.v * log(v);
    return p;
}


// (a+da)^(b+db) => a^b + b * a^(b-1) da + a^b log(a) * db
template <typename T, int N> 
inline DualVector<T, N> pow(DualVector<T, N> v, DualVector<T, N> const & p)
{
    T pow_p_1 = pow(v.v, p.v-T(1.0)),
      pow_p   = v.v * pow_p_1;
    v.d = p.v * pow_p_1 * v.d + pow_p * log(v.v) * p.d;
    v.v = pow_p;
    return v;
}

using vigra::min;
template <class T, int N>
inline DualVector<T, N> min(DualVector<T, N> const & v1, DualVector<T, N> const & v2)
{
    return v1.v < v2.v
               ? v1
               : v2;
}

template <class T, int N>
inline DualVector<T, N> min(T v1, DualVector<T, N> const & v2)
{
    return v1 < v2.v
               ? DualVector<T, N>(v1)
               : v2;
}

template <class T, int N>
inline DualVector<T, N> min(DualVector<T, N> const & v1, T v2)
{
    return v1.v < v2
               ? v1
               : DualVector<T, N>(v2);
}

using vigra::max;
template <class T, int N>
inline DualVector<T, N> max(DualVector<T, N> const & v1, DualVector<T, N> const & v2)
{
    return v1.v > v2.v
               ? v1
               : v2;
}

template <class T, int N>
inline DualVector<T, N> max(T v1, DualVector<T, N> const & v2)
{
    return v1 > v2.v
               ? DualVector<T, N>(v1)
               : v2;
}

template <class T, int N>
inline DualVector<T, N> max(DualVector<T, N> const & v1, T v2)
{
    return v1.v > v2
               ? v1
               : DualVector<T, N>(v2);
}

template <class T, int N>
inline bool 
operator==(DualVector<T, N> const & v1, DualVector<T, N> const & v2)
{
    return v1.v == v2.v && v1.d == v2.d;
}

template <class T, int N>
inline bool 
operator!=(DualVector<T, N> const & v1, DualVector<T, N> const & v2)
{
    return v1.v != v2.v || v1.d != v2.d;
}

#define VIGRA_DUALVECTOR_RELATIONAL_OPERATORS(op) \
template <class T, int N> \
inline bool  \
operator op(DualVector<T, N> const & v1, DualVector<T, N> const & v2) \
{ \
    return v1.v op v2.v; \
} \
 \
template <class T, int N> \
inline bool  \
operator op(T v1, DualVector<T, N> const & v2) \
{ \
    return v1 op v2.v; \
} \
 \
template <class T, int N> \
inline bool  \
operator op(DualVector<T, N> const & v1, T v2) \
{ \
    return v1.v op v2; \
}

VIGRA_DUALVECTOR_RELATIONAL_OPERATORS(<)
VIGRA_DUALVECTOR_RELATIONAL_OPERATORS(<=)
VIGRA_DUALVECTOR_RELATIONAL_OPERATORS(>)
VIGRA_DUALVECTOR_RELATIONAL_OPERATORS(>=)

#undef VIGRA_DUALVECTOR_RELATIONAL_OPERATORS

template <class T, int N>
inline bool 
closeAtTolerance(DualVector<T, N> const & v1, DualVector<T, N> const & v2, 
                 T epsilon = NumericTraits<T>::epsilon())
{
    return vigra::closeAtTolerance(v1.v, v2.v, epsilon) && vigra::closeAtTolerance(v1.d, v2.d, epsilon);
}

} // namespace autodiff

} // namespace vigra

namespace std {

   /// stream output
template <class T, int N>
ostream &
operator<<(ostream & out, vigra::autodiff::DualVector<T, N> const & l)
{
    out << l.v << " " << l.d;
    return out;
}

} // namespace std

#endif // VIGRA_AUTODIFF_HXX
