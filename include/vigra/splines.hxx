/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2004 by Ullrich Koethe                  */
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

#ifndef VIGRA_SPLINES_HXX
#define VIGRA_SPLINES_HXX

#include <cmath>
#include "config.hxx"
#include "mathutil.hxx"
#include "array_vector.hxx"
#include "fixedpoint.hxx"

namespace vigra {

namespace autodiff {

template <class T, int N>
class DualVector;

} // namespace autodiff

/** \addtogroup MathFunctions Mathematical Functions
*/
//@{
/* B-Splines of arbitrary order and interpolating Catmull/Rom splines.

    <b>\#include</b> \<vigra/splines.hxx\><br>
    Namespace: vigra
*/
#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

/** Basic interface of the spline functors.

    Implements the spline functions defined by the recursion

    \f[ B_0(x) = \left\{ \begin{array}{ll}
                                  1 & -\frac{1}{2} \leq x < \frac{1}{2} \\
                                  0 & \mbox{otherwise}
                        \end{array}\right.
    \f]

    and

    \f[ B_n(x) = B_0(x) * B_{n-1}(x)
    \f]

    where * denotes convolution, and <i>n</i> is the spline order given by the template
    parameter <tt>ORDER</tt> with <tt>ORDER < 18</tt>. These spline classes can be used as
    unary and binary functors, as kernels for \ref resamplingConvolveImage(),
    and as arguments for \ref vigra::SplineImageView. Note that the spline order
    is given as a template argument.

    <b>\#include</b> \<vigra/splines.hxx\><br>
    Namespace: vigra
*/
template <int ORDER, class T = double>
class BSplineBase
{
  public:

    static_assert (ORDER < 18 , "BSpline: ORDER must be less than 18." );

        /** the value type if used as a kernel in \ref resamplingConvolveImage().
        */
    typedef T            value_type;
        /** the functor's unary argument type
        */
    typedef T            argument_type;
        /** the functor's first binary argument type
        */
    typedef T            first_argument_type;
        /** the functor's second binary argument type
        */
    typedef unsigned int second_argument_type;
        /** the functor's result type (unary and binary)
        */
    typedef T            result_type;
        /** the spline order
        */
    enum StaticOrder { order = ORDER };

        /** Create functor for given derivative of the spline. The spline's order
            is specified spline by the template argument <TT>ORDER</tt>.
        */
    explicit BSplineBase(unsigned int derivativeOrder = 0)
    : s1_(derivativeOrder)
    {}

        /** Unary function call.
            Returns the value of the spline with the derivative order given in the
            constructor. Note that only derivatives up to <tt>ORDER-1</tt> are
            continuous, and derivatives above <tt>ORDER+1</tt> are zero.
        */
    result_type operator()(argument_type x) const
    {
        return exec(x, derivativeOrder());
    }

        /** Binary function call.
            The given derivative order is added to the derivative order
            specified in the constructor. Note that only derivatives up to <tt>ORDER-1</tt> are
            continuous, and derivatives above <tt>ORDER+1</tt> are zero.
        */
    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
         return exec(x, derivativeOrder() + derivative_order);
    }

        /** Index operator. Same as unary function call.
        */
    value_type operator[](value_type x) const
        { return operator()(x); }

        /** Get the required filter radius for a discrete approximation of the
            spline. Always equal to <tt>(ORDER + 1) / 2.0</tt>.
        */
    double radius() const
        { return (ORDER + 1) * 0.5; }

        /** Get the derivative order of the Gaussian.
        */
    unsigned int derivativeOrder() const
        { return s1_.derivativeOrder(); }

        /** Get the prefilter coefficients required for interpolation.
            To interpolate with a B-spline, \ref resamplingConvolveImage()
            can be used. However, the image to be interpolated must be
            pre-filtered using \ref recursiveFilterX() and \ref recursiveFilterY()
            with the filter coefficients given by this function. The length of the array
            corresponds to how many times the above recursive filtering
            has to be applied (zero length means no prefiltering necessary).
        */
    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

    typedef ArrayVector<ArrayVector<T> > WeightMatrix;

        /** Get the coefficients to transform spline coefficients into
            the coefficients of the corresponding polynomial.
            Currently internally used in SplineImageView; needs more
            documentation ???
        */
    static WeightMatrix const & weights()
    {
        return weightMatrix_;
    }

    static ArrayVector<double> getPrefilterCoefficients();

  protected:
    result_type exec(first_argument_type x, second_argument_type derivative_order) const;

        // factory function for the weight matrix
    static WeightMatrix calculateWeightMatrix();

    BSplineBase<ORDER-1, T> s1_;
    static ArrayVector<double> prefilterCoefficients_;
    static WeightMatrix weightMatrix_;
};

template <int ORDER, class T>
ArrayVector<double> BSplineBase<ORDER, T>::prefilterCoefficients_(getPrefilterCoefficients());

template <int ORDER, class T>
typename BSplineBase<ORDER, T>::WeightMatrix BSplineBase<ORDER, T>::weightMatrix_(calculateWeightMatrix());

template <int ORDER, class T>
typename BSplineBase<ORDER, T>::result_type
BSplineBase<ORDER, T>::exec(first_argument_type x, second_argument_type derivative_order) const
{
    if(derivative_order == 0)
    {
        T n12 = (ORDER + 1.0) / 2.0;
        return ((n12 + x) * s1_(x + 0.5) + (n12 - x) * s1_(x - 0.5)) / ORDER;
    }
    else
    {
        --derivative_order;
        return s1_(x + 0.5, derivative_order) - s1_(x - 0.5, derivative_order);
    }
}

template <int ORDER, class T>
ArrayVector<double>
BSplineBase<ORDER, T>::getPrefilterCoefficients()
{
   static const double coeffs[18][8] = {
        { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.17157287525380971, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.26794919243112281, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.36134122590022018, -0.01372542929733912, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.43057534709997379, -0.04309628820326465, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.48829458930304398, -0.081679271076237972, -0.0014141518083258175, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.53528043079643672, -0.1225546151923274, -0.0091486948096082786, 0.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.57468690924876631, -0.16303526929728085, -0.023632294694844857, -0.00015382131064169087, 0.0, 0.0, 0.0, 0.0 },
        { -0.60799738916862989, -0.20175052019315337, -0.043222608540481752, -0.0021213069031808186, 0.0, 0.0, 0.0, 0.0 },
        { -0.63655066396942439, -0.2381827983775629, -0.065727033228308585, -0.0075281946755486927, -1.6982762823274658e-5, 0.0, 0.0, 0.0 },
        { -0.66126606890072925, -0.27218034929478602, -0.089759599793713341, -0.016669627366234657, -0.00051055753444650205, 0.0, 0.0, 0.0 },
        { -0.68286488419772362, -0.30378079328825425, -0.11435052002713579, -0.028836190198663809, -0.0025161662172613372, -1.8833056450639017e-6, 0.0, 0.0 },
        { -0.70189425181681642, -0.33310723293062366, -0.13890111319431958, -0.043213866740363663, -0.0067380314152449142, -0.00012510011321441875, 0.0, 0.0 },
        { -0.71878378723997516, -0.3603190719169625, -0.1630335147992984, -0.059089482194831018, -0.013246756734847919, -0.00086402404095333829, -2.0913096775275374e-7, 0.0 },
        { -0.73387257168487741, -0.38558573427843323, -0.18652010845096478, -0.075907592047668185, -0.02175206579654047, -0.0028011514820764556, -3.093568045147443e-5, 0.0 },
        { -0.747432387772212103, -0.409073604757528353, -0.29228719338953817, -9.32547189803214355e-2 -3.18677061204386616e-2, -6.25840678512839046e-3, -3.01565363306955866e-4, -2.32324863642097035e-8 },
        { -0.75968322407189071, -0.43093965318039579, -0.23108984359927232, -0.1108289933162471, -0.043213911456684129, -0.011258183689471605, -0.0011859331251521767, -7.6875625812546846e-6 }
    };
    return ArrayVector<double>(coeffs[ORDER], coeffs[ORDER]+ORDER/2);
}

template <int ORDER, class T>
typename BSplineBase<ORDER, T>::WeightMatrix
BSplineBase<ORDER, T>::calculateWeightMatrix()
{
    WeightMatrix res(ORDER+1, ArrayVector<T>(ORDER+1));
    double faculty = 1.0;
    for(int d = 0; d <= ORDER; ++d)
    {
        if(d > 1)
            faculty *= d;
        double x = ORDER / 2; // (note: integer division)
        BSplineBase spline;
        for(int i = 0; i <= ORDER; ++i, --x)
            res[d][i] = spline(x, d) / faculty;
    }
    return res;
}

/********************************************************/
/*                                                      */
/*                     BSpline<N, T>                    */
/*                                                      */
/********************************************************/

/** Spline functors for arbitrary orders.

    Provides the interface of \ref vigra::BSplineBase with a more convenient
    name -- see there for more documentation.
*/
template <int ORDER, class T = double>
class BSpline
: public BSplineBase<ORDER, T>
{
  public:
        /** Constructor forwarded to the base class constructor..
        */
    explicit BSpline(unsigned int derivativeOrder = 0)
    : BSplineBase<ORDER, T>(derivativeOrder)
    {}
};

/********************************************************/
/*                                                      */
/*                     BSpline<0, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSplineBase<0, T>
{
  public:

    typedef T            value_type;
    typedef T            argument_type;
    typedef T            first_argument_type;
    typedef unsigned int second_argument_type;
    typedef T            result_type;
    enum StaticOrder { order = 0 };

    explicit BSplineBase(unsigned int derivativeOrder = 0)
    : derivativeOrder_(derivativeOrder)
    {}

    result_type operator()(argument_type x) const
    {
         return exec(x, derivativeOrder_);
    }

    template <unsigned int IntBits, unsigned int FracBits>
    FixedPoint<IntBits, FracBits> operator()(FixedPoint<IntBits, FracBits> x) const
    {
        typedef FixedPoint<IntBits, FracBits> Value;
        return x.value < Value::ONE_HALF && -Value::ONE_HALF <= x.value
                   ? Value(Value::ONE, FPNoShift)
                   : Value(0, FPNoShift);
    }

    template <class U, int N>
    autodiff::DualVector<U, N> operator()(autodiff::DualVector<U, N> const & x) const
    {
        return x < 0.5 && -0.5 <= x
                   ? autodiff::DualVector<U, N>(1.0)
                   : autodiff::DualVector<U, N>(0.0);
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
         return exec(x, derivativeOrder_ + derivative_order);
    }

    value_type operator[](value_type x) const
        { return operator()(x); }

    double radius() const
        { return 0.5; }

    unsigned int derivativeOrder() const
        { return derivativeOrder_; }

    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

    typedef T WeightMatrix[1][1];

    static WeightMatrix const & weights()
    {
        return weightMatrix_;
    }

  protected:
    result_type exec(first_argument_type x, second_argument_type derivative_order) const
    {
        if(derivative_order == 0)
            return x < 0.5 && -0.5 <= x ?
                     1.0
                   : 0.0;
        else
            return 0.0;
    }

    unsigned int derivativeOrder_;
    static ArrayVector<double> prefilterCoefficients_;
    static WeightMatrix weightMatrix_;
};

template <class T>
ArrayVector<double> BSplineBase<0, T>::prefilterCoefficients_;

template <class T>
typename BSplineBase<0, T>::WeightMatrix BSplineBase<0, T>::weightMatrix_ = {{ 1.0 }};

/********************************************************/
/*                                                      */
/*                     BSpline<1, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<1, T>
{
  public:

    typedef T            value_type;
    typedef T            argument_type;
    typedef T            first_argument_type;
    typedef unsigned int second_argument_type;
    typedef T            result_type;
    enum  StaticOrder { order = 1 };

    explicit BSpline(unsigned int derivativeOrder = 0)
    : derivativeOrder_(derivativeOrder)
    {}

    result_type operator()(argument_type x) const
    {
        return exec(x, derivativeOrder_);
    }

    template <unsigned int IntBits, unsigned int FracBits>
    FixedPoint<IntBits, FracBits> operator()(FixedPoint<IntBits, FracBits> x) const
    {
        typedef FixedPoint<IntBits, FracBits> Value;
        int v = abs(x.value);
        return v < Value::ONE ?
                Value(Value::ONE - v, FPNoShift)
                : Value(0, FPNoShift);
    }

    template <class U, int N>
    autodiff::DualVector<U, N> operator()(autodiff::DualVector<U, N> x) const
    {
        x = abs(x);
        return x < 1.0
                    ? 1.0 - x
                    : autodiff::DualVector<U, N>(0.0);
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
         return exec(x, derivativeOrder_ + derivative_order);
    }

    value_type operator[](value_type x) const
        { return operator()(x); }

    double radius() const
        { return 1.0; }

    unsigned int derivativeOrder() const
        { return derivativeOrder_; }

    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

    typedef T WeightMatrix[2][2];

    static WeightMatrix const & weights()
    {
        return weightMatrix_;
    }

  protected:
    T exec(T x, unsigned int derivative_order) const;

    unsigned int derivativeOrder_;
    static ArrayVector<double> prefilterCoefficients_;
    static WeightMatrix weightMatrix_;
};

template <class T>
ArrayVector<double> BSpline<1, T>::prefilterCoefficients_;

template <class T>
typename BSpline<1, T>::WeightMatrix BSpline<1, T>::weightMatrix_ = {{ 1.0, 0.0}, {-1.0, 1.0}};

template <class T>
T BSpline<1, T>::exec(T x, unsigned int derivative_order) const
{
    switch(derivative_order)
    {
        case 0:
        {
            x = VIGRA_CSTD::fabs(x);
            return x < 1.0 ?
                    1.0 - x
                    : 0.0;
        }
        case 1:
        {
            return x < 0.0 ?
                     -1.0 <= x ?
                          1.0
                     : 0.0
                   : x < 1.0 ?
                       -1.0
                     : 0.0;
        }
        default:
            return 0.0;
    }
}

/********************************************************/
/*                                                      */
/*                     BSpline<2, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<2, T>
{
  public:

    typedef T            value_type;
    typedef T            argument_type;
    typedef T            first_argument_type;
    typedef unsigned int second_argument_type;
    typedef T            result_type;
    enum StaticOrder { order = 2 };

    explicit BSpline(unsigned int derivativeOrder = 0)
    : derivativeOrder_(derivativeOrder)
    {}

    result_type operator()(argument_type x) const
    {
        return exec(x, derivativeOrder_);
    }

    template <unsigned int IntBits, unsigned int FracBits>
    FixedPoint<IntBits, FracBits> operator()(FixedPoint<IntBits, FracBits> x) const
    {
        typedef FixedPoint<IntBits, FracBits> Value;
        enum { ONE_HALF = Value::ONE_HALF, THREE_HALVES = ONE_HALF * 3, THREE_QUARTERS = THREE_HALVES / 2,
               PREMULTIPLY_SHIFT1 = FracBits <= 16 ? 0 : FracBits - 16,
               PREMULTIPLY_SHIFT2 = FracBits - 1 <= 16 ? 0 : FracBits - 17,
               POSTMULTIPLY_SHIFT1 = FracBits - 2*PREMULTIPLY_SHIFT1,
               POSTMULTIPLY_SHIFT2 = FracBits - 2*PREMULTIPLY_SHIFT2  };
        int v = abs(x.value);
        return v == ONE_HALF
                   ? Value(ONE_HALF, FPNoShift)
                   : v <= ONE_HALF
                       ? Value(THREE_QUARTERS -
                               (int)(sq((unsigned)v >> PREMULTIPLY_SHIFT2) >> POSTMULTIPLY_SHIFT2), FPNoShift)
                       : v < THREE_HALVES
                            ? Value((int)(sq((unsigned)(THREE_HALVES-v) >> PREMULTIPLY_SHIFT1) >> (POSTMULTIPLY_SHIFT1 + 1)), FPNoShift)
                            : Value(0, FPNoShift);
    }

    template <class U, int N>
    autodiff::DualVector<U, N> operator()(autodiff::DualVector<U, N> x) const
    {
        x = abs(x);
        return x < 0.5
                   ? 0.75 - x*x
                   : x < 1.5
                        ? 0.5 * sq(1.5 - x)
                        : autodiff::DualVector<U, N>(0.0);
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
         return exec(x, derivativeOrder_ + derivative_order);
    }

    result_type dx(argument_type x) const
        { return operator()(x, 1); }

    value_type operator[](value_type x) const
        { return operator()(x); }

    double radius() const
        { return 1.5; }

    unsigned int derivativeOrder() const
        { return derivativeOrder_; }

    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

    typedef T WeightMatrix[3][3];

    static WeightMatrix const & weights()
    {
        return weightMatrix_;
    }

  protected:
    result_type exec(first_argument_type x, second_argument_type derivative_order) const;

    unsigned int derivativeOrder_;
    static ArrayVector<double> prefilterCoefficients_;
    static WeightMatrix weightMatrix_;
};

template <class T>
ArrayVector<double> BSpline<2, T>::prefilterCoefficients_(1, 2.0*M_SQRT2 - 3.0);

template <class T>
typename BSpline<2, T>::WeightMatrix BSpline<2, T>::weightMatrix_ =
                           {{ 0.125, 0.75, 0.125},
                            {-0.5, 0.0, 0.5},
                            { 0.5, -1.0, 0.5}};

template <class T>
typename BSpline<2, T>::result_type
BSpline<2, T>::exec(first_argument_type x, second_argument_type derivative_order) const
{
    switch(derivative_order)
    {
        case 0:
        {
            x = VIGRA_CSTD::fabs(x);
            return x < 0.5 ?
                    0.75 - x*x
                    : x < 1.5 ?
                        0.5 * sq(1.5 - x)
                    : 0.0;
        }
        case 1:
        {
            return x >= -0.5 ?
                     x <= 0.5 ?
                       -2.0 * x
                     : x < 1.5 ?
                         x - 1.5
                       : 0.0
                   : x > -1.5 ?
                       x + 1.5
                     : 0.0;
        }
        case 2:
        {
            return x >= -0.5 ?
                     x < 0.5 ?
                         -2.0
                     : x < 1.5 ?
                         1.0
                       : 0.0
                   : x >= -1.5 ?
                       1.0
                     : 0.0;
        }
        default:
            return 0.0;
    }
}

/********************************************************/
/*                                                      */
/*                     BSpline<3, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<3, T>
{
  public:

    typedef T            value_type;
    typedef T            argument_type;
    typedef T            first_argument_type;
    typedef unsigned int second_argument_type;
    typedef T            result_type;
    enum StaticOrder { order = 3 };

    explicit BSpline(unsigned int derivativeOrder = 0)
    : derivativeOrder_(derivativeOrder)
    {}

    result_type operator()(argument_type x) const
    {
        return exec(x, derivativeOrder_);
    }

    template <unsigned int IntBits, unsigned int FracBits>
    FixedPoint<IntBits, FracBits> operator()(FixedPoint<IntBits, FracBits> x) const
    {
        typedef FixedPoint<IntBits, FracBits> Value;
        enum { ONE = Value::ONE, TWO = 2 * ONE, TWO_THIRDS = TWO / 3, ONE_SIXTH = ONE / 6,
               PREMULTIPLY_SHIFT = FracBits <= 16 ? 0 : FracBits - 16,
               POSTMULTIPLY_SHIFT = FracBits - 2*PREMULTIPLY_SHIFT };
        int v = abs(x.value);
        return v == ONE
                   ? Value(ONE_SIXTH, FPNoShift)
                   : v < ONE
                       ? Value(TWO_THIRDS +
                               (((int)(sq((unsigned)v >> PREMULTIPLY_SHIFT) >> (POSTMULTIPLY_SHIFT + PREMULTIPLY_SHIFT))
                                       * (((v >> 1) - ONE) >> PREMULTIPLY_SHIFT)) >> POSTMULTIPLY_SHIFT), FPNoShift)
                       : v < TWO
                            ? Value((int)((sq((unsigned)(TWO-v) >> PREMULTIPLY_SHIFT) >> (POSTMULTIPLY_SHIFT + PREMULTIPLY_SHIFT))
                                      * ((unsigned)(TWO-v) >> PREMULTIPLY_SHIFT) / 6) >> POSTMULTIPLY_SHIFT, FPNoShift)
                            : Value(0, FPNoShift);
    }

    template <class U, int N>
    autodiff::DualVector<U, N> operator()(autodiff::DualVector<U, N> x) const
    {
        x = abs(x);
        if(x < 1.0)
        {
            return 2.0/3.0 + x*x*(-1.0 + 0.5*x);
        }
        else if(x < 2.0)
        {
            x = 2.0 - x;
            return x*x*x/6.0;
        }
        else
            return autodiff::DualVector<U, N>(0.0);
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
         return exec(x, derivativeOrder_ + derivative_order);
    }

    result_type dx(argument_type x) const
        { return operator()(x, 1); }

    result_type dxx(argument_type x) const
        { return operator()(x, 2); }

    value_type operator[](value_type x) const
        { return operator()(x); }

    double radius() const
        { return 2.0; }

    unsigned int derivativeOrder() const
        { return derivativeOrder_; }

    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

    typedef T WeightMatrix[4][4];

    static WeightMatrix const & weights()
    {
        return weightMatrix_;
    }

  protected:
    result_type exec(first_argument_type x, second_argument_type derivative_order) const;

    unsigned int derivativeOrder_;
    static ArrayVector<double> prefilterCoefficients_;
    static WeightMatrix weightMatrix_;
};

template <class T>
ArrayVector<double> BSpline<3, T>::prefilterCoefficients_(1, VIGRA_CSTD::sqrt(3.0) - 2.0);

template <class T>
typename BSpline<3, T>::WeightMatrix BSpline<3, T>::weightMatrix_ =
                           {{ 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0, 0.0},
                            {-0.5, 0.0, 0.5, 0.0},
                            { 0.5, -1.0, 0.5, 0.0},
                            {-1.0 / 6.0, 0.5, -0.5, 1.0 / 6.0}};

template <class T>
typename BSpline<3, T>::result_type
BSpline<3, T>::exec(first_argument_type x, second_argument_type derivative_order) const
{
    switch(derivative_order)
    {
        case 0:
        {
            x = VIGRA_CSTD::fabs(x);
            if(x < 1.0)
            {
                return 2.0/3.0 + x*x*(-1.0 + 0.5*x);
            }
            else if(x < 2.0)
            {
                x = 2.0 - x;
                return x*x*x/6.0;
            }
            else
                return 0.0;
        }
        case 1:
        {
            double s = x < 0.0 ?
                         -1.0
                       :  1.0;
            x = VIGRA_CSTD::fabs(x);
            return x < 1.0 ?
                     s*x*(-2.0 + 1.5*x)
                   : x < 2.0 ?
                       -0.5*s*sq(2.0 - x)
                     : 0.0;
        }
        case 2:
        {
            x = VIGRA_CSTD::fabs(x);
            return x < 1.0 ?
                     3.0*x - 2.0
                   : x < 2.0 ?
                       2.0 - x
                     : 0.0;
        }
        case 3:
        {
            return x < 0.0 ?
                     x < -1.0 ?
                       x < -2.0 ?
                         0.0
                       : 1.0
                     : -3.0
                   : x < 1.0 ?
                       3.0
                     : x < 2.0 ?
                         -1.0
                       : 0.0;
        }
        default:
            return 0.0;
    }
}

typedef BSpline<3, double> CubicBSplineKernel;

/********************************************************/
/*                                                      */
/*                     BSpline<4, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<4, T>
{
  public:

    typedef T            value_type;
    typedef T            argument_type;
    typedef T            first_argument_type;
    typedef unsigned int second_argument_type;
    typedef T            result_type;
    enum StaticOrder { order = 4 };

    explicit BSpline(unsigned int derivativeOrder = 0)
    : derivativeOrder_(derivativeOrder)
    {}

    result_type operator()(argument_type x) const
    {
        return exec(x, derivativeOrder_);
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
         return exec(x, derivativeOrder_ + derivative_order);
    }

    template <class U, int N>
    autodiff::DualVector<U, N> operator()(autodiff::DualVector<U, N> x) const
    {
        x = abs(x);
        if(x <= 0.5)
        {
            return 115.0/192.0 + x*x*(-0.625 + x*x*0.25);
        }
        else if(x < 1.5)
        {
            return (55.0/16.0 + x*(1.25 + x*(-7.5 + x*(5.0 - x)))) / 6.0;
        }
        else if(x < 2.5)
        {
            x = 2.5 - x;
            return sq(x*x) / 24.0;
        }
        else
            return autodiff::DualVector<U, N>(0.0);
    }

    result_type dx(argument_type x) const
        { return operator()(x, 1); }

    result_type dxx(argument_type x) const
        { return operator()(x, 2); }

    result_type dx3(argument_type x) const
        { return operator()(x, 3); }

    value_type operator[](value_type x) const
        { return operator()(x); }

    double radius() const
        { return 2.5; }

    unsigned int derivativeOrder() const
        { return derivativeOrder_; }

    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

    typedef T WeightMatrix[5][5];

    static WeightMatrix const & weights()
    {
        return weightMatrix_;
    }

  protected:
    result_type exec(first_argument_type x, second_argument_type derivative_order) const;

    unsigned int derivativeOrder_;
    static ArrayVector<double> prefilterCoefficients_;
    static WeightMatrix weightMatrix_;
};

template <class T>
ArrayVector<double> BSpline<4, T>::prefilterCoefficients_(BSplineBase<4, T>::getPrefilterCoefficients());

template <class T>
typename BSpline<4, T>::WeightMatrix BSpline<4, T>::weightMatrix_ =
                           {{ 1.0/384.0, 19.0/96.0, 115.0/192.0, 19.0/96.0, 1.0/384.0},
                            {-1.0/48.0, -11.0/24.0, 0.0, 11.0/24.0, 1.0/48.0},
                            { 1.0/16.0, 1.0/4.0, -5.0/8.0, 1.0/4.0, 1.0/16.0},
                            {-1.0/12.0, 1.0/6.0, 0.0, -1.0/6.0, 1.0/12.0},
                            { 1.0/24.0, -1.0/6.0, 0.25, -1.0/6.0, 1.0/24.0}};

template <class T>
typename BSpline<4, T>::result_type
BSpline<4, T>::exec(first_argument_type x, second_argument_type derivative_order) const
{
    switch(derivative_order)
    {
        case 0:
        {
            x = VIGRA_CSTD::fabs(x);
            if(x <= 0.5)
            {
                return 115.0/192.0 + x*x*(-0.625 + x*x*0.25);
            }
            else if(x < 1.5)
            {
                return (55.0/16.0 + x*(1.25 + x*(-7.5 + x*(5.0 - x)))) / 6.0;
            }
            else if(x < 2.5)
            {
                x = 2.5 - x;
                return sq(x*x) / 24.0;
            }
            else
                return 0.0;
        }
        case 1:
        {
            double s = x < 0.0 ?
                          -1.0 :
                           1.0;
            x = VIGRA_CSTD::fabs(x);
            if(x <= 0.5)
            {
                return s*x*(-1.25 + x*x);
            }
            else if(x < 1.5)
            {
                return s*(5.0 + x*(-60.0 + x*(60.0 - 16.0*x))) / 24.0;
            }
            else if(x < 2.5)
            {
                x = 2.5 - x;
                return s*x*x*x / -6.0;
            }
            else
                return 0.0;
        }
        case 2:
        {
            x = VIGRA_CSTD::fabs(x);
            if(x <= 0.5)
            {
                return -1.25 + 3.0*x*x;
            }
            else if(x < 1.5)
            {
                return -2.5 + x*(5.0 - 2.0*x);
            }
            else if(x < 2.5)
            {
                x = 2.5 - x;
                return x*x / 2.0;
            }
            else
                return 0.0;
        }
        case 3:
        {
            double s = x < 0.0 ?
                          -1.0 :
                           1.0;
            x = VIGRA_CSTD::fabs(x);
            if(x <= 0.5)
            {
                return s*x*6.0;
            }
            else if(x < 1.5)
            {
                return s*(5.0 - 4.0*x);
            }
            else if(x < 2.5)
            {
                return s*(x - 2.5);
            }
            else
                return 0.0;
        }
        case 4:
        {
            return x < 0.0
                     ? x < -2.5
                         ? 0.0
                         : x < -1.5
                             ? 1.0
                             : x < -0.5
                                 ? -4.0
                                 : 6.0
                     : x < 0.5
                         ? 6.0
                         : x < 1.5
                             ? -4.0
                             : x < 2.5
                                 ? 1.0
                                 : 0.0;
        }
        default:
            return 0.0;
    }
}

typedef BSpline<4, double> QuarticBSplineKernel;

/********************************************************/
/*                                                      */
/*                     BSpline<5, T>                    */
/*                                                      */
/********************************************************/

template <class T>
class BSpline<5, T>
{
  public:

    typedef T            value_type;
    typedef T            argument_type;
    typedef T            first_argument_type;
    typedef unsigned int second_argument_type;
    typedef T            result_type;
    enum StaticOrder { order = 5 };

    explicit BSpline(unsigned int derivativeOrder = 0)
    : derivativeOrder_(derivativeOrder)
    {}

    result_type operator()(argument_type x) const
    {
        return exec(x, derivativeOrder_);
    }

    result_type operator()(first_argument_type x, second_argument_type derivative_order) const
    {
         return exec(x, derivativeOrder_ + derivative_order);
    }

    template <class U, int N>
    autodiff::DualVector<U, N> operator()(autodiff::DualVector<U, N> x) const
    {
        x = abs(x);
        if(x <= 1.0)
        {
            return 0.55 + x*x*(-0.5 + x*x*(0.25 - x/12.0));
        }
        else if(x < 2.0)
        {
            return 17.0/40.0 + x*(0.625 + x*(-1.75 + x*(1.25 + x*(-0.375 + x/24.0))));
        }
        else if(x < 3.0)
        {
            x = 3.0 - x;
            return x*sq(x*x) / 120.0;
        }
        else
            return autodiff::DualVector<U, N>(0.0);
    }

    result_type dx(argument_type x) const
        { return operator()(x, 1); }

    result_type dxx(argument_type x) const
        { return operator()(x, 2); }

    result_type dx3(argument_type x) const
        { return operator()(x, 3); }

    result_type dx4(argument_type x) const
        { return operator()(x, 4); }

    value_type operator[](value_type x) const
        { return operator()(x); }

    double radius() const
        { return 3.0; }

    unsigned int derivativeOrder() const
        { return derivativeOrder_; }

    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

    typedef T WeightMatrix[6][6];

    static WeightMatrix const & weights()
    {
        return weightMatrix_;
    }

  protected:
    result_type exec(first_argument_type x, second_argument_type derivative_order) const;

    unsigned int derivativeOrder_;
    static ArrayVector<double> prefilterCoefficients_;
    static WeightMatrix weightMatrix_;
};

template <class T>
ArrayVector<double> BSpline<5, T>::prefilterCoefficients_(BSplineBase<5, T>::getPrefilterCoefficients());

template <class T>
typename BSpline<5, T>::WeightMatrix BSpline<5, T>::weightMatrix_ =
                           {{ 1.0/120.0, 13.0/60.0, 11.0/20.0, 13.0/60.0, 1.0/120.0, 0.0},
                            {-1.0/24.0, -5.0/12.0, 0.0, 5.0/12.0, 1.0/24.0, 0.0},
                            { 1.0/12.0, 1.0/6.0, -0.5, 1.0/6.0, 1.0/12.0, 0.0},
                            {-1.0/12.0, 1.0/6.0, 0.0, -1.0/6.0, 1.0/12.0, 0.0},
                            { 1.0/24.0, -1.0/6.0, 0.25, -1.0/6.0, 1.0/24.0, 0.0},
                            {-1.0/120.0, 1.0/24.0, -1.0/12.0, 1.0/12.0, -1.0/24.0, 1.0/120.0}};

template <class T>
typename BSpline<5, T>::result_type
BSpline<5, T>::exec(first_argument_type x, second_argument_type derivative_order) const
{
    switch(derivative_order)
    {
        case 0:
        {
            x = VIGRA_CSTD::fabs(x);
            if(x <= 1.0)
            {
                return 0.55 + x*x*(-0.5 + x*x*(0.25 - x/12.0));
            }
            else if(x < 2.0)
            {
                return 17.0/40.0 + x*(0.625 + x*(-1.75 + x*(1.25 + x*(-0.375 + x/24.0))));
            }
            else if(x < 3.0)
            {
                x = 3.0 - x;
                return x*sq(x*x) / 120.0;
            }
            else
                return 0.0;
        }
        case 1:
        {
            double s = x < 0.0 ?
                          -1.0 :
                           1.0;
            x = VIGRA_CSTD::fabs(x);
            if(x <= 1.0)
            {
                return s*x*(-1.0 + x*x*(1.0 - 5.0/12.0*x));
            }
            else if(x < 2.0)
            {
                return s*(0.625 + x*(-3.5 + x*(3.75 + x*(-1.5 + 5.0/24.0*x))));
            }
            else if(x < 3.0)
            {
                x = 3.0 - x;
                return s*sq(x*x) / -24.0;
            }
            else
                return 0.0;
        }
        case 2:
        {
            x = VIGRA_CSTD::fabs(x);
            if(x <= 1.0)
            {
                return -1.0 + x*x*(3.0 -5.0/3.0*x);
            }
            else if(x < 2.0)
            {
                return -3.5 + x*(7.5 + x*(-4.5 + 5.0/6.0*x));
            }
            else if(x < 3.0)
            {
                x = 3.0 - x;
                return x*x*x / 6.0;
            }
            else
                return 0.0;
        }
        case 3:
        {
            double s = x < 0.0 ?
                          -1.0 :
                           1.0;
            x = VIGRA_CSTD::fabs(x);
            if(x <= 1.0)
            {
                return s*x*(6.0 - 5.0*x);
            }
            else if(x < 2.0)
            {
                return s*(7.5 + x*(-9.0 + 2.5*x));
            }
            else if(x < 3.0)
            {
                x = 3.0 - x;
                return -0.5*s*x*x;
            }
            else
                return 0.0;
        }
        case 4:
        {
            x = VIGRA_CSTD::fabs(x);
            if(x <= 1.0)
            {
                return 6.0 - 10.0*x;
            }
            else if(x < 2.0)
            {
                return -9.0 + 5.0*x;
            }
            else if(x < 3.0)
            {
                return 3.0 - x;
            }
            else
                return 0.0;
        }
        case 5:
        {
            return x < 0.0 ?
                     x < -2.0 ?
                       x < -3.0 ?
                         0.0
                       : 1.0
                     : x < -1.0 ?
                         -5.0
                       : 10.0
                   : x < 2.0 ?
                       x < 1.0 ?
                         -10.0
                       : 5.0
                     : x < 3.0 ?
                         -1.0
                       : 0.0;
        }
        default:
            return 0.0;
    }
}

typedef BSpline<5, double> QuinticBSplineKernel;

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

/********************************************************/
/*                                                      */
/*                      CatmullRomSpline                */
/*                                                      */
/********************************************************/

/** Interpolating 3-rd order splines.

    Implements the Catmull/Rom cardinal function

    \f[ f(x) = \left\{ \begin{array}{ll}
                                  \frac{3}{2}x^3 - \frac{5}{2}x^2 + 1 & |x| \leq 1 \\
                                  -\frac{1}{2}x^3 + \frac{5}{2}x^2 -4x + 2 & |x| \leq 2 \\
                                  0 & \mbox{otherwise}
                        \end{array}\right.
    \f]

    It can be used as a functor, and as a kernel for
    \ref resamplingConvolveImage() to create a differentiable interpolant
    of an image. However, it should be noted that a twice differentiable
    interpolant can be created with only slightly more effort by recursive
    prefiltering followed by convolution with a 3rd order B-spline.

    <b>\#include</b> \<vigra/splines.hxx\><br>
    Namespace: vigra
*/
template <class T = double>
class CatmullRomSpline
{
public:
        /** the kernel's value type
        */
    typedef T value_type;
        /** the unary functor's argument type
        */
    typedef T argument_type;
        /** the unary functor's result type
        */
    typedef T result_type;
        /** the splines polynomial order
        */
    enum StaticOrder { order = 3 };

        /** function (functor) call
        */
    result_type operator()(argument_type x) const;

        /** index operator -- same as operator()
        */
    T operator[] (T x) const
        { return operator()(x); }

        /** Radius of the function's support.
            Needed for  \ref resamplingConvolveImage(), always 2.
        */
    int radius() const
        {return 2;}

        /** Derivative order of the function: always 0.
        */
    unsigned int derivativeOrder() const
        { return 0; }

        /** Prefilter coefficients for compatibility with \ref vigra::BSpline.
            (array has zero length, since prefiltering is not necessary).
        */
    ArrayVector<double> const & prefilterCoefficients() const
    {
        return prefilterCoefficients_;
    }

  protected:
    static ArrayVector<double> prefilterCoefficients_;
};

template <class T>
ArrayVector<double> CatmullRomSpline<T>::prefilterCoefficients_;

template <class T>
typename CatmullRomSpline<T>::result_type
CatmullRomSpline<T>::operator()(argument_type x) const
{
    x = VIGRA_CSTD::fabs(x);
    if (x <= 1.0)
    {
        return 1.0 + x * x * (-2.5 + 1.5 * x);
    }
    else if (x >= 2.0)
    {
        return 0.0;
    }
    else
    {
        return 2.0 + x * (-4.0 + x * (2.5 - 0.5 * x));
    }
}


//@}

} // namespace vigra


#endif /* VIGRA_SPLINES_HXX */
