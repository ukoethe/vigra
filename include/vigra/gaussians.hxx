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

#ifndef VIGRA_GAUSSIANS_HXX
#define VIGRA_GAUSSIANS_HXX

#include <cmath>
#include "config.hxx"
#include "mathutil.hxx"
#include "array_vector.hxx"
#include "error.hxx"

namespace vigra {

#if 0
/** \addtogroup MathFunctions Mathematical Functions
*/
//@{
#endif
/*! The Gaussian function and its derivatives.

    Implemented as a unary functor. Since it supports the <tt>radius()</tt> function
    it can also be used as a kernel in \ref resamplingConvolveImage().

    <b>\#include</b> \<<a href="gaussians_8hxx-source.html">vigra/gaussians.hxx</a>\><br>
    Namespace: vigra

    \ingroup MathFunctions
*/
template <class T = double>
class Gaussian
{
  public:

        /** the value type if used as a kernel in \ref resamplingConvolveImage().
        */
    typedef T            value_type;
        /** the functor's argument type
        */
    typedef T            argument_type;
        /** the functor's result type
        */
    typedef T            result_type;

        /** Create functor for the given standard deviation <tt>sigma</tt> and
            derivative order <i>n</i>. The functor then realizes the function

            \f[ f_{\sigma,n}(x)=\frac{\partial^n}{\partial x^n}
                 \frac{1}{\sqrt{2\pi}\sigma}e^{-\frac{x^2}{2\sigma^2}}
            \f]

            Precondition:
            \code
            sigma > 0.0
            \endcode
        */
    explicit Gaussian(T sigma = 1.0, unsigned int derivativeOrder = 0)
    : sigma_(sigma),
      sigma2_(T(-0.5 / sigma / sigma)),
      norm_(0.0),
      order_(derivativeOrder),
      hermitePolynomial_(derivativeOrder / 2 + 1)
    {
        vigra_precondition(sigma_ > 0.0,
            "Gaussian::Gaussian(): sigma > 0 required.");
        switch(order_)
        {
            case 1:
            case 2:
                norm_ = T(-1.0 / (VIGRA_CSTD::sqrt(2.0 * M_PI) * sq(sigma) * sigma));
                break;
            case 3:
                norm_ = T(1.0 / (VIGRA_CSTD::sqrt(2.0 * M_PI) * sq(sigma) * sq(sigma) * sigma));
                break;
            default:
                norm_ = T(1.0 / VIGRA_CSTD::sqrt(2.0 * M_PI) / sigma);
        }
        calculateHermitePolynomial();
    }

        /** Function (functor) call.
        */
    result_type operator()(argument_type x) const;

        /** Get the standard deviation of the Gaussian.
        */
    value_type sigma() const
        { return sigma_; }

        /** Get the derivative order of the Gaussian.
        */
    unsigned int derivativeOrder() const
        { return order_; }

        /** Get the required filter radius for a discrete approximation of the Gaussian.
            The radius is given as a multiple of the Gaussian's standard deviation
            (default: <tt>sigma * (3 + 1/2 * derivativeOrder()</tt> -- the second term
            accounts for the fact that the derivatives of the Gaussian become wider
            with increasing order). The result is rounded to the next higher integer.
        */
    double radius(double sigmaMultiple = 3.0) const
        { return VIGRA_CSTD::ceil(sigma_ * (sigmaMultiple + 0.5 * derivativeOrder())); }

  private:
    void calculateHermitePolynomial();
    T horner(T x) const;

    T sigma_, sigma2_, norm_;
    unsigned int order_;
    ArrayVector<T> hermitePolynomial_;
};

template <class T>
typename Gaussian<T>::result_type
Gaussian<T>::operator()(argument_type x) const
{
    T x2 = x * x;
    T g  = norm_ * VIGRA_CSTD::exp(x2 * sigma2_);
    switch(order_)
    {
        case 0:
            return detail::RequiresExplicitCast<result_type>::cast(g);
        case 1:
            return detail::RequiresExplicitCast<result_type>::cast(x * g);
        case 2:
            return detail::RequiresExplicitCast<result_type>::cast((1.0 - sq(x / sigma_)) * g);
        case 3:
            return detail::RequiresExplicitCast<result_type>::cast((3.0 - sq(x / sigma_)) * x * g);
        default:
            return order_ % 2 == 0 ?
                       detail::RequiresExplicitCast<result_type>::cast(g * horner(x2))
                     : detail::RequiresExplicitCast<result_type>::cast(x * g * horner(x2));
    }
}

template <class T>
T Gaussian<T>::horner(T x) const
{
    int i = order_ / 2;
    T res = hermitePolynomial_[i];
    for(--i; i >= 0; --i)
        res = x * res + hermitePolynomial_[i];
    return res;
}

template <class T>
void Gaussian<T>::calculateHermitePolynomial()
{
    if(order_ == 0)
    {
        hermitePolynomial_[0] = 1.0;
    }
    else if(order_ == 1)
    {
        hermitePolynomial_[0] = T(-1.0 / sigma_ / sigma_);
    }
    else
    {
        // calculate Hermite polynomial for requested derivative
        // recursively according to
        //     (0)
        //    h   (x) = 1
        //
        //     (1)
        //    h   (x) = -x / s^2
        //
        //     (n+1)                        (n)           (n-1)
        //    h     (x) = -1 / s^2 * [ x * h   (x) + n * h     (x) ]
        //
        T s2 = T(-1.0 / sigma_ / sigma_);
        ArrayVector<T> hn(3*order_+3, 0.0);
        typename ArrayVector<T>::iterator hn0 = hn.begin(),
                                          hn1 = hn0 + order_+1,
                                          hn2 = hn1 + order_+1,
                                          ht;
        hn2[0] = 1.0;
        hn1[1] = s2;
        for(unsigned int i = 2; i <= order_; ++i)
        {
            hn0[0] = s2 * (i-1) * hn2[0];
            for(unsigned int j = 1; j <= i; ++j)
                hn0[j] = s2 * (hn1[j-1] + (i-1) * hn2[j]);
            ht = hn2;
            hn2 = hn1;
            hn1 = hn0;
            hn0 = ht;
        }
        // keep only non-zero coefficients of the polynomial
        for(unsigned int i = 0; i < hermitePolynomial_.size(); ++i)
            hermitePolynomial_[i] = order_ % 2 == 0 ?
                                         hn1[2*i]
                                       : hn1[2*i+1];
    }
}


////@}

} // namespace vigra


#endif /* VIGRA_GAUSSIANS_HXX */
