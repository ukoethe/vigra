/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2004 by Ullrich Koethe                  */
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

#ifndef VIGRA_GAUSSIANS_HXX
#define VIGRA_GAUSSIANS_HXX

#include <cmath>
#include "vigra/config.hxx"
#include "vigra/mathutil.hxx"
#include "vigra/polynomial.hxx"

namespace vigra {

/** \addtogroup MathFunctions Mathematical Functions
*/
//@{
/*! 

    <b>\#include</b> "<a href="gaussians_8hxx-source.html">vigra/gaussians.hxx</a>"<br>
    Namespace: vigra
*/
template <class T = double>
class Gaussian
{
  public:
  
    typedef T            value_type;  
    typedef T            argument_type;  
    typedef T            result_type; 
    
    explicit Gaussian(T sigma = 1.0, unsigned int derivativeOrder = 0)
    : sigma_(sigma),
      sigma2_(-0.5 / sigma / sigma),
      norm_(derivativeOrder == 1 ?
               -1.0 / (VIGRA_CSTD::sqrt(2.0 * M_PI) * sigma * sigma * sigma)
             :  1.0 / VIGRA_CSTD::sqrt(2.0 * M_PI) / sigma),
      order_(derivativeOrder),
      hermitePolynomial_(derivativeOrder / 2)
    {
        vigra_precondition(sigma_ > 0.0,
            "Gaussian::Gaussian(): sigma > 0 required.");
        calculateHermitePolynomial();
    }

    result_type operator()(argument_type x) const
    {
        T x2 = x * x;
        if(order_ == 0)
            return norm_ * VIGRA_CSTD::exp(x2 * sigma2_);
        else if(order_ == 1)
            return x * norm_ * VIGRA_CSTD::exp(x2 * sigma2_);
        else if(order_ % 2 == 0)
            return hermitePolynomial_(x2) * norm_ * VIGRA_CSTD::exp(x2 * sigma2_);
        else
            return x * hermitePolynomial_(x2) * norm_ * VIGRA_CSTD::exp(x2 * sigma2_);
    }

    value_type sigma() const
        { return sigma_; }
    
    unsigned int derivativeOrder() const
        { return order_; }
    
    double radius(double sigmaMultiple = 3.0) const
        { return sigmaMultiple * sigma_ + 0.5 * derivativeOrder(); }

  private:
    void calculateHermitePolynomial();
    
    T sigma_, sigma2_, norm_;
    unsigned int order_;
    Polynomial<T> hermitePolynomial_;
};

template <class T>
void Gaussian<T>::calculateHermitePolynomial()
{
    if(order_ == 0)
    {
        hermitePolynomial_[0] = 1.0;
    }
    else if(order_ == 1)
    {
        hermitePolynomial_[0] = 2.0 * sigma2_;
    }
    else
    {
        T s2 = 2.0 * sigma2_;
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
        for(unsigned int i = 0; i <= hermitePolynomial_.order(); ++i)
            hermitePolynomial_[i] = order_ % 2 == 0 ?
                                         hn1[2*i]
                                       : hn1[2*i+1];
    }
}


//@}

} // namespace vigra


#endif /* VIGRA_GAUSSIANS_HXX */
