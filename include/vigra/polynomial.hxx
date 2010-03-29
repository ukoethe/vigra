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


#ifndef VIGRA_POLYNOMIAL_HXX
#define VIGRA_POLYNOMIAL_HXX

#include <cmath>
#include <complex>
#include <algorithm>
#include <iosfwd>
#include "error.hxx"
#include "mathutil.hxx"
#include "numerictraits.hxx"
#include "array_vector.hxx"

namespace vigra {

template <class T> class Polynomial;
template <unsigned int MAXORDER, class T> class StaticPolynomial;

/*****************************************************************/
/*                                                               */
/*                         PolynomialView                        */
/*                                                               */
/*****************************************************************/

/** Polynomial interface for an externally managed array.

    The coefficient type <tt>T</tt> can be either a scalar or complex 
    (compatible to <tt>std::complex</tt>) type.
    
    \see vigra::Polynomial, vigra::StaticPolynomial, polynomialRoots()

    <b>\#include</b> \<<a href="polynomial_8hxx-source.html">vigra/polynomial.hxx</a>\><br>
    Namespace: vigra
    
    \ingroup Polynomials
*/
template <class T>
class PolynomialView
{
  public:
    
        /** Coefficient type of the polynomial
        */
    typedef T         value_type;

        /** Promote type of <tt>value_type</tt>
            used for polynomial calculations
        */
    typedef typename NumericTraits<T>::RealPromote RealPromote;

        /** Scalar type associated with <tt>RealPromote</tt>
        */
    typedef typename NumericTraits<RealPromote>::ValueType Real;

        /** Complex type associated with <tt>RealPromote</tt>
        */
    typedef typename NumericTraits<RealPromote>::ComplexPromote Complex;

        /** Iterator for the coefficient sequence
        */
    typedef T *       iterator;

        /** Const iterator for the coefficient sequence
        */
    typedef T const * const_iterator;
    
    typedef Polynomial<Real> RealPolynomial;
    typedef Polynomial<Complex> ComplexPolynomial;
    
    
        /** Construct from a coefficient array of given <tt>order</tt>.

            The externally managed array must have length <tt>order+1</tt>
            and is interpreted as representing the polynomial:
            
            \code
            coeffs[0] + x * coeffs[1] + x * x * coeffs[2] + ...
            \endcode
            
            The coefficients are not copied, we only store a pointer to the 
            array.<tt>epsilon</tt> (default: 1.0e-14) determines the precision 
            of subsequent algorithms (especially root finding) performed on the
            polynomial.
        */
    PolynomialView(T * coeffs, unsigned int order, double epsilon = 1.0e-14)
    : coeffs_(coeffs),
      order_(order),
      epsilon_(epsilon)
    {}
    
        /// Access the coefficient of x^i
    T & operator[](unsigned int i)
        { return coeffs_[i]; }
    
        /// Access the coefficient of x^i
    T const & operator[](unsigned int i) const
        { return coeffs_[i]; }
    
        /** Evaluate the polynomial at the point <tt>v</tt> 
        
            Multiplication must be defined between the types
            <tt>V</tt> and <tt>PromoteTraits<T, V>::Promote</tt>.
            If both <tt>V</tt> and <tt>V</tt> are scalar, the result will
            be a scalar, otherwise it will be complex.
        */
    template <class V>
    typename PromoteTraits<T, V>::Promote
    operator()(V v) const;
    
        /** Differentiate the polynomial <tt>n</tt> times.
        */
    void differentiate(unsigned int n = 1);
    
        /** Deflate the polynomial at the root <tt>r</tt> with 
            the given <tt>multiplicity</tt>.
            
            The behavior of this function is undefined if <tt>r</tt>
            is not a root with at least the given multiplicity.
            This function calls forwardBackwardDeflate().
        */
    void deflate(T const & r, unsigned int multiplicity = 1);
    
        /** Forward-deflate the polynomial at the root <tt>r</tt>.
            
            The behavior of this function is undefined if <tt>r</tt>
            is not a root. Forward deflation is best if <tt>r</tt> is
            the biggest root (by magnitude).
        */
    void forwardDeflate(T const & v);
    
        /** Forward/backward eflate the polynomial at the root <tt>r</tt>.
            
            The behavior of this function is undefined if <tt>r</tt>
            is not a root. Combined forward/backward deflation is best 
            if <tt>r</tt> is an ontermediate root or we don't know
            <tt>r</tt>'s relation to the other roots of the polynomial.
        */
    void forwardBackwardDeflate(T v);
    
        /** Backward-deflate the polynomial at the root <tt>r</tt>.
            
            The behavior of this function is undefined if <tt>r</tt>
            is not a root. Backward deflation is best if <tt>r</tt> is
            the snallest root (by magnitude).
        */
    void backwardDeflate(T v);
    
        /** Deflate the polynomial with the complex conjugate roots 
            <tt>r</tt> and <tt>conj(r)</tt>.
            
            The behavior of this function is undefined if these are not
            roots.
        */
    void deflateConjugatePair(Complex const & v);
    
        /** Adjust the polynomial's order if the highest coefficients are near zero.
            The order is reduced as long as the absolute value does not exceed 
            the given \a epsilon.
        */
    void minimizeOrder(double epsilon = 0.0);    
    
        /** Normalize the polynomial, i.e. dived by the highest coefficient.
        */
    void normalize();
     
    void balance();
        
        /** Get iterator for the coefficient sequence.
        */
    iterator begin()
        { return coeffs_; }
    
        /** Get end iterator for the coefficient sequence.
        */
    iterator end()
        { return begin() + size(); }
    
        /** Get const_iterator for the coefficient sequence.
        */
    const_iterator begin() const
        { return coeffs_; }
    
        /** Get end const_iterator for the coefficient sequence.
        */
    const_iterator end() const
        { return begin() + size(); }
    
        /** Get length of the coefficient sequence (<tt>order() + 1</tt>).
        */
    unsigned int size() const
        { return order_ + 1; }
        
        /** Get order of the polynomial.
        */
    unsigned int order() const
        { return order_; }
        
        /** Get requested precision for polynomial algorithms 
            (especially root finding).
        */
    double epsilon() const
        { return epsilon_; }
        
        /** Set requested precision for polynomial algorithms 
            (especially root finding).
        */
    void setEpsilon(double eps)
        { epsilon_ = eps; }

  protected:
    PolynomialView(double epsilon = 1e-14)
    : coeffs_(0),
      order_(0),
      epsilon_(epsilon)
    {}
    
    void setCoeffs(T * coeffs, unsigned int order)
    {
        coeffs_ = coeffs;
        order_ = order;
    }
  
    T * coeffs_;
    unsigned int order_;
    double epsilon_;
};

template <class T>
template <class U>
typename PromoteTraits<T, U>::Promote
PolynomialView<T>::operator()(U v) const
{
    typename PromoteTraits<T, U>::Promote p(coeffs_[order_]);
    for(int i = order_ - 1; i >= 0; --i)
    {
       p = v * p + coeffs_[i];
    }
    return p;
}

/*
template <class T>
typename PolynomialView<T>::Complex 
PolynomialView<T>::operator()(Complex const & v) const
{
    Complex p(coeffs_[order_]);
    for(int i = order_ - 1; i >= 0; --i)
    {
       p = v * p + coeffs_[i];
    }
    return p;
}
*/

template <class T>
void
PolynomialView<T>::differentiate(unsigned int n)
{
    if(n == 0)
        return;
    if(order_ == 0)
    {
        coeffs_[0] = 0.0;
        return;
    }
    for(unsigned int i = 1; i <= order_; ++i)
    {
        coeffs_[i-1] = double(i)*coeffs_[i];
    }
    --order_;
    if(n > 1)
        differentiate(n-1);
}

template <class T>
void
PolynomialView<T>::deflate(T const & v, unsigned int multiplicity)
{
    vigra_precondition(order_ > 0,
        "PolynomialView<T>::deflate(): cannot deflate 0th order polynomial.");
    if(v == 0.0)
    {
        ++coeffs_;
        --order_;
    }
    else
    {
        // we use combined forward/backward deflation because
        // our initial guess seems to favour convergence to 
        // a root with magnitude near the median among all roots
        forwardBackwardDeflate(v);
    }
    if(multiplicity > 1)
        deflate(v, multiplicity-1);
}

template <class T>
void
PolynomialView<T>::forwardDeflate(T const & v)
{
    for(int i = order_-1; i > 0; --i)
    {
        coeffs_[i] += v * coeffs_[i+1];
    }
    ++coeffs_;
    --order_;
}

template <class T>
void
PolynomialView<T>::forwardBackwardDeflate(T v)
{
    unsigned int order2 = order_ / 2;
    T tmp = coeffs_[order_];
    for(unsigned int i = order_-1; i >= order2; --i)
    {
        T tmp1 = coeffs_[i];
        coeffs_[i] = tmp;
        tmp = tmp1 + v * tmp;
    }
    v = -1.0 / v;
    coeffs_[0] *= v;
    for(unsigned int i = 1; i < order2; ++i)
    {
        coeffs_[i] = v * (coeffs_[i] - coeffs_[i-1]);
    }
    --order_;
}

template <class T>
void
PolynomialView<T>::backwardDeflate(T v)
{
    v = -1.0 / v;
    coeffs_[0] *= v;
    for(unsigned int i = 1; i < order_; ++i)
    {
        coeffs_[i] = v * (coeffs_[i] - coeffs_[i-1]);
    }
    --order_;
}

template <class T>
void
PolynomialView<T>::deflateConjugatePair(Complex const & v)
{
    vigra_precondition(order_ > 1,
        "PolynomialView<T>::deflateConjugatePair(): cannot deflate 2 roots "
        "from 1st order polynomial.");
    Real a = 2.0*v.real();
    Real b = -sq(v.real()) - sq(v.imag());
    coeffs_[order_-1] += a * coeffs_[order_];
    for(int i = order_-2; i > 1; --i)
    {
        coeffs_[i] += a * coeffs_[i+1] + b*coeffs_[i+2];
    }
    coeffs_ += 2;
    order_ -= 2;
}
    
template <class T>
void 
PolynomialView<T>::minimizeOrder(double epsilon)
{
    while(std::abs(coeffs_[order_]) <= epsilon && order_ > 0)
            --order_;
}

template <class T>
void 
PolynomialView<T>::normalize()
{
    for(unsigned int i = 0; i<order_; ++i)
        coeffs_[i] /= coeffs_[order_];
    coeffs_[order_] = T(1.0);
}

template <class T>
void 
PolynomialView<T>::balance()
{
    Real p0 = abs(coeffs_[0]), po = abs(coeffs_[order_]);
    Real norm = (p0 > 0.0)
                    ? VIGRA_CSTD::sqrt(p0*po) 
                    : po;
    for(unsigned int i = 0; i<=order_; ++i)
        coeffs_[i] /= norm;
}

/*****************************************************************/
/*                                                               */
/*                           Polynomial                          */
/*                                                               */
/*****************************************************************/

/** Polynomial with internally managed array.

    Most interesting functionality is inherited from \ref vigra::PolynomialView.

    \see vigra::PolynomialView, vigra::StaticPolynomial, polynomialRoots()

    <b>\#include</b> \<<a href="polynomial_8hxx-source.html">vigra/polynomial.hxx</a>\><br>
    Namespace: vigra
    
    \ingroup Polynomials
*/
template <class T>
class Polynomial
: public PolynomialView<T>
{
    typedef PolynomialView<T> BaseType;
  public:
    typedef typename BaseType::Real    Real;
    typedef typename BaseType::Complex Complex;
    typedef Polynomial<Real>           RealPolynomial;
    typedef Polynomial<Complex>        ComplexPolynomial;

    typedef T         value_type;
    typedef T *       iterator;
    typedef T const * const_iterator;    
    
        /** Construct polynomial with given <tt>order</tt> and all coefficients
            set to zero (they can be set later using <tt>operator[]</tt>
            or the iterators). <tt>epsilon</tt> (default: 1.0e-14) determines 
            the precision of subsequent algorithms (especially root finding) 
            performed on the polynomial.
        */
    Polynomial(unsigned int order = 0, double epsilon = 1.0e-14)
    : BaseType(epsilon),
      polynomial_(order + 1, T())
    {
        this->setCoeffs(&polynomial_[0], order);
    }
    
        /** Copy constructor
        */
    Polynomial(Polynomial const & p)
    : BaseType(p.epsilon()),
      polynomial_(p.begin(), p.end())
    {
        this->setCoeffs(&polynomial_[0], p.order());
    }

        /** Construct polynomial by copying the given coefficient sequence.
        */
    template <class ITER>
    Polynomial(ITER i, unsigned int order)
    : BaseType(),
      polynomial_(i, i + order + 1)
    {
        this->setCoeffs(&polynomial_[0], order);
    }
    
        /** Construct polynomial by copying the given coefficient sequence.
            Set <tt>epsilon</tt> (default: 1.0e-14) as 
            the precision of subsequent algorithms (especially root finding) 
            performed on the polynomial.
        */
    template <class ITER>
    Polynomial(ITER i, unsigned int order, double epsilon)
    : BaseType(epsilon),
      polynomial_(i, i + order + 1)
    {
        this->setCoeffs(&polynomial_[0], order);
    }
    
        /** Assigment
        */
    Polynomial & operator=(Polynomial const & p)
    {
        if(this == &p)
            return *this;
        ArrayVector<T> tmp(p.begin(), p.end());
        polynomial_.swap(tmp);
        this->setCoeffs(&polynomial_[0], p.order());
        this->epsilon_ = p.epsilon_;
        return *this;
    }
    
        /** Construct new polynomial representing the derivative of this
            polynomial.
        */
    Polynomial<T> getDerivative(unsigned int n = 1) const
    {
        Polynomial<T> res(*this);
        res.differentiate(n);
        return res;
    }

        /** Construct new polynomial representing this polynomial after
            deflation at the real root <tt>r</tt>.
        */
    Polynomial<T> 
    getDeflated(Real r) const
    {
        Polynomial<T> res(*this);
        res.deflate(r);
        return res;
    }

        /** Construct new polynomial representing this polynomial after
            deflation at the complex root <tt>r</tt>. The resulting
            polynomial will have complex coefficients, even if this
            polynomial had real ones.
        */
    Polynomial<Complex> 
    getDeflated(Complex const & r) const
    {
        Polynomial<Complex> res(this->begin(), this->order(), this->epsilon());
        res.deflate(r);
        return res;
    }

  protected:
    ArrayVector<T> polynomial_;
};

/*****************************************************************/
/*                                                               */
/*                        StaticPolynomial                       */
/*                                                               */
/*****************************************************************/

/** Polynomial with internally managed array of static length.

    Most interesting functionality is inherited from \ref vigra::PolynomialView.
    This class differs from \ref vigra::Polynomial in that it allocates
    its memory statically which is much faster. Therefore, <tt>StaticPolynomial</tt>
    can only represent polynomials up to the given <tt>MAXORDER</tt>.

    \see vigra::PolynomialView, vigra::Polynomial, polynomialRoots()

    <b>\#include</b> \<<a href="polynomial_8hxx-source.html">vigra/polynomial.hxx</a>\><br>
    Namespace: vigra
    
    \ingroup Polynomials
*/
template <unsigned int MAXORDER, class T>
class StaticPolynomial
: public PolynomialView<T>
{
    typedef PolynomialView<T> BaseType;
    
  public:
    typedef typename BaseType::Real    Real;
    typedef typename BaseType::Complex Complex;
    typedef StaticPolynomial<MAXORDER, Real> RealPolynomial;
    typedef StaticPolynomial<MAXORDER, Complex> ComplexPolynomial;

    typedef T         value_type;
    typedef T *       iterator;
    typedef T const * const_iterator;
    
    
        /** Construct polynomial with given <tt>order <= MAXORDER</tt> and all 
            coefficients set to zero (they can be set later using <tt>operator[]</tt>
            or the iterators). <tt>epsilon</tt> (default: 1.0e-14) determines 
            the precision of subsequent algorithms (especially root finding) 
            performed on the polynomial.
        */
    StaticPolynomial(unsigned int order = 0, double epsilon = 1.0e-14)
    : BaseType(epsilon)
    {
        vigra_precondition(order <= MAXORDER,
            "StaticPolynomial(): order exceeds MAXORDER.");
        std::fill_n(polynomial_, order+1, T());
        this->setCoeffs(polynomial_, order);
    }
    
        /** Copy constructor
        */
    StaticPolynomial(StaticPolynomial const & p)
    : BaseType(p.epsilon())
    {
        std::copy(p.begin(), p.end(), polynomial_);
        this->setCoeffs(polynomial_, p.order());
    }

        /** Construct polynomial by copying the given coefficient sequence.
            <tt>order <= MAXORDER</tt> is required.
        */
    template <class ITER>
    StaticPolynomial(ITER i, unsigned int order)
    : BaseType()
    {
        vigra_precondition(order <= MAXORDER,
            "StaticPolynomial(): order exceeds MAXORDER.");
        std::copy(i, i + order + 1, polynomial_);
        this->setCoeffs(polynomial_, order);
    }
    
        /** Construct polynomial by copying the given coefficient sequence.
            <tt>order <= MAXORDER</tt> is required. Set <tt>epsilon</tt> (default: 1.0e-14) as 
            the precision of subsequent algorithms (especially root finding) 
            performed on the polynomial.
        */
    template <class ITER>
    StaticPolynomial(ITER i, unsigned int order, double epsilon)
    : BaseType(epsilon)
    {
        vigra_precondition(order <= MAXORDER,
            "StaticPolynomial(): order exceeds MAXORDER.");
        std::copy(i, i + order + 1, polynomial_);
        this->setCoeffs(polynomial_, order);
    }
    
        /** Assigment.
        */
    StaticPolynomial & operator=(StaticPolynomial const & p)
    {
        if(this == &p)
            return *this;
        std::copy(p.begin(), p.end(), polynomial_);
        this->setCoeffs(polynomial_, p.order());
        this->epsilon_ = p.epsilon_;
        return *this;
    }
    
        /** Construct new polynomial representing the derivative of this
            polynomial.
        */
    StaticPolynomial getDerivative(unsigned int n = 1) const
    {
        StaticPolynomial res(*this);
        res.differentiate(n);
        return res;
    }

        /** Construct new polynomial representing this polynomial after
            deflation at the real root <tt>r</tt>.
        */
    StaticPolynomial 
    getDeflated(Real r) const
    {
        StaticPolynomial res(*this);
        res.deflate(r);
        return res;
    }

        /** Construct new polynomial representing this polynomial after
            deflation at the complex root <tt>r</tt>. The resulting
            polynomial will have complex coefficients, even if this
            polynomial had real ones.
        */
    StaticPolynomial<MAXORDER, Complex> 
    getDeflated(Complex const & r) const
    {
        StaticPolynomial<MAXORDER, Complex>  res(this->begin(), this->order(), this->epsilon());
        res.deflate(r);
        return res;
    }
    
    void setOrder(unsigned int order)
    {
        vigra_precondition(order <= MAXORDER,
            "taticPolynomial::setOrder(): order exceeds MAXORDER.");
        this->order_ = order;
    }

  protected:
    T polynomial_[MAXORDER+1];
};

/************************************************************/

namespace detail {

// replacement for complex division (some compilers have numerically
// less stable implementations); code from python complexobject.c
template <class T>
std::complex<T> complexDiv(std::complex<T> const & a, std::complex<T> const & b)
{
 	 const double abs_breal = b.real() < 0 ? -b.real() : b.real();
	 const double abs_bimag = b.imag() < 0 ? -b.imag() : b.imag();

	 if (abs_breal >= abs_bimag) 
	 {
 		/* divide tops and bottom by b.real() */
	 	if (abs_breal == 0.0) 
	 	{
	 		return std::complex<T>(a.real() / abs_breal, a.imag() / abs_breal);
	 	}
	 	else 
	 	{
	 		const double ratio = b.imag() / b.real();
	 		const double denom = b.real() + b.imag() * ratio;
	 		return std::complex<T>((a.real() + a.imag() * ratio) / denom,
	 		                       (a.imag() - a.real() * ratio) / denom);
	 	}
	}
	else 
	{
		/* divide tops and bottom by b.imag() */
		const double ratio = b.real() / b.imag();
		const double denom = b.real() * ratio + b.imag();
		return std::complex<T>((a.real() * ratio + a.imag()) / denom,
                               (a.imag() * ratio - a.real()) / denom);
	}
}

template <class T>
std::complex<T> deleteBelowEpsilon(std::complex<T> const & x, double eps)
{
    return std::abs(x.imag()) <= 2.0*eps*std::abs(x.real()) 
                ? std::complex<T>(x.real())
                : std::abs(x.real()) <= 2.0*eps*std::abs(x.imag()) 
                    ? std::complex<T>(NumericTraits<T>::zero(), x.imag())
                    :  x;
}

template <class POLYNOMIAL>
typename POLYNOMIAL::value_type
laguerreStartingGuess(POLYNOMIAL const & p)
{
    double N = p.order();
    typename POLYNOMIAL::value_type centroid = -p[p.order()-1] / N / p[p.order()];
    double dist = VIGRA_CSTD::pow(std::abs(p(centroid) / p[p.order()]), 1.0 / N);
    return centroid + dist;
}

template <class POLYNOMIAL, class Complex>
int laguerre1Root(POLYNOMIAL const & p, Complex & x, unsigned int multiplicity)
{
    typedef typename NumericTraits<Complex>::ValueType Real;
    
    static double frac[] = {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};
    int maxiter = 80, 
        count;
    double N = p.order();
    double eps  = p.epsilon(),
           eps2 = VIGRA_CSTD::sqrt(eps);
        
    if(multiplicity == 0)
        x = laguerreStartingGuess(p);
        
    bool mayTryDerivative = true;  // try derivative for multiple roots
    
    for(count = 0; count < maxiter; ++count)
    {
        // Horner's algorithm to calculate values of polynomial and its
        // first two derivatives and estimate error for current x
        Complex p0(p[p.order()]);
        Complex p1(0.0);
        Complex p2(0.0);
        Real ax    = std::abs(x);
        Real err = std::abs(p0);
        for(int i = p.order()-1; i >= 0; --i)
        {
            p2  = p2  * x  + p1;
            p1  = p1  * x  + p0;
            p0  = p0  * x  + p[i];
            err = err * ax + std::abs(p0);
        }
        p2 *= 2.0;
        err *= eps;
        Real ap0 = std::abs(p0);
        if(ap0 <= err)
        {
            break;  // converged
        }
        Complex g = complexDiv(p1, p0);
        Complex g2 = g * g;
        Complex h = g2 - complexDiv(p2, p0);
        // estimate root multiplicity according to Tien Chen
        if(g2 != 0.0)
        {
            multiplicity = (unsigned int)VIGRA_CSTD::floor(N / 
                                (std::abs(N * complexDiv(h, g2) - 1.0) + 1.0) + 0.5);
            if(multiplicity < 1)
                multiplicity = 1;
        }
        // improve accuracy of multiple roots on the derivative, as suggested by C. Bond
        // (do this only if we are already near the root, otherwise we may converge to 
        //  a different root of the derivative polynomial)
        if(mayTryDerivative && multiplicity > 1 && ap0 < eps2)
        {
            Complex x1 = x;
            int derivativeMultiplicity = laguerre1Root(p.getDerivative(), x1, multiplicity-1);
            if(derivativeMultiplicity && std::abs(p(x1)) < std::abs(p(x)))
            {
                // successful search on derivative
                x = x1;
                return derivativeMultiplicity + 1;
            }
            else
            {
                // unsuccessful search on derivative => don't do it again
                mayTryDerivative = false;
            }
        }
        Complex sq = VIGRA_CSTD::sqrt((N - 1.0) * (N * h - g2));
        Complex gp = g + sq;
        Complex gm = g - sq;
        if(std::abs(gp) < std::abs(gm))
            gp = gm;
        Complex dx;
        if(gp != 0.0)
        {
            dx = complexDiv(Complex(N) , gp);
        }
        else
        {
            // re-initialisation trick due to Numerical Recipes
            dx = (1.0 + ax) * Complex(VIGRA_CSTD::cos(double(count)), VIGRA_CSTD::sin(double(count)));
        }
        Complex x1 = x - dx;

        if(x1 - x == 0.0)
        {
            break;  // converged
        }
        if((count + 1) % 10)
            x = x1;
        else
            // cycle breaking trick according to Numerical Recipes
            x = x - frac[(count+1)/10] * dx;
    }
    return count < maxiter ? 
        multiplicity : 
        0;
}

template <class Real>
struct PolynomialRootCompare
{
    Real epsilon;

    PolynomialRootCompare(Real eps)
    : epsilon(eps)
    {}
    
    template <class T>
    bool operator()(T const & l, T const & r)
    {
        return closeAtTolerance(l.real(), r.real(), epsilon)
                     ? l.imag() < r.imag()
                     : l.real() < r.real();
    }
};

} // namespace detail 

/** \addtogroup Polynomials Polynomials and root determination

    Classes to represent polynomials and functions to find polynomial roots.
*/
//@{

/*****************************************************************/
/*                                                               */
/*                         polynomialRoots                       */
/*                                                               */
/*****************************************************************/

/** Determine the roots of the polynomial <tt>poriginal</tt>.

    The roots are appended to the vector <tt>roots</tt>, with optional root
    polishing as specified by <tt>polishRoots</tt> (default: do polishing). The function uses an 
    improved version of Laguerre's algorithm. The improvements are as follows:
    
    <ul>
    <li>It uses a clever initial guess for the iteration, according to a proposal by Tien Chen</li>
    <li>It estimates each root's multiplicity, again according to Tien Chen, and reduces multiplicity
        by switching to the polynomial's derivative (which has the same root, with multiplicity
        reduced by one), as proposed by C. Bond.</li>
    </ul>
    
    The algorithm has been successfully used for polynomials up to order 80.
    The function stops and returns <tt>false</tt> if an iteration fails to converge within 
    80 steps. The type <tt>POLYNOMIAL</tt> must be compatible to 
    \ref vigra::PolynomialView, <tt>VECTOR</tt> must be compatible to <tt>std::vector</tt>
    with a <tt>value_type</tt> compatible to the type <tt>POLYNOMIAL::Complex</tt>.

    <b> Declaration:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class POLYNOMIAL, class VECTOR>
        bool 
        polynomialRoots(POLYNOMIAL const & poriginal, VECTOR & roots, bool polishRoots = true);
    }
    \endcode


    <b> Usage:</b>

        <b>\#include</b> \<<a href="polynomial_8hxx-source.html">vigra/polynomial.hxx</a>\><br>
        Namespace: vigra

    \code
    // encode the polynomial  x^4 - 1
    Polynomial<double> poly(4);
    poly[0] = -1.0;
    poly[4] =  1.0;

    ArrayVector<std::complex<double> > roots;
    polynomialRoots(poly, roots);
    \endcode

    \see polynomialRootsEigenvalueMethod()
*/
template <class POLYNOMIAL, class VECTOR>
bool polynomialRoots(POLYNOMIAL const & poriginal, VECTOR & roots, bool polishRoots)
{
    typedef typename POLYNOMIAL::value_type T;
    typedef typename POLYNOMIAL::Real    Real;
    typedef typename POLYNOMIAL::Complex Complex;
    typedef typename POLYNOMIAL::ComplexPolynomial WorkPolynomial;
    
    double eps  = poriginal.epsilon();

    WorkPolynomial p(poriginal.begin(), poriginal.order(), eps);
    p.minimizeOrder();
    if(p.order() == 0)
        return true;

    Complex x = detail::laguerreStartingGuess(p);
    
    unsigned int multiplicity = 1;
    bool triedConjugate = false;
    
    // handle the high order cases
    while(p.order() > 2)
    {
        p.balance();
        
        // find root estimate using Laguerre's method on deflated polynomial p;
        // zero return indicates failure to converge
        multiplicity = detail::laguerre1Root(p, x, multiplicity);
    
        if(multiplicity == 0)
            return false;
        // polish root on original polynomial poriginal;
        // zero return indicates failure to converge
        if(polishRoots && !detail::laguerre1Root(poriginal, x, multiplicity))
            return false;
        x = detail::deleteBelowEpsilon(x, eps);
        roots.push_back(x);
        p.deflate(x);
        // determine the next starting guess
        if(multiplicity > 1)
        {
            // probably multiple root => keep current root as starting guess
            --multiplicity;
            triedConjugate = false;
        }
        else
        {
            // need a new starting guess
            if(x.imag() != 0.0 && !triedConjugate)
            {
                // if the root is complex and we don't already have 
                // the conjugate root => try the conjugate as starting guess
                triedConjugate = true;
                x = conj(x);
            }
            else
            {
                // otherwise generate new starting guess
                triedConjugate = false;
                x = detail::laguerreStartingGuess(p);
            }
        }
    }
    
    // handle the low order cases
    if(p.order() == 2)
    {
        Complex a = p[2];
        Complex b = p[1];
        Complex c = p[0];
        Complex b2 = std::sqrt(b*b - 4.0*a*c);
        Complex q;
        if((conj(b)*b2).real() >= 0.0)
            q = -0.5 * (b + b2);
        else
            q = -0.5 * (b - b2);
        x = detail::complexDiv(q, a);
        if(polishRoots)
            detail::laguerre1Root(poriginal, x, 1);
        roots.push_back(detail::deleteBelowEpsilon(x, eps));
        x = detail::complexDiv(c, q);
        if(polishRoots)
            detail::laguerre1Root(poriginal, x, 1);
        roots.push_back(detail::deleteBelowEpsilon(x, eps));
    }
    else if(p.order() == 1)
    {
        x = detail::complexDiv(-p[0], p[1]);
        if(polishRoots)
            detail::laguerre1Root(poriginal, x, 1);
        roots.push_back(detail::deleteBelowEpsilon(x, eps));
    }
    std::sort(roots.begin(), roots.end(), detail::PolynomialRootCompare<Real>(eps));
    return true;
}

template <class POLYNOMIAL, class VECTOR>
inline bool 
polynomialRoots(POLYNOMIAL const & poriginal, VECTOR & roots)
{
    return polynomialRoots(poriginal, roots, true);
}

/** Determine the real roots of the polynomial <tt>p</tt>.

    This function simply calls \ref polynomialRoots() and than throws away all complex roots.
    Accordingly, <tt>VECTOR</tt> must be compatible to <tt>std::vector</tt>
    with a <tt>value_type</tt> compatible to the type <tt>POLYNOMIAL::Real</tt>.

    <b> Declaration:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class POLYNOMIAL, class VECTOR>
        bool 
        polynomialRealRoots(POLYNOMIAL const & p, VECTOR & roots, bool polishRoots = true);
    }
    \endcode


    <b> Usage:</b>

        <b>\#include</b> \<<a href="polynomial_8hxx-source.html">vigra/polynomial.hxx</a>\><br>
        Namespace: vigra

    \code
    // encode the polynomial  x^4 - 1
    Polynomial<double> poly(4);
    poly[0] = -1.0;
    poly[4] =  1.0;

    ArrayVector<double> roots;
    polynomialRealRoots(poly, roots);
    \endcode

    \see polynomialRealRootsEigenvalueMethod()
*/
template <class POLYNOMIAL, class VECTOR>
bool polynomialRealRoots(POLYNOMIAL const & p, VECTOR & roots, bool polishRoots)
{
    typedef typename NumericTraits<typename VECTOR::value_type>::ComplexPromote Complex;
    ArrayVector<Complex> croots;
    if(!polynomialRoots(p, croots, polishRoots))
        return false;
    for(unsigned int i = 0; i < croots.size(); ++i)
        if(croots[i].imag() == 0.0)
            roots.push_back(croots[i].real());
    return true;
}

template <class POLYNOMIAL, class VECTOR>
inline bool 
polynomialRealRoots(POLYNOMIAL const & poriginal, VECTOR & roots)
{
    return polynomialRealRoots(poriginal, roots, true);
}

//@}

} // namespace vigra

namespace std {

template <class T>
ostream & operator<<(ostream & o, vigra::PolynomialView<T> const & p)
{
    for(unsigned int k=0; k < p.order(); ++k)
        o << p[k] << " ";
    o << p[p.order()];
    return o;
}

} // namespace std 

#endif // VIGRA_POLYNOMIAL_HXX
