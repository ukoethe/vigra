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


#ifndef VIGRA_POLYNOMIAL_HXX
#define VIGRA_POLYNOMIAL_HXX

#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include "vigra/numerictraits.hxx"

namespace vigra {

namespace detail {

template <class T>
struct ComplexTraits
{
    typedef T Real;
    typedef std::complex<T> Complex;
    static T zero() { return 0.0; }
};

template <class T>
struct ComplexTraits<std::complex<T> >
{
    typedef T Real;
    typedef std::complex<T> Complex;
    static Complex zero() { return Complex(0.0); }
};

} // namespace detail

template <class T = double>
class Polynomial
{
  public:
    typedef typename detail::ComplexTraits<T>::Real    Real;
    typedef typename detail::ComplexTraits<T>::Complex Complex;
    typedef T value_type;
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    
    Polynomial(unsigned int ord = 0, double epsilon = 1.0e-14)
    : coeffs_(ord + 1, T()),
      epsilon_(epsilon)
    {}

    template <class ITER>
    Polynomial(ITER i, ITER end)
    : coeffs_(i, end),
      epsilon_(1.0e-14)
    {}
    
    template <class ITER>
    Polynomial(ITER i, ITER end, double epsilon)
    : coeffs_(i, end),
      epsilon_(epsilon)
    {}
    
    T & operator[](unsigned int i)
        { return coeffs_[i]; }
    
    T const & operator[](unsigned int i) const
        { return coeffs_[i]; }
    
    template <class U>
    typename PromoteTraits<T, U>::Promote
    operator()(U const & v) const;
    
    Polynomial<T> derivative() const;

    template <class U>
    Polynomial<typename PromoteTraits<T, U>::Promote > 
    deflate(U const & v) const;

    template <class U>
    Polynomial<typename PromoteTraits<T, U>::Promote > 
    forwardDeflate(U const & v) const;

    template <class U>
    Polynomial<typename PromoteTraits<T, U>::Promote > 
    backwardDeflate(U v) const;

    void minimizeOrder();
    
    void normalize();
    
    iterator begin()
        { return coeffs_.begin(); }
    
    iterator end()
        { return coeffs_.end(); }
    
    const_iterator begin() const
        { return coeffs_.begin(); }
    
    const_iterator end() const
        { return coeffs_.end(); }
    
    unsigned int size() const
        { return coeffs_.size(); }
        
    unsigned int order() const
        { return coeffs_.size() - 1; }
        
    double epsilon() const
        { return epsilon_; }
        
    void setEpsilon(double eps)
        { epsilon_ = eps; }

  protected:
    std::vector<T> coeffs_;
    double epsilon_;
};

template <class T>
template <class U>
typename PromoteTraits<T, U>::Promote
Polynomial<T>::operator()(U const & v) const
{
    typename PromoteTraits<T, U>::Promote p(coeffs_[order()]);
    for(int i = order() - 1; i >= 0; --i)
    {
       p = v * p + coeffs_[i];
    }
    return p;
}

template <class T>
Polynomial<T> 
Polynomial<T>::derivative() const
{
    Polynomial<T> p(order() - 1, epsilon_);
    for(unsigned int i = 1; i <= order(); ++i)
    {
        p[i-1] = double(i)*coeffs_[i];
    }
    p.minimizeOrder();
    return p;
}

template <class T>
template <class U>
Polynomial<typename PromoteTraits<T, U>::Promote > 
Polynomial<T>::deflate(U const & v) const
{
    if(std::abs(v) <= 1.0)
        return forwardDeflate(v);
    else
        return backwardDeflate(v);
}

template <class T>
template <class U>
Polynomial<typename PromoteTraits<T, U>::Promote > 
Polynomial<T>::forwardDeflate(U const & v) const
{
    typedef typename PromoteTraits<T, U>::Promote Promote;
    Polynomial<Promote> p(order() - 1, epsilon_);
    Promote tmp = coeffs_[order()];
    for(int i = order()-1; i >= 0; --i)
    {
        p[i] = tmp;
        tmp = coeffs_[i] + v * tmp;
    }
    return p;
}

template <class T>
template <class U>
Polynomial<typename PromoteTraits<T, U>::Promote > 
Polynomial<T>::backwardDeflate(U v) const
{
    typedef typename PromoteTraits<T, U>::Promote Promote;
    Polynomial<Promote> p(order() - 1, epsilon_);
    Promote tmp = coeffs_[0];
    v = 1.0 / v;
    for(unsigned int i = 0; i < order(); ++i)
    {
        p[i] = tmp;
        tmp = coeffs_[i+1] + v * tmp;
    }
    return p;
}

template <class T>
void 
Polynomial<T>::minimizeOrder()
{
    for(int i = order(); i >= 0; --i)
        if(std::abs(coeffs_[i]) < epsilon_)
            coeffs_.erase(coeffs_.begin() + i);
        else
            break;
}

template <class T>
void 
Polynomial<T>::normalize()
{
    for(unsigned int i = 0; i<order(); ++i)
        coeffs_[i] = coeffs_[order()];
    coeffs_[order()] = T(1.0);
}

namespace detail {

// replacement for complex division
// code form python complexobject.c
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
std::complex<T> deleteImaginaryBelowEpsilon(std::complex<T> const & x, double eps)
{
    return std::abs(x.imag()) <= 2.0*eps*std::abs(x.real()) ?
              std::complex<T>(x.real())
           :  x;
}

template <class T, class W, class U>
bool polynomialRoots(Polynomial<T> p, 
                     std::vector<std::complex<W> > & roots, 
                     U const & initial, int multiplicity)
{
    typedef typename detail::ComplexTraits<T>::Real    Real;
    typedef typename detail::ComplexTraits<T>::Complex Complex;
    
    p.minimizeOrder();
    
    double eps  = p.epsilon(),
           eps2 = VIGRA_CSTD::sqrt(eps);
    
    // handle the easy cases
    if(p.order() == 0)
        return true;
    if(std::abs(p[0]) < eps)
    {
        roots.push_back(Complex(0.0));
        polynomialRoots(Polynomial<T>(p.begin()+1, p.end()), roots, Complex(0.0), 0);
        return true;
    }
    if(p.order() == 1)
    {
        roots.push_back(deleteImaginaryBelowEpsilon(
                    complexDiv(Complex(-p[0]), Complex(p[1])), eps));
        return true;
    }
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
        roots.push_back(deleteImaginaryBelowEpsilon(complexDiv(q, a), eps));
        roots.push_back(deleteImaginaryBelowEpsilon(complexDiv(c, q), eps));
        return true;
    }

    // now for the standard case
    //    * initialze with the given initial value or 
    //      a value based on the root distribution
    //    * iteratively find root using Laguerre's method
    //    * if estimated multiplicity is > 1, improve root accuracy on derivative polynomial
    static double frac[] = {0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0};
    int maxiter = 80, 
        count, 
        work_with_derivative = 0;
    double N = p.order();
    Polynomial<T> polysave;
    Complex x, xsave;
        
    if(multiplicity > 0)
    {
        x = initial;
    }
    else
    {
        T centroid = -p[p.order()-1] / N / p[p.order()];
        x = centroid + VIGRA_CSTD::pow(Complex(-p(centroid) / p[p.order()]), 1.0 / N);
        multiplicity = 1;
    }
    
    for(count = 0; count < maxiter; ++count)
    {
        Complex p0(p[p.order()]);
        Complex p1 = Complex(0.0);
        Complex p2 = Complex(0.0);
        Real ax    = std::abs(x);
        double err = std::abs(p0);
        for(int i = p.order()-1; i >= 0; --i)
        {
            p2  = p2  * x  + p1;
            p1  = p1  * x  + p0;
            p0  = p0  * x  + p[i];
            err = err * ax + std::abs(p0);
        }
        p2 *= 2.0;
        err *= eps;
//        std::cerr << count << " " << multiplicity << " " << 
//                     x << " " << p0 << " " << p1 << " " << p2 << "\n";
        Real ap0 = std::abs(p0);
        if(ap0 <= err)
            break;  // converged
        Complex g = complexDiv(p1, p0);
        Complex g2 = g * g;
        Complex h = g2 - complexDiv(p2, p0);
        // estimate of root multiplicity according to Tien Chen
        if(g2 != 0.0)
            multiplicity = (unsigned int)VIGRA_CSTD::floor(N / 
                        (std::abs(N * complexDiv(h, g2) - 1.0) + 1.0) + 0.5);
        // improve accuracy of multiple roots on the derivative, as suggested by C. Bond
        // (do this only if we are already near the root, otherwise we may converge to 
        //  a different root of the derivative polynomial)
        if(multiplicity > 1 && ap0 < eps2)
        {
            if (work_with_derivative == 0)
            {
                polysave = p;
                xsave = x;  // remember for possible recovery from convergence to wrong root
            }
            ++work_with_derivative;
            --multiplicity;
            p = p.derivative();
            N = p.order();
            continue;
        }
        Complex sq = VIGRA_CSTD::sqrt((N - 1.0) * (N * h - g2));
        Complex gp = g + sq;
        Complex gm = g - sq;
        Real abp = std::abs(gp);
        Real abm = std::abs(gm);
        Complex dx;
        if(abp < abm)
            gp = gm;
        if(gp != 0.0)
        {
            dx = complexDiv(Complex(N) , gp);
        }
        else
        {
            dx = (1.0 + ax) * Complex(VIGRA_CSTD::cos(double(count)), VIGRA_CSTD::sin(double(count)));
        }
        Complex x1 = x - dx;

        if(x1 - x == 0.0)
            break; // convergence
        if((count + 1) % 10)
            x = x1;
        else
            x = x - frac[(count+1)/10] * dx;
    }
    if(count == maxiter)
    {
        return false;  // no convergence
    }
    if(work_with_derivative && std::abs(polysave(xsave)) < std::abs(polysave(x)))
    {
            // unsuccessful derivative iteration => restore xsave
            x = deleteImaginaryBelowEpsilon(xsave, eps);
    }
    else
    {
        x = deleteImaginaryBelowEpsilon(x, eps);
        if(work_with_derivative + 1 > multiplicity)
            multiplicity = work_with_derivative + 1;
    }
    roots.push_back(x);
    if(multiplicity > 1)
    {
        // multiple root => try same x again
        if(work_with_derivative)
            return polynomialRoots(polysave.deflate(x), roots, x, multiplicity - 1);
        else
            return polynomialRoots(p.deflate(x), roots, x, multiplicity - 1);
    }
    if(x.imag() != 0.0 && (roots.size() == 1 || roots[roots.size()-2] != conj(x)))
    {
        // new complex root => use conjugate as initial of next iteration
        return polynomialRoots(p.deflate(x), roots, conj(x), 1);
    }
    return polynomialRoots(p.deflate(x), roots, Complex(0.0), 0);
}

} // namespace detail 

template <class T, class U>
bool polynomialRoots(Polynomial<T> const & p, 
                     std::vector<std::complex<U> > & roots)
{
    typedef typename detail::ComplexTraits<T>::Complex Complex;
    return detail::polynomialRoots(p, roots, Complex(0.0), 0);
}

template <class T, class U>
bool realPolynomialRoots(Polynomial<T> const & p, std::vector<U> & roots)
{
    std::vector<std::complex<T> > croots;
    if(!polynomialRoots(p, croots))
        return false;
    for(unsigned int i = 0; i < croots.size(); ++i)
        if(croots[i].imag() == 0.0)
            roots.push_back(croots[i].real());
    return true;
}

} // namespace vigra

#endif // VIGRA_POLYNOMIAL_HXX
