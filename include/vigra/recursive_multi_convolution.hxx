//-- -*- c++ -*-
/************************************************************************/
/*                                                                      */
/*               Copyright 2016 by Sven Peter and Ullrich Koethe        */
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

#ifndef VIGRA_RECURSIVE_MULTI_CONVOLUTION_H
#define VIGRA_RECURSIVE_MULTI_CONVOLUTION_H

// TODO: cleanup includes
#include "kerneltraits.hxx"
#include "separableconvolution.hxx"
#include "array_vector.hxx"
#include "multi_array.hxx"
#include "accessor.hxx"
#include "numerictraits.hxx"
#include "navigator.hxx"
#include "metaprogramming.hxx"
#include "multi_pointoperators.hxx"
#include "multi_math.hxx"
#include "convolution.hxx"
#include "multi_convolution.hxx"
#include "functorexpression.hxx"
#include "tinyvector.hxx"
#include "algorithm.hxx"
#include "matrix.hxx"
#include "linear_algebra.hxx"
#include "linear_solve.hxx"

#include <xmmintrin.h>
#include <emmintrin.h>
#include <complex>
#include <type_traits>

namespace vigra {

namespace detail {
// Rachid Deriche. Recursively implementating the Gaussian and its derivatives.
// [Research Report] RR-1893, 1993, pp.24. <inria-00074778>

struct DerichePrecomputed
{
    const std::complex<double> alpha0, alpha1, alpha2, alpha3;
    const std::complex<double> lambda0, lambda1, lambda2, lambda3;

    DerichePrecomputed(const double a0, const double a1, const double a2, const double a3, const double l0, const double l1, const double l2, const double l3) :
    alpha0(std::complex<double>(a0 / 2.0, a1 / 2.0)), alpha1(std::complex<double>(a0 / 2.0, -a1 / 2.0)), alpha2(std::complex<double>(a2 / 2.0, a3 / 2.0)), alpha3(std::complex<double>(a2 / 2.0, -a3 / 2.0)),
    lambda0(std::complex<double>(l0, l1)), lambda1(std::complex<double>(l0, -l1)), lambda2(std::complex<double>(l2, l3)), lambda3(std::complex<double>(l2, -l3))
    {
    }

    DerichePrecomputed(const double a0, const double a1, const double a2, const double l0, const double l1, const double l2) :
    alpha0(std::complex<double>(a0 / 2.0, a1 / 2.0)), alpha1(std::complex<double>(a0 / 2.0, -a1 / 2.0)), alpha2(std::complex<double>(a2, 0)), alpha3(std::complex<double>(a2, 0)),
    lambda0(std::complex<double>(l0, l1)), lambda1(std::complex<double>(l0, -l1)), lambda2(std::complex<double>(l2, 0)), lambda3(std::complex<double>(l2, 0))
    {
    }

    DerichePrecomputed(const double a0, const double a1, const double l0, const double l1) :
    alpha0(std::complex<double>(a0 / 2.0, a1 / 2.0)), alpha1(std::complex<double>(a0 / 2.0, -a1 / 2.0)), alpha2(std::complex<double>(0, 0)), alpha3(std::complex<double>(0, 0)),
    lambda0(std::complex<double>(l0, l1)), lambda1(std::complex<double>(l0, -l1)), lambda2(std::complex<double>(0, 0)), lambda3(std::complex<double>(0, 0))
    {
    }

    const std::complex<double> get_alpha(const unsigned int idx) const {
        switch (idx) {
            case 0: return alpha0;
            case 1: return alpha1;
            case 2: return alpha2;
            case 3: return alpha3;
        }

        vigra_fail("DerichePrecomputed::get_alpha: index out of bounds.");
        return 0;
    }

    const std::complex<double> get_lambda(const unsigned int idx) const {
        switch (idx) {
            case 0: return lambda0;
            case 1: return lambda1;
            case 2: return lambda2;
            case 3: return lambda3;
        }

        vigra_fail("DerichePrecomputed::get_lambda: index out of bounds.");
        return 0;
    }
};

// eq. 35 - 46
// 4th order coefficients from The Mines Java Toolkit
// https://github.com/dhale/jtk/blob/master/core/src/main/java/edu/mines/jtk/dsp/RecursiveGaussianFilter.java

static const DerichePrecomputed deriche_precomputed_coefs[3][3] = {
    { // smoothing
        DerichePrecomputed(0.9629, 1.942, 1.26, 0.8448),
        DerichePrecomputed(-0.8929, 1.021, 1.898, 1.512, 1.475, 1.556),
        DerichePrecomputed(1.6797292232361107, 3.7348298269103580, -0.6802783501806897, -0.2598300478959625, 1.7831906544515104, 0.6318113174569493, 1.7228297663338028, 1.9969276832487770)
    },
    { // first derivative
        DerichePrecomputed(0.1051, -1.898, 0.9338, 0.9459),
        DerichePrecomputed(1.666, -0.4557, -1.682, 1.266, 1.558, 1.244),
        DerichePrecomputed(0.6494024008440620, 0.9557370760729773, -0.6472105276644291, -4.5306923044570760, 1.5159726670750566, 2.0718953658782650, 1.5267608734791140, 0.6719055957689513)
    },
    { // second derivative
        DerichePrecomputed(-1.228, 0.2976, 0.6134, 1.332),
        DerichePrecomputed(-1.773, -1.129, 0.8309, 1.032, 1.664, 0.8459),
        DerichePrecomputed(0.3224570510072559, -1.7382843963561239, -1.3312275593739595, 3.6607035671974897, 1.3138054926516880, 2.1656041357418863, 1.2402181393295362, 0.7479888745408682)
    }
};

// Van Vliet, Lucas J., Ian T. Young, and Piet W. Verbeek. "Recursive Gaussian
// derivative filters." Pattern Recognition, 1998. Proceedings. Fourteenth
// International Conference on. Vol. 1. IEEE, 1998.
struct VYVPrecomputed
{
    const double pole0, pole1, pole2, pole3, pole4;

    VYVPrecomputed(const double p0, const double p1, const double p2, const double p3 = 0, const double p4 = 0) :
    pole0(p0), pole1(p1), pole2(p2), pole3(p3), pole4(p4)
    {
    }

    double get_pole(const unsigned int idx) const {
        switch (idx) {
            case 0: return pole0;
            case 1: return pole1;
            case 2: return pole2;
            case 3: return pole3;
            case 4: return pole4;
        }

        return 0;
    }

    std::complex<double> get_complex_pole(const unsigned int idx) const {
        double real, imag;
        unsigned int offset = idx & ~1;

        real = get_pole(offset);
        imag = get_pole(offset + 1);

        if (idx % 2)
            imag = -imag;

        return std::complex<double>(real, imag);
    }
};


// Table 1 and 2
static const VYVPrecomputed vyv_poles[3][3] = {
#if 1
    { // smoothing, L_2
        VYVPrecomputed(1.4165, 1.00829, 1.86543),
        VYVPrecomputed(1.13228, 1.28114, 1.78534, 0.46763),
        VYVPrecomputed(0.8643, 1.45389, 1.61433, 0.83134, 1.87504)
    },
#else
    { // smoothing, L_infinity
        VYVPrecomputed(1.40098, 1.00236, 1.85132),
        VYVPrecomputed(1.12075, 1.27788, 1.76952, 0.46611),
        VYVPrecomputed(0.85480, 1.43749, 1.161231, 0.82053, 1.87415)
    },
#endif
    { // 1st deriv, L_infinity
        VYVPrecomputed(1.31553, 0.97057, 1.77635),
        VYVPrecomputed(1.04185, 1.24034, 1.69747, 0.44790),
        VYVPrecomputed(0.77934, 1.41423, 1.50941, 0.80828, 1.77181)

    },
    { // 2nd deriv, L_infinity
        VYVPrecomputed(1.22886, 0.93058, 1.70493),
        VYVPrecomputed(0.94570, 1.21064, 1.60161, 0.42647),
        VYVPrecomputed(0.69843, 1.37655, 1.42631, 0.77399, 1.69668)
    }
};

} // namespace detail

namespace detail
{

struct deriche_tag {};
struct vyv_tag {};

struct deriche_4_tag : public deriche_tag { static const unsigned int order = 4; };
struct deriche_3_tag : public deriche_tag { static const unsigned int order = 3; };
struct deriche_2_tag : public deriche_tag { static const unsigned int order = 2; };

struct vyv_3_tag : public vyv_tag { static const unsigned int order = 3; };
struct vyv_4_tag : public vyv_tag { static const unsigned int order = 4; };
struct vyv_5_tag : public vyv_tag { static const unsigned int order = 5; };


template<typename X>
struct is_deriche_kernel : public std::is_base_of<deriche_tag, typename X::vigra_recursive_kernel_type>
{
};

template<typename X>
struct is_vyv_kernel : public std::is_base_of<vyv_tag, typename X::vigra_recursive_kernel_type>
{
};


// based on code Copyright (c) 2012-2013, Pascal Getreuer
// <getreuer@cmla.ens-cachan.fr>
// licensed under the terms of the simplified BSD license.
template<typename ARITHTYPE, typename kernel_type, unsigned order>
class RecursiveConvolutionKernelDericheMembers
{
public:
    ARITHTYPE n_causal[order];
    ARITHTYPE n_anticausal[order];
    ARITHTYPE d[order];

    static const unsigned int vigra_recursive_kernel_order = order;
    typedef kernel_type vigra_recursive_kernel_type;

    RecursiveConvolutionKernelDericheMembers() {
        for (unsigned int i = 0; i < vigra_recursive_kernel_order; ++i) {
            d[i] = 0;
            n_causal[i] = 0;
            n_anticausal[i] = 0;
        }
    }

    void scale(double factor) {
        for (unsigned int i = 0; i < vigra_recursive_kernel_order; ++i) {
            n_causal[i] *= factor;
            n_anticausal[i] *= factor;
        }
    }

protected:
    void compute_coefs(unsigned deriv_order, double sigma) {
        std::complex<ARITHTYPE> alpha[order];
        std::complex<ARITHTYPE> lambda[order];
        std::complex<ARITHTYPE> beta[order];

        for (unsigned int i = 0; i < vigra_recursive_kernel_order; ++i) {
            alpha[i] = detail::deriche_precomputed_coefs[deriv_order][vigra_recursive_kernel_order - 2].get_alpha(i);
            lambda[i] = detail::deriche_precomputed_coefs[deriv_order][vigra_recursive_kernel_order - 2].get_lambda(i) / sigma;
            beta[i] = std::complex<ARITHTYPE>(
                -exp(-lambda[i].real()) * cos(lambda[i].imag()),
                exp(-lambda[i].real()) * sin(lambda[i].imag()));
        }

        std::complex<ARITHTYPE> a[vigra_recursive_kernel_order + 1];
        std::complex<ARITHTYPE> b[vigra_recursive_kernel_order];

        b[0] = alpha[0];
        a[0] = std::complex<ARITHTYPE>(1, 0);
        a[1] = beta[0];

        for (unsigned int k = 1; k < vigra_recursive_kernel_order; ++k) {
            b[k] = beta[k] * b[k - 1];

            for (unsigned int j = k - 1; j > 0; --j)
                b[j] += beta[k] * b[j - 1];

            for (unsigned int j = 0; j <= k; ++j)
                b[j] += alpha[k] * a[j];

            a[k + 1] = beta[k] * a[k];

            for (unsigned int j = k; j > 0; --j)
                a[j] += beta[k] * a[j - 1];
        }

        for (unsigned int i = 0; i < vigra_recursive_kernel_order; ++i) {
            n_causal[i] = b[i].real() / (M_SQRT2PI * pow(sigma, deriv_order + 1));
            d[i] = a[i + 1].real();
        }

        // eq. 31 for symmetrical filters and eq. 32 for antisymmetrical
        int sign = 1;
        if (deriv_order == 1)
            sign = -1;

        for (unsigned int i = 0; i < vigra_recursive_kernel_order - 1; ++i)
            n_anticausal[i] = sign * (n_causal[i + 1] - n_causal[0] * d[i]);
        n_anticausal[vigra_recursive_kernel_order - 1] = sign * (-1.0) * n_causal[0] * d[vigra_recursive_kernel_order - 1];
    }
};


// based on code Copyright (c) 2012-2013, Pascal Getreuer
// <getreuer@cmla.ens-cachan.fr>
// licensed under the terms of the simplified BSD license.
template<typename ARITHTYPE, typename kernel_type, unsigned order>
class RecursiveConvolutionKernelVYVMembers
{
public:
    static const unsigned int vigra_recursive_kernel_order = order;
    typedef kernel_type vigra_recursive_kernel_type;

    ARITHTYPE coefs[order + 1];
    Matrix<ARITHTYPE> M;
    double norm_;

    RecursiveConvolutionKernelVYVMembers() {
        for (unsigned int i = 0; i <= order; ++i)
            coefs[i] = 0;
        norm_ = 1.0;
        M = Matrix<ARITHTYPE>(Shape2(order, order));
    }

    void scale(double factor) {
        norm_ *= factor;
    }

protected:
    void compute_coefs(unsigned deriv_order, double sigma) {
        for (unsigned int i = 0; i < order; ++i)
            poles[i] = (std::complex<ARITHTYPE>)detail::vyv_poles[deriv_order][order - 3].get_complex_pole(i);

        ARITHTYPE q = compute_q(sigma);

        for (unsigned int i = 0; i < order; ++i)
            poles[i] = pow(poles[i], 1 / q);

        expand_pole_product();
        compute_border_mtx();
    }

private:
    std::complex<ARITHTYPE> poles[order];

    ARITHTYPE compute_q(ARITHTYPE sigma) {
        ARITHTYPE q = sigma / 2;
        ARITHTYPE sigma2 = sigma * sigma;

        for (unsigned int i = 0; i < 20; ++i)
            q -= (variance(q) - sigma2) / dq_variance(q);
        return q;
    }

    ARITHTYPE variance(ARITHTYPE q) {
        const std::complex<ARITHTYPE> one(1, 0);
        std::complex<ARITHTYPE> sum(0, 0);
        std::complex<ARITHTYPE> z;

        for (unsigned int i = 0; i < order; ++i) {
            z = pow(poles[i], 1 / q);
            sum += z / (std::complex<ARITHTYPE>)pow(z - one, 2);
        }

        return 2.0 * sum.real();
    }

    ARITHTYPE dq_variance(ARITHTYPE q) {
        const std::complex<ARITHTYPE> one(1, 0);
        std::complex<ARITHTYPE> sum(0, 0);

        for (unsigned int i = 0; i < order; ++i) {
            std::complex<ARITHTYPE> z = pow(poles[i], 1 / q);

            sum += z * log(z) * (z + one) / (std::complex<ARITHTYPE>)pow(z - one, 3);
        }

        return (2.0 / q) * sum.real();
    }

    void expand_pole_product() {
        std::complex<ARITHTYPE> denom[order + 1];

        denom[0] = poles[0];
        denom[1] = std::complex<ARITHTYPE>(-1, 0);

        for (unsigned int i = 1; i < order; ++i) {
            denom[i + 1] = std::complex<ARITHTYPE>(-1, 0) * denom[i];

            for (unsigned int j = i; j > 0; --j)
                denom[j] = (denom[j] * poles[i]) - denom[j - 1];

            denom[0] *= poles[i];
        }

        for (unsigned int i = 1; i <= order; ++i)
            coefs[i] = (denom[i] / denom[0]).real();

        coefs[0] = 1;
        for (unsigned int i = 1; i <= order; ++i)
            coefs[0] += coefs[i];
    }

    void compute_border_mtx() {
        Matrix<ARITHTYPE> A(Shape2(order, order));
        Matrix<ARITHTYPE> Apow(Shape2(order, order));
        Matrix<ARITHTYPE> Asum(Shape2(order, order));
        Matrix<ARITHTYPE> Asuminv(Shape2(order, order));

        A.init(0);
        identityMatrix(Asum);

        for (unsigned int i = 0; i < order; ++i)
            A(0,i) = coefs[1+i];

        for (unsigned int i = 0; i < order - 1; ++i)
            A(1+i,i) = 1.0;

        Apow = A;
        for (unsigned int i = 0; i < order; ++i) {
            Asum -= coefs[1+i] * Apow;
            Apow = Apow * A;
        }

        linalg::inverse(Asum, Asuminv);        
    
        M.rowVector(0) = Asuminv.rowVector(0);
        Apow = A;
        for (unsigned int i = 1; i < order; ++i){
            M.rowVector(i) = M.rowVector(0) * Apow;
            Apow = Apow * A;
        }
    }
};



template <typename SrcIterator, typename SrcAccessor>
class IIRBorderSrcAccessorZero
{
public:
    typedef typename SrcAccessor::value_type value_type;

    SrcAccessor sa_;
    SrcIterator is_;
    SrcIterator iend_;
    int start_, stop_;

    IIRBorderSrcAccessorZero(SrcIterator is, SrcIterator iend, SrcAccessor sa, int start, int stop) : sa_(sa), is_(is), iend_(iend), start_(start), stop_(stop) { }

    inline value_type operator()( SrcIterator i )
    {
            return NumericTraits<value_type>::zero();
    }
};

template <typename SrcIterator, typename SrcAccessor>
class IIRBorderSrcAccessorRepeatLeft
{
public:
    typedef typename SrcAccessor::value_type value_type;

    SrcAccessor sa_;
    SrcIterator is_;
    SrcIterator iend_;
    int start_, stop_;

    IIRBorderSrcAccessorRepeatLeft(SrcIterator is, SrcIterator iend, SrcAccessor sa, int start, int stop) : sa_(sa), is_(is), iend_(iend), start_(start), stop_(stop) { }

    inline value_type operator()( SrcIterator i )
    {
            return this->sa_(this->is_ + this->start_);
    }
};

template <typename SrcIterator, typename SrcAccessor>
class IIRBorderSrcAccessorRepeatRight
{
public:
    typedef typename SrcAccessor::value_type value_type;

    SrcAccessor sa_;
    SrcIterator is_;
    SrcIterator iend_;
    int start_, stop_;

    IIRBorderSrcAccessorRepeatRight(SrcIterator is, SrcIterator iend, SrcAccessor sa, int start, int stop) : sa_(sa), is_(is), iend_(iend), start_(start), stop_(stop) { }

    inline value_type operator()( SrcIterator i )
    {
            return this->sa_(this->is_ + this->stop_ - 1);
    }
};

template <typename SrcIterator, typename SrcAccessor>
class IIRBorderSrcAccessorReflectLeft
{
public:
    typedef typename SrcAccessor::value_type value_type;

    SrcAccessor sa_;
    SrcIterator is_;
    SrcIterator iend_;
    int start_, stop_;

    IIRBorderSrcAccessorReflectLeft(SrcIterator is, SrcIterator iend, SrcAccessor sa, int start, int stop) : sa_(sa), is_(is), iend_(iend), start_(start), stop_(stop) { }

    inline value_type operator()( SrcIterator i )
    {
            return this->sa_(this->is_ + this->start_ + std::distance(i, this->is_) + this->start_);
    }
};

template <typename SrcIterator, typename SrcAccessor>
class IIRBorderSrcAccessorReflectRight
{
public:
    typedef typename SrcAccessor::value_type value_type;

    SrcAccessor sa_;
    SrcIterator is_;
    SrcIterator iend_;
    int start_, stop_;

    IIRBorderSrcAccessorReflectRight(SrcIterator is, SrcIterator iend, SrcAccessor sa, int start, int stop) : sa_(sa), is_(is), iend_(iend), start_(start), stop_(stop) { }

    inline value_type operator()( SrcIterator i )
    {
            return this->sa_(this->iend_ - std::distance(this->iend_, i) - 2);
    }
};

template <typename SrcIterator, typename SrcAccessor>
class IIRBorderSrcAccessorWrapLeft
{
public:
    typedef typename SrcAccessor::value_type value_type;

    SrcAccessor sa_;
    SrcIterator is_;
    SrcIterator iend_;
    int start_, stop_;

    IIRBorderSrcAccessorWrapLeft(SrcIterator is, SrcIterator iend, SrcAccessor sa, int start, int stop) : sa_(sa), is_(is), iend_(iend), start_(start), stop_(stop) { }

    inline value_type operator()( SrcIterator i )
    {
            return this->sa_(this->iend_ - std::distance(i, this->is_) + this->start_);
    }
};

template <typename SrcIterator, typename SrcAccessor>
class IIRBorderSrcAccessorWrapRight
{
public:
    typedef typename SrcAccessor::value_type value_type;

    SrcAccessor sa_;
    SrcIterator is_;
    SrcIterator iend_;
    int start_, stop_;

    IIRBorderSrcAccessorWrapRight(SrcIterator is, SrcIterator iend, SrcAccessor sa, int start, int stop) : sa_(sa), is_(is), iend_(iend), start_(start), stop_(stop) { }

    inline value_type operator()( SrcIterator i )
    {
            return this->sa_(this->is_ + this->start_ + std::distance(this->iend_, i));
    }
};


template <typename DestIterator, typename DestAccessor>
class DericheBorderDstAccessor
{
public:
    typedef typename DestAccessor::value_type value_type;

    inline value_type operator()( DestIterator i )
    {
        return NumericTraits<value_type>::zero();
    }

    inline void set(value_type dummy, DestIterator dummy2) { }
};


template <typename DestIterator, typename DestAccessor>
class VYVBorderDstAccessor
{
public:
    typedef typename DestAccessor::value_type value_type;
    value_type tmp;

    inline value_type operator()( DestIterator i )
    {
        return tmp;
    }

    inline void set(value_type value, DestIterator dummy2)
    {
        tmp = value;
    }
};

template <class SrcIterator, class SrcAccessor, class DestIterator,
          class DestAccessor, class RecursiveConvolutionKernel, class SumType>
void dericheApplyCausal(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                        DestIterator id, DestAccessor da,
                        RecursiveConvolutionKernel kernel, SumType xtmp[],
                        SumType ytmp[], int start, int stop) {
    SumType xi, dyi;

    for (int x = start; x < stop; ++x) {
        xi = sa(is + x);

        xtmp[0] = xi;

        dyi = NumericTraits<SumType>::zero();
        for (unsigned int i = 0; i < kernel.order; ++i)
            dyi += kernel.n_causal[i] * xtmp[i];
        for (unsigned int i = 0; i < kernel.order; ++i)
            dyi -= kernel.d[i] * ytmp[i];

        for (unsigned int i = kernel.order - 1; i > 0; --i) {
            xtmp[i] = xtmp[i - 1];
            ytmp[i] = ytmp[i - 1];
        }

        ytmp[0] = dyi;
        da.set(detail::RequiresExplicitCast<
                   typename DestAccessor::value_type>::cast(dyi),
               id + x);
    }
}

template <class SrcIterator, class SrcAccessor, class DestIterator,
          class DestAccessor, class RecursiveConvolutionKernel, class SumType>
void dericheApplyAntiCausal(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                            DestIterator id, DestAccessor da,
                            RecursiveConvolutionKernel kernel, SumType xtmp[],
                            SumType ytmp[], int start, int stop) {
    SumType xi, yi, dyi;

    for (int x = stop - 1; x >= start; --x) {
        xi = sa(is + x);
        yi = da(id + x);

        dyi = NumericTraits<SumType>::zero();
        for (unsigned int i = 0; i < kernel.order; ++i)
            dyi += kernel.n_anticausal[i] * xtmp[i];
        for (unsigned int i = 0; i < kernel.order; ++i)
            dyi -= kernel.d[i] * ytmp[i];

        for (unsigned int i = kernel.order - 1; i > 0; --i) {
            xtmp[i] = xtmp[i - 1];
            ytmp[i] = ytmp[i - 1];
        }

        ytmp[0] = dyi;
        xtmp[0] = xi;

        da.set(detail::RequiresExplicitCast<
                   typename DestAccessor::value_type>::cast(yi + dyi),
               id + x);
    }
}

template <class BorderTreatmentLeft, class BorderTreatmentRight,
          class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel>
inline void recursiveConvolveLineDericheBorder(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  DestIterator id, DestAccessor da,
                  RecursiveConvolutionKernel kernel,
                  int start = 0, int stop = 0)
{
    typedef typename PromoteTraits<
        typename SrcAccessor::value_type,
        typename RecursiveConvolutionKernel::value_type>::Promote SumType;
    SumType xtmp[kernel.order];
    SumType ytmp[kernel.order];

    BorderTreatmentLeft sa_causal = BorderTreatmentLeft(is, iend, sa, start, stop);
    BorderTreatmentRight sa_anticausal = BorderTreatmentRight(is, iend, sa, start, stop);
    auto da_fake = DericheBorderDstAccessor<DestIterator, DestAccessor>();

    for (unsigned int i = 0; i < kernel.order; ++i)
        xtmp[i] = ytmp[i] = NumericTraits<SumType>::zero();
    dericheApplyCausal(is, iend, sa_causal, id, da_fake, kernel, xtmp, ytmp, start + kernel.left(), start);
    dericheApplyCausal(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);

    for (unsigned int i = 0; i < kernel.order; ++i)
        xtmp[i] = ytmp[i] = NumericTraits<SumType>::zero();
    dericheApplyAntiCausal(is, iend, sa_anticausal, id, da_fake, kernel, xtmp, ytmp, stop, stop + kernel.right());
    dericheApplyAntiCausal(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);

}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel>
inline void recursiveConvolveLineDericheZeroBorder(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  DestIterator id, DestAccessor da,
                  RecursiveConvolutionKernel kernel,
                  int start = 0, int stop = 0)
{
    typedef typename PromoteTraits<
        typename SrcAccessor::value_type,
        typename RecursiveConvolutionKernel::value_type>::Promote SumType;
    SumType xtmp[kernel.order];
    SumType ytmp[kernel.order];

    for (unsigned int i = 0; i < kernel.order; ++i)
        xtmp[i] = ytmp[i] = NumericTraits<SumType>::zero();
    dericheApplyCausal(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);

    for (unsigned int i = 0; i < kernel.order; ++i)
        xtmp[i] = ytmp[i] = NumericTraits<SumType>::zero();
    dericheApplyAntiCausal(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);

}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel>
inline typename std::enable_if<is_deriche_kernel<RecursiveConvolutionKernel>::value>::type
recursiveConvolveLine(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  DestIterator id, DestAccessor da,
                  RecursiveConvolutionKernel kernel,
                  BorderTreatmentMode border,
                  int start = 0, int stop = 0)
{ 
    switch(kernel.border_treatment) {
        case BORDER_TREATMENT_ZEROPAD:
            recursiveConvolveLineDericheZeroBorder(is, iend, sa, id, da, kernel, start, stop);
            break;

        case BORDER_TREATMENT_REPEAT:
            recursiveConvolveLineDericheBorder<detail::IIRBorderSrcAccessorRepeatLeft<SrcIterator, SrcAccessor>, detail::IIRBorderSrcAccessorRepeatRight<SrcIterator, SrcAccessor>>(is, iend, sa, id, da, kernel, start, stop);
            break;

        case BORDER_TREATMENT_REFLECT:
            recursiveConvolveLineDericheBorder<detail::IIRBorderSrcAccessorReflectLeft<SrcIterator, SrcAccessor>, detail::IIRBorderSrcAccessorReflectRight<SrcIterator, SrcAccessor>>(is, iend, sa, id, da, kernel, start, stop);
            break;

        case BORDER_TREATMENT_WRAP:
            recursiveConvolveLineDericheBorder<detail::IIRBorderSrcAccessorWrapLeft<SrcIterator, SrcAccessor>, detail::IIRBorderSrcAccessorWrapRight<SrcIterator, SrcAccessor>>(is, iend, sa, id, da, kernel, start, stop);
            break;

        case BORDER_TREATMENT_AVOID:
        case BORDER_TREATMENT_CLIP:
        default:
            vigra_precondition(false, "recursiveConvolveLine: unsupported border treatment for Deriche filters.");
            break;
    }
}



template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel,
          class SumType>
void vyvApplyCausal(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  DestIterator id, DestAccessor da,
                  RecursiveConvolutionKernel kernel,
                  SumType xtmp, SumType ytmp[],
                  int start, int stop)
{
    SumType dyi;
    for (int x = start; x < stop; ++x) {

        if (kernel.deriv_order == 1)
            dyi = (sa(is + x + 1) - xtmp) * kernel.coefs[0] / 2.0;
        else if (kernel.deriv_order == 2)
            dyi = (sa(is + x) - xtmp) * kernel.coefs[0];
        else
            dyi = sa(is + x) * kernel.coefs[0];

        if (kernel.deriv_order > 0)
            xtmp = sa(is + x);

        for (unsigned int i = 0; i < kernel.order; ++i)
            dyi -= kernel.coefs[i + 1] * ytmp[i];

        for (unsigned int i = kernel.order - 1; i > 0; --i)
            ytmp[i] = ytmp[i - 1];

        ytmp[0] = dyi;
        da.set(detail::RequiresExplicitCast<
                   typename DestAccessor::value_type>::cast(dyi),
               id + x);
    }
}

template <class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel,
          class SumType>
void vyvApplyAntiCausal(DestIterator id, DestAccessor da,
                  RecursiveConvolutionKernel kernel,
                  SumType xtmp, SumType ytmp[],
                  int start, int stop)
{
    SumType dyi;
    for (int x = stop - 1; x >= start; --x) {
        if (kernel.deriv_order == 2) {
            dyi = (xtmp - da(id + x)) * kernel.coefs[0];
            xtmp = da(id + x);
        } else
            dyi = da(id + x) * kernel.coefs[0];

        for (unsigned int i = 0; i < kernel.order; ++i)
            dyi -= kernel.coefs[i + 1] * ytmp[i];

        for (unsigned int i = kernel.order - 1; i > 0; --i)
            ytmp[i] = ytmp[i - 1];
        ytmp[0] = dyi;

        da.set(detail::RequiresExplicitCast<
                   typename DestAccessor::value_type>::cast(dyi * kernel.norm_),
               id + x);
    }
}




template <class BorderTreatmentLeft, class BorderTreatmentRight,
          class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel,
          class SumType>
inline void recursiveConvolveLineVYVBorder(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  DestIterator id, DestAccessor da,
                  RecursiveConvolutionKernel kernel,
                  SumType xtmp, SumType ytmp[],
                  int start = 0, int stop = 0)
{
    BorderTreatmentLeft sa_left = BorderTreatmentLeft(is, iend, sa, start, stop);
    BorderTreatmentRight sa_right = BorderTreatmentRight(is, iend, sa, start, stop);
    auto da_fake = VYVBorderDstAccessor<DestIterator, DestAccessor>();


    typedef typename AccessorTraits<SumType>::default_accessor SumAccessor;
    SumAccessor va;
    ArrayVector<SumType> v;

    v.reserve(kernel.right());

    for (unsigned int i = 0; i < kernel.order; ++i)
        ytmp[i] = NumericTraits<SumType>::zero();
    xtmp = NumericTraits<SumType>::zero();

    // causal filter
    vyvApplyCausal(is, iend, sa_left, id, da_fake, kernel, xtmp, ytmp, start + kernel.left(), start);
    vyvApplyCausal(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);

    // handle right border
    vyvApplyCausal(is + stop, iend, sa_right, v.begin(), va, kernel, xtmp, ytmp, 0, kernel.right());

    TinyVector<SumType, kernel.order> boundary;
    for (unsigned int i = 0; i < kernel.order; ++i)
        boundary[i] = ytmp[i];

    for (unsigned int i = 0; i < kernel.order; ++i) {
        ytmp[i] = NumericTraits<SumType>::zero();
        for (unsigned int j = 0; j < kernel.order; ++j)
            ytmp[i] += kernel.M(i,j) * boundary[j];
    }

    xtmp = NumericTraits<SumType>::zero();

    vyvApplyAntiCausal(v.begin(), va, kernel, xtmp, ytmp, 0, kernel.right() - kernel.order);

    // anticausal filter
    xtmp = NumericTraits<SumType>::zero();
    vyvApplyAntiCausal(id, da, kernel, xtmp, ytmp, start, stop);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel>
inline typename std::enable_if<is_vyv_kernel<RecursiveConvolutionKernel>::value>::type
recursiveConvolveLine(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  DestIterator id, DestAccessor da,
                  RecursiveConvolutionKernel kernel,
                  BorderTreatmentMode border,
                  int start = 0, int stop = 0)
{
    typedef typename PromoteTraits<
        typename SrcAccessor::value_type,
        typename RecursiveConvolutionKernel::value_type>::Promote SumType;

    SumType ytmp[kernel.order];
    SumType xtmp = NumericTraits<SumType>::zero();

    switch(kernel.border_treatment) {
        case BORDER_TREATMENT_ZEROPAD:
            recursiveConvolveLineVYVBorder<detail::IIRBorderSrcAccessorZero<SrcIterator, SrcAccessor>, detail::IIRBorderSrcAccessorZero<SrcIterator, SrcAccessor>>(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
            break;

        case BORDER_TREATMENT_REPEAT:
            recursiveConvolveLineVYVBorder<detail::IIRBorderSrcAccessorRepeatLeft<SrcIterator, SrcAccessor>, detail::IIRBorderSrcAccessorRepeatRight<SrcIterator, SrcAccessor>>(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
            break;

        case BORDER_TREATMENT_REFLECT:
            recursiveConvolveLineVYVBorder<detail::IIRBorderSrcAccessorReflectLeft<SrcIterator, SrcAccessor>, detail::IIRBorderSrcAccessorReflectRight<SrcIterator, SrcAccessor>>(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
            break;

        case BORDER_TREATMENT_WRAP:
            recursiveConvolveLineVYVBorder<detail::IIRBorderSrcAccessorWrapLeft<SrcIterator, SrcAccessor>, detail::IIRBorderSrcAccessorWrapRight<SrcIterator, SrcAccessor>>(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
            break;

        case BORDER_TREATMENT_AVOID:
        case BORDER_TREATMENT_CLIP:
        default:
            vigra_precondition(false, "recursiveConvolveLine: unsupported border treatment for VYV filters.");
            break;
    }
}



template<typename ARITHTYPE, typename kernel_type>
class RecursiveConvolutionKernelMembers {
};

template<typename ARITHTYPE>
class RecursiveConvolutionKernelMembers<ARITHTYPE, deriche_2_tag> : public RecursiveConvolutionKernelDericheMembers<ARITHTYPE, deriche_2_tag, 2>
{
};

template<typename ARITHTYPE>
class RecursiveConvolutionKernelMembers<ARITHTYPE, deriche_3_tag> : public RecursiveConvolutionKernelDericheMembers<ARITHTYPE, deriche_3_tag, 3>
{
};

template<typename ARITHTYPE>
class RecursiveConvolutionKernelMembers<ARITHTYPE, deriche_4_tag> : public RecursiveConvolutionKernelDericheMembers<ARITHTYPE, deriche_4_tag, 4>
{
};

template<typename ARITHTYPE>
class RecursiveConvolutionKernelMembers<ARITHTYPE, vyv_3_tag> : public RecursiveConvolutionKernelVYVMembers<ARITHTYPE, vyv_3_tag, 3>
{
};

template<typename ARITHTYPE>
class RecursiveConvolutionKernelMembers<ARITHTYPE, vyv_4_tag> : public RecursiveConvolutionKernelVYVMembers<ARITHTYPE, vyv_4_tag, 4>
{
};

template<typename ARITHTYPE>
class RecursiveConvolutionKernelMembers<ARITHTYPE, vyv_5_tag> : public RecursiveConvolutionKernelVYVMembers<ARITHTYPE, vyv_5_tag, 5>
{
};

} // namespace detail

template <class ARITHTYPE, typename kernel_tag>
class RecursiveConvolutionKernel : public detail::RecursiveConvolutionKernelMembers<ARITHTYPE, kernel_tag>
{
public:
    static const unsigned int order = kernel_tag::order;

    typedef detail::iir_kernel1d_tag vigra_kernel_category;

    typedef ARITHTYPE value_type;

    unsigned int deriv_order;
    double kernel_radius;
    BorderTreatmentMode border_treatment;

    RecursiveConvolutionKernel() : deriv_order(0), kernel_radius(0.0), border_treatment(BORDER_TREATMENT_REFLECT) {};

    void initGaussian(double std_dev, value_type norm, double windowRatio = 0.0)
    {
        initGaussianDerivative(std_dev, 0, norm, windowRatio);
    }

    void initGaussian(double std_dev)
    {
        initGaussian(std_dev, 1.0);
    }

    void initGaussianDerivative(double std_dev, int order)
    {
        initGaussianDerivative(std_dev, order, 1.0);
    }

    void initGaussianDerivative(double std_dev, int order, value_type norm, double windowRatio = 0.0)
    {
        vigra_precondition(
            order >= 0 && order <= 2,
            "RecursiveConvolutionKernel::initGaussianDerivative(): order must be between 0 and 2.");

        deriv_order = order;
        this->compute_coefs(order, std_dev);
        this->scale((double)norm);

        if (windowRatio == 0.0)
            windowRatio = 3.0 + 0.5*order;

        kernel_radius = windowRatio * std_dev;
    }

    int left() {
        return rational_cast<int>(-kernel_radius);
    }

    int right() {
        return rational_cast<int>(kernel_radius);
    }


        /** current border treatment mode
        */
    BorderTreatmentMode borderTreatment() const
    { return border_treatment; }

        /** Set border treatment mode.
        */
    void setBorderTreatment( BorderTreatmentMode new_mode)
    { border_treatment = new_mode; }
};

template<typename RecursiveConvolutionKernel>
inline typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value, RecursiveConvolutionKernel>::type
kernel1d(RecursiveConvolutionKernel k)
{
    return k;
}


template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel>
typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value>::type
convolveLine(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  DestIterator id, DestAccessor da,
                  RecursiveConvolutionKernel ik,
                  BorderTreatmentMode border,
                  int start = 0, int stop = 0)
{
    if (stop == 0)
        stop = start + std::distance(is, iend);

    detail::recursiveConvolveLine(is, iend, sa, id, da, ik, border, start, stop);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel>
inline typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value>::type 
convolveLine(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  DestIterator id, DestAccessor da,
                  RecursiveConvolutionKernel ik,
                  int start = 0, int stop = 0)
{
    convolveLine(is, iend, sa, id, da, ik, ik.border_treatment, start, stop);
}


template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel>
inline typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value>::type
convolveLine(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIterator, DestAccessor> dest,
                  RecursiveConvolutionKernel kernel,
                  int start = 0, int stop = 0)
{
    convolveLine(src.first, src.second, src.third, dest.first, dest.second, kernel, kernel.border_treatment, start, stop);
}


} // namespace vigra

#endif
// VIGRA_RECURSIVE_MULTI_CONVOLUTION_H
