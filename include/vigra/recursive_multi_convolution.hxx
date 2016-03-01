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
#include "linear_solve.hxx"

#include <xmmintrin.h>
#include <emmintrin.h>
#include <complex>
#include <type_traits>

namespace vigra {

namespace detail {
// Rachid Deriche. Recursively implementating the Gaussian and its derivatives.
// [Research Report]
// RR-1893, 1993, pp.24. <inria-00074778>

#define DERICHE_LAMBDA_FROM_BW(b0, w0) std::complex<double>(b0, w0)
#define DERICHE_LAMBDA_CONJ_FROM_BW(b0, w0) std::complex<double>(b0, -w0)
#define DERICHE_ALPHA_FROM_AC(ac0, ac1)                                        \
    std::complex<double>(ac0 / 2.0, ac1 / 2.0)
#define DERICHE_ALPHA_CONJ_FROM_AC(ac0, ac1)                                   \
    std::complex<double>(ac0 / 2.0, -ac1 / 2.0)

#define DERICHE_LAMBDAS_FROM_BW(b0, w0)                                        \
    DERICHE_LAMBDA_FROM_BW(b0, w0), DERICHE_LAMBDA_CONJ_FROM_BW(b0, w0)
#define DERICHE_ALPHAS_FROM_AC(ac0, ac1)                                       \
    DERICHE_ALPHA_FROM_AC(ac0, ac1), DERICHE_ALPHA_CONJ_FROM_AC(ac0, ac1)

typedef struct deriche_precomputed {
    const std::complex<double> alpha[4];
    const std::complex<double> lambda[4];
} deriche_precomputed;

// eq. 35 - 46
// 4th order coefficients from The Mines Java Toolkit
// https://github.com/dhale/jtk/blob/master/core/src/main/java/edu/mines/jtk/dsp/RecursiveGaussianFilter.java
static const deriche_precomputed deriche_precomputed_filters[3][3] = {
    // smoothing
    {// 2nd order
     {.alpha = {DERICHE_ALPHAS_FROM_AC(0.9629, 1.942),
                DERICHE_ALPHAS_FROM_AC(0, 0)},
      .lambda = {DERICHE_LAMBDAS_FROM_BW(1.26, 0.8448),
                 DERICHE_LAMBDAS_FROM_BW(0, 0)}},
     // 3rd order
     {.alpha = {DERICHE_ALPHAS_FROM_AC(-0.8929, 1.021),
                std::complex<double>(1.898, 0), std::complex<double>(0, 0)},
      .lambda = {DERICHE_LAMBDAS_FROM_BW(1.512, 1.475),
                 std::complex<double>(1.556, 0), std::complex<double>(0, 0)}},
     // 4th order
     {.alpha = {
          DERICHE_ALPHAS_FROM_AC(1.6797292232361107, 3.7348298269103580),
          DERICHE_ALPHAS_FROM_AC(-0.6802783501806897, -0.2598300478959625)},
      .lambda = {
          DERICHE_LAMBDAS_FROM_BW(1.7831906544515104, 0.6318113174569493),
          DERICHE_LAMBDAS_FROM_BW(1.7228297663338028, 1.9969276832487770)}}},
    // first derivative
    {// 2nd order
     {.alpha = {DERICHE_ALPHAS_FROM_AC(0.1051, -1.898),
                DERICHE_ALPHAS_FROM_AC(0, 0)},
      .lambda = {DERICHE_LAMBDAS_FROM_BW(0.9338, 0.9459),
                 DERICHE_LAMBDAS_FROM_BW(0, 0)}},
     // 3rd order
     {.alpha = {DERICHE_ALPHAS_FROM_AC(1.666, -0.4557),
                std::complex<double>(-1.682, 0), std::complex<double>(0, 0)},

      .lambda = {DERICHE_LAMBDAS_FROM_BW(1.266, 1.558),
                 std::complex<double>(1.244, 0), std::complex<double>(0, 0)}},
     // 4th order
     {.alpha = {
          DERICHE_ALPHAS_FROM_AC(0.6494024008440620, 0.9557370760729773),
          DERICHE_ALPHAS_FROM_AC(-0.6472105276644291, -4.5306923044570760)},
      .lambda = {
          DERICHE_LAMBDAS_FROM_BW(1.5159726670750566, 2.0718953658782650),
          DERICHE_LAMBDAS_FROM_BW(1.5267608734791140, 0.6719055957689513)}}},
    // second derivative
    {// 2nd order
     {.alpha = {DERICHE_ALPHAS_FROM_AC(-1.228, 0.2976),
                DERICHE_ALPHAS_FROM_AC(0, 0)},
      .lambda = {DERICHE_LAMBDAS_FROM_BW(0.6134, 1.332),
                 DERICHE_LAMBDAS_FROM_BW(0, 0)}},
     // 3rd order
     {.alpha = {DERICHE_ALPHAS_FROM_AC(-1.773, -1.129),
                std::complex<double>(0.8309, 0), std::complex<double>(0, 0)},

      .lambda = {DERICHE_LAMBDAS_FROM_BW(1.032, 1.664),
                 std::complex<double>(0.8459, 0), std::complex<double>(0, 0)}},
     // 4th order
     {.alpha = {
          DERICHE_ALPHAS_FROM_AC(0.3224570510072559, -1.7382843963561239),
          DERICHE_ALPHAS_FROM_AC(-1.3312275593739595, 3.6607035671974897)},
      .lambda = {
          DERICHE_LAMBDAS_FROM_BW(1.3138054926516880, 2.1656041357418863),
          DERICHE_LAMBDAS_FROM_BW(1.2402181393295362, 0.7479888745408682)}}}};

// Van Vliet, Lucas J., Ian T. Young, and Piet W. Verbeek. "Recursive Gaussian
// derivative filters." Pattern Recognition, 1998. Proceedings. Fourteenth
// International Conference on. Vol. 1. IEEE, 1998.
// Table 1 and 2
static const double vyv_poles[3][3][5] = {
#if 1
    // smoothing, L_2
    {// 3rd order
     {1.4165, 1.00829, 1.86543, 0, 0},
     // 4th order
     {1.13228, 1.28114, 1.78534, 0.46763, 0},
     // 5th order
     {0.8643, 1.45389, 1.61433, 0.83134, 1.87504}},
#else
    // smoothing, L_infinity
    {// 3rd order
     {1.40098, 1.00236, 1.85132, 0, 0},
     // 4th order
     {1.12075, 1.27788, 1.76952, 0.46611, 0},
     // 5th order
     {0.85480, 1.43749, 1.161231, 0.82053, 1.87415}},
#endif
    // 1st deriv, L_infinity
    {// 3rd order
     {1.31553, 0.97057, 1.77635, 0, 0},
     // 4th order
     {1.04185, 1.24034, 1.69747, 0.44790},
     // 5th order
     {0.77934, 1.41423, 1.50941, 0.80828, 1.77181}},
    // 2nd deriv, L_infinity
    {// 3rd order
     {1.22886, 0.93058, 1.70493, 0, 0},
     // 4th order
     {0.94570, 1.21064, 1.60161, 0.42647},
     // 5th order
     {0.69843, 1.37655, 1.42631, 0.77399, 1.69668}}};

} // namespace detail

namespace detail
{

struct deriche_tag {};
struct vyv_tag {};

template<typename X, typename Y>
struct is_specific_iir_kernel
{
  static const bool value = std::is_same<typename X::vigra_recursive_kernel_type, Y>::value;
};

template<typename X>
using is_deriche_kernel = is_specific_iir_kernel<X, deriche_tag>;
template<typename X>
using is_vyv_kernel = is_specific_iir_kernel<X, vyv_tag>;


// based on code Copyright (c) 2012-2013, Pascal Getreuer
// <getreuer@cmla.ens-cachan.fr>
// licensed under the terms of the simplified BSD license.
template<unsigned order, typename ARITHTYPE>
class RecursiveConvolutionKernelDericheMembers
{
public:
    typedef deriche_tag vigra_recursive_kernel_type;

    ARITHTYPE n_causal[order];
    ARITHTYPE n_anticausal[order];
    ARITHTYPE d[order];

    RecursiveConvolutionKernelDericheMembers() {
        for (unsigned int i = 0; i < order; ++i) {
            d[i] = 0;
            n_causal[i] = 0;
            n_anticausal[i] = 0;
        }
    }

    void scale(double factor) {
        for (unsigned int i = 0; i < order; ++i) {
            n_causal[i] *= factor;
            n_anticausal[i] *= factor;
        }
    }

protected:
    void compute_coefs(unsigned deriv_order, double sigma) {
        std::complex<ARITHTYPE> alpha[order];
        std::complex<ARITHTYPE> lambda[order];
        std::complex<ARITHTYPE> beta[order];

        for (unsigned int i = 0; i < order; ++i) {
            alpha[i] = detail::deriche_precomputed_filters[deriv_order][order - 2]
                           .alpha[i];
            lambda[i] = detail::deriche_precomputed_filters[deriv_order][order - 2]
                            .lambda[i];
        }

        for (unsigned int i = 0; i < order; ++i) {
            lambda[i] /= sigma;
            beta[i] = std::complex<ARITHTYPE>(
                -exp(-lambda[i].real()) * cos(lambda[i].imag()),
                exp(-lambda[i].real()) * sin(lambda[i].imag()));
        }

        std::complex<ARITHTYPE> a[order + 1];
        std::complex<ARITHTYPE> b[order];

        b[0] = alpha[0];
        a[0] = std::complex<ARITHTYPE>(1, 0);
        a[1] = beta[0];

        for (unsigned int k = 1; k < order; ++k) {
            b[k] = beta[k] * b[k - 1];

            for (unsigned int j = k - 1; j > 0; --j)
                b[j] += beta[k] * b[j - 1];

            for (unsigned int j = 0; j <= k; ++j)
                b[j] += alpha[k] * a[j];

            a[k + 1] = beta[k] * a[k];

            for (unsigned int j = k; j > 0; --j)
                a[j] += beta[k] * a[j - 1];
        }

        for (unsigned int i = 0; i < order; ++i) {
            n_causal[i] = b[i].real() / (M_SQRT2PI * pow(sigma, deriv_order + 1));
            d[i] = a[i + 1].real();
        }

        // eq. 31 for symmetrical filters and eq. 32 for antisymmetrical
        int sign = 1;
        if (deriv_order == 1)
            sign = -1;

        for (unsigned int i = 0; i < order - 1; ++i)
            n_anticausal[i] = sign * (n_causal[i + 1] - n_causal[0] * d[i]);
        n_anticausal[order - 1] = sign * (-1.0) * n_causal[0] * d[order - 1];
    }
};

// based on code Copyright (c) 2012-2013, Pascal Getreuer
// <getreuer@cmla.ens-cachan.fr>
// licensed under the terms of the simplified BSD license.
template<unsigned order, typename ARITHTYPE>
class RecursiveConvolutionKernelVYVMembers
{
public:
    typedef vyv_tag vigra_recursive_kernel_type;

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
        for (unsigned int i = 0; i < order; ++i) {
            unsigned int offset = i & ~1;
            ARITHTYPE real;
            ARITHTYPE imag;

            real = detail::vyv_poles[deriv_order][order - 3][offset];

            if (offset + 1 == order)
                imag = 0;
            else
                imag = detail::vyv_poles[deriv_order][order - 3][offset + 1];

            if (i % 2)
                imag = -imag;

            poles[i] = std::complex<ARITHTYPE>(real, imag);
        }

        ARITHTYPE q;

        q = compute_q(sigma);

        for (unsigned int i = 0; i < order; ++i)
            poles[i] = pow(poles[i], 1 / q);

        expand_pole_product();
        compute_border_mtx();
    }

private:
    std::complex<ARITHTYPE> poles[order];

    double compute_q(double sigma) {
        double q = sigma / 2;
        double sigma2 = sigma * sigma;

        for (unsigned int i = 0; i < 20; ++i)
            q -= (variance(q) - sigma2) / dq_variance(q);
        return q;
    }

    double variance(double q) {
        const std::complex<ARITHTYPE> one(1, 0);
        std::complex<ARITHTYPE> sum(0, 0);
        std::complex<ARITHTYPE> z;

        for (unsigned int i = 0; i < order; ++i) {
            z = pow(poles[i], 1 / q);
            sum += z / (pow(z - one, 2));
        }

        return 2.0 * sum.real();
    }

    double dq_variance(double q) {
        const std::complex<ARITHTYPE> one(1, 0);
        std::complex<ARITHTYPE> sum(0, 0);
        std::complex<ARITHTYPE> z;

        for (unsigned int i = 0; i < order; ++i) {
            z = pow(poles[i], 1 / q);

            sum += z * log(z) * (z + one) / (pow(z - one, 3));
        }

        return (2.0 / q) * sum.real();
    }

    void expand_pole_product(void) {
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

    void compute_border_mtx(void) {
        Matrix<ARITHTYPE> A(Shape2(order, order));
        Matrix<ARITHTYPE> Asum(Shape2(order, order));
        Matrix<ARITHTYPE> Asuminv(Shape2(order, order));

        A.init(0);
        identityMatrix(Asum);

        for (unsigned int i = 0; i < order; ++i)
            A(0,i) = coefs[1+i];

        for (unsigned int i = 0; i < order - 1; ++i)
            A(1+i,i) = 1.0;

        for (unsigned int i = 0; i < order; ++i)
            Asum -= coefs[1+i] * pow(A, 1+i);

        linalg::inverse(Asum, Asuminv);        
    
        M.rowVector(0) = Asuminv.rowVector(0);
        for (unsigned int i = 1; i < order; ++i)
            M.rowVector(i) = M.rowVector(0) * pow(A, i);

        for (unsigned int i = 0; i < order; ++i) {
            for (unsigned int j = 0; j < order; ++j)
                std::cout << M(i,j) << " ";
            std::cout << std::endl;
        }
    }
};

template <unsigned BorderTreatmentType, bool left_border, typename SrcIterator, typename SrcAccessor>
class IIRBorderSrcAccessor
{
public:
    typedef typename SrcAccessor::value_type value_type;

    SrcAccessor sa_;
    SrcIterator is_;
    SrcIterator iend_;
    BorderTreatmentMode border_;
    bool left_border_;
    int start_, stop_;

    IIRBorderSrcAccessor(SrcIterator is, SrcIterator iend, SrcAccessor sa, int start, int stop) : sa_(sa), is_(is), iend_(iend), start_(start), stop_(stop)
    {  
    }

    template<unsigned T = BorderTreatmentType>
    typename std::enable_if<T == BORDER_TREATMENT_ZEROPAD, value_type>::type
    inline operator()( SrcIterator i )
    {
            return NumericTraits<value_type>::zero();
    }

    template<unsigned T = BorderTreatmentType>
    typename std::enable_if<T == BORDER_TREATMENT_REPEAT && left_border, value_type>::type
    inline operator()( SrcIterator i )
    {
            return sa_(is_ + start_);
    }

    template<unsigned T = BorderTreatmentType>
    typename std::enable_if<T == BORDER_TREATMENT_REPEAT && !left_border, value_type>::type
    inline operator()( SrcIterator i )
    {
            return sa_(is_ + stop_ - 1);
    }

    template<unsigned T = BorderTreatmentType>
    typename std::enable_if<T == BORDER_TREATMENT_REFLECT && left_border, value_type>::type
    inline operator()( SrcIterator i )
    {
        return sa_(is_ + start_ + std::distance(i, is_) + start_);
    }

    template<unsigned T = BorderTreatmentType>
    typename std::enable_if<T == BORDER_TREATMENT_REFLECT && !left_border, value_type>::type
    inline operator()( SrcIterator i )
    {
        return sa_(iend_ - std::distance(iend_, i) - 2);
    }

    template<unsigned T = BorderTreatmentType>
    typename std::enable_if<T == BORDER_TREATMENT_WRAP && left_border, value_type>::type
    inline operator()( SrcIterator i )
    {
        return sa_(iend_ - std::distance(i, is_) + start_);
    }

    template<unsigned T = BorderTreatmentType>
    typename std::enable_if<T == BORDER_TREATMENT_WRAP && !left_border, value_type>::type
    inline operator()( SrcIterator i )
    {
        return sa_(is_ + start_ + std::distance(iend_, i));
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
          class DestAccessor, class RecursiveConvolutionKernel, class SumType,
          typename std::enable_if<is_deriche_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr>
void dericheApplyCausal(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                        DestIterator id, DestAccessor da,
                        RecursiveConvolutionKernel kernel, SumType xtmp[],
                        SumType ytmp[], int start, int stop) {
    const unsigned int order = kernel.order;

    SumType xi, yi, dyi;

    for (int x = start; x < stop; ++x) {
        xi = sa(is + x);

        xtmp[0] = xi;

        dyi = NumericTraits<SumType>::zero();
        for (unsigned int i = 0; i < order; ++i)
            dyi += kernel.n_causal[i] * xtmp[i];
        for (unsigned int i = 0; i < order; ++i)
            dyi -= kernel.d[i] * ytmp[i];

        for (unsigned int i = order - 1; i > 0; --i) {
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
          class DestAccessor, class RecursiveConvolutionKernel, class SumType,
          typename std::enable_if<is_deriche_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr>
void dericheApplyAntiCausal(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                            DestIterator id, DestAccessor da,
                            RecursiveConvolutionKernel kernel, SumType xtmp[],
                            SumType ytmp[], int start, int stop) {
    const unsigned int order = kernel.order;

    SumType xi, yi, dyi;

    for (int x = stop - 1; x >= start; --x) {
        xi = sa(is + x);
        yi = da(id + x);

        dyi = NumericTraits<SumType>::zero();
        for (unsigned int i = 0; i < order; ++i)
            dyi += kernel.n_anticausal[i] * xtmp[i];
        for (unsigned int i = 0; i < order; ++i)
            dyi -= kernel.d[i] * ytmp[i];

        for (unsigned int i = order - 1; i > 0; --i) {
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

// SSE
template <class SrcIterator, class SrcAccessor, class DestIterator,
          class DestAccessor, class RecursiveConvolutionKernel, class SumType,
          typename std::enable_if<is_deriche_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr,
          typename std::enable_if<RecursiveConvolutionKernel::sse>::type * = nullptr>
void dericheApplyCausalSIMD(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                        DestIterator id, DestAccessor da,
                        RecursiveConvolutionKernel kernel, SumType xtmp[],
                        SumType ytmp[], int start, int stop)
{
    dericheApplyCausal(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
}

template <class SrcIterator, class SrcAccessor, class DestIterator,
          class DestAccessor, class RecursiveConvolutionKernel, class SumType,
          typename std::enable_if<is_deriche_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr,
          typename std::enable_if<RecursiveConvolutionKernel::sse>::type * = nullptr>
void dericheApplyAntiCausalSIMD(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                        DestIterator id, DestAccessor da,
                        RecursiveConvolutionKernel kernel, SumType xtmp[],
                        SumType ytmp[], int start, int stop)
{
    dericheApplyAntiCausal(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
}

// AVX
template <class SrcIterator, class SrcAccessor, class DestIterator,
          class DestAccessor, class RecursiveConvolutionKernel, class SumType,
          typename std::enable_if<is_deriche_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr,
          typename std::enable_if<RecursiveConvolutionKernel::avx>::type * = nullptr>
void dericheApplyCausalSIMD(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                        DestIterator id, DestAccessor da,
                        RecursiveConvolutionKernel kernel, SumType xtmp[],
                        SumType ytmp[], int start, int stop)
{
    dericheApplyCausal(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
}

template <class SrcIterator, class SrcAccessor, class DestIterator,
          class DestAccessor, class RecursiveConvolutionKernel, class SumType,
          typename std::enable_if<is_deriche_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr,
          typename std::enable_if<RecursiveConvolutionKernel::avx>::type * = nullptr>
void dericheApplyAntiCausalSIMD(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                        DestIterator id, DestAccessor da,
                        RecursiveConvolutionKernel kernel, SumType xtmp[],
                        SumType ytmp[], int start, int stop)
{
    dericheApplyAntiCausal(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
}

// no SSE, no AVX
template <class SrcIterator, class SrcAccessor, class DestIterator,
          class DestAccessor, class RecursiveConvolutionKernel, class SumType,
          typename std::enable_if<is_deriche_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr,
          typename std::enable_if<!RecursiveConvolutionKernel::sse && !RecursiveConvolutionKernel::avx>::type * = nullptr>
inline void dericheApplyCausalSIMD(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                        DestIterator id, DestAccessor da,
                        RecursiveConvolutionKernel kernel, SumType xtmp[],
                        SumType ytmp[], int start, int stop)
{
    dericheApplyCausal(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
}

template <class SrcIterator, class SrcAccessor, class DestIterator,
          class DestAccessor, class RecursiveConvolutionKernel, class SumType,
          typename std::enable_if<is_deriche_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr,
          typename std::enable_if<!RecursiveConvolutionKernel::sse && !RecursiveConvolutionKernel::avx>::type * = nullptr>
inline void dericheApplyAntiCausalSIMD(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                        DestIterator id, DestAccessor da,
                        RecursiveConvolutionKernel kernel, SumType xtmp[],
                        SumType ytmp[], int start, int stop)
{
    dericheApplyAntiCausal(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
}

template <unsigned BorderTreatmentType,
          class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel,
          typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr,
          typename std::enable_if<is_deriche_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr>
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

    auto sa_causal = IIRBorderSrcAccessor<BorderTreatmentType, true, SrcIterator, SrcAccessor>(is, iend, sa, start, stop);
    auto sa_anticausal = IIRBorderSrcAccessor<BorderTreatmentType, false, SrcIterator, SrcAccessor>(is, iend, sa, start, stop);
    auto da_fake = DericheBorderDstAccessor<DestIterator, DestAccessor>();

    for (unsigned int i = 0; i < kernel.order; ++i)
        xtmp[i] = ytmp[i] = NumericTraits<SumType>::zero();
    dericheApplyCausal(is, iend, sa_causal, id, da_fake, kernel, xtmp, ytmp, start + kernel.left(), start);
    dericheApplyCausalSIMD(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);

    for (unsigned int i = 0; i < kernel.order; ++i)
        xtmp[i] = ytmp[i] = NumericTraits<SumType>::zero();
    dericheApplyAntiCausal(is, iend, sa_anticausal, id, da_fake, kernel, xtmp, ytmp, stop, stop + kernel.right());
    dericheApplyAntiCausalSIMD(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);

}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel,
          typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr,
          typename std::enable_if<is_deriche_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr>
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
    dericheApplyCausalSIMD(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);

    for (unsigned int i = 0; i < kernel.order; ++i)
        xtmp[i] = ytmp[i] = NumericTraits<SumType>::zero();
    dericheApplyAntiCausalSIMD(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);

}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel,
          typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr,
          typename std::enable_if<is_deriche_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr>
inline void recursiveConvolveLine(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  DestIterator id, DestAccessor da,
                  RecursiveConvolutionKernel kernel,
                  BorderTreatmentMode border,
                  int start = 0, int stop = 0)
{ 
    switch(kernel.border_treatment) {
        case BORDER_TREATMENT_ZEROPAD:
            recursiveConvolveLineDericheZeroBorder(is, iend, sa, id, da, kernel, start, stop);
            //recursiveConvolveLineDericheBorder<BORDER_TREATMENT_ZEROPAD>(is, iend, sa, id, da, kernel, start, stop);
            break;

        case BORDER_TREATMENT_REPEAT:
            recursiveConvolveLineDericheBorder<BORDER_TREATMENT_REPEAT>(is, iend, sa, id, da, kernel, start, stop);
            break;

        case BORDER_TREATMENT_REFLECT:
            recursiveConvolveLineDericheBorder<BORDER_TREATMENT_REFLECT>(is, iend, sa, id, da, kernel, start, stop);
            break;

        case BORDER_TREATMENT_WRAP:
            recursiveConvolveLineDericheBorder<BORDER_TREATMENT_WRAP>(is, iend, sa, id, da, kernel, start, stop);
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
          class SumType,
          typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr,
          typename std::enable_if<is_vyv_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr>
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
          class SumType,
          typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr,
          typename std::enable_if<is_vyv_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr>
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


template <class T, class A, int N, int M>
TinyVector<TinyVector<T, N> , M>
operator*(const Matrix<T, A> &a, const TinyVector<TinyVector<T, N>, M> &b)
{
    vigra_precondition(N == rowCount(a) && N == columnCount(a),
         "operator*(Matrix, TinyVector): Shape mismatch.");

    typedef TinyVector<T, N> InnerVector;

    TinyVector<InnerVector, M> res = TinyVectorView<InnerVector, M>();
    for(MultiArrayIndex i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j)
            res[j] += TinyVectorView<T, N>(&a(0,i)) * b[i][j];
    }
    return res;
}


template <unsigned BorderTreatmentType, class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel,
          class SumType,
          typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr,
          typename std::enable_if<is_vyv_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr>
inline void recursiveConvolveLineVYVBorder(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  DestIterator id, DestAccessor da,
                  RecursiveConvolutionKernel kernel,
                  SumType xtmp, SumType ytmp[],
                  int start = 0, int stop = 0)
{
    auto sa_left = IIRBorderSrcAccessor<BorderTreatmentType, true, SrcIterator, SrcAccessor>(is, iend, sa, start, stop);
    auto sa_right = IIRBorderSrcAccessor<BorderTreatmentType, false, SrcIterator, SrcAccessor>(is, iend, sa, start, stop);
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

    boundary = kernel.M * boundary;

    for (unsigned int i = 0; i < kernel.order; ++i) {
        std::cout << "HUH " << ytmp[i] << " " << boundary[i] << std::endl;
        ytmp[i] = boundary[i];
        v[kernel.right() - 5 + i] = boundary[i];
    }

    xtmp = NumericTraits<SumType>::zero();


    vyvApplyAntiCausal(v.begin(), va, kernel, xtmp, ytmp, 0, kernel.right() - 4);

    // anticausal filter
    vyvApplyAntiCausal(id, da, kernel, xtmp, ytmp, start, stop);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel,
          typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr,
          typename std::enable_if<is_vyv_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr>
inline void recursiveConvolveLine(SrcIterator is, SrcIterator iend, SrcAccessor sa,
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
            //recursiveConvolveLineVYVBorderZero(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
            recursiveConvolveLineVYVBorder<BORDER_TREATMENT_ZEROPAD>(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
            break;

        case BORDER_TREATMENT_REPEAT:
            recursiveConvolveLineVYVBorder<BORDER_TREATMENT_REPEAT>(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
            break;

        case BORDER_TREATMENT_REFLECT:
            recursiveConvolveLineVYVBorder<BORDER_TREATMENT_REFLECT>(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
            break;

        case BORDER_TREATMENT_WRAP:
            recursiveConvolveLineVYVBorder<BORDER_TREATMENT_WRAP>(is, iend, sa, id, da, kernel, xtmp, ytmp, start, stop);
            break;

        case BORDER_TREATMENT_AVOID:
        case BORDER_TREATMENT_CLIP:
        default:
            vigra_precondition(false, "recursiveConvolveLine: unsupported border treatment for VYV filters.");
            break;
    }
}

} // namespace detail


template <unsigned t_order, bool is_vyv, bool use_sse, bool use_avx, class ARITHTYPE = double>
class RecursiveConvolutionKernel :
public std::conditional<is_vyv, detail::RecursiveConvolutionKernelVYVMembers<t_order,ARITHTYPE>, detail::RecursiveConvolutionKernelDericheMembers<t_order,ARITHTYPE>>::type
{
#ifndef __SSE2__
    static_assert(!use_sse, "SSE requested but compiled without SSE support!");
#endif
#ifndef __AVX__
    static_assert(!use_avx, "AVX requested but compiled without AVX support!");
#endif
    static_assert(!(use_avx && is_vyv), "VYV filters do not support AVX");
    static_assert(!(use_sse && is_vyv), "VYV filters do not support SSE");
    static_assert(!(use_sse || use_avx) || std::is_same<ARITHTYPE, float>::value, "SSE/AVX require single precision arithmetic.");
    static_assert(!(use_sse || use_avx) || (use_sse ^ use_avx), "SSE and AVX cannot both be enabled at the same time.");
    static_assert(is_vyv || (t_order >= 2 && t_order <= 4), "Deriche order must be between two and four.");
    static_assert(!is_vyv || (t_order >= 3 && t_order <= 5), "VYV order must be between three and five.");

public:
    static const bool sse = use_sse;
    static const bool avx = use_avx;
    static const unsigned int order = t_order;

    typedef detail::iir_kernel_tag vigra_kernel_category;

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
        return round(-kernel_radius);
    }

    int right() {
        return round(kernel_radius);
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


#if 0

namespace detail {

static inline __m128
    __attribute__((__gnu_inline__, __always_inline__, __artificial__))
    _mm_hsum_ps_vec(__m128 a) {
#ifdef __SSE3__
    a = _mm_hadd_ps(a, a);
    a = _mm_hadd_ps(a, a);
#else
    __m128 tmp = _mm_setzero_ps();
    tmp = _mm_movehl_ps(tmp, a);
    a = _mm_add_ps(tmp, a);
    tmp = _mm_shuffle_ps(a, a, _MM_SHUFFLE(1, 1, 1, 1));
    a = _mm_add_ps(tmp, a);
#endif

    return a;
}

static inline float
    __attribute__((__gnu_inline__, __always_inline__, __artificial__))
    _mm_hsum_ps(__m128 a) {
    return _mm_cvtss_f32(_mm_hsum_ps_vec(a));
}

// SSE2 implementation for 4th order filter with single precision
template <class SrcIterator, class SrcAccessor, class DestIterator,
          class DestAccessor, class ConvolutionCoefficients, class SumType>
void dericheApplyCasual(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                        DestIterator id, DestAccessor da,
                        ConvolutionCoefficients coef, SumType xtmp[],
                        SumType ytmp[], int start, int stop,
                        std::true_type use_sse_float,
                        std::false_type use_avx_double) {
    unsigned int aligned_stop;
    float dyi;

    float xloadtmp[4];
    float ystoretmp[4];

    __m128 ps_x_vec[4];
    __m128 x_carry, y_carry;
    __m128 ps_x, ps_y;
    __m128 ps_n, ps_d;
    __m128 tmp;

    ps_d = _mm_set_ps(-coef.d[0], -coef.d[1], -coef.d[2], -coef.d[3]);
    ps_n = _mm_set_ps(coef.n_causal[0], coef.n_causal[1], coef.n_causal[2],
                      coef.n_causal[3]);

    x_carry = _mm_setr_ps(xtmp[0], xtmp[1], xtmp[2], xtmp[3]);
    y_carry = _mm_setr_ps(ytmp[0], ytmp[1], ytmp[2], ytmp[3]);

    aligned_stop = stop & ~3;

    for (int x = start; x < aligned_stop; x += 4) {
        for (unsigned int i = 0; i < 4; ++i)
            xloadtmp[i] = sa(is + x + i);

        ps_x_vec[3] = _mm_loadu_ps(xloadtmp); // 7 6 5 4
        ps_x_vec[1] = _mm_shuffle_ps(x_carry, ps_x_vec[3],
                                     _MM_SHUFFLE(1, 0, 3, 2)); // 5 4 3 2

        ps_x_vec[0] = _mm_shuffle_ps(x_carry, x_carry,
                                     _MM_SHUFFLE(0, 3, 2, 1)); // 0 3 2 1
        ps_x_vec[2] = _mm_shuffle_ps(ps_x_vec[3], ps_x_vec[3],
                                     _MM_SHUFFLE(0, 0, 1, 2)); // 4 4 5 6

        tmp = _mm_shuffle_ps(ps_x_vec[1], ps_x_vec[1],
                             _MM_SHUFFLE(0, 3, 2, 1)); // 4 3 2 5
        ps_x_vec[2] = _mm_shuffle_ps(tmp, ps_x_vec[2],
                                     _MM_SHUFFLE(0, 1, 1, 0)); // 6 5 4 3
        ps_x_vec[0] = _mm_shuffle_ps(ps_x_vec[0], tmp,
                                     _MM_SHUFFLE(1, 0, 1, 0)); // 4 3 2 1
        x_carry = ps_x_vec[3];

        for (unsigned int i = 0; i < 4; ++i) {
            ps_x = _mm_mul_ps(ps_n, ps_x_vec[i]);
            ps_y = _mm_mul_ps(y_carry, ps_d);

            tmp = _mm_hsum_ps_vec(_mm_add_ps(ps_x, ps_y));
            tmp = _mm_shuffle_ps(y_carry, tmp, _MM_SHUFFLE(0, 0, 3, 3));
            y_carry = _mm_shuffle_ps(y_carry, tmp, _MM_SHUFFLE(3, 0, 2, 1));
        }

        _mm_store_ps(ystoretmp, y_carry);

        for (unsigned int i = 0; i < 4; ++i)
            da.set(detail::RequiresExplicitCast<
                       typename DestAccessor::value_type>::cast(ystoretmp[i]),
                   id + x + i);
    }

    if (aligned_stop != stop) {
        unsigned int left = stop - aligned_stop;
        for (int x = aligned_stop; x < stop; x += 4) {
            for (unsigned int i = 0; i < left; ++i)
                xloadtmp[i] = sa(is + x + i);

            ps_x_vec[3] = _mm_loadu_ps(xloadtmp); // 7 6 5 4
            ps_x_vec[1] = _mm_shuffle_ps(x_carry, ps_x_vec[3],
                                         _MM_SHUFFLE(1, 0, 3, 2)); // 5 4 3 2

            ps_x_vec[0] = _mm_shuffle_ps(x_carry, x_carry,
                                         _MM_SHUFFLE(0, 3, 2, 1)); // 0 3 2 1
            ps_x_vec[2] = _mm_shuffle_ps(ps_x_vec[3], ps_x_vec[3],
                                         _MM_SHUFFLE(0, 0, 1, 2)); // 4 4 5 6

            tmp = _mm_shuffle_ps(ps_x_vec[1], ps_x_vec[1],
                                 _MM_SHUFFLE(0, 3, 2, 1)); // 4 3 2 5
            ps_x_vec[2] = _mm_shuffle_ps(tmp, ps_x_vec[2],
                                         _MM_SHUFFLE(0, 1, 1, 0)); // 6 5 4 3
            ps_x_vec[0] = _mm_shuffle_ps(ps_x_vec[0], tmp,
                                         _MM_SHUFFLE(1, 0, 1, 0)); // 4 3 2 1
            x_carry = ps_x_vec[3];

            for (unsigned int i = 0; i < 4; ++i) {
                ps_x = _mm_mul_ps(ps_n, ps_x_vec[i]);
                ps_y = _mm_mul_ps(y_carry, ps_d);

                tmp = _mm_hsum_ps_vec(_mm_add_ps(ps_x, ps_y));
                tmp = _mm_shuffle_ps(y_carry, tmp, _MM_SHUFFLE(0, 0, 3, 3));
                y_carry = _mm_shuffle_ps(y_carry, tmp, _MM_SHUFFLE(3, 0, 2, 1));
            }

            _mm_store_ps(ystoretmp, y_carry);

            for (unsigned int i = 0; i < left; ++i)
                da.set(
                    detail::RequiresExplicitCast<
                        typename DestAccessor::value_type>::cast(ystoretmp[i]),
                    id + x + i);
        }
    }
}

template <class SrcIterator, class SrcAccessor, class DestIterator,
          class DestAccessor, class ConvolutionCoefficients, class SumType>
void dericheApplyAntiCasual(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                            DestIterator id, DestAccessor da,
                            ConvolutionCoefficients coef, SumType xtmp[],
                            SumType ytmp[], int start, int stop,
                            std::true_type use_sse_float,
                            std::false_type use_avx_double) {

    __m128 ps_n = _mm_setr_ps(coef.n_anticausal[0], coef.n_anticausal[1],
                              coef.n_anticausal[2], coef.n_anticausal[3]);
    __m128 ps_d = _mm_setr_ps(-coef.d[0], -coef.d[1], -coef.d[2], -coef.d[3]);

    float xi[4] __attribute__((aligned(16)));
    float yi[4] __attribute__((aligned(16)));
    float dyi;

    for (unsigned int i = 0; i < 4; ++i) {
        yi[i] = ytmp[i];
        xi[i] = xtmp[i];
    }

    for (int x = stop - 1; x >= start; --x) {
        __m128 ps_x = _mm_load_ps(xi);
        __m128 ps_y = _mm_load_ps(yi);

        ps_x = _mm_mul_ps(ps_x, ps_n);
        ps_y = _mm_mul_ps(ps_y, ps_d);
        dyi = _mm_hsum_ps(_mm_add_ps(ps_x, ps_y));

        da.set(detail::RequiresExplicitCast<
                   typename DestAccessor::value_type>::cast(da(id + x) + dyi),
               id + x);

        for (unsigned int i = 3; i > 0; --i)
            yi[i] = yi[i - 1];
        for (unsigned int i = 3; i > 0; --i)
            xi[i] = xi[i - 1];
        yi[0] = dyi;
        xi[0] = sa(is + x);
    }
}

#endif



template<typename RecursiveConvolutionKernel, typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr>
inline RecursiveConvolutionKernel kernel1d(RecursiveConvolutionKernel k)
{
    return k;
}


template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel,
          typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr>
void convolveLine(SrcIterator is, SrcIterator iend, SrcAccessor sa,
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
          class RecursiveConvolutionKernel,
          typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr>
inline void convolveLine(SrcIterator is, SrcIterator iend, SrcAccessor sa,
                  DestIterator id, DestAccessor da,
                  RecursiveConvolutionKernel ik,
                  int start = 0, int stop = 0)
{
    convolveLine(is, iend, sa, id, da, ik, ik.border_treatment, start, stop);
}


template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class RecursiveConvolutionKernel,
          typename std::enable_if<detail::is_iir_kernel<RecursiveConvolutionKernel>::value>::type * = nullptr>
inline void convolveLine(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIterator, DestAccessor> dest,
                  RecursiveConvolutionKernel kernel,
                  int start = 0, int stop = 0)
{
    convolveLine(src.first, src.second, src.third, dest.first, dest.second, kernel, kernel.border_treatment, start, stop);
}


} // namespace vigra

#endif
// VIGRA_RECURSIVE_MULTI_CONVOLUTION_H
