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

#ifndef VIGRA_KERNEL_H
#define VIGRA_KERNEL_H

#include <type_traits>

namespace vigra {
namespace detail {

struct fir_kernel1d_tag {};
struct fir_kernel2d_tag {};
struct iir_kernel1d_tag {};
struct gaussian_kernel1d_tag {};

template<typename X>
struct is_kernel
{
	template <typename Y>
	static char test(int, typename Y::vigra_kernel_category*);

	template <typename Y>
	static int test(int, ...);

	static const bool value = (sizeof(char) == sizeof(test<X>(0,0)));
};

template<typename X, bool = is_kernel<X>::value >
struct is_fir_kernel1d : public std::is_same<typename X::vigra_kernel_category, fir_kernel1d_tag>
{
};

template<typename X, bool = is_kernel<X>::value >
struct is_iir_kernel1d : public std::is_same<typename X::vigra_kernel_category, iir_kernel1d_tag>
{
};

template<typename X, bool = is_kernel<X>::value >
struct is_gaussian_kernel1d : public std::is_same<typename X::vigra_kernel_category, gaussian_kernel1d_tag>
{
};

template<typename X>
struct is_fir_kernel1d<X, false>
{
	static const bool value = false;
};

template<typename X>
struct is_iir_kernel1d<X, false>
{
	static const bool value = false;
};

template<typename X>
struct is_gaussian_kernel1d<X, false>
{
	static const bool value = false;
};


} // namespace detail

} // namespace vigra

#endif