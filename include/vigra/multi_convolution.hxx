//-- -*- c++ -*-
/************************************************************************/
/*                                                                      */
/*               Copyright 2003 by Christian-Dennis Rahn                */
/*                        and Ullrich Koethe                            */
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

#ifndef VIGRA_MULTI_CONVOLUTION_H
#define VIGRA_MULTI_CONVOLUTION_H

#include "separableconvolution.hxx"
#include "array_vector.hxx"
#include "multi_array.hxx"
#include "accessor.hxx"
#include "numerictraits.hxx"
#include "navigator.hxx"
#include "metaprogramming.hxx"
#include "multi_pointoperators.hxx"
#include "functorexpression.hxx"
#include "tinyvector.hxx"

namespace vigra
{

namespace detail
{

struct DoubleYielder
{
    const double value;
    DoubleYielder(double v, unsigned, const char *const) : value(v) {}
    DoubleYielder(double v)                              : value(v) {}
    void operator++() {}
    double operator*() const { return value; }
};

template <typename X>
struct IteratorDoubleYielder
{
    X it;
    IteratorDoubleYielder(X i, unsigned, const char *const) : it(i) {}
    IteratorDoubleYielder(X i)                              : it(i) {}
    void operator++() { ++it; }
    double operator*() const { return *it; }
};

template <typename X>
struct SequenceDoubleYielder
{
    typename X::const_iterator it;
    SequenceDoubleYielder(const X & seq, unsigned dim,
                          const char *const function_name = "SequenceDoubleYielder")
        : it(seq.begin())
    {
        if (seq.size() == dim)
            return;
        std::string msg = "(): Parameter number be equal to the number of spatial dimensions.";
        vigra_precondition(false, function_name + msg);
    }
    void operator++() { ++it; }
    double operator*() const { return *it; }
};

template <typename X>
struct WrapDoubleIterator
{
    typedef
    typename IfBool< IsConvertibleTo<X, double>::value,
        DoubleYielder,
        typename IfBool< IsIterator<X>::value || IsArray<X>::value,
            IteratorDoubleYielder<X>,
            SequenceDoubleYielder<X>
        >::type
    >::type type;
};

template <class Param1, class Param2, class Param3>
struct WrapDoubleIteratorTriple
{
    typename WrapDoubleIterator<Param1>::type sigma_eff_it;
    typename WrapDoubleIterator<Param2>::type sigma_d_it;
    typename WrapDoubleIterator<Param3>::type step_size_it;
    WrapDoubleIteratorTriple(Param1 sigma_eff, Param2 sigma_d, Param3 step_size)
        : sigma_eff_it(sigma_eff), sigma_d_it(sigma_d), step_size_it(step_size) {}
    void operator++()
    {
        ++sigma_eff_it;
        ++sigma_d_it;
        ++step_size_it;
    }
    double sigma_eff() const { return *sigma_eff_it; }
    double sigma_d() const { return *sigma_d_it; }
    double step_size() const { return *step_size_it; }
    static double sqr(double x) { return x * x; }
    static void sigma_precondition(double sigma, const char *const function_name)
    {
        if (sigma < 0.0)
        {
             std::string msg = "(): Scale must be positive.";
             vigra_precondition(false, function_name + msg);
        }
    }
    double sigma_scaled(const char *const function_name = "unknown function ") const
    {
        sigma_precondition(sigma_eff(), function_name);
        sigma_precondition(sigma_d(), function_name);
        double sigma_squared = sqr(sigma_eff()) - sqr(sigma_d());
        if (sigma_squared > 0.0)
        {
            return std::sqrt(sigma_squared) / step_size();
        }
        else
        {
             std::string msg = "(): Scale would be imaginary or zero.";
             vigra_precondition(false, function_name + msg);
             return 0;
        }
    }
};

template <unsigned dim>
struct multiArrayScaleParam
{
    typedef TinyVector<double, dim> p_vector;
    typedef typename p_vector::const_iterator return_type;
    p_vector vec;

    template <class Param>
    multiArrayScaleParam(Param val, const char *const function_name = "multiArrayScaleParam")
    {
        typename WrapDoubleIterator<Param>::type in(val, dim, function_name);
        for (unsigned i = 0; i != dim; ++i, ++in)
            vec[i] = *in;
    }
    return_type operator()() const
    {
        return vec.begin();
    }
    static void precondition(unsigned n_par, const char *const function_name = "multiArrayScaleParam")
    {
        char n[3] = "0.";
        n[0] += dim;
        std::string msg = "(): dimension parameter must be ";
        vigra_precondition(dim == n_par, function_name + msg + n);
    }
    multiArrayScaleParam(double v0, double v1, const char *const function_name = "multiArrayScaleParam")
    {
        precondition(2, function_name);
        vec = p_vector(v0, v1);
    }
    multiArrayScaleParam(double v0, double v1, double v2, const char *const function_name = "multiArrayScaleParam")
    {
        precondition(3, function_name);
        vec = p_vector(v0, v1, v2);
    }
    multiArrayScaleParam(double v0, double v1, double v2,  double v3, const char *const function_name = "multiArrayScaleParam")
    {
        precondition(4, function_name);
        vec = p_vector(v0, v1, v2, v3);
    }
    multiArrayScaleParam(double v0, double v1, double v2,  double v3, double v4, const char *const function_name = "multiArrayScaleParam")
    {
        precondition(5, function_name);
        vec = p_vector(v0, v1, v2, v3, v4);
    }
};

} // namespace detail

#define VIGRA_CONVOLUTION_OPTIONS(function_name, default_value, member_name) \
    template <class Param> \
    ConvolutionOptions & function_name(const Param & val) \
    { \
        member_name = ParamVec(val, "ConvolutionOptions::" #function_name); \
        return *this; \
    } \
    ConvolutionOptions & function_name() \
    { \
        member_name = ParamVec(default_value, "ConvolutionOptions::" #function_name); \
        return *this; \
    } \
    ConvolutionOptions & function_name(double v0, double v1) \
    { \
        member_name = ParamVec(v0, v1, "ConvolutionOptions::" #function_name); \
        return *this; \
    } \
    ConvolutionOptions & function_name(double v0, double v1, double v2) \
    { \
        member_name = ParamVec(v0, v1, v2, "ConvolutionOptions::" #function_name); \
        return *this; \
    } \
    ConvolutionOptions & function_name(double v0, double v1, double v2, double v3) \
    { \
        member_name = ParamVec(v0, v1, v2, v3, "ConvolutionOptions::" #function_name); \
        return *this; \
    } \
    ConvolutionOptions & function_name(double v0, double v1, double v2, double v3, double v4) \
    { \
        member_name = ParamVec(v0, v1, v2, v3, v4, "ConvolutionOptions::" #function_name); \
        return *this; \
    }


/** \brief  Options class template for convolutions.
 
  <b>\#include</b> \<vigra/multi_convolution.hxx\>
  
  This class enables the calculation of scale space convolutions
  such as \ref gaussianGradientMultiArray() on data with anistropic
  discretization. For these, the result of the ordinary caluculation
  has to be multiplied by factors of \f$1/w^{n}\f$ for each dimension,
  where \f$w\f$ is the step size of the grid in said dimension and
  \f$n\f$ is the differiantial order of the convolution, e.g., 1 for
  gaussianGradientMultiArray(), and 0 for gaussianSmoothMultiArray(),
  respectively. Also for each dimension in turn, the convolution's scale
  parameter \f$\sigma\f$ has to be replaced by
  \f$\sqrt{\sigma_\mathrm{eff}^2 - \sigma_\mathrm{D}^2}\Big/w\f$,
  where \f$\sigma_\mathrm{eff}\f$ is the resulting effective filtering
  scale. The data is assumed to be already filtered by a 
  gaussian smoothing with the scale parameter \f$\sigma_\mathrm{D}\f$
  (such as by measuring equipment). All of the above changes are
  automatically employed by the convolution functions for <tt>MultiArray</tt>s
  if a correspondig options object is provided.

  The <tt>ConvolutionOptions</tt> class must be parameterized by the dimension
  <tt>dim</tt>
  of the <tt>MultiArray</tt>s on which it is used. The actual per-axis
  options are set by (overloaded) member functions explained below,
  or else default to neutral values corresponding to the absence of the
  particular option.
  
  All member functions set <tt>dim</tt> values of the respective convolution
  option, one for each dimension. They may be set explicitly by multiple
  arguments for up to five dimensions, or by a single argument to the same
  value for all dimensions. For the general case, a single argument that is
  either a C-syle array, an iterator, or a C++ standard library style
  sequence (such as <tt>std::vector</tt>, with member functions <tt>begin()</tt>
  and <tt>size()</tt>) supplies the option values for any number of dimensions.
  
  Note that the return value of all member functions is <tt>*this</tt>, which
  provides the mechanism for concatenating member function calls as shown below.

  <b>usage with explicit parameters:</b>

  \code
  ConvolutionOptions<2> opt = ConvolutionOptions<2>().stepSize(1, 2.3);
  \endcode
 
  <b>usage with arrays:</b>
 
  \code
  const double step_size[3] = { x_scale, y_scale, z_scale };
  ConvolutionOptions<3> opt = ConvolutionOptions<3>().stepSize(step_size);
  \endcode

  <b>usage with C++ standard library style sequences:</b>
 
  \code
  TinyVector<double, 4> step_size(1, 1, 2.0, 1.5);
  TinyVector<double, 4>  r_sigmas(1, 1, 2.3, 3.2);
  ConvolutionOptions<4> opt = ConvolutionOptions<4>().stepSize(step_size).resolutionStdDev(r_sigmas);
  \endcode

  <b>usage with iterators:</b>

  \code
  ArrayVector<double> step_size;
  step_size.push_back(0);
  step_size.push_back(3);
  step_size.push_back(4);
  ArrayVector<double>::iterator i = step_size.begin();
  ++i;
  ConvolutionOptions<2> opt = ConvolutionOptions<2>().stepSize(i);
  \endcode

  <b>general usage in a convolution function call:</b>

  \code
  MultiArray<3, double> test_image;
  MultiArray<3, double> out_image;
  gaussianSmoothMultiArray(srcMultiArrayRange(test_image),
                           destMultiArray(out_image),
                           5.0,
                           ConvolutionOptions<3>()
                              .stepSize        (1, 1, 3.2)
                              .resolutionStdDev(1, 1, 4)
                          );
  \endcode
 
*/
template <unsigned dim>
class ConvolutionOptions
{
     typedef detail::multiArrayScaleParam<dim> ParamVec;
     typedef typename ParamVec::return_type    ParamIt;
     ParamVec sigma_eff;
     ParamVec sigma_d;
     ParamVec step_size;
     ParamVec outer_scale;
  public:
    ConvolutionOptions()
    : sigma_eff(0.0),
      sigma_d(0.0),
      step_size(1.0),
      outer_scale(0.0)
    {}

    typedef typename detail::WrapDoubleIteratorTriple<ParamIt, ParamIt, ParamIt>
        ScaleIterator;
    typedef typename detail::WrapDoubleIterator<ParamIt>::type
        StepIterator;

    ScaleIterator scaleParams() const
    {
        return ScaleIterator(sigma_eff(), sigma_d(), step_size());
    }
    StepIterator stepParams() const
    {
        return StepIterator(step_size());
    }

    ConvolutionOptions outerOptions() const
    {
        ConvolutionOptions outer = *this;
        // backward-compatible values:
        return outer.stdDev(outer_scale()).resolutionStdDev(0.0);
    }

    // Step size per axis.
    // Default: dim values of 1.0
    VIGRA_CONVOLUTION_OPTIONS(stepSize, 1.0, step_size)
#ifdef DOXYGEN
        /** Step size(s) per axis, i.e., the distance between two
            adjacent pixels. Required for <tt>MultiArray</tt>
            containig anisotropic data.
 
            Note that a convolution containing a derivative operator
            of order <tt>n</tt> results in a multiplication by 
            \f${\rm stepSize}^{-n}\f$ for each axis.
            Also, the above standard deviations
            are scaled according to the step size of each axis.
            Default value for the options object if this member funtion is not
            used: A value of 1.0 for each dimension.
        */
    ConvolutionOptions<dim> & stepSize(...);
#endif

    // Resolution standard deviation per axis.
    // Default: dim values of 0.0
    VIGRA_CONVOLUTION_OPTIONS(resolutionStdDev, 0.0, sigma_d)
#ifdef DOXYGEN
        /** Resolution standard deviation(s) per axis, i.e., a supposed
            pre-existing gaussian filtering by this value.
       
            The standard deviation actually used by the convolution operators
            is \f$\sqrt{{\rm sigma}^{2} - {\rm resolutionStdDev}^{2}}\f$ for each
            axis.
            Default value for the options object if this member funtion is not
            used: A value of 0.0 for each dimension.
        */
    ConvolutionOptions<dim> & resolutionStdDev(...);
#endif

    // Standard deviation of scale space operators.
    // Default: dim values of 0.0
    VIGRA_CONVOLUTION_OPTIONS(stdDev, 0.0, sigma_eff)
    VIGRA_CONVOLUTION_OPTIONS(innerScale, 0.0, sigma_eff)

#ifdef DOXYGEN
        /** Standard deviation(s) of scale space operators, or inner scale(s) for \ref structureTensorMultiArray().
        
            Usually not
            needed, since a single value for all axes may be specified as a parameter
            <tt>sigma</tt> to the call of
            an convolution operator such as \ref gaussianGradientMultiArray(), and
            anisotropic data requiring the use of the stepSize() member function.
            Default value for the options object if this member funtion is not
            used: A value of 0.0 for each dimension.
        */
    ConvolutionOptions<dim> & stdDev(...);

        /** Standard deviation(s) of scale space operators, or inner scale(s) for \ref structureTensorMultiArray().
        
            Usually not
            needed, since a single value for all axes may be specified as a parameter
            <tt>sigma</tt> to the call of
            an convolution operator such as \ref gaussianGradientMultiArray(), and
            anisotropic data requiring the use of the stepSize() member function.
            Default value for the options object if this member funtion is not
            used: A value of 0.0 for each dimension.
        */
    ConvolutionOptions<dim> & innerScale(...);
#endif

    // Outer scale, for structure tensor.
    // Default: dim values of 0.0
    VIGRA_CONVOLUTION_OPTIONS(outerScale, 0.0, outer_scale)
#ifdef DOXYGEN
        /** Standard deviation(s) of the second convolution of the
            structure tensor. 

            Usually not needed, since a single value for
            all axes may be specified as a parameter <tt>outerScale</tt> to
            the call of \ref structureTensorMultiArray(), and
            anisotropic data requiring the use of the stepSize() member
            function.
            Default value for the options object if this member funtion is not
            used: A value of 0.0 for each dimension.
        */
    ConvolutionOptions<dim> & outerScale(...);
#endif


};

namespace detail
{

/********************************************************/
/*                                                      */
/*        internalSeparableConvolveMultiArray           */
/*                                                      */
/********************************************************/

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class KernelIterator>
void
internalSeparableConvolveMultiArrayTmp(
                      SrcIterator si, SrcShape const & shape, SrcAccessor src,
                      DestIterator di, DestAccessor dest, KernelIterator kit)
{
    enum { N = 1 + SrcIterator::level };

    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;

    // temporay array to hold the current line to enable in-place operation
    ArrayVector<TmpType> tmp( shape[0] );

    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;

    { // only operate on first dimension here
        SNavigator snav( si, shape, 0 );
        DNavigator dnav( di, shape, 0 );

        for( ; snav.hasMore(); snav++, dnav++ )
        {
             // first copy source to temp for maximum cache efficiency
             copyLine( snav.begin(), snav.end(), src,
                       tmp.begin(), typename AccessorTraits<TmpType>::default_accessor() );

             convolveLine( srcIterRange(tmp.begin(), tmp.end(),
                                        typename AccessorTraits<TmpType>::default_const_accessor()),
                           destIter( dnav.begin(), dest ),
                           kernel1d( *kit ) );
        }
        ++kit;
    }

    // operate on further dimensions
    for( int d = 1; d < N; ++d, ++kit )
    {
        DNavigator dnav( di, shape, d );

        tmp.resize( shape[d] );

        for( ; dnav.hasMore(); dnav++ )
        {
             // first copy source to temp for maximum cache efficiency
             copyLine( dnav.begin(), dnav.end(), dest,
                       tmp.begin(), typename AccessorTraits<TmpType>::default_accessor() );

             convolveLine( srcIterRange(tmp.begin(), tmp.end(),
                                        typename AccessorTraits<TmpType>::default_const_accessor()),
                           destIter( dnav.begin(), dest ),
                           kernel1d( *kit ) );
        }
    }
}

template <class K>
void 
scaleKernel(K & kernel, double a)
{
    for(int i = kernel.left(); i <= kernel.right(); ++i)
        kernel[i] = detail::RequiresExplicitCast<typename K::value_type>::cast(kernel[i] * a);
}


} // namespace detail

/** \addtogroup MultiArrayConvolutionFilters Convolution filters for multi-dimensional arrays.

    These functions realize a separable convolution on an arbitrary dimensional
    array that is specified by iterators (compatible to \ref MultiIteratorPage)
    and shape objects. It can therefore be applied to a wide range of data structures
    (\ref vigra::MultiArrayView, \ref vigra::MultiArray etc.).
*/
//@{

/********************************************************/
/*                                                      */
/*             separableConvolveMultiArray              */
/*                                                      */
/********************************************************/

/** \brief Separated convolution on multi-dimensional arrays.

    This function computes a separated convolution on all dimensions
    of the given multi-dimensional array. Both source and destination
    arrays are represented by iterators, shape objects and accessors.
    The destination array is required to already have the correct size.

    There are two variants of this functions: one takes a single kernel
    of type \ref vigra::Kernel1D which is then applied to all dimensions,
    whereas the other requires an iterator referencing a sequence of
    \ref vigra::Kernel1D objects, one for every dimension of the data.
    Then the first kernel in this sequence is applied to the innermost
    dimension (e.g. the x-dimension of an image), while the last is applied to the
    outermost dimension (e.g. the z-dimension in a 3D image).

    This function may work in-place, which means that <tt>siter == diter</tt> is allowed.
    A full-sized internal array is only allocated if working on the destination
    array directly would cause round-off errors (i.e. if
    <tt>typeid(typename NumericTraits<typename DestAccessor::value_type>::RealPromote)
    != typeid(typename DestAccessor::value_type)</tt>.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        // apply the same kernel to all dimensions
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class T>
        void
        separableConvolveMultiArray(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest,
                                    Kernel1D<T> const & kernel);

        // apply each kernel from the sequence 'kernels' in turn
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class KernelIterator>
        void
        separableConvolveMultiArray(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest,
                                    KernelIterator kernels);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        // apply the same kernel to all dimensions
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class T>
        void
        separableConvolveMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest,
                                    Kernel1D<T> const & kernel);

        // apply each kernel from the sequence 'kernels' in turn
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class KernelIterator>
        void
        separableConvolveMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest,
                                    KernelIterator kernels);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\>

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, float> dest(shape);
    ...
    Kernel1D<float> gauss;
    gauss.initGaussian(sigma);
    // create 3 Gauss kernels, one for each dimension
    ArrayVector<Kernel1D<float> > kernels(3, gauss);

    // perform Gaussian smoothing on all dimensions
    separableConvolveMultiArray(srcMultiArrayRange(source), destMultiArray(dest), 
                                kernels.begin());
    \endcode

    <b> Required Interface:</b>

    see \ref separableConvolveMultiArray(), in addition:

    \code
    int dimension = 0;
    VectorElementAccessor<DestAccessor> elementAccessor(0, dest);
    \endcode

    \see vigra::Kernel1D, convolveLine()
*/
doxygen_overloaded_function(template <...> void separableConvolveMultiArray)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class KernelIterator>
void
separableConvolveMultiArray( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                             DestIterator d, DestAccessor dest, KernelIterator kernels )
{
    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;

    if(!IsSameType<TmpType, typename DestAccessor::value_type>::boolResult)
    {
        // need a temporary array to avoid rounding errors
        MultiArray<SrcShape::static_size, TmpType> tmpArray(shape);
        detail::internalSeparableConvolveMultiArrayTmp( s, shape, src,
             tmpArray.traverser_begin(), typename AccessorTraits<TmpType>::default_accessor(), kernels );
        copyMultiArray(srcMultiArrayRange(tmpArray), destIter(d, dest));
    }
    else
    {
        // work directly on the destination array
        detail::internalSeparableConvolveMultiArrayTmp( s, shape, src, d, dest, kernels );
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class KernelIterator>
inline
void separableConvolveMultiArray(
    triple<SrcIterator, SrcShape, SrcAccessor> const & source,
    pair<DestIterator, DestAccessor> const & dest, KernelIterator kit )
{
    separableConvolveMultiArray( source.first, source.second, source.third,
                                 dest.first, dest.second, kit );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class T>
inline void
separableConvolveMultiArray( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                             DestIterator d, DestAccessor dest,
                             Kernel1D<T> const & kernel )
{
    ArrayVector<Kernel1D<T> > kernels(shape.size(), kernel);

    separableConvolveMultiArray( s, shape, src, d, dest, kernels.begin() );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class T>
inline void
separableConvolveMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                            pair<DestIterator, DestAccessor> const & dest,
                            Kernel1D<T> const & kernel )
{
    ArrayVector<Kernel1D<T> > kernels(source.second.size(), kernel);

    separableConvolveMultiArray( source.first, source.second, source.third,
                                 dest.first, dest.second, kernels.begin() );
}

/********************************************************/
/*                                                      */
/*            convolveMultiArrayOneDimension            */
/*                                                      */
/********************************************************/

/** \brief Convolution along a single dimension of a multi-dimensional arrays.

    This function computes a convolution along one dimension (specified by
    the parameter <tt>dim</tt> of the given multi-dimensional array with the given
    <tt>kernel</tt>. Both source and destination arrays are represented by
    iterators, shape objects and accessors. The destination array is required to
    already have the correct size.

    This function may work in-place, which means that <tt>siter == diter</tt> is allowed.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class T>
        void
        convolveMultiArrayOneDimension(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                       DestIterator diter, DestAccessor dest,
                                       unsigned int dim, vigra::Kernel1D<T> const & kernel);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor, class T>
        void
        convolveMultiArrayOneDimension(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                       pair<DestIterator, DestAccessor> const & dest,
                                       unsigned int dim, vigra::Kernel1D<T> const & kernel);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\>

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, float> dest(shape);
    ...
    Kernel1D<float> gauss;
    gauss.initGaussian(sigma);

    // perform Gaussian smoothing along dimensions 1 (height)
    convolveMultiArrayOneDimension(srcMultiArrayRange(source), destMultiArray(dest), 1, gauss);
    \endcode

    \see separableConvolveMultiArray()
*/
doxygen_overloaded_function(template <...> void convolveMultiArrayOneDimension)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class T>
void
convolveMultiArrayOneDimension(SrcIterator s, SrcShape const & shape, SrcAccessor src,
                               DestIterator d, DestAccessor dest,
                               unsigned int dim, vigra::Kernel1D<T> const & kernel )
{
    enum { N = 1 + SrcIterator::level };
    vigra_precondition( dim < N,
                        "convolveMultiArrayOneDimension(): The dimension number to convolve must be smaller "
                        "than the data dimensionality" );

    typedef typename NumericTraits<typename DestAccessor::value_type>::RealPromote TmpType;
    ArrayVector<TmpType> tmp( shape[dim] );

    typedef MultiArrayNavigator<SrcIterator, N> SNavigator;
    typedef MultiArrayNavigator<DestIterator, N> DNavigator;

    SNavigator snav( s, shape, dim );
    DNavigator dnav( d, shape, dim );

    for( ; snav.hasMore(); snav++, dnav++ )
    {
         // first copy source to temp for maximum cache efficiency
         copyLine( snav.begin(), snav.end(), src,
           tmp.begin(), typename AccessorTraits<TmpType>::default_accessor() );

         convolveLine( srcIterRange( tmp.begin(), tmp.end(), typename AccessorTraits<TmpType>::default_const_accessor()),
                       destIter( dnav.begin(), dest ),
                       kernel1d( kernel ) );
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor, class T>
inline void
convolveMultiArrayOneDimension(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                               pair<DestIterator, DestAccessor> const & dest,
                               unsigned int dim, vigra::Kernel1D<T> const & kernel )
{
    convolveMultiArrayOneDimension( source.first, source.second, source.third,
                                   dest.first, dest.second, dim, kernel );
}

/********************************************************/
/*                                                      */
/*             gaussianSmoothMultiArray                 */
/*                                                      */
/********************************************************/

/** \brief Isotropic Gaussian smoothing of a multi-dimensional arrays.

    This function computes an isotropic convolution of the given N-dimensional
    array with a Gaussian filter at the given standard deviation <tt>sigma</tt>.
    Both source and destination arrays are represented by
    iterators, shape objects and accessors. The destination array is required to
    already have the correct size. This function may work in-place, which means
    that <tt>siter == diter</tt> is allowed. It is implemented by a call to
    \ref separableConvolveMultiArray() with the appropriate kernel.

    Anisotropic data should be passed with appropiate
    \ref ConvolutionOptions, the parameter <tt>opt</tt> is otherwise optional
    unless the parameter <tt>sigma</tt> is left out.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        gaussianSmoothMultiArray(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                 DestIterator diter, DestAccessor dest,
                                 double sigma, const ConvolutionOptions<N> & opt);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        gaussianSmoothMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                 pair<DestIterator, DestAccessor> const & dest,
                                 double sigma, const ConvolutionOptions<N> & opt);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\>

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, float> dest(shape);
    ...
    // perform isotropic Gaussian smoothing at scale 'sigma'
    gaussianSmoothMultiArray(srcMultiArrayRange(source), destMultiArray(dest), sigma);
    \endcode

    <b> Usage with anisotropic data:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\>

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, float> dest(shape);
    TinyVector<float, 3> step_size;
    TinyVector<float, 3> resolution_sigmas;
    ...
    // perform anisotropic Gaussian smoothing at scale 'sigma'
    gaussianSmoothMultiArray(srcMultiArrayRange(source), destMultiArray(dest), sigma,
                             ConvolutionOptions<3>().stepSize(step_size).resolutionStdDev(resolution_sigmas));
    \endcode

    \see separableConvolveMultiArray()
*/
doxygen_overloaded_function(template <...> void gaussianSmoothMultiArray)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
gaussianSmoothMultiArray( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                   DestIterator d, DestAccessor dest,
                   const ConvolutionOptions<SrcShape::static_size> & opt,
                   const char *const function_name = "gaussianSmoothMultiArray" )
{
    typedef typename DestAccessor::value_type DestType;

    static const int N = SrcShape::static_size;

    typename ConvolutionOptions<N>::ScaleIterator params = opt.scaleParams();
    ArrayVector<Kernel1D<double> > kernels(N);

    for (int dim = 0; dim < N; ++dim, ++params)
        kernels[dim].initGaussian(params.sigma_scaled(function_name));

    separableConvolveMultiArray(s, shape, src, d, dest, kernels.begin());
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
gaussianSmoothMultiArray( SrcIterator s, SrcShape const & shape, SrcAccessor src,
                   DestIterator d, DestAccessor dest, double sigma,
                   const ConvolutionOptions<SrcShape::static_size> & opt = ConvolutionOptions<SrcShape::static_size>())
{
    ConvolutionOptions<SrcShape::static_size> par = opt;
    gaussianSmoothMultiArray(s, shape, src, d,  dest, par.stdDev(sigma));
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
gaussianSmoothMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                  pair<DestIterator, DestAccessor> const & dest,
                  const ConvolutionOptions<SrcShape::static_size> & opt)
{
    gaussianSmoothMultiArray( source.first, source.second, source.third,
                              dest.first, dest.second, opt );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
gaussianSmoothMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                  pair<DestIterator, DestAccessor> const & dest, double sigma,
                  const ConvolutionOptions<SrcShape::static_size> & opt = ConvolutionOptions<SrcShape::static_size>())
{
    gaussianSmoothMultiArray( source.first, source.second, source.third,
                              dest.first, dest.second, sigma, opt );
}


/********************************************************/
/*                                                      */
/*             gaussianGradientMultiArray               */
/*                                                      */
/********************************************************/

/** \brief Calculate Gaussian gradient of a multi-dimensional arrays.

    This function computes the Gaussian gradient of the given N-dimensional
    array with a sequence of first-derivative-of-Gaussian filters at the given
    standard deviation <tt>sigma</tt> (differentiation is applied to each dimension
    in turn, starting with the innermost dimension). Both source and destination arrays
    are represented by iterators, shape objects and accessors. The destination array is
    required to have a vector valued pixel type with as many elements as the number of
    dimensions. This function is implemented by calls to
    \ref separableConvolveMultiArray() with the appropriate kernels.

    Anisotropic data should be passed with appropiate
    \ref ConvolutionOptions, the parameter <tt>opt</tt> is otherwise optional
    unless the parameter <tt>sigma</tt> is left out.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        gaussianGradientMultiArray(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                   DestIterator diter, DestAccessor dest,
                                   double sigma, const ConvolutionOptions<N> & opt);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        gaussianGradientMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                   pair<DestIterator, DestAccessor> const & dest,
                                   double sigma, const ConvolutionOptions<N> & opt);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\>

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, TinyVector<float, 3> > dest(shape);
    ...
    // compute Gaussian gradient at scale sigma
    gaussianGradientMultiArray(srcMultiArrayRange(source), destMultiArray(dest), sigma);
    \endcode

    <b> Usage with anisotropic data:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\>

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, TinyVector<float, 3> > dest(shape);
    TinyVector<float, 3> step_size;
    TinyVector<float, 3> resolution_sigmas;
    ...
    // compute Gaussian gradient at scale sigma
    gaussianGradientMultiArray(srcMultiArrayRange(source), destMultiArray(dest), sigma,
                               ConvolutionOptions<3>().stepSize(step_size).resolutionStdDev(resolution_sigmas));
    \endcode

    <b> Required Interface:</b>

    see \ref separableConvolveMultiArray(), in addition:

    \code
    int dimension = 0;
    VectorElementAccessor<DestAccessor> elementAccessor(0, dest);
    \endcode

    \see separableConvolveMultiArray()
*/
doxygen_overloaded_function(template <...> void gaussianGradientMultiArray)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
gaussianGradientMultiArray(SrcIterator si, SrcShape const & shape, SrcAccessor src,
                           DestIterator di, DestAccessor dest,
                           ConvolutionOptions<SrcShape::static_size> const & opt,
                           const char *const function_name = "gaussianGradientMultiArray")
{
    typedef typename DestAccessor::value_type DestType;
    typedef typename DestType::value_type     DestValueType;
    typedef typename NumericTraits<DestValueType>::RealPromote KernelType;
   
    static const int N = SrcShape::static_size;
    typedef typename ConvolutionOptions<N>::ScaleIterator ParamType;

    for(int k=0; k<N; ++k)
        if(shape[k] <=0)
            return;

    vigra_precondition(N == (int)dest.size(di),
        "gaussianGradientMultiArray(): Wrong number of channels in output array.");

    ParamType params = opt.scaleParams();
    ParamType params2(params);

    ArrayVector<Kernel1D<KernelType> > plain_kernels(N);
    for (int dim = 0; dim < N; ++dim, ++params)
    {
        double sigma = params.sigma_scaled(function_name);
        plain_kernels[dim].initGaussian(sigma);
    }

    typedef VectorElementAccessor<DestAccessor> ElementAccessor;

    // compute gradient components
    for (int dim = 0; dim < N; ++dim, ++params2)
    {
        ArrayVector<Kernel1D<KernelType> > kernels(plain_kernels);
        kernels[dim].initGaussianDerivative(params2.sigma_scaled(), 1);
        detail::scaleKernel(kernels[dim], 1 / params2.step_size());
        separableConvolveMultiArray(si, shape, src, di, ElementAccessor(dim, dest), kernels.begin());
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
gaussianGradientMultiArray(SrcIterator si, SrcShape const & shape, SrcAccessor src,
                           DestIterator di, DestAccessor dest, double sigma,
                           const ConvolutionOptions<SrcShape::static_size> & opt = ConvolutionOptions<SrcShape::static_size>())
{
    ConvolutionOptions<SrcShape::static_size> par = opt;
    gaussianGradientMultiArray(si, shape, src, di, dest, par.stdDev(sigma));
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
gaussianGradientMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                           pair<DestIterator, DestAccessor> const & dest,
                           ConvolutionOptions<SrcShape::static_size> const & opt )
{
    gaussianGradientMultiArray( source.first, source.second, source.third,
                                dest.first, dest.second, opt );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
gaussianGradientMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                           pair<DestIterator, DestAccessor> const & dest,
                           double sigma,
                           const ConvolutionOptions<SrcShape::static_size> & opt = ConvolutionOptions<SrcShape::static_size>())
{
    gaussianGradientMultiArray( source.first, source.second, source.third,
                                dest.first, dest.second, sigma, opt );
}

/********************************************************/
/*                                                      */
/*             symmetricGradientMultiArray              */
/*                                                      */
/********************************************************/

/** \brief Calculate gradient of a multi-dimensional arrays using symmetric difference filters.

    This function computes the gradient of the given N-dimensional
    array with a sequence of symmetric difference filters a (differentiation is applied
    to each dimension in turn, starting with the innermost dimension). Both source and
    destination arrays are represented by iterators, shape objects and accessors.
    The destination array is required to have a vector valued pixel type with as many
    elements as the number of dimensions. This function is implemented by calls to
    \ref convolveMultiArrayOneDimension() with the symmetric difference kernel.

    Anisotropic data should be passed with appropiate
    \ref ConvolutionOptions, the parameter <tt>opt</tt> is optional
    otherwise.
    
    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        symmetricGradientMultiArray(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest,
                                    const ConvolutionOptions<N> & opt);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        symmetricGradientMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest,
                                    const ConvolutionOptions<N> & opt);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\>

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, TinyVector<float, 3> > dest(shape);
    ...
    // compute gradient
    symmetricGradientMultiArray(srcMultiArrayRange(source), destMultiArray(dest));
    \endcode

    <b> Usage with anisotropic data:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\>

    \code
    MultiArray<3, unsigned char>::size_type shape(width, height, depth);
    MultiArray<3, unsigned char> source(shape);
    MultiArray<3, TinyVector<float, 3> > dest(shape);
    TinyVector<float, 3> step_size;
    ...
    // compute gradient
    symmetricGradientMultiArray(srcMultiArrayRange(source), destMultiArray(dest),
                                ConvolutionOptions<3>().stepSize(step_size));
    \endcode

    <b> Required Interface:</b>

    see \ref convolveMultiArrayOneDimension(), in addition:

    \code
    int dimension = 0;
    VectorElementAccessor<DestAccessor> elementAccessor(0, dest);
    \endcode

    \see convolveMultiArrayOneDimension()
*/
doxygen_overloaded_function(template <...> void symmetricGradientMultiArray)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
symmetricGradientMultiArray(SrcIterator si, SrcShape const & shape, SrcAccessor src,
                            DestIterator di, DestAccessor dest,
                            const ConvolutionOptions<SrcShape::static_size> & opt = ConvolutionOptions<SrcShape::static_size>())
{
    typedef typename DestAccessor::value_type DestType;
    typedef typename DestType::value_type     DestValueType;
    typedef typename NumericTraits<DestValueType>::RealPromote KernelType;

    static const int N = SrcShape::static_size;
    typedef typename ConvolutionOptions<N>::StepIterator StepType;

    for(int k=0; k<N; ++k)
        if(shape[k] <=0)
            return;

    vigra_precondition(N == (int)dest.size(di),
        "symmetricGradientMultiArray(): Wrong number of channels in output array.");

    Kernel1D<KernelType> filter;
    filter.initSymmetricGradient();

    StepType step_size_it = opt.stepParams();

    typedef VectorElementAccessor<DestAccessor> ElementAccessor;

    // compute gradient components
    for (int d = 0; d < N; ++d, ++step_size_it)
    {
        Kernel1D<KernelType> symmetric(filter);
        detail::scaleKernel(symmetric, 1 / *step_size_it);
        convolveMultiArrayOneDimension(si, shape, src,
                                       di, ElementAccessor(d, dest),
                                       d, symmetric);
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
symmetricGradientMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                            pair<DestIterator, DestAccessor> const & dest,
                            const ConvolutionOptions<SrcShape::static_size> & opt = ConvolutionOptions<SrcShape::static_size>())
{
    symmetricGradientMultiArray(source.first, source.second, source.third,
                                dest.first, dest.second, opt);
}

/********************************************************/
/*                                                      */
/*            laplacianOfGaussianMultiArray             */
/*                                                      */
/********************************************************/

/** \brief Calculate Laplacian of a N-dimensional arrays using Gaussian derivative filters.

    This function computes the Laplacian the given N-dimensional
    array with a sequence of second-derivative-of-Gaussian filters at the given
    standard deviation <tt>sigma</tt>. Both source and destination arrays
    are represented by iterators, shape objects and accessors. Both source and destination 
    arrays must have scalar value_type. This function is implemented by calls to
    \ref separableConvolveMultiArray() with the appropriate kernels, followed by summation.

    Anisotropic data should be passed with appropiate
    \ref ConvolutionOptions, the parameter <tt>opt</tt> is otherwise optional
    unless the parameter <tt>sigma</tt> is left out.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        laplacianOfGaussianMultiArray(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                      DestIterator diter, DestAccessor dest,
                                      double sigma, const ConvolutionOptions<N> & opt);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        laplacianOfGaussianMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                      pair<DestIterator, DestAccessor> const & dest,
                                      double sigma, const ConvolutionOptions<N> & opt);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\>

    \code
    MultiArray<3, float> source(shape);
    MultiArray<3, float> laplacian(shape);
    ...
    // compute Laplacian at scale sigma
    laplacianOfGaussianMultiArray(srcMultiArrayRange(source), destMultiArray(laplacian), sigma);
    \endcode

    <b> Usage with anisotropic data:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\>

    \code
    MultiArray<3, float> source(shape);
    MultiArray<3, float> laplacian(shape);
    TinyVector<float, 3> step_size;
    TinyVector<float, 3> resolution_sigmas;
    ...
    // compute Laplacian at scale sigma
    laplacianOfGaussianMultiArray(srcMultiArrayRange(source), destMultiArray(laplacian), sigma,
                                  ConvolutionOptions<3>().stepSize(step_size).resolutionStdDev(resolution_sigmas));
    \endcode

    <b> Required Interface:</b>

    see \ref separableConvolveMultiArray(), in addition:

    \code
    int dimension = 0;
    VectorElementAccessor<DestAccessor> elementAccessor(0, dest);
    \endcode

    \see separableConvolveMultiArray()
*/
doxygen_overloaded_function(template <...> void laplacianOfGaussianMultiArray)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
laplacianOfGaussianMultiArray(SrcIterator si, SrcShape const & shape, SrcAccessor src,
                              DestIterator di, DestAccessor dest,
                              ConvolutionOptions<SrcShape::static_size> const & opt )
{ 
    using namespace functor;
    
    typedef typename DestAccessor::value_type DestType;
    typedef typename NumericTraits<DestType>::RealPromote KernelType;
    typedef typename AccessorTraits<KernelType>::default_accessor DerivativeAccessor;

    static const int N = SrcShape::static_size;
    typedef typename ConvolutionOptions<N>::ScaleIterator ParamType;
    
    ParamType params = opt.scaleParams();
    ParamType params2(params);

    ArrayVector<Kernel1D<KernelType> > plain_kernels(N);
    for (int dim = 0; dim < N; ++dim, ++params)
    {
        double sigma = params.sigma_scaled("laplacianOfGaussianMultiArray");
        plain_kernels[dim].initGaussian(sigma);
    }
    
    MultiArray<N, KernelType> derivative(shape);

    // compute 2nd derivatives and sum them up
    for (int dim = 0; dim < N; ++dim, ++params2)
    {
        ArrayVector<Kernel1D<KernelType> > kernels(plain_kernels);
        kernels[dim].initGaussianDerivative(params2.sigma_scaled(), 2);
        detail::scaleKernel(kernels[dim], 1 / (params2.step_size() * params2.step_size()));

        if (dim == 0)
        {
            separableConvolveMultiArray( si, shape, src, 
                                         di, dest, kernels.begin());
        }
        else
        {
            separableConvolveMultiArray( si, shape, src, 
                                         derivative.traverser_begin(), DerivativeAccessor(), 
                                         kernels.begin());
            combineTwoMultiArrays(di, shape, dest, derivative.traverser_begin(), DerivativeAccessor(), 
                                  di, dest, Arg1() + Arg2() );
        }
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
laplacianOfGaussianMultiArray(SrcIterator si, SrcShape const & shape, SrcAccessor src,
                              DestIterator di, DestAccessor dest, double sigma,
                              const ConvolutionOptions<SrcShape::static_size> & opt = ConvolutionOptions<SrcShape::static_size>())
{
    ConvolutionOptions<SrcShape::static_size> par = opt;
    laplacianOfGaussianMultiArray(si, shape, src, di, dest, par.stdDev(sigma));
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
laplacianOfGaussianMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                              pair<DestIterator, DestAccessor> const & dest,
                              ConvolutionOptions<SrcShape::static_size> const & opt )
{
    laplacianOfGaussianMultiArray( source.first, source.second, source.third,
                                   dest.first, dest.second, opt );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
laplacianOfGaussianMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                              pair<DestIterator, DestAccessor> const & dest,
                              double sigma,
                              const ConvolutionOptions<SrcShape::static_size> & opt = ConvolutionOptions<SrcShape::static_size>())
{
    laplacianOfGaussianMultiArray( source.first, source.second, source.third,
                                   dest.first, dest.second,  sigma, opt );
}

/********************************************************/
/*                                                      */
/*              hessianOfGaussianMultiArray             */
/*                                                      */
/********************************************************/

/** \brief Calculate Hessian matrix of a N-dimensional arrays using Gaussian derivative filters.

    This function computes the Hessian matrix the given scalar N-dimensional
    array with a sequence of second-derivative-of-Gaussian filters at the given
    standard deviation <tt>sigma</tt>. Both source and destination arrays
    are represented by iterators, shape objects and accessors. The destination array must 
    have a vector valued element type with N*(N+1)/2 elements (it represents the
    upper triangular part of the symmetric Hessian matrix). This function is implemented by calls to
    \ref separableConvolveMultiArray() with the appropriate kernels.

    Anisotropic data should be passed with appropiate
    \ref ConvolutionOptions, the parameter <tt>opt</tt> is otherwise optional
    unless the parameter <tt>sigma</tt> is left out.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        hessianOfGaussianMultiArray(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                    DestIterator diter, DestAccessor dest,
                                    double sigma, const ConvolutionOptions<N> & opt);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        hessianOfGaussianMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                    pair<DestIterator, DestAccessor> const & dest,
                                    double sigma, const ConvolutionOptions<N> & opt);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\>

    \code
    MultiArray<3, float> source(shape);
    MultiArray<3, TinyVector<float, 6> > dest(shape);
    ...
    // compute Hessian at scale sigma
    hessianOfGaussianMultiArray(srcMultiArrayRange(source), destMultiArray(dest), sigma);
    \endcode

    <b> Usage with anisotropic data:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\>

    \code
    MultiArray<3, float> source(shape);
    MultiArray<3, TinyVector<float, 6> > dest(shape);
    TinyVector<float, 3> step_size;
    TinyVector<float, 3> resolution_sigmas;
    ...
    // compute Hessian at scale sigma
    hessianOfGaussianMultiArray(srcMultiArrayRange(source), destMultiArray(dest), sigma,
                                ConvolutionOptions<3>().stepSize(step_size).resolutionStdDev(resolution_sigmas));
    \endcode

    <b> Required Interface:</b>

    see \ref separableConvolveMultiArray(), in addition:

    \code
    int dimension = 0;
    VectorElementAccessor<DestAccessor> elementAccessor(0, dest);
    \endcode

    \see separableConvolveMultiArray(), vectorToTensorMultiArray()
*/
doxygen_overloaded_function(template <...> void hessianOfGaussianMultiArray)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
hessianOfGaussianMultiArray(SrcIterator si, SrcShape const & shape, SrcAccessor src,
                            DestIterator di, DestAccessor dest,
                            ConvolutionOptions<SrcShape::static_size> const & opt )
{ 
    typedef typename DestAccessor::value_type DestType;
    typedef typename DestType::value_type     DestValueType;
    typedef typename NumericTraits<DestValueType>::RealPromote KernelType;

    static const int N = SrcShape::static_size;
    static const int M = N*(N+1)/2;
    typedef typename ConvolutionOptions<N>::ScaleIterator ParamType;
    
    for(int k=0; k<N; ++k)
        if(shape[k] <=0)
            return;

    vigra_precondition(M == (int)dest.size(di),
        "hessianOfGaussianMultiArray(): Wrong number of channels in output array.");

    ParamType params_init = opt.scaleParams();

    ArrayVector<Kernel1D<KernelType> > plain_kernels(N);
    ParamType params(params_init);
    for (int dim = 0; dim < N; ++dim, ++params)
    {
        double sigma = params.sigma_scaled("hessianOfGaussianMultiArray");
        plain_kernels[dim].initGaussian(sigma);
    }

    typedef VectorElementAccessor<DestAccessor> ElementAccessor;

    // compute elements of the Hessian matrix
    ParamType params_i(params_init);
    for (int b=0, i=0; i<N; ++i, ++params_i)
    {
        ParamType params_j(params_i);
        for (int j=i; j<N; ++j, ++b, ++params_j)
        {
            ArrayVector<Kernel1D<KernelType> > kernels(plain_kernels);
            if(i == j)
            {
                kernels[i].initGaussianDerivative(params_i.sigma_scaled(), 2);
            }
            else
            {
                kernels[i].initGaussianDerivative(params_i.sigma_scaled(), 1);
                kernels[j].initGaussianDerivative(params_j.sigma_scaled(), 1);
            }
            detail::scaleKernel(kernels[i], 1 / params_i.step_size());
            detail::scaleKernel(kernels[j], 1 / params_j.step_size());
            separableConvolveMultiArray(si, shape, src, di, ElementAccessor(b, dest),
                                        kernels.begin());
        }
    }
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
hessianOfGaussianMultiArray(SrcIterator si, SrcShape const & shape, SrcAccessor src,
                            DestIterator di, DestAccessor dest, double sigma,
                            const ConvolutionOptions<SrcShape::static_size> & opt = ConvolutionOptions<SrcShape::static_size>())
{
    ConvolutionOptions<SrcShape::static_size> par = opt;
    hessianOfGaussianMultiArray(si, shape, src, di, dest, par.stdDev(sigma));
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
hessianOfGaussianMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                            pair<DestIterator, DestAccessor> const & dest,
                            ConvolutionOptions<SrcShape::static_size> const & opt )
{
    hessianOfGaussianMultiArray( source.first, source.second, source.third,
                                 dest.first, dest.second, opt );
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
hessianOfGaussianMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                            pair<DestIterator, DestAccessor> const & dest,
                            double sigma,
                            const ConvolutionOptions<SrcShape::static_size> & opt = ConvolutionOptions<SrcShape::static_size>())
{
    hessianOfGaussianMultiArray( source.first, source.second, source.third,
                                 dest.first, dest.second, sigma, opt );
}

namespace detail {

template<int N, class VectorType>
struct StructurTensorFunctor
{
    typedef VectorType result_type;
    typedef typename VectorType::value_type ValueType;
    
    template <class T>
    VectorType operator()(T const & in) const
    {
        VectorType res;
        for(int b=0, i=0; i<N; ++i)
        {
            for(int j=i; j<N; ++j, ++b)
            {
                res[b] = detail::RequiresExplicitCast<ValueType>::cast(in[i]*in[j]);
            }
        }
        return res;
    }
};

} // namespace detail

/********************************************************/
/*                                                      */
/*               structureTensorMultiArray              */
/*                                                      */
/********************************************************/

/** \brief Calculate th structure tensor of a multi-dimensional arrays.

    This function computes the gradient (outer product) tensor for each element
    of the given N-dimensional array with first-derivative-of-Gaussian filters at 
    the given <tt>innerScale</tt>, followed by Gaussian smoothing at <tt>outerScale</tt>.
    Both source and destination arrays are represented by iterators, shape objects and 
    accessors. The destination array must have a vector valued pixel type with 
    N*(N+1)/2 elements (it represents the upper triangular part of the symmetric 
    structure tensor matrix). If the source array is also vector valued, the 
    resulting structure tensor is the sum of the individual tensors for each channel.
    This function is implemented by calls to
    \ref separableConvolveMultiArray() with the appropriate kernels.

    Anisotropic data should be passed with appropiate
    \ref ConvolutionOptions, the parameter <tt>opt</tt> is otherwise optional
    unless the parameters <tt>innerScale</tt> and <tt>outerScale</tt> are
    both left out.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        structureTensorMultiArray(SrcIterator siter, SrcShape const & shape, SrcAccessor src,
                                  DestIterator diter, DestAccessor dest,
                                  double innerScale, double outerScale,
                                  const ConvolutionOptions<N> & opt);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcShape, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        structureTensorMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                                  pair<DestIterator, DestAccessor> const & dest,
                                  double innerScale, double outerScale,
                                  const ConvolutionOptions<N> & opt);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\>

    \code
    MultiArray<3, RGBValue<float> > source(shape);
    MultiArray<3, TinyVector<float, 6> > dest(shape);
    ...
    // compute structure tensor at scales innerScale and outerScale
    structureTensorMultiArray(srcMultiArrayRange(source), destMultiArray(dest), innerScale, outerScale);
    \endcode

    <b> Usage with anisotropic data:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\>

    \code
    MultiArray<3, RGBValue<float> > source(shape);
    MultiArray<3, TinyVector<float, 6> > dest(shape);
    TinyVector<float, 3> step_size;
    TinyVector<float, 3> resolution_sigmas;
    ...
    // compute structure tensor at scales innerScale and outerScale
    structureTensorMultiArray(srcMultiArrayRange(source), destMultiArray(dest), innerScale, outerScale,
                              ConvolutionOptions<3>().stepSize(step_size).resolutionStdDev(resolution_sigmas));
    \endcode

    <b> Required Interface:</b>

    see \ref separableConvolveMultiArray(), in addition:

    \code
    int dimension = 0;
    VectorElementAccessor<DestAccessor> elementAccessor(0, dest);
    \endcode

    \see separableConvolveMultiArray(), vectorToTensorMultiArray()
*/
doxygen_overloaded_function(template <...> void structureTensorMultiArray)

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
structureTensorMultiArray(SrcIterator si, SrcShape const & shape, SrcAccessor src,
                          DestIterator di, DestAccessor dest, 
                          ConvolutionOptions<SrcShape::static_size> const & opt)
{
    static const int N = SrcShape::static_size;
    static const int M = N*(N+1)/2;
    
    typedef typename DestAccessor::value_type DestType;
    typedef typename DestType::value_type     DestValueType;
    typedef typename NumericTraits<DestValueType>::RealPromote KernelType;
    typedef TinyVector<KernelType, N> GradientVector;
    typedef typename AccessorTraits<GradientVector>::default_accessor GradientAccessor;

    for(int k=0; k<N; ++k)
        if(shape[k] <=0)
            return;

    vigra_precondition(M == (int)dest.size(di),
        "structureTensorMultiArray(): Wrong number of channels in output array.");

    MultiArray<N, GradientVector> gradient(shape);
    gaussianGradientMultiArray(si, shape, src, 
                               gradient.traverser_begin(), GradientAccessor(), 
                               opt,
                               "structureTensorMultiArray");

    transformMultiArray(gradient.traverser_begin(), shape, GradientAccessor(), 
                        di, dest, 
                        detail::StructurTensorFunctor<N, DestType>());

    gaussianSmoothMultiArray(di, shape, dest, di, dest, opt.outerOptions(),
                             "structureTensorMultiArray");
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
structureTensorMultiArray(SrcIterator si, SrcShape const & shape, SrcAccessor src,
                          DestIterator di, DestAccessor dest,
                          double innerScale, double outerScale,
                          const ConvolutionOptions<SrcShape::static_size> & opt = ConvolutionOptions<SrcShape::static_size>())
{
    ConvolutionOptions<SrcShape::static_size> par = opt;
    structureTensorMultiArray(si, shape, src, di, dest,
                              par.stdDev(innerScale).outerScale(outerScale));
}

template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
structureTensorMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                          pair<DestIterator, DestAccessor> const & dest, 
                          ConvolutionOptions<SrcShape::static_size> const & opt )
{
    structureTensorMultiArray( source.first, source.second, source.third,
                               dest.first, dest.second, opt );
}


template <class SrcIterator, class SrcShape, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
structureTensorMultiArray(triple<SrcIterator, SrcShape, SrcAccessor> const & source,
                          pair<DestIterator, DestAccessor> const & dest,
                          double innerScale, double outerScale,
                          const ConvolutionOptions<SrcShape::static_size> & opt = ConvolutionOptions<SrcShape::static_size>())
{
    structureTensorMultiArray( source.first, source.second, source.third,
                               dest.first, dest.second,
                               innerScale, outerScale, opt);
}

//@}

} //-- namespace vigra


#endif        //-- VIGRA_MULTI_CONVOLUTION_H
