/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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


#ifndef VIGRA_CONVOLUTION_HXX
#define VIGRA_CONVOLUTION_HXX

#include <functional>
#include "stdconvolution.hxx"
#include "separableconvolution.hxx"
#include "recursiveconvolution.hxx"
#include "nonlineardiffusion.hxx"
#include "combineimages.hxx"
#include "multi_shape.hxx"


/** \page Convolution Functions to Convolve Images and Signals

    1D and 2D filters, including separable and recursive convolution, and non-linear diffusion

    <b>\#include</b> \<vigra/convolution.hxx\><br>
    Namespace: vigra

    <UL style="list-style-image:url(documents/bullet.gif)">
    <LI> \ref CommonConvolutionFilters
         <BR>&nbsp;&nbsp;&nbsp;<em>Short-hands for many common 2D convolution filters (including normalized convolution)</em>
    <LI> \ref MultiArrayConvolutionFilters
         <BR>&nbsp;&nbsp;&nbsp;<em>Convolution filters for arbitrary dimensional arrays (MultiArray etc.)</em>
    <LI> \ref ResamplingConvolutionFilters
         <BR>&nbsp;&nbsp;&nbsp;<em>Resampling convolution filters</em>
    <LI> \ref vigra::Kernel2D
         <BR>&nbsp;&nbsp;&nbsp;<em>Generic 2-dimensional discrete convolution kernel </em>
    <LI> \ref SeparableConvolution
         <BR>&nbsp;&nbsp;&nbsp;<em>1D convolution and separable filters in 2 dimensions </em>
    <LI> \ref vigra::Kernel1D
         <BR>&nbsp;&nbsp;&nbsp;<em>Generic 1-dimensional discrete convolution kernel </em>
    <LI> \ref RecursiveConvolution
         <BR>&nbsp;&nbsp;&nbsp;<em>Recursive filters (1st and 2nd order)</em>
    <LI> \ref NonLinearDiffusion
         <BR>&nbsp;&nbsp;&nbsp;<em>Edge-preserving smoothing </em>
    <LI> \ref BorderTreatmentMode
         <BR>&nbsp;&nbsp;&nbsp;<em>Choose between different border treatment modes </em>
    <LI> \ref KernelArgumentObjectFactories
         <BR>&nbsp;&nbsp;&nbsp;<em>Factory functions to create argument objects to simplify passing kernels</em>
    </UL>
*/

/** \page KernelArgumentObjectFactories Kernel Argument Object Factories

    These factory functions allow to create argument objects for 1D
    and 2D convolution kernel analogously to
    \ref ArgumentObjectFactories for images.

    \section Kernel1dFactory kernel1d()

        Pass a \ref vigra::Kernel1D to a 1D or separable convolution algorithm.

        These factories can be used to create argument objects when we
        are given instances or subclasses of \ref vigra::Kernel1D
        (analogous to the \ref ArgumentObjectFactories for images).
        These factory functions access <TT>kernel.center()</TT>,
        <TT>kernel.left()</TT>, <TT>kernel.right()</TT>, <TT>kernel.accessor()</TT>,
        and  <TT>kernel.borderTreatment()</TT> to obtain the necessary
        information. The following factory functions are provided:

        <table>
        <tr><th bgcolor="#f0e0c0" colspan=2 align=left>
            <TT>\ref vigra::Kernel1D "vigra::Kernel1D<SomeType>" kernel;</TT>
            </th>
        </tr>
        <tr><td>
        <TT>kernel1d(kernel)</TT>
        </td><td>
            create argument object from information provided by
            kernel

        </td></tr>
        <tr><td>
        <TT>kernel1d(kernel, vigra::BORDER_TREATMENT_CLIP)</TT>
        </td><td>
            create argument object from information provided by
            kernel, but use given border treatment mode

        </td></tr>
        <tr><td>
        <TT>kernel1d(kerneliterator, kernelaccessor,</TT><br>
        <TT>                kernelleft, kernelright,</TT><br>
        <TT>                vigra::BORDER_TREATMENT_CLIP)</TT>
        </td><td>
            create argument object from explicitly given iterator
            (pointing to the center of th kernel), accessor,
            left and right boundaries, and border treatment mode

        </table>

        For usage examples see \ref convolveImage().

    \section Kernel2dFactory kernel2d()

        Pass a \ref vigra::Kernel2D to a 2D (non-separable) convolution algorithm.

        These factories can be used to create argument objects when we
        are given instances or subclasses of \ref vigra::Kernel2D
        (analogous to the \ref ArgumentObjectFactories for images).
        These factory functions access <TT>kernel.center()</TT>,
        <TT>kernel.upperLeft()</TT>, <TT>kernel.lowerRight()</TT>, <TT>kernel.accessor()</TT>,
        and  <TT>kernel.borderTreatment()</TT> to obtain the necessary
        information. The following factory functions are provided:

        <table>
        <tr><th bgcolor="#f0e0c0" colspan=2 align=left>
            <TT>\ref vigra::Kernel2D "vigra::Kernel2D<SomeType>" kernel;</TT>
            </th>
        </tr>
        <tr><td>
        <TT>kernel2d(kernel)</TT>
        </td><td>
            create argument object from information provided by
            kernel

        </td></tr>
        <tr><td>
        <TT>kernel2d(kernel, vigra::BORDER_TREATMENT_CLIP)</TT>
        </td><td>
            create argument object from information provided by
            kernel, but use given border treatment mode

        </td></tr>
        <tr><td>
        <TT>kernel2d(kerneliterator, kernelaccessor,</TT>
        <TT>                upperleft, lowerright,</TT>
        <TT>                vigra::BORDER_TREATMENT_CLIP)</TT>
        </td><td>
            create argument object from explicitly given iterator
            (pointing to the center of th kernel), accessor,
            upper left and lower right corners, and border treatment mode

        </table>

        For usage examples see \ref convolveImage().
*/

namespace vigra {



/********************************************************/
/*                                                      */
/*             Common convolution filters               */
/*                                                      */
/********************************************************/

/** \addtogroup CommonConvolutionFilters Common Filters

    These functions calculate common filters by appropriate sequences of calls 
    to \ref separableConvolveX() and \ref separableConvolveY() or explicit 2-dimensional
    convolution.
*/
//@{

/** \brief Convolve an image with the given kernel(s).

    If you pass \ref vigra::Kernel2D to this function, it will perform an explicit 2-dimensional 
    convolution. If you pass a single \ref vigra::Kernel1D, it performs a separable convolution,
    i.e. it concatenates two 1D convolutions (along the x-axis and along the y-axis) with the same
    kernel via internal calls to \ref separableConvolveX() and \ref separableConvolveY(). If two
    1D kernels are specified, separable convolution uses different kernels for the x- and y-axis.

    All \ref BorderTreatmentMode "border treatment modes" are supported.

    The input pixel type <tt>T1</tt> must be a \ref LinearSpace "linear space" over 
    the kernel's value_type <tt>T</tt>, i.e. addition of source values, multiplication with kernel values,
    and NumericTraits must be defined. The kernel's value_type must be an \ref AlgebraicField "algebraic field",
    i.e. the arithmetic operations (+, -, *, /) and NumericTraits must be defined. Typically, you will use 
    <tt>double</tt> for the kernel type.
    
    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        // use the same 1D kernel for all axes
        template <class T1, class S1,
                  class T2, class S2,
                  class T>
        void
        convolveImage(MultiArrayView<2, T1, S1> const & src,
                      MultiArrayView<2, T2, S2> dest,
                      Kernel1D<T> const & k);

        // use a different kernel for each axis
        template <class T1, class S1,
                  class T2, class S2,
                  class T>
        void
        convolveImage(MultiArrayView<2, T1, S1> const & src,
                      MultiArrayView<2, T2, S2> dest,
                      Kernel1D<T> const & kx, Kernel1D<T> const & ky);
                      
        // use a non-separable 2D kernel
        template <class T1, class S1,
                  class T2, class S2,
                  class T3>
        void
        convolveImage(MultiArrayView<2, T1, S1> const & src,
                      MultiArrayView<2, T2, S2> dest,
                      Kernel2D<T3> const & kernel);
    }
    \endcode

    \deprecatedAPI{convolveImage}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        // use a different kernel for each axis
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class T>
        void convolveImage(SrcIterator supperleft,
                           SrcIterator slowerright, SrcAccessor sa,
                           DestIterator dupperleft, DestAccessor da,
                           Kernel1D<T> const & kx, Kernel1D<T> const & ky);

        // use a non-separable 2D kernel
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class KernelIterator, class KernelAccessor>
        void convolveImage(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                           DestIterator dest_ul, DestAccessor dest_acc,
                           KernelIterator ki, KernelAccessor ak,
                           Diff2D kul, Diff2D klr, BorderTreatmentMode border);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        // use a different kernel for each axis
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class T>
        void
        convolveImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                      pair<DestIterator, DestAccessor> dest,
                      Kernel1D<T> const & kx, Kernel1D<T> const & ky);

        // use a non-separable 2D kernel
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor,
                  class KernelIterator, class KernelAccessor>
        void convolveImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                           pair<DestIterator, DestAccessor> dest,
                           tuple5<KernelIterator, KernelAccessor, Diff2D, Diff2D,
                           BorderTreatmentMode> kernel);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/convolution.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<2, float> src(w,h), dest1(w,h), dest2(w,h);
    ...

    // create horizontal sobel filter (symmetric difference in x-direction, smoothing in y direction)
    Kernel1D<double> kx, ky;
    kx.initSymmetricDifference();
    ky.initBinomial(1);
    
    // calls separable convolution with the two 1D kernels
    convolveImage(src, dest1, kx, ky);
    
    // create a 3x3 Laplacian filter
    Kernel2D<double> laplace;
    laplace.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =
            0.375,  0.25, 0.375,
            0.25,  -2.5,  0.25,
            0.375,  0.25, 0.375;
    
    // calls 2D convolution
    convolveImage(src, dest2, laplace);
    \endcode

    \deprecatedUsage{convolveImage}
    \code
    vigra::FImage src(w,h), dest(w,h);
    ...

    // create horizontal sobel filter (symmetric difference in x-direction, smoothing in y direction)
    Kernel1D<double> kx, ky;
    kx.initSymmetricDifference();
    ky.initBinomial(1);
    
    // calls separable convolution with the two 1D kernels
    vigra::convolveImage(srcImageRange(src), destImage(dest), kx, ky);
    
    // create a 3x3 Laplacian filter
    Kernel2D<double> laplace;
    laplace.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) =
            0.375,  0.25, 0.375,
            0.25,  -2.5,  0.25,
            0.375,  0.25, 0.375;

    // calls 2D convolution
    vigra::convolveImage(srcImageRange(src), destImage(dest), kernel2d(laplace));
    \endcode
    \deprecatedEnd

    <b> Preconditions:</b>

    The image must be larger than the kernel radius. 
    <ul>
    <li>For 1D kernels, <tt>w > std::max(xkernel.right(), -xkernel.keft())</tt> and 
         <tt>h > std::max(ykernel.right(), -ykernel.left())</tt> are required.
    <li>For 2D kernels, <tt>w > std::max(kernel.lowerRight().x, -kernel.upperLeft().x)</tt> and 
         <tt>h > std::max(kernel.lowerRight().y, -kernel.upperLeft().y)</tt> are required.
    </ul>
    If <tt>BORDER_TREATMENT_CLIP</tt> is requested: the sum of kernel elements must be != 0.
*/
doxygen_overloaded_function(template <...> void convolveImage)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class T>
void convolveImage(SrcIterator supperleft,
                   SrcIterator slowerright, SrcAccessor sa,
                   DestIterator dupperleft, DestAccessor da,
                   Kernel1D<T> const & kx, Kernel1D<T> const & ky)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    BasicImage<TmpType> tmp(slowerright - supperleft, SkipInitialization);

    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(kx));
    separableConvolveY(srcImageRange(tmp),
                       destIter(dupperleft, da), kernel1d(ky));
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor,
          class T>
inline void
convolveImage(triple<SrcIterator, SrcIterator, SrcAccessor> src,
              pair<DestIterator, DestAccessor> dest,
              Kernel1D<T> const & kx, Kernel1D<T> const & ky)
{
    convolveImage(src.first, src.second, src.third,
                  dest.first, dest.second, kx, ky);
}

template <class T1, class S1,
          class T2, class S2,
          class T>
inline void
convolveImage(MultiArrayView<2, T1, S1> const & src,
              MultiArrayView<2, T2, S2> dest,
              Kernel1D<T> const & k)
{
    vigra_precondition(src.shape() == dest.shape(),
        "convolveImage(): shape mismatch between input and output.");
    convolveImage(srcImageRange(src),
                  destImage(dest), k, k);
}

template <class T1, class S1,
          class T2, class S2,
          class T>
inline void
convolveImage(MultiArrayView<2, T1, S1> const & src,
              MultiArrayView<2, T2, S2> dest,
              Kernel1D<T> const & kx, Kernel1D<T> const & ky)
{
    vigra_precondition(src.shape() == dest.shape(),
        "convolveImage(): shape mismatch between input and output.");
    convolveImage(srcImageRange(src),
                  destImage(dest), kx, ky);
}

/********************************************************/
/*                                                      */
/*                    simpleSharpening                  */
/*                                                      */
/********************************************************/

/** \brief Perform simple sharpening function.

    This function uses \ref convolveImage() with the following 3x3 filter:
    
    \code
    -sharpening_factor/16.0,    -sharpening_factor/8.0,    -sharpening_factor/16.0,
    -sharpening_factor/8.0,   1.0+sharpening_factor*0.75,  -sharpening_factor/8.0,
    -sharpening_factor/16.0,    -sharpening_factor/8.0,    -sharpening_factor/16.0;    
    \endcode
    
    and uses <TT>BORDER_TREATMENT_REFLECT</TT> as border treatment mode.

    <b> Preconditions:</b>
    \code  
    1. sharpening_factor >= 0
    2. scale >= 0
    \endcode

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1, 
                  class T2, class S2>
        void
        simpleSharpening(MultiArrayView<2, T1, S1> const & src,
                         MultiArrayView<2, T2, S2> dest, 
                         double sharpening_factor);
    }
    \endcode

    \deprecatedAPI{simpleSharpening}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void simpleSharpening(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                              DestIterator dest_ul, DestAccessor dest_acc, double sharpening_factor);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor, 
                  class DestIterator, class DestAccessor>
        void simpleSharpening(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                              pair<DestIterator, DestAccessor> dest, double sharpening_factor);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/convolution.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<2, float> src(w,h), dest(w,h);
    ...

    // sharpening with sharpening_factor = 0.1
    vigra::simpleSharpening(src, dest, 0.1);
    \endcode

    \deprecatedUsage{simpleSharpening}
    \code
    vigra::FImage src(w,h), dest(w,h);
    ...

    // sharpening with sharpening_factor = 0.1
    vigra::simpleSharpening(srcImageRange(src), destImage(dest), 0.1);

    \endcode
    \deprecatedEnd
*/    
doxygen_overloaded_function(template <...> void simpleSharpening)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void simpleSharpening(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                    DestIterator dest_ul, DestAccessor dest_acc, double sharpening_factor)
{

    vigra_precondition(sharpening_factor >= 0.0,
                       "simpleSharpening(): amount of sharpening must be >= 0.");

    Kernel2D<double> kernel;

    kernel.initExplicitly(Diff2D(-1,-1), Diff2D(1,1)) = -sharpening_factor/16.0,    -sharpening_factor/8.0,    -sharpening_factor/16.0,
                                                        -sharpening_factor/8.0,   1.0+sharpening_factor*0.75,  -sharpening_factor/8.0,
                                                        -sharpening_factor/16.0,    -sharpening_factor/8.0,    -sharpening_factor/16.0;

    convolveImage(src_ul, src_lr, src_acc, dest_ul, dest_acc, 
                  kernel.center(), kernel.accessor(), 
                  kernel.upperLeft(), kernel.lowerRight() , BORDER_TREATMENT_REFLECT );
}

template <class SrcIterator, class SrcAccessor, 
          class DestIterator, class DestAccessor>
inline
void simpleSharpening(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                      pair<DestIterator, DestAccessor> dest, double sharpening_factor)
{
    simpleSharpening(src.first, src.second, src.third,
                     dest.first, dest.second, sharpening_factor);
}

template <class T1, class S1, 
          class T2, class S2>
inline void
simpleSharpening(MultiArrayView<2, T1, S1> const & src,
                 MultiArrayView<2, T2, S2> dest, 
                 double sharpening_factor)
{
    vigra_precondition(src.shape() == dest.shape(),
        "simpleSharpening(): shape mismatch between input and output.");
    simpleSharpening(srcImageRange(src),
                     destImage(dest), sharpening_factor);
}


/********************************************************/
/*                                                      */
/*                    gaussianSharpening                */
/*                                                      */
/********************************************************/

/** \brief Perform sharpening function with gaussian filter.


    This function uses \ref gaussianSmoothing() at the given scale to create a
    temporary image 'smooth' and than blends the original and smoothed image 
    according to the formula    

    \code
    dest = (1 + sharpening_factor)*src - sharpening_factor*smooth
    \endcode

    <b> Preconditions:</b>
    \code  
    1. sharpening_factor >= 0
    2. scale >= 0
    \endcode

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        gaussianSharpening(MultiArrayView<2, T1, S1> const & src,
                           MultiArrayView<2, T2, S2> dest, 
                           double sharpening_factor, 
                           double scale);
    }
    \endcode

    \deprecatedAPI{gaussianSharpening}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
      template <class SrcIterator, class SrcAccessor,
                class DestIterator, class DestAccessor>
      void gaussianSharpening(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                              DestIterator dest_ul, DestAccessor dest_acc, 
                              double sharpening_factor, double scale)
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
      template <class SrcIterator, class SrcAccessor,
                class DestIterator, class DestAccessor>
      void gaussianSharpening(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                               pair<DestIterator, DestAccessor> dest, 
                               double sharpening_factor, double scale)
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/convolution.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<2, float> src(w,h), dest(w,h);
    ...

    // sharpening with sharpening_factor = 3.0
    // smoothing with scale = 0.5
    gaussianSharpening(src, dest, 3.0, 0.5);
    \endcode

    \deprecatedUsage{gaussianSharpening}
    \code
    vigra::FImage src(w,h), dest(w,h);
    ...

    // sharpening with sharpening_factor = 3.0
    // smoothing with scale = 0.5
    vigra::gaussianSharpening(srcImageRange(src), destImage(dest), 3.0, 0.5);
    \endcode
    \deprecatedEnd
*/    
doxygen_overloaded_function(template <...> void gaussianSharpening)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void gaussianSharpening(SrcIterator src_ul, SrcIterator src_lr, SrcAccessor src_acc,
                        DestIterator dest_ul, DestAccessor dest_acc, double sharpening_factor, 
                        double scale)
{
    vigra_precondition(sharpening_factor >= 0.0,
                       "gaussianSharpening(): amount of sharpening must be >= 0");
    vigra_precondition(scale >= 0.0,
                       "gaussianSharpening(): scale parameter should be >= 0.");

    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote ValueType;

    BasicImage<ValueType> tmp(src_lr - src_ul, SkipInitialization);

    gaussianSmoothing(src_ul, src_lr, src_acc, tmp.upperLeft(), tmp.accessor(), scale);

    SrcIterator i_src = src_ul;
    DestIterator i_dest = dest_ul;
    typename BasicImage<ValueType>::traverser tmp_ul = tmp.upperLeft();
    typename BasicImage<ValueType>::traverser i_tmp = tmp_ul;
    typename BasicImage<ValueType>::Accessor tmp_acc = tmp.accessor();

    for(; i_src.y != src_lr.y ; i_src.y++, i_dest.y++, i_tmp.y++ )
    {
        for (;i_src.x != src_lr.x ; i_src.x++, i_dest.x++, i_tmp.x++ )
        {
            dest_acc.set((1.0 + sharpening_factor)*src_acc(i_src) - sharpening_factor*tmp_acc(i_tmp), i_dest);
        }
        i_src.x = src_ul.x;
        i_dest.x = dest_ul.x;
        i_tmp.x = tmp_ul.x;
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
gaussianSharpening(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                   pair<DestIterator, DestAccessor> dest, double sharpening_factor, 
                   double scale)
{
    gaussianSharpening(src.first, src.second, src.third,
                       dest.first, dest.second,
                       sharpening_factor, scale);
}

template <class T1, class S1,
          class T2, class S2>
inline void
gaussianSharpening(MultiArrayView<2, T1, S1> const & src,
                   MultiArrayView<2, T2, S2> dest, 
                   double sharpening_factor, 
                   double scale)
{
    vigra_precondition(src.shape() == dest.shape(),
        "gaussianSharpening(): shape mismatch between input and output.");
    gaussianSharpening(srcImageRange(src),
                       destImage(dest),
                       sharpening_factor, scale);
}



/********************************************************/
/*                                                      */
/*                    gaussianSmoothing                 */
/*                                                      */
/********************************************************/

/** \brief Perform isotropic Gaussian convolution.

    This function is a shorthand for the concatenation of a call to
    \ref separableConvolveX() and \ref separableConvolveY() with a
    Gaussian kernel of the given scale. If two scales are provided, 
    smoothing in x and y direction will have different strength. 
    The function uses <TT>BORDER_TREATMENT_REFLECT</TT>. 
    
    Function \ref gaussianSmoothMultiArray() performs the same filter operation
    on arbitrary dimensional arrays.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        gaussianSmoothing(MultiArrayView<2, T1, S1> const & src,
                          MultiArrayView<2, T2, S2> dest,
                          double scale_x, double scale_y = scale_x);
    }
    \endcode

    \deprecatedAPI{gaussianSmoothing}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void gaussianSmoothing(SrcIterator supperleft,
                                SrcIterator slowerright, SrcAccessor sa,
                                DestIterator dupperleft, DestAccessor da,
                                double scale_x, double scale_y = scale_x);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        gaussianSmoothing(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIterator, DestAccessor> dest,
                          double scale_x, double scale_y = scale_x);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/convolution.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<2, float> src(w,h), dest(w,h);
    ...

    // smooth with scale = 3.0
    gaussianSmoothing(src, dest, 3.0);
    \endcode

    \deprecatedUsage{gaussianSmoothing}
    \code
    vigra::FImage src(w,h), dest(w,h);
    ...

    // smooth with scale = 3.0
    vigra::gaussianSmoothing(srcImageRange(src), destImage(dest), 3.0);
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void gaussianSmoothing)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void 
gaussianSmoothing(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor sa,
                  DestIterator dupperleft, DestAccessor da,
                  double scale_x, double scale_y)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    BasicImage<TmpType> tmp(slowerright - supperleft, SkipInitialization);

    Kernel1D<double> smooth_x, smooth_y;
    smooth_x.initGaussian(scale_x);
    smooth_x.setBorderTreatment(BORDER_TREATMENT_REFLECT);
    smooth_y.initGaussian(scale_y);
    smooth_y.setBorderTreatment(BORDER_TREATMENT_REFLECT);

    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(smooth_x));
    separableConvolveY(srcImageRange(tmp),
                       destIter(dupperleft, da), kernel1d(smooth_y));
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void 
gaussianSmoothing(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor sa,
                  DestIterator dupperleft, DestAccessor da,
                  double scale)
{
    gaussianSmoothing(supperleft, slowerright, sa,
                      dupperleft, da,
                      scale, scale);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
gaussianSmoothing(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIterator, DestAccessor> dest,
                  double scale_x, double scale_y)
{
    gaussianSmoothing(src.first, src.second, src.third,
                 dest.first, dest.second, scale_x, scale_y);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
gaussianSmoothing(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIterator, DestAccessor> dest,
                  double scale)
{
    gaussianSmoothing(src.first, src.second, src.third,
                      dest.first, dest.second, scale, scale);
}

template <class T1, class S1,
          class T2, class S2>
inline void
gaussianSmoothing(MultiArrayView<2, T1, S1> const & src,
                  MultiArrayView<2, T2, S2> dest,
                  double scale_x, double scale_y)
{
    vigra_precondition(src.shape() == dest.shape(),
        "gaussianSmoothing(): shape mismatch between input and output.");
    gaussianSmoothing(srcImageRange(src),
                      destImage(dest), scale_x, scale_y);
}

template <class T1, class S1,
          class T2, class S2>
inline void
gaussianSmoothing(MultiArrayView<2, T1, S1> const & src,
                  MultiArrayView<2, T2, S2> dest,
                  double scale)
{
    vigra_precondition(src.shape() == dest.shape(),
        "gaussianSmoothing(): shape mismatch between input and output.");
    gaussianSmoothing(srcImageRange(src),
                      destImage(dest), scale, scale);
}

/********************************************************/
/*                                                      */
/*                     gaussianGradient                 */
/*                                                      */
/********************************************************/

/** \brief Calculate the gradient vector by means of a 1st derivatives of
    Gaussian filter.

    This function is a shorthand for the concatenation of a call to
    \ref separableConvolveX() and \ref separableConvolveY() with the
    appropriate kernels at the given scale. Note that this function can either produce
    two separate result images for the x- and y-components of the gradient, or write
    into a vector valued image (with at least two components).

    Function \ref gaussianGradientMultiArray() performs the same filter operation
    on arbitrary dimensional arrays.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        // write x and y component of the gradient into separate images
        template <class T1, class S1,
                  class T2X, class S2X,
                  class T2Y, class S2Y>
        void
        gaussianGradient(MultiArrayView<2, T1, S1> const & src,
                         MultiArrayView<2, T2X, S2X> destx,
                         MultiArrayView<2, T2Y, S2Y> desty,
                         double scale);

       // write x and y component of the gradient into a vector-valued image
        template <class T1, class S1,
                  class T2, class S2>
        void
        gaussianGradient(MultiArrayView<2, T1, S1> const & src,
                         MultiArrayView<2, TinyVector<T2, 2>, S2> dest,
                         double scale);
    }
    \endcode

    \deprecatedAPI{gaussianGradient}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        // write x and y component of the gradient into separate images
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorX, class DestAccessorX,
                  class DestIteratorY, class DestAccessorY>
        void gaussianGradient(SrcIterator supperleft,
                              SrcIterator slowerright, SrcAccessor sa,
                              DestIteratorX dupperleftx, DestAccessorX dax,
                              DestIteratorY dupperlefty, DestAccessorY day,
                              double scale);

        // write x and y component of the gradient into a vector-valued image
        template <class SrcIterator, class SrcAccessor,
                 class DestIterator, class DestAccessor>
        void gaussianGradient(SrcIterator supperleft,
                              SrcIterator slowerright, SrcAccessor src,
                              DestIterator dupperleft, DestAccessor dest,
                              double scale);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        // write x and y component of the gradient into separate images
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorX, class DestAccessorX,
                  class DestIteratorY, class DestAccessorY>
        void
        gaussianGradient(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                         pair<DestIteratorX, DestAccessorX> destx,
                         pair<DestIteratorY, DestAccessorY> desty,
                         double scale);

        // write x and y component of the gradient into a vector-valued image
        template <class SrcIterator, class SrcAccessor,
                 class DestIterator, class DestAccessor>
        void
        gaussianGradient(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                         pair<DestIterator, DestAccessor> dest,
                         double scale);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/convolution.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<2, float> src(w,h), gradx(w,h), grady(w,h);
    ...

    // calculate gradient vector at scale = 3.0
    gaussianGradient(src, gradx, grady, 3.0);
    
    // likewise, but use a vector image to store the gradient
    MultiArray<2, TinyVector<float, 2> > dest(w,h);
    gaussianGradient(src, dest, 3.0);
    \endcode

    \deprecatedUsage{gaussianGradient}
    \code
    vigra::FImage src(w,h), gradx(w,h), grady(w,h);
    ...

    // calculate gradient vector at scale = 3.0
    vigra::gaussianGradient(srcImageRange(src),
                             destImage(gradx), destImage(grady), 3.0);
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void gaussianGradient)

template <class SrcIterator, class SrcAccessor,
          class DestIteratorX, class DestAccessorX,
          class DestIteratorY, class DestAccessorY>
void gaussianGradient(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor sa,
                        DestIteratorX dupperleftx, DestAccessorX dax,
                        DestIteratorY dupperlefty, DestAccessorY day,
                        double scale)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    BasicImage<TmpType> tmp(slowerright - supperleft, SkipInitialization);

    Kernel1D<double> smooth, grad;
    smooth.initGaussian(scale);
    grad.initGaussianDerivative(scale, 1);

    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(grad));
    separableConvolveY(srcImageRange(tmp),
                       destIter(dupperleftx, dax), kernel1d(smooth));
    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(smooth));
    separableConvolveY(srcImageRange(tmp),
                       destIter(dupperlefty, day), kernel1d(grad));
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void gaussianGradient(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor src,
                        DestIterator dupperleft, DestAccessor dest,
                        double scale)
{
    VectorElementAccessor<DestAccessor> gradx(0, dest), grady(1, dest);
    gaussianGradient(supperleft, slowerright, src, 
                     dupperleft, gradx, dupperleft, grady, scale);
}

template <class SrcIterator, class SrcAccessor,
          class DestIteratorX, class DestAccessorX,
          class DestIteratorY, class DestAccessorY>
inline void
gaussianGradient(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                 pair<DestIteratorX, DestAccessorX> destx,
                 pair<DestIteratorY, DestAccessorY> desty,
                 double scale)
{
    gaussianGradient(src.first, src.second, src.third,
                 destx.first, destx.second, desty.first, desty.second, scale);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
gaussianGradient(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                 pair<DestIterator, DestAccessor> dest,
                 double scale)
{
    gaussianGradient(src.first, src.second, src.third,
                     dest.first, dest.second, scale);
}

template <class T1, class S1,
          class T2X, class S2X,
          class T2Y, class S2Y>
inline void
gaussianGradient(MultiArrayView<2, T1, S1> const & src,
                 MultiArrayView<2, T2X, S2X> destx,
                 MultiArrayView<2, T2Y, S2Y> desty,
                 double scale)
{
    vigra_precondition(src.shape() == destx.shape(),
        "gaussianGradient(): shape mismatch between input and output.");
    gaussianGradient(srcImageRange(src),
                     destImage(destx), destImage(desty), scale);
}

template <class T1, class S1,
          class T2, class S2>
inline void
gaussianGradient(MultiArrayView<2, T1, S1> const & src,
                 MultiArrayView<2, TinyVector<T2, 2>, S2> dest,
                 double scale)
{
    vigra_precondition(src.shape() == dest.shape(),
        "gaussianGradient(): shape mismatch between input and output.");
    gaussianGradient(srcImageRange(src),
                     destImage(dest), scale);
}

/** \brief Calculate the gradient magnitude by means of a 1st derivatives of
    Gaussian filter.

    This function calls gaussianGradient() and returns the pixel-wise magnitude of
    the resulting gradient vectors. If the original image has multiple bands,
    the squared gradient magnitude is computed for each band separately, and the
    return value is the square root of the sum of these squared magnitudes.

    <b> Declarations:</b>

    use arbitrary-dimensional arrays:
    \code
    namespace vigra {
        // pass filter scale explicitly
        template <unsigned int N, class T1, class S1,
                                  class T2, class S2>
        void 
        gaussianGradientMagnitude(MultiArrayView<N, T1, S1> const & src,
                                  MultiArrayView<N, T2, S2> dest,
                                  double sigma,
                                  ConvolutionOptions<N> opt = ConvolutionOptions<N>());

        template <unsigned int N, class MT, class S1,
                                  class T2, class S2>
        void 
        gaussianGradientMagnitude(MultiArrayView<N+1, Multiband<MT>, S1> const & src,
                                  MultiArrayView<N,   T2, S2> dest,
                                  double sigma,
                                  ConvolutionOptions<N> opt = ConvolutionOptions<N>());
                                  
        // pass filter scale(s) in option object
        template <unsigned int N, class T1, class S1,
                                  class T2, class S2>
        void 
        gaussianGradientMagnitude(MultiArrayView<N, T1, S1> const & src,
                                  MultiArrayView<N, T2, S2> dest,
                                  ConvolutionOptions<N> const & opt);

        template <unsigned int N, class MT, class S1,
                                  class T2, class S2>
        void 
        gaussianGradientMagnitude(MultiArrayView<N+1, Multiband<MT>, S1> const & src,
                                  MultiArrayView<N,   T2, S2> dest,
                                  ConvolutionOptions<N> const & opt);
    }
    \endcode
    Here, the input element types <tt>T1</tt> and <tt>MT</tt> can be arbitrary scalar types, and <tt>T1</tt> 
    may also be <tt>TinyVector</tt> or <tt>RGBValue</tt>. The output element type <tt>T2</tt> should 
    be the corresponding norm type (see \ref NormTraits "NormTraits"). In the <tt>Multiband<MT></tt>-version,
    the input array's right-most dimension is interpreted as a channel axis, therefore it must 
    have one dimension more than the output array.

    \deprecatedAPI{gaussianGradientMagnitude}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void gaussianGradientMagnitude(SrcIterator sul,
                                       SrcIterator slr, SrcAccessor src,
                                       DestIterator dupperleft, DestAccessor dest,
                                       double scale);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        gaussianGradientMagnitude(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                  pair<DestIterator, DestAccessor> dest,
                                  double scale);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/multi_convolution.hxx\> (or \<vigra/convolution.hxx\> to use the old API)<br/>
    Namespace: vigra

    \code
    // example 1
    {
        // use a 3-dimensional float array
        MultiArray<3, float> volume(Shape3(w, h, d)), grad(volume.shape());
        ...

        // calculate gradient magnitude at scale = 3.0
        gaussianGradientMagnitude(volume, grad, 3.0);
    }
    
    // example 2
    {
        // use a 2-dimensional RGB array
        MultiArray<2, RGBValue<float> > rgb(Shape2(w, h));
        MultiArray<2, float> grad(rgb.shape());
        ...

        // calculate the color gradient magnitude at scale = 3.0
        gaussianGradientMagnitude(rgb, grad, 3.0);
    }
    
    // example 3
    {
        // use a 3-dimensional array whose right-most axis is interpreted as 
        // a multi-spectral axis with arbitrary many channels
        MultiArray<3, Multiband<float> > spectral(Shape3(w, h, channelCount));
        MultiArray<2, float> grad(Shape2(w, h));
        ...

        // calculate the multi-channel gradient magnitude at scale = 3.0
        // (note that the template parameter N (number of spatial dimensions)
        //  must be provided explicitly as gaussianGradientMagnitude<2>(...) )
        MultiArrayView<3, Multiband<float> > view(spectral);
        gaussianGradientMagnitude<2>(view, grad, 3.0);
    }
    \endcode

    \deprecatedUsage{gaussianGradientMagnitude}
    \code
    // use a traditional float or RGB image
    FImage image(w, h), grad(w, h);
    FRGBImage rgb(w, h);
    ...

    // calculate gradient magnitude at scale = 3.0
    gaussianGradientMagnitude(srcImageRange(image), destImage(grad), 3.0);
    
    // calculate color gradient magnitude at scale = 3.0
    gaussianGradientMagnitude(srcImageRange(rgb), destImage(grad), 3.0);
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void gaussianGradientMagnitude)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void gaussianGradientMagnitude(SrcIterator sul,
                               SrcIterator slr, SrcAccessor src,
                               DestIterator dupperleft, DestAccessor dest,
                               double scale)
{
    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    BasicImage<TmpType> gradx(slr-sul, SkipInitialization), grady(slr-sul, SkipInitialization);

    gaussianGradient(srcIterRange(sul, slr, src),
                     destImage(gradx), destImage(grady), scale);
    combineTwoImages(srcImageRange(gradx), srcImage(grady), destIter(dupperleft, dest),
                     MagnitudeFunctor<TmpType>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
gaussianGradientMagnitude(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIterator, DestAccessor> dest,
                          double scale)
{
    gaussianGradientMagnitude(src.first, src.second, src.third,
                              dest.first, dest.second, scale);
}

/********************************************************/
/*                                                      */
/*                 laplacianOfGaussian                  */
/*                                                      */
/********************************************************/

/** \brief Filter image with the Laplacian of Gaussian operator
    at the given scale.

    This function calls \ref separableConvolveX() and \ref separableConvolveY() with the appropriate 2nd derivative
    of Gaussian kernels in x- and y-direction and then sums the results
    to get the Laplacian.
    
    Function \ref laplacianOfGaussianMultiArray() performs the same filter operation
    on arbitrary dimensional arrays.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        laplacianOfGaussian(MultiArrayView<2, T1, S1> const & src,
                            MultiArrayView<2, T2, S2> dest,
                            double scale);
    }
    \endcode

    \deprecatedAPI{laplacianOfGaussian}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void laplacianOfGaussian(SrcIterator supperleft,
                                SrcIterator slowerright, SrcAccessor sa,
                                DestIterator dupperleft, DestAccessor da,
                                double scale);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        laplacianOfGaussian(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIterator, DestAccessor> dest,
                          double scale);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/convolution.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<2, float> src(w,h), dest(w,h);
    ...

    // calculate Laplacian of Gaussian at scale = 3.0
    laplacianOfGaussian(src, dest, 3.0);
    \endcode

    \deprecatedUsage{laplacianOfGaussian}
    \code
    vigra::FImage src(w,h), dest(w,h);
    ...

    // calculate Laplacian of Gaussian at scale = 3.0
    vigra::laplacianOfGaussian(srcImageRange(src), destImage(dest), 3.0);
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void laplacianOfGaussian)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void laplacianOfGaussian(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor sa,
                        DestIterator dupperleft, DestAccessor da,
                        double scale)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    BasicImage<TmpType> tmp(slowerright - supperleft, SkipInitialization),
                        tmpx(slowerright - supperleft, SkipInitialization),
                        tmpy(slowerright - supperleft, SkipInitialization);

    Kernel1D<double> smooth, deriv;
    smooth.initGaussian(scale);
    deriv.initGaussianDerivative(scale, 2);

    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(deriv));
    separableConvolveY(srcImageRange(tmp),
                       destImage(tmpx), kernel1d(smooth));
    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(smooth));
    separableConvolveY(srcImageRange(tmp),
                       destImage(tmpy), kernel1d(deriv));
    combineTwoImages(srcImageRange(tmpx), srcImage(tmpy),
                       destIter(dupperleft, da), std::plus<TmpType>());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
laplacianOfGaussian(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<DestIterator, DestAccessor> dest,
                    double scale)
{
    laplacianOfGaussian(src.first, src.second, src.third,
                        dest.first, dest.second, scale);
}

template <class T1, class S1,
          class T2, class S2>
inline void
laplacianOfGaussian(MultiArrayView<2, T1, S1> const & src,
                    MultiArrayView<2, T2, S2> dest,
                    double scale)
{
    vigra_precondition(src.shape() == dest.shape(),
        "laplacianOfGaussian(): shape mismatch between input and output.");
    laplacianOfGaussian(srcImageRange(src),
                        destImage(dest), scale);
}

/********************************************************/
/*                                                      */
/*               hessianMatrixOfGaussian                */
/*                                                      */
/********************************************************/

/** \brief Filter image with the 2nd derivatives of the Gaussian
    at the given scale to get the Hessian matrix.

    The Hessian matrix is a symmetric matrix defined as:

    \f[
        \mbox{\rm Hessian}(I) = \left(
        \begin{array}{cc}
        G_{xx} \ast I & G_{xy} \ast I \\
        G_{xy} \ast I & G_{yy} \ast I
        \end{array} \right)
    \f]

    where \f$G_{xx}, G_{xy}, G_{yy}\f$ denote 2nd derivatives of Gaussians
    at the given scale, and
    \f$\ast\f$ is the convolution symbol. This function calls
    \ref separableConvolveX() and \ref separableConvolveY()
    with the appropriate 2nd derivative
    of Gaussian kernels and puts the results in
    the three destination images. The first destination image will
    contain the second derivative in x-direction, the second one the mixed
    derivative, and the third one holds the derivative in y-direction.
    
    Function \ref hessianOfGaussianMultiArray() performs the same filter operation
    on arbitrary dimensional arrays.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class T1, class S1,
                  class T2, class S2>
        void
        hessianMatrixOfGaussian(MultiArrayView<2, T1, S1> const & src,
                                MultiArrayView<2, TinyVector<T2, 3>, S2> dest,
                                double scale);
    }
    \endcode

    \deprecatedAPI{hessianMatrixOfGaussian}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorX, class DestAccessorX,
                  class DestIteratorXY, class DestAccessorXY,
                  class DestIteratorY, class DestAccessorY>
        void hessianMatrixOfGaussian(SrcIterator supperleft,
                                SrcIterator slowerright, SrcAccessor sa,
                                DestIteratorX dupperleftx, DestAccessorX dax,
                                DestIteratorXY dupperleftxy, DestAccessorXY daxy,
                                DestIteratorY dupperlefty, DestAccessorY day,
                                double scale);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorX, class DestAccessorX,
                  class DestIteratorXY, class DestAccessorXY,
                  class DestIteratorY, class DestAccessorY>
        void
        hessianMatrixOfGaussian(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIteratorX, DestAccessorX> destx,
                          pair<DestIteratorXY, DestAccessorXY> destxy,
                          pair<DestIteratorY, DestAccessorY> desty,
                          double scale);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/convolution.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<2, float>                  src(w,h);
    MultiArray<2, TinyVector<float, 3> >  hessian(w,h);  // will hold the three components of the Hessian
    ...

    // calculate Hessian of Gaussian at scale = 3.0, use a 3-band output image
    hessianMatrixOfGaussian(src, hessian, 3.0);
    \endcode

    \deprecatedUsage{hessianMatrixOfGaussian}
    \code
    vigra::FImage src(w,h), 
                  hxx(w,h), hxy(w,h), hyy(w,h); // use a separate image for each component of the Hessian
    ...

    // calculate Hessian of Gaussian at scale = 3.0, use 3 single.band output images
    vigra::hessianMatrixOfGaussian(srcImageRange(src),
                                   destImage(hxx), destImage(hxy), destImage(hyy), 3.0);
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void hessianMatrixOfGaussian)

template <class SrcIterator, class SrcAccessor,
          class DestIteratorX, class DestAccessorX,
          class DestIteratorXY, class DestAccessorXY,
          class DestIteratorY, class DestAccessorY>
void hessianMatrixOfGaussian(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor sa,
                        DestIteratorX dupperleftx, DestAccessorX dax,
                        DestIteratorXY dupperleftxy, DestAccessorXY daxy,
                        DestIteratorY dupperlefty, DestAccessorY day,
                        double scale)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    BasicImage<TmpType> tmp(slowerright - supperleft, SkipInitialization);

    Kernel1D<double> smooth, deriv1, deriv2;
    smooth.initGaussian(scale);
    deriv1.initGaussianDerivative(scale, 1);
    deriv2.initGaussianDerivative(scale, 2);

    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(deriv2));
    separableConvolveY(srcImageRange(tmp),
                       destIter(dupperleftx, dax), kernel1d(smooth));
    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(smooth));
    separableConvolveY(srcImageRange(tmp),
                       destIter(dupperlefty, day), kernel1d(deriv2));
    separableConvolveX(srcIterRange(supperleft, slowerright, sa),
                       destImage(tmp), kernel1d(deriv1));
    separableConvolveY(srcImageRange(tmp),
                       destIter(dupperleftxy, daxy), kernel1d(deriv1));
}

template <class SrcIterator, class SrcAccessor,
          class DestIteratorX, class DestAccessorX,
          class DestIteratorXY, class DestAccessorXY,
          class DestIteratorY, class DestAccessorY>
inline void
hessianMatrixOfGaussian(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                        pair<DestIteratorX, DestAccessorX> destx,
                        pair<DestIteratorXY, DestAccessorXY> destxy,
                        pair<DestIteratorY, DestAccessorY> desty,
                        double scale)
{
    hessianMatrixOfGaussian(src.first, src.second, src.third,
                            destx.first, destx.second,
                            destxy.first, destxy.second,
                            desty.first, desty.second,
                            scale);
}

template <class T1, class S1,
          class T2X, class S2X,
          class T2XY, class S2XY,
          class T2Y, class S2Y>
inline void
hessianMatrixOfGaussian(MultiArrayView<2, T1, S1> const & src,
                        MultiArrayView<2, T2X, S2X> destx,
                        MultiArrayView<2, T2XY, S2XY> destxy,
                        MultiArrayView<2, T2Y, S2Y> desty,
                        double scale)
{
    vigra_precondition(src.shape() == destx.shape() && src.shape() == destxy.shape() && src.shape() == desty.shape(),
        "hessianMatrixOfGaussian(): shape mismatch between input and output.");
    hessianMatrixOfGaussian(srcImageRange(src),
                            destImage(destx),
                            destImage(destxy),
                            destImage(desty),
                            scale);
}

template <class T1, class S1,
          class T2, class S2>
inline void
hessianMatrixOfGaussian(MultiArrayView<2, T1, S1> const & src,
                        MultiArrayView<2, TinyVector<T2, 3>, S2> dest,
                        double scale)
{
    vigra_precondition(src.shape() == dest.shape(),
        "hessianMatrixOfGaussian(): shape mismatch between input and output.");
        
    MultiArrayView<3, T2> expanded(dest.expandElements(0));
    MultiArrayView<2, T2> dxx(expanded.template bind<0>(0));
    MultiArrayView<2, T2> dxy(expanded.template bind<0>(1));
    MultiArrayView<2, T2> dyy(expanded.template bind<0>(2));
    
    hessianMatrixOfGaussian(srcImageRange(src),
                            destImage(dxx),
                            destImage(dxy),
                            destImage(dyy),
                            scale);
}

/********************************************************/
/*                                                      */
/*                   structureTensor                    */
/*                                                      */
/********************************************************/

/** \brief Calculate the Structure Tensor for each pixel of
 and image, using Gaussian (derivative) filters.

    The Structure Tensor is is a smoothed version of the Euclidean product
    of the gradient vector with itself. I.e. it's a symmetric matrix defined as:

    \f[
        \mbox{\rm StructurTensor}(I) = \left(
        \begin{array}{cc}
        G \ast (I_x I_x) & G \ast (I_x I_y) \\
        G \ast (I_x I_y) & G \ast (I_y I_y)
        \end{array} \right) = \left(
        \begin{array}{cc}
        A & C \\
        C & B
        \end{array} \right)
    \f]

    where \f$G\f$ denotes Gaussian smoothing at the <i>outer scale</i>,
    \f$I_x, I_y\f$ are the gradient components taken at the <i>inner scale</i>,
    \f$\ast\f$ is the convolution symbol, and \f$I_x I_x\f$ etc. are pixelwise
    products of the 1st derivative images. This function calls
    \ref separableConvolveX() and \ref separableConvolveY() with the
    appropriate Gaussian kernels and puts the results in
    the three separate destination images (where the first one will
    contain \f$G \ast (I_x I_x)\f$, the second one \f$G \ast (I_x I_y)\f$, and the
    third one holds \f$G \ast (I_y I_y)\f$), or into a single 3-band image (where the bands
    hold the result in the same order as above). The latter form is also applicable when
    the source image is a multi-band image (e.g. RGB). In this case, tensors are
    first computed for each band separately, and then summed up to get a single result tensor.
    
    Function \ref structureTensorMultiArray() performs the same filter operation
    on arbitrary dimensional arrays.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        // create three separate destination images
        template <class T, class S,
                  class TX, class SX,
                  class TXY, class SXY,
                  class TY, class SY>
        void
        structureTensor(MultiArrayView<2, S, T> const & src,
                        MultiArrayView<2, TX, SX>       destx,
                        MultiArrayView<2, TXY, SXY>     destxy,
                        MultiArrayView<2, TY, SY>       desty,
                        double inner_scale, double outer_scale);

        // create a single 3-band destination image
        template <class T1, class S1,
                  class T2, class S2>
        void
        structureTensor(MultiArrayView<2, T1, S1> const & src,
                        MultiArrayView<2, TinyVector<T2, 3>, S2> dest,
                        double inner_scale, double outer_scale);
    }
    \endcode

    \deprecatedAPI{structureTensor}
    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        // create three separate destination images
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorX, class DestAccessorX,
                  class DestIteratorXY, class DestAccessorXY,
                  class DestIteratorY, class DestAccessorY>
        void structureTensor(SrcIterator supperleft,
                                SrcIterator slowerright, SrcAccessor sa,
                                DestIteratorX dupperleftx, DestAccessorX dax,
                                DestIteratorXY dupperleftxy, DestAccessorXY daxy,
                                DestIteratorY dupperlefty, DestAccessorY day,
                                double inner_scale, double outer_scale);

        // create a single 3-band destination image
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void structureTensor(SrcIterator supperleft,
                                SrcIterator slowerright, SrcAccessor sa,
                                DestIterator dupperleft, DestAccessor da,
                                double inner_scale, double outer_scale);
    }
    \endcode
    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        // create three separate destination images
        template <class SrcIterator, class SrcAccessor,
                  class DestIteratorX, class DestAccessorX,
                  class DestIteratorXY, class DestAccessorXY,
                  class DestIteratorY, class DestAccessorY>
        void
        structureTensor(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIteratorX, DestAccessorX> destx,
                          pair<DestIteratorXY, DestAccessorXY> destxy,
                          pair<DestIteratorY, DestAccessorY> desty,
                          double nner_scale, double outer_scale);

        // create a single 3-band destination image
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void
        structureTensor(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                          pair<DestIterator, DestAccessor> dest,
                          double nner_scale, double outer_scale);
    }
    \endcode
    \deprecatedEnd

    <b> Usage:</b>

    <b>\#include</b> \<vigra/convolution.hxx\><br/>
    Namespace: vigra

    \code
    MultiArray<2, flost> src(w,h), 
                         stxx(w,h), stxy(w,h), styy(w,h);  // use a separate image for each component
    ...

    // calculate Structure Tensor at inner scale = 1.0 and outer scale = 3.0
    structureTensor(src, stxx, stxy, styy, 1.0, 3.0);

    // likwise with a single 3-band destination image
    MultiArray<2, TinyVector<float, 3> > st(w,h);
    structureTensor(src, st, 1.0, 3.0);
    \endcode

    \deprecatedUsage{structureTensor}
    \code
    vigra::FImage src(w,h), 
                  stxx(w,h), stxy(w,h), styy(w,h);
    vigra::BasicImage<TinyVector<float, 3> > st(w,h);
    ...

    vigra::structureTensor(srcImageRange(src),
                           destImage(stxx), destImage(stxy), destImage(styy), 1.0, 3.0);

    vigra::structureTensor(srcImageRange(src), destImage(st), 1.0, 3.0);
    \endcode
    \deprecatedEnd
*/
doxygen_overloaded_function(template <...> void structureTensor)

template <class SrcIterator, class SrcAccessor,
          class DestIteratorX, class DestAccessorX,
          class DestIteratorXY, class DestAccessorXY,
          class DestIteratorY, class DestAccessorY>
void structureTensor(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor sa,
                        DestIteratorX dupperleftx, DestAccessorX dax,
                        DestIteratorXY dupperleftxy, DestAccessorXY daxy,
                        DestIteratorY dupperlefty, DestAccessorY day,
                        double inner_scale, double outer_scale)
{
    typedef typename
        NumericTraits<typename SrcAccessor::value_type>::RealPromote
        TmpType;
    BasicImage<TmpType> tmp(slowerright - supperleft, SkipInitialization),
                        tmpx(slowerright - supperleft, SkipInitialization),
                        tmpy(slowerright - supperleft, SkipInitialization);

    gaussianGradient(srcIterRange(supperleft, slowerright, sa),
           destImage(tmpx), destImage(tmpy), inner_scale);
    combineTwoImages(srcImageRange(tmpx), srcImage(tmpx),
                     destImage(tmp), std::multiplies<TmpType>());
    gaussianSmoothing(srcImageRange(tmp),
                      destIter(dupperleftx, dax), outer_scale);
    combineTwoImages(srcImageRange(tmpy), srcImage(tmpy),
                     destImage(tmp), std::multiplies<TmpType>());
    gaussianSmoothing(srcImageRange(tmp),
                      destIter(dupperlefty, day), outer_scale);
    combineTwoImages(srcImageRange(tmpx), srcImage(tmpy),
                     destImage(tmp), std::multiplies<TmpType>());
    gaussianSmoothing(srcImageRange(tmp),
                      destIter(dupperleftxy, daxy), outer_scale);
}

template <class SrcIterator, class SrcAccessor,
          class DestIteratorX, class DestAccessorX,
          class DestIteratorXY, class DestAccessorXY,
          class DestIteratorY, class DestAccessorY>
inline void
structureTensor(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                  pair<DestIteratorX, DestAccessorX> destx,
                  pair<DestIteratorXY, DestAccessorXY> destxy,
                  pair<DestIteratorY, DestAccessorY> desty,
                  double inner_scale, double outer_scale)
{
    structureTensor(src.first, src.second, src.third,
                 destx.first, destx.second,
                 destxy.first, destxy.second,
                 desty.first, desty.second,
                 inner_scale, outer_scale);
}

template <class T, class S,
          class TX, class SX,
          class TXY, class SXY,
          class TY, class SY>
inline void
structureTensor(MultiArrayView<2, S, T> const & src,
                MultiArrayView<2, TX, SX> destx,
                MultiArrayView<2, TXY, SXY> destxy,
                MultiArrayView<2, TY, SY> desty,
                double inner_scale, double outer_scale)
{
    vigra_precondition(src.shape() == destx.shape(),
        "structureTensor(): shape mismatch between input and output.");
    structureTensor(srcImageRange(src),
                    destImage(destx), destImage(destxy), destImage(desty),
                    inner_scale, outer_scale);
}

namespace detail {

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void structureTensor(SrcIterator supperleft,
                     SrcIterator slowerright, SrcAccessor src,
                     DestIterator dupperleft, DestAccessor dest,
                     double inner_scale, double outer_scale,
                     VigraTrueType /* isScalar */)
{
    typedef VectorElementAccessor<DestAccessor> DA;
    structureTensor(supperleft, slowerright, src,
                    dupperleft, DA(0, dest),
                    dupperleft, DA(1, dest),
                    dupperleft, DA(2, dest),
                    inner_scale, outer_scale);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void structureTensor(SrcIterator supperleft,
                     SrcIterator slowerright, SrcAccessor src,
                     DestIterator dupperleft, DestAccessor dest,
                     double inner_scale, double outer_scale,
                     VigraFalseType /* isScalar */)
{
    int bands = src.size(supperleft);
    typedef VectorElementAccessor<SrcAccessor> SA;
    
    structureTensor(supperleft, slowerright, SA(0, src),
                    dupperleft, dest,
                    inner_scale, outer_scale,
                    VigraTrueType() /* isScalar */);
                    
    BasicImage<typename DestAccessor::value_type> st(slowerright - supperleft, SkipInitialization);
    for(int k=1; k < bands; ++k)
    {
        structureTensor(supperleft, slowerright, SA(k, src),
                        st.upperLeft(), st.accessor(),
                        inner_scale, outer_scale,
                        VigraTrueType() /* isScalar */);
        combineTwoImages(srcImageRange(st), srcIter(dupperleft, dest), destIter(dupperleft, dest),
                         std::plus<typename DestAccessor::value_type>());
    }
}

} // namespace detail

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void structureTensor(SrcIterator supperleft,
                        SrcIterator slowerright, SrcAccessor src,
                        DestIterator dupperleft, DestAccessor dest,
                        double inner_scale, double outer_scale)
{
    typedef typename 
        NumericTraits<typename SrcAccessor::value_type>::isScalar isScalar;
    detail::structureTensor(supperleft, slowerright, src,
                            dupperleft, dest, inner_scale, outer_scale, isScalar());
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline void
structureTensor(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                pair<DestIterator, DestAccessor> dest,
                double inner_scale, double outer_scale)
{
    structureTensor(src.first, src.second, src.third,
                    dest.first, dest.second,
                    inner_scale, outer_scale);
}

template <class T1, class S1,
          class T2, class S2>
inline void
structureTensor(MultiArrayView<2, T1, S1> const & src,
                MultiArrayView<2, TinyVector<T2, 3>, S2> dest,
                double inner_scale, double outer_scale)
{
    vigra_precondition(src.shape() == dest.shape(),
        "structureTensor(): shape mismatch between input and output.");
    structureTensor(srcImageRange(src),
                    destImage(dest),
                    inner_scale, outer_scale);
}

//@}

} // namespace vigra

#endif // VIGRA_CONVOLUTION_HXX
