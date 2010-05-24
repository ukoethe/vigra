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


#ifndef VIGRA_BOUNDARYTENSOR_HXX
#define VIGRA_BOUNDARYTENSOR_HXX

#include <cmath>
#include <functional>
#include "utilities.hxx"
#include "array_vector.hxx"
#include "basicimage.hxx"
#include "combineimages.hxx"
#include "numerictraits.hxx"
#include "convolution.hxx"

namespace vigra {

namespace detail {

/***********************************************************************/

typedef ArrayVector<Kernel1D<double> > KernelArray;

template <class KernelArray>
void
initGaussianPolarFilters1(double std_dev, KernelArray & k)
{
    typedef typename KernelArray::value_type Kernel;
    typedef typename Kernel::iterator iterator;

    vigra_precondition(std_dev >= 0.0,
              "initGaussianPolarFilter1(): "
              "Standard deviation must be >= 0.");

    k.resize(4);

    int radius = (int)(4.0*std_dev + 0.5);
    std_dev *= 1.08179074376;
    double f = 1.0 / VIGRA_CSTD::sqrt(2.0 * M_PI) / std_dev;  // norm
    double a = 0.558868151788 / VIGRA_CSTD::pow(std_dev, 5);
    double b = -2.04251639729 / VIGRA_CSTD::pow(std_dev, 3);
    double sigma22 = -0.5 / std_dev / std_dev;


    for(unsigned int i=0; i<k.size(); ++i)
    {
        k[i].initExplicitly(-radius, radius);
        k[i].setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }

    int ix;
    iterator c = k[0].center();
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * VIGRA_CSTD::exp(sigma22 * x * x);
    }

    c = k[1].center();
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * x * VIGRA_CSTD::exp(sigma22 * x * x);
    }

    c = k[2].center();
    double b2 = b / 3.0;
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * (b2 + a * x * x) * VIGRA_CSTD::exp(sigma22 * x * x);
    }

    c = k[3].center();
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * x * (b + a * x * x) * VIGRA_CSTD::exp(sigma22 * x * x);
    }
}

template <class KernelArray>
void
initGaussianPolarFilters2(double std_dev, KernelArray & k)
{
    typedef typename KernelArray::value_type Kernel;
    typedef typename Kernel::iterator iterator;

    vigra_precondition(std_dev >= 0.0,
              "initGaussianPolarFilter2(): "
              "Standard deviation must be >= 0.");

    k.resize(3);

    int radius = (int)(4.0*std_dev + 0.5);
    double f = 1.0 / VIGRA_CSTD::sqrt(2.0 * M_PI) / std_dev;  // norm
    double sigma2 = std_dev*std_dev;
    double sigma22 = -0.5 / sigma2;

    for(unsigned int i=0; i<k.size(); ++i)
    {
        k[i].initExplicitly(-radius, radius);
        k[i].setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }

    int ix;
    iterator c = k[0].center();
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * VIGRA_CSTD::exp(sigma22 * x * x);
    }

    c = k[1].center();
    double f1 = f / sigma2;
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f1 * x * VIGRA_CSTD::exp(sigma22 * x * x);
    }

    c = k[2].center();
    double f2 = f / (sigma2 * sigma2);
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f2 * (x * x - sigma2) * VIGRA_CSTD::exp(sigma22 * x * x);
    }
}

template <class KernelArray>
void
initGaussianPolarFilters3(double std_dev, KernelArray & k)
{
    typedef typename KernelArray::value_type Kernel;
    typedef typename Kernel::iterator iterator;

    vigra_precondition(std_dev >= 0.0,
              "initGaussianPolarFilter3(): "
              "Standard deviation must be >= 0.");

    k.resize(4);

    int radius = (int)(4.0*std_dev + 0.5);
    std_dev *= 1.15470053838;
    double sigma22 = -0.5 / std_dev / std_dev;
    double f = 1.0 / VIGRA_CSTD::sqrt(2.0 * M_PI) / std_dev;  // norm
    double a = 0.883887052922 / VIGRA_CSTD::pow(std_dev, 5);

    for(unsigned int i=0; i<k.size(); ++i)
    {
        k[i].initExplicitly(-radius, radius);
        k[i].setBorderTreatment(BORDER_TREATMENT_REFLECT);
    }

    //double b = -1.3786348292 / VIGRA_CSTD::pow(std_dev, 3);

    int ix;
    iterator c = k[0].center();
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * VIGRA_CSTD::exp(sigma22 * x * x);
    }

    c = k[1].center();
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * x * VIGRA_CSTD::exp(sigma22 * x * x);
    }

    c = k[2].center();
    double a2 = 3.0 * a;
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * a2 * x * x * VIGRA_CSTD::exp(sigma22 * x * x);
    }

    c = k[3].center();
    for(ix=-radius; ix<=radius; ++ix)
    {
        double x = (double)ix;
        c[ix] = f * a * x * x * x * VIGRA_CSTD::exp(sigma22 * x * x);
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
evenPolarFilters(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor src,
                 DestIterator dupperleft, DestAccessor dest,
                 double scale, bool noLaplacian)
{
    vigra_precondition(dest.size(dupperleft) == 3,
                       "evenPolarFilters(): image for even output must have 3 bands.");

    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;

    typedef typename
       NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    typedef BasicImage<TinyVector<TmpType, 3> > TmpImage;
    typedef typename TmpImage::traverser TmpTraverser;
    TmpImage t(w, h);

    KernelArray k2;
    initGaussianPolarFilters2(scale, k2);

    // calculate filter responses for even filters
    VectorElementAccessor<typename TmpImage::Accessor> tmpBand(0, t.accessor());
    convolveImage(srcIterRange(supperleft, slowerright, src),
                  destImage(t, tmpBand), k2[2], k2[0]);
    tmpBand.setIndex(1);
    convolveImage(srcIterRange(supperleft, slowerright, src),
                  destImage(t, tmpBand), k2[1], k2[1]);
    tmpBand.setIndex(2);
    convolveImage(srcIterRange(supperleft, slowerright, src),
                  destImage(t, tmpBand), k2[0], k2[2]);

    // create even tensor from filter responses
    TmpTraverser tul(t.upperLeft());
    TmpTraverser tlr(t.lowerRight());
    for(; tul.y != tlr.y; ++tul.y, ++dupperleft.y)
    {
        typename TmpTraverser::row_iterator tr = tul.rowIterator();
        typename TmpTraverser::row_iterator trend = tr + w;
        typename DestIterator::row_iterator d = dupperleft.rowIterator();
        if(noLaplacian)
        {
            for(; tr != trend; ++tr, ++d)
            {
                TmpType v = detail::RequiresExplicitCast<TmpType>::cast(0.5*sq((*tr)[0]-(*tr)[2]) + 2.0*sq((*tr)[1]));
                dest.setComponent(v, d, 0);
                dest.setComponent(0, d, 1);
                dest.setComponent(v, d, 2);
            }
        }
        else
        {
            for(; tr != trend; ++tr, ++d)
            {
                dest.setComponent(sq((*tr)[0]) + sq((*tr)[1]), d, 0);
                dest.setComponent(-(*tr)[1] * ((*tr)[0] + (*tr)[2]), d, 1);
                dest.setComponent(sq((*tr)[1]) + sq((*tr)[2]), d, 2);
            }
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void
oddPolarFilters(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor src,
                DestIterator dupperleft, DestAccessor dest,
                double scale, bool addResult)
{
    vigra_precondition(dest.size(dupperleft) == 3,
                       "oddPolarFilters(): image for odd output must have 3 bands.");

    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;

    typedef typename
       NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    typedef BasicImage<TinyVector<TmpType, 4> > TmpImage;
    typedef typename TmpImage::traverser TmpTraverser;
    TmpImage t(w, h);

    detail::KernelArray k1;
    detail::initGaussianPolarFilters1(scale, k1);

    // calculate filter responses for odd filters
    VectorElementAccessor<typename TmpImage::Accessor> tmpBand(0, t.accessor());
    convolveImage(srcIterRange(supperleft, slowerright, src),
                  destImage(t, tmpBand), k1[3], k1[0]);
    tmpBand.setIndex(1);
    convolveImage(srcIterRange(supperleft, slowerright, src),
                  destImage(t, tmpBand), k1[2], k1[1]);
    tmpBand.setIndex(2);
    convolveImage(srcIterRange(supperleft, slowerright, src),
                  destImage(t, tmpBand), k1[1], k1[2]);
    tmpBand.setIndex(3);
    convolveImage(srcIterRange(supperleft, slowerright, src),
                  destImage(t, tmpBand), k1[0], k1[3]);

    // create odd tensor from filter responses
    TmpTraverser tul(t.upperLeft());
    TmpTraverser tlr(t.lowerRight());
    for(; tul.y != tlr.y; ++tul.y, ++dupperleft.y)
    {
        typename TmpTraverser::row_iterator tr = tul.rowIterator();
        typename TmpTraverser::row_iterator trend = tr + w;
        typename DestIterator::row_iterator d = dupperleft.rowIterator();
        if(addResult)
        {
            for(; tr != trend; ++tr, ++d)
            {
                TmpType d0 = (*tr)[0] + (*tr)[2];
                TmpType d1 = -(*tr)[1] - (*tr)[3];

                dest.setComponent(dest.getComponent(d, 0) + sq(d0), d, 0);
                dest.setComponent(dest.getComponent(d, 1) + d0 * d1, d, 1);
                dest.setComponent(dest.getComponent(d, 2) + sq(d1), d, 2);
            }
        }
        else
        {
            for(; tr != trend; ++tr, ++d)
            {
                TmpType d0 = (*tr)[0] + (*tr)[2];
                TmpType d1 = -(*tr)[1] - (*tr)[3];

                dest.setComponent(sq(d0), d, 0);
                dest.setComponent(d0 * d1, d, 1);
                dest.setComponent(sq(d1), d, 2);
            }
        }
    }
}

} // namespace detail

/** \addtogroup CommonConvolutionFilters Common Filters
*/
//@{

/********************************************************/
/*                                                      */
/*                   rieszTransformOfLOG                */
/*                                                      */
/********************************************************/

/** \brief Calculate Riesz transforms of the Laplacian of Gaussian.

    The Riesz transforms of the Laplacian of Gaussian have the following transfer
    functions (defined in a polar coordinate representation of the frequency domain):

    \f[
        F_{\sigma}(r, \phi)=(i \cos \phi)^n (i \sin \phi)^m r^2 e^{-r^2 \sigma^2 / 2}
    \f]

    where <i>n</i> = <tt>xorder</tt> and <i>m</i> = <tt>yorder</tt> determine th e
    order of the transform, and <tt>sigma > 0</tt> is the scale of the Laplacian
    of Gaussian. This function computes a good spatial domain approximation of
    these transforms for <tt>xorder + yorder <= 2</tt>. The filter responses may be used
    to calculate the monogenic signal or the boundary tensor.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                class DestIterator, class DestAccessor>
        void rieszTransformOfLOG(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor src,
                                 DestIterator dupperleft, DestAccessor dest,
                                 double scale, unsigned int xorder, unsigned int yorder);
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                class DestIterator, class DestAccessor>
        void rieszTransformOfLOG(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                                 pair<DestIterator, DestAccessor> dest,
                                 double scale, unsigned int xorder, unsigned int yorder);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<<a href="boundarytensor_8hxx-source.html">vigra/boundarytensor.hxx</a>\>

    \code
    FImage impulse(17,17), res(17, 17);
    impulse(8,8) = 1.0;

    // calculate the impulse response of the first order Riesz transform in x-direction
    rieszTransformOfLOG(srcImageRange(impulse), destImage(res), 2.0, 1, 0);
    \endcode

*/
doxygen_overloaded_function(template <...> void rieszTransformOfLOG)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void rieszTransformOfLOG(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor src,
                         DestIterator dupperleft, DestAccessor dest,
                         double scale, unsigned int xorder, unsigned int yorder)
{
    unsigned int order = xorder + yorder;

    vigra_precondition(order <= 2,
            "rieszTransformOfLOG(): can only compute Riesz transforms up to order 2.");
    vigra_precondition(scale > 0.0,
            "rieszTransformOfLOG(): scale must be positive.");

    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;

    typedef typename NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    typedef BasicImage<TmpType> TmpImage;

    switch(order)
    {
        case 0:
        {
            detail::KernelArray k2;
            detail::initGaussianPolarFilters2(scale, k2);

            TmpImage t1(w, h), t2(w, h);

            convolveImage(srcIterRange(supperleft, slowerright, src),
                          destImage(t1), k2[2], k2[0]);
            convolveImage(srcIterRange(supperleft, slowerright, src),
                          destImage(t2), k2[0], k2[2]);
            combineTwoImages(srcImageRange(t1), srcImage(t2),
                             destIter(dupperleft, dest), std::plus<TmpType>());
            break;
        }
        case 1:
        {
            detail::KernelArray k1;
            detail::initGaussianPolarFilters1(scale, k1);

            TmpImage t1(w, h), t2(w, h);

            if(xorder == 1)
            {
                convolveImage(srcIterRange(supperleft, slowerright, src),
                            destImage(t1), k1[3], k1[0]);
                convolveImage(srcIterRange(supperleft, slowerright, src),
                            destImage(t2), k1[1], k1[2]);
            }
            else
            {
                convolveImage(srcIterRange(supperleft, slowerright, src),
                            destImage(t1), k1[0], k1[3]);
                convolveImage(srcIterRange(supperleft, slowerright, src),
                            destImage(t2), k1[2], k1[1]);
            }
            combineTwoImages(srcImageRange(t1), srcImage(t2),
                             destIter(dupperleft, dest), std::plus<TmpType>());
            break;
        }
        case 2:
        {
            detail::KernelArray k2;
            detail::initGaussianPolarFilters2(scale, k2);

            convolveImage(srcIterRange(supperleft, slowerright, src),
                          destIter(dupperleft, dest), k2[xorder], k2[yorder]);
            break;
        }
        /* for test purposes only: compute 3rd order polar filters */
        case 3:
        {
            detail::KernelArray k3;
            detail::initGaussianPolarFilters3(scale, k3);

            TmpImage t1(w, h), t2(w, h);

            if(xorder == 3)
            {
                convolveImage(srcIterRange(supperleft, slowerright, src),
                            destImage(t1), k3[3], k3[0]);
                convolveImage(srcIterRange(supperleft, slowerright, src),
                            destImage(t2), k3[1], k3[2]);
            }
            else
            {
                convolveImage(srcIterRange(supperleft, slowerright, src),
                            destImage(t1), k3[0], k3[3]);
                convolveImage(srcIterRange(supperleft, slowerright, src),
                            destImage(t2), k3[2], k3[1]);
            }
            combineTwoImages(srcImageRange(t1), srcImage(t2),
                             destIter(dupperleft, dest), std::minus<TmpType>());
            break;
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void rieszTransformOfLOG(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                         pair<DestIterator, DestAccessor> dest,
                         double scale, unsigned int xorder, unsigned int yorder)
{
    rieszTransformOfLOG(src.first, src.second, src.third, dest.first, dest.second,
                        scale, xorder, yorder);
}
//@}

/** \addtogroup TensorImaging Tensor Image Processing
*/
//@{

/********************************************************/
/*                                                      */
/*                     boundaryTensor                   */
/*                                                      */
/********************************************************/

/** \brief Calculate the boundary tensor for a scalar valued image.

    These functions calculate a spatial domain approximation of
    the boundary tensor as described in

    U. K&ouml;the: <a href="http://hci.iwr.uni-heidelberg.de/people/ukoethe/papers/index.php#cite_polarfilters">
    <i>"Integrated Edge and Junction Detection with the Boundary Tensor"</i></a>,
     in: ICCV 03, Proc. of 9th Intl. Conf. on Computer Vision, Nice 2003, vol. 1,
     pp. 424-431, Los Alamitos: IEEE Computer Society, 2003

    with the Laplacian of Gaussian as the underlying bandpass filter (see
    \ref rieszTransformOfLOG()). The output image must have 3 bands which will hold the
    tensor components in the order t11, t12 (== t21), t22. The function
    \ref boundaryTensor1() with the same interface implements a variant of the
    boundary tensor where the 0th-order Riesz transform has been dropped, so that the
    tensor is no longer sensitive to blobs.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void boundaryTensor(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor src,
                            DestIterator dupperleft, DestAccessor dest,
                            double scale);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void boundaryTensor(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                            pair<DestIterator, DestAccessor> dest,
                            double scale);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<<a href="boundarytensor_8hxx-source.html">vigra/boundarytensor.hxx</a>\>

    \code
    FImage img(w,h);
    FVector3Image bt(w,h);
    ...
    boundaryTensor(srcImageRange(img), destImage(bt), 2.0);
    \endcode

*/
doxygen_overloaded_function(template <...> void boundaryTensor)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void boundaryTensor(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor src,
                    DestIterator dupperleft, DestAccessor dest,
                    double scale)
{
    vigra_precondition(dest.size(dupperleft) == 3,
                       "boundaryTensor(): image for even output must have 3 bands.");
    vigra_precondition(scale > 0.0,
                       "boundaryTensor(): scale must be positive.");

    detail::evenPolarFilters(supperleft, slowerright, src,
                             dupperleft, dest, scale, false);
    detail::oddPolarFilters(supperleft, slowerright, src,
                             dupperleft, dest, scale, true);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void boundaryTensor(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                    pair<DestIterator, DestAccessor> dest,
                    double scale)
{
    boundaryTensor(src.first, src.second, src.third,
                   dest.first, dest.second, scale);
}

/** \brief Boundary tensor variant.

    This function implements a variant of the boundary tensor where the 
    0th-order Riesz transform has been dropped, so that the tensor is no 
    longer sensitive to blobs. See \ref boundaryTensor() for more detailed 
    documentation.

    <b> Declarations:</b>

    <b>\#include</b> \<<a href="boundarytensor_8hxx-source.html">vigra/boundarytensor.hxx</a>\>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void boundaryTensor1(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor src,
                             DestIterator dupperleft, DestAccessor dest,
                             double scale);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcIterator, class SrcAccessor,
                  class DestIterator, class DestAccessor>
        void boundaryTensor1(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                             pair<DestIterator, DestAccessor> dest,
                             double scale);
    }
    \endcode
*/
doxygen_overloaded_function(template <...> void boundaryTensor1)

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
void boundaryTensor1(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor src,
                    DestIterator dupperleft, DestAccessor dest,
                    double scale)
{
    vigra_precondition(dest.size(dupperleft) == 3,
                       "boundaryTensor1(): image for even output must have 3 bands.");
    vigra_precondition(scale > 0.0,
                       "boundaryTensor1(): scale must be positive.");

    detail::evenPolarFilters(supperleft, slowerright, src,
                             dupperleft, dest, scale, true);
    detail::oddPolarFilters(supperleft, slowerright, src,
                             dupperleft, dest, scale, true);
}

template <class SrcIterator, class SrcAccessor,
          class DestIterator, class DestAccessor>
inline
void boundaryTensor1(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                     pair<DestIterator, DestAccessor> dest,
                     double scale)
{
    boundaryTensor1(src.first, src.second, src.third,
                    dest.first, dest.second, scale);
}

/********************************************************/
/*                                                      */
/*                    boundaryTensor3                   */
/*                                                      */
/********************************************************/

/*  Add 3rd order Riesz transform to boundary tensor
    ??? Does not work -- bug or too coarse approximation for 3rd order ???
*/
template <class SrcIterator, class SrcAccessor,
          class DestIteratorEven, class DestAccessorEven,
          class DestIteratorOdd, class DestAccessorOdd>
void boundaryTensor3(SrcIterator supperleft, SrcIterator slowerright, SrcAccessor sa,
                     DestIteratorEven dupperleft_even, DestAccessorEven even,
                     DestIteratorOdd dupperleft_odd, DestAccessorOdd odd,
                     double scale)
{
    vigra_precondition(even.size(dupperleft_even) == 3,
                       "boundaryTensor3(): image for even output must have 3 bands.");
    vigra_precondition(odd.size(dupperleft_odd) == 3,
                       "boundaryTensor3(): image for odd output must have 3 bands.");

    detail::evenPolarFilters(supperleft, slowerright, sa,
                             dupperleft_even, even, scale, false);

    int w = slowerright.x - supperleft.x;
    int h = slowerright.y - supperleft.y;

    typedef typename
       NumericTraits<typename SrcAccessor::value_type>::RealPromote TmpType;
    typedef BasicImage<TinyVector<TmpType, 4> > TmpImage;
    TmpImage t1(w, h), t2(w, h);

    detail::KernelArray k1, k3;
    detail::initGaussianPolarFilters1(scale, k1);
    detail::initGaussianPolarFilters3(scale, k3);

    // calculate filter responses for odd filters
    VectorElementAccessor<typename TmpImage::Accessor> tmpBand(0, t1.accessor());
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t1, tmpBand), k1[3], k1[0]);
    tmpBand.setIndex(1);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t1, tmpBand), k1[1], k1[2]);
    tmpBand.setIndex(2);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t1, tmpBand), k3[3], k3[0]);
    tmpBand.setIndex(3);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t1, tmpBand), k3[1], k3[2]);

    tmpBand.setIndex(0);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t2, tmpBand), k1[0], k1[3]);
    tmpBand.setIndex(1);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t2, tmpBand), k1[2], k1[1]);
    tmpBand.setIndex(2);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t2, tmpBand), k3[0], k3[3]);
    tmpBand.setIndex(3);
    convolveImage(srcIterRange(supperleft, slowerright, sa),
                  destImage(t2, tmpBand), k3[2], k3[1]);

    // create odd tensor from filter responses
    typedef typename TmpImage::traverser TmpTraverser;
    TmpTraverser tul1(t1.upperLeft());
    TmpTraverser tlr1(t1.lowerRight());
    TmpTraverser tul2(t2.upperLeft());
    for(; tul1.y != tlr1.y; ++tul1.y, ++tul2.y, ++dupperleft_odd.y)
    {
        typename TmpTraverser::row_iterator tr1 = tul1.rowIterator();
        typename TmpTraverser::row_iterator trend1 = tr1 + w;
        typename TmpTraverser::row_iterator tr2 = tul2.rowIterator();
        typename DestIteratorOdd::row_iterator o = dupperleft_odd.rowIterator();
        for(; tr1 != trend1; ++tr1, ++tr2, ++o)
        {
            TmpType d11 =  (*tr1)[0] + (*tr1)[2];
            TmpType d12 = -(*tr1)[1] - (*tr1)[3];
            TmpType d31 =  (*tr2)[0] - (*tr2)[2];
            TmpType d32 =  (*tr2)[1] - (*tr2)[3];
            TmpType d111 = 0.75 * d11 + 0.25 * d31;
            TmpType d112 = 0.25 * (d12 + d32);
            TmpType d122 = 0.25 * (d11 - d31);
            TmpType d222 = 0.75 * d12 - 0.25 * d32;
            TmpType d2 = sq(d112);
            TmpType d3 = sq(d122);

            odd.setComponent(0.25 * (sq(d111) + 2.0*d2 + d3), o, 0);
            odd.setComponent(0.25 * (d111*d112 + 2.0*d112*d122 + d122*d222), o, 1);
            odd.setComponent(0.25 * (d2 + 2.0*d3 + sq(d222)), o, 2);
        }
    }
}

template <class SrcIterator, class SrcAccessor,
          class DestIteratorEven, class DestAccessorEven,
          class DestIteratorOdd, class DestAccessorOdd>
inline
void boundaryTensor3(triple<SrcIterator, SrcIterator, SrcAccessor> src,
                     pair<DestIteratorEven, DestAccessorEven> even,
                     pair<DestIteratorOdd, DestAccessorOdd> odd,
                     double scale)
{
    boundaryTensor3(src.first, src.second, src.third,
                    even.first, even.second, odd.first, odd.second, scale);
}

//@}

} // namespace vigra

#endif // VIGRA_BOUNDARYTENSOR_HXX
