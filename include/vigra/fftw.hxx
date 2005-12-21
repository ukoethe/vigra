/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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

#ifndef VIGRA_FFTW_HXX
#define VIGRA_FFTW_HXX

#include <cmath>
#include <functional>
#include "vigra/stdimage.hxx"
#include "vigra/copyimage.hxx"
#include "vigra/transformimage.hxx"
#include "vigra/combineimages.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/imagecontainer.hxx"
#include <fftw.h>

namespace vigra {

/********************************************************/
/*                                                      */
/*                    FFTWComplex                       */
/*                                                      */
/********************************************************/

/* documentation: see fftw3.hxx
*/
class FFTWComplex
: public fftw_complex
{
  public:
        /** The complex' component type, as defined in '<TT>fftw.h</TT>'
        */
    typedef fftw_real value_type;

        /** reference type (result of operator[])
        */
    typedef fftw_real & reference;

        /** const reference type (result of operator[] const)
        */
    typedef fftw_real const & const_reference;

        /** iterator type (result of begin() )
        */
    typedef fftw_real * iterator;

        /** const iterator type (result of begin() const)
        */
    typedef fftw_real const * const_iterator;

        /** The norm type (result of magnitde())
        */
    typedef fftw_real NormType;

        /** The squared norm type (result of squaredMagnitde())
        */
    typedef fftw_real SquaredNormType;

        /** Construct from real and imaginary part.
            Default: 0.
        */
    FFTWComplex(value_type const & ire = 0.0, value_type const & iim = 0.0)
    {
        re = ire;
        im = iim;
    }

        /** Copy constructor.
        */
    FFTWComplex(FFTWComplex const & o)
    : fftw_complex(o)
    {}

        /** Construct from plain <TT>fftw_complex</TT>.
        */
    FFTWComplex(fftw_complex const & o)
    : fftw_complex(o)
    {}

        /** Construct from TinyVector.
        */
    template <class T>
    FFTWComplex(TinyVector<T, 2> const & o)
    {
        re = o[0];
        im = o[1];
    }

        /** Assignment.
        */
    FFTWComplex& operator=(FFTWComplex const & o)
    {
        re = o.re;
        im = o.im;
        return *this;
    }

        /** Assignment.
        */
    FFTWComplex& operator=(fftw_complex const & o)
    {
        re = o.re;
        im = o.im;
        return *this;
    }

        /** Assignment.
        */
    FFTWComplex& operator=(fftw_real const & o)
    {
        re = o;
        im = 0.0;
        return *this;
    }

        /** Assignment.
        */
    template <class T>
    FFTWComplex& operator=(TinyVector<T, 2> const & o)
    {
        re = o[0];
        im = o[1];
        return *this;
    }

        /** Unary negation.
        */
    FFTWComplex operator-() const
        { return FFTWComplex(-re, -im); }

        /** Squared magnitude x*conj(x)
        */
    SquaredNormType squaredMagnitude() const
        { return c_re(*this)*c_re(*this)+c_im(*this)*c_im(*this); }

        /** Magnitude (length of radius vector).
        */
    NormType magnitude() const
        { return VIGRA_CSTD::sqrt(squaredMagnitude()); }

        /** Phase angle.
        */
    value_type phase() const
        { return VIGRA_CSTD::atan2(c_im(*this),c_re(*this)); }

        /** Access components as if number were a vector.
        */
    reference operator[](int i)
        { return (&re)[i]; }

        /** Read components as if number were a vector.
        */
    const_reference operator[](int i) const
        { return (&re)[i]; }

        /** Length of complex number (always 2).
        */
    int size() const
        { return 2; }

    iterator begin()
        { return &re; }

    iterator end()
        { return &re + 2; }

    const_iterator begin() const
        { return &re; }

    const_iterator end() const
        { return &re + 2; }
};

/********************************************************/
/*                                                      */
/*                    FFTWComplex Traits                */
/*                                                      */
/********************************************************/

/* documentation: see fftw3.hxx
*/
template<>
struct NumericTraits<fftw_complex>
{
    typedef fftw_complex Type;
    typedef fftw_complex Promote;
    typedef fftw_complex RealPromote;
    typedef fftw_complex ComplexPromote;
    typedef fftw_real    ValueType;

    typedef VigraFalseType isIntegral;
    typedef VigraFalseType isScalar;
    typedef NumericTraits<fftw_real>::isSigned isSigned;
    typedef VigraFalseType isOrdered;
    typedef VigraTrueType  isComplex;

    static FFTWComplex zero() { return FFTWComplex(0.0, 0.0); }
    static FFTWComplex one() { return FFTWComplex(1.0, 0.0); }
    static FFTWComplex nonZero() { return one(); }

    static const Promote & toPromote(const Type & v) { return v; }
    static const RealPromote & toRealPromote(const Type & v) { return v; }
    static const Type & fromPromote(const Promote & v) { return v; }
    static const Type & fromRealPromote(const RealPromote & v) { return v; }
};

template<>
struct NumericTraits<FFTWComplex>
: public NumericTraits<fftw_complex>
{
    typedef FFTWComplex Type;
    typedef FFTWComplex Promote;
    typedef FFTWComplex RealPromote;
    typedef FFTWComplex ComplexPromote;
    typedef fftw_real   ValueType;

    typedef VigraFalseType isIntegral;
    typedef VigraFalseType isScalar;
    typedef NumericTraits<fftw_real>::isSigned isSigned;
    typedef VigraFalseType isOrdered;
    typedef VigraTrueType  isComplex;
};

template<>
struct NormTraits<fftw_complex>
{
    typedef fftw_complex Type;
    typedef fftw_real    SquaredNormType;
    typedef fftw_real    NormType;
};

template<>
struct NormTraits<FFTWComplex>
{
    typedef FFTWComplex Type;
    typedef fftw_real   SquaredNormType;
    typedef fftw_real   NormType;
};


template <>
struct PromoteTraits<fftw_complex, fftw_complex>
{
    typedef fftw_complex Promote;
};

template <>
struct PromoteTraits<fftw_complex, double>
{
    typedef fftw_complex Promote;
};

template <>
struct PromoteTraits<double, fftw_complex>
{
    typedef fftw_complex Promote;
};

template <>
struct PromoteTraits<FFTWComplex, FFTWComplex>
{
    typedef FFTWComplex Promote;
};

template <>
struct PromoteTraits<FFTWComplex, double>
{
    typedef FFTWComplex Promote;
};

template <>
struct PromoteTraits<double, FFTWComplex>
{
    typedef FFTWComplex Promote;
};


/********************************************************/
/*                                                      */
/*                    FFTWComplex Operations            */
/*                                                      */
/********************************************************/

/* documentation: see fftw3.hxx
*/
inline bool operator ==(FFTWComplex const &a, const FFTWComplex &b) {
    return c_re(a) == c_re(b) && c_im(a) == c_im(b);
}

inline bool operator !=(FFTWComplex const &a, const FFTWComplex &b) {
    return c_re(a) != c_re(b) || c_im(a) != c_im(b);
}

inline FFTWComplex &operator +=(FFTWComplex &a, const FFTWComplex &b) {
    c_re(a) += c_re(b);
    c_im(a) += c_im(b);
    return a;
}

inline FFTWComplex &operator -=(FFTWComplex &a, const FFTWComplex &b) {
    c_re(a) -= c_re(b);
    c_im(a) -= c_im(b);
    return a;
}

inline FFTWComplex &operator *=(FFTWComplex &a, const FFTWComplex &b) {
    FFTWComplex::value_type t = c_re(a)*c_re(b)-c_im(a)*c_im(b);
    c_im(a) = c_re(a)*c_im(b)+c_im(a)*c_re(b);
    c_re(a) = t;
    return a;
}

inline FFTWComplex &operator /=(FFTWComplex &a, const FFTWComplex &b) {
    FFTWComplex::value_type sm = b.squaredMagnitude();
    FFTWComplex::value_type t = (c_re(a)*c_re(b)+c_im(a)*c_im(b))/sm;
    c_im(a) = (c_re(b)*c_im(a)-c_re(a)*c_im(b))/sm;
    c_re(a) = t;
    return a;
}

inline FFTWComplex &operator *=(FFTWComplex &a, const double &b) {
    c_re(a) *= b;
    c_im(a) *= b;
    return a;
}

inline FFTWComplex &operator /=(FFTWComplex &a, const double &b) {
    c_re(a) /= b;
    c_im(a) /= b;
    return a;
}

inline FFTWComplex operator +(FFTWComplex a, const FFTWComplex &b) {
    a += b;
    return a;
}

inline FFTWComplex operator -(FFTWComplex a, const FFTWComplex &b) {
    a -= b;
    return a;
}

inline FFTWComplex operator *(FFTWComplex a, const FFTWComplex &b) {
    a *= b;
    return a;
}

inline FFTWComplex operator *(FFTWComplex a, const double &b) {
    a *= b;
    return a;
}

inline FFTWComplex operator *(const double &a, FFTWComplex b) {
    b *= a;
    return b;
}

inline FFTWComplex operator /(FFTWComplex a, const FFTWComplex &b) {
    a /= b;
    return a;
}

inline FFTWComplex operator /(FFTWComplex a, const double &b) {
    a /= b;
    return a;
}

using VIGRA_CSTD::abs;

inline FFTWComplex::value_type abs(const FFTWComplex &a)
{
    return a.magnitude();
}

inline FFTWComplex conj(const FFTWComplex &a)
{
    return FFTWComplex(a.re, -a.im);
}

inline FFTWComplex::NormType norm(const FFTWComplex &a)
{
    return a.magnitude();
}

inline FFTWComplex::SquaredNormType squaredNorm(const FFTWComplex &a)
{
    return a.squaredMagnitude();
}

/********************************************************/
/*                                                      */
/*                      FFTWRealImage                   */
/*                                                      */
/********************************************************/

/* documentation: see fftw3.hxx
*/
typedef BasicImage<fftw_real> FFTWRealImage;

/********************************************************/
/*                                                      */
/*                     FFTWComplexImage                 */
/*                                                      */
/********************************************************/

template<>
struct IteratorTraits<
        BasicImageIterator<FFTWComplex, FFTWComplex **> >
{
    typedef BasicImageIterator<FFTWComplex, FFTWComplex **>  Iterator;
    typedef Iterator                             iterator;
    typedef BasicImageIterator<FFTWComplex, FFTWComplex **>         mutable_iterator;
    typedef ConstBasicImageIterator<FFTWComplex, FFTWComplex **>    const_iterator;
    typedef iterator::iterator_category          iterator_category;
    typedef iterator::value_type                 value_type;
    typedef iterator::reference                  reference;
    typedef iterator::index_reference            index_reference;
    typedef iterator::pointer                    pointer;
    typedef iterator::difference_type            difference_type;
    typedef iterator::row_iterator               row_iterator;
    typedef iterator::column_iterator            column_iterator;
    typedef VectorAccessor<FFTWComplex>          default_accessor;
    typedef VectorAccessor<FFTWComplex>          DefaultAccessor;
    typedef VigraTrueType                        hasConstantStrides;
};

template<>
struct IteratorTraits<
        ConstBasicImageIterator<FFTWComplex, FFTWComplex **> >
{
    typedef ConstBasicImageIterator<FFTWComplex, FFTWComplex **>    Iterator;
    typedef Iterator                             iterator;
    typedef BasicImageIterator<FFTWComplex, FFTWComplex **>         mutable_iterator;
    typedef ConstBasicImageIterator<FFTWComplex, FFTWComplex **>    const_iterator;
    typedef iterator::iterator_category          iterator_category;
    typedef iterator::value_type                 value_type;
    typedef iterator::reference                  reference;
    typedef iterator::index_reference            index_reference;
    typedef iterator::pointer                    pointer;
    typedef iterator::difference_type            difference_type;
    typedef iterator::row_iterator               row_iterator;
    typedef iterator::column_iterator            column_iterator;
    typedef VectorAccessor<FFTWComplex>          default_accessor;
    typedef VectorAccessor<FFTWComplex>          DefaultAccessor;
    typedef VigraTrueType                        hasConstantStrides;
};

/* documentation: see fftw3.hxx
*/
typedef BasicImage<FFTWComplex> FFTWComplexImage;

/********************************************************/
/*                                                      */
/*                  FFTWComplex-Accessors               */
/*                                                      */
/********************************************************/

/* documentation: see fftw3.hxx
*/
class FFTWRealAccessor
{
  public:

        /// The accessor's value type.
    typedef fftw_real value_type;

        /// Read real part at iterator position.
    template <class ITERATOR>
    value_type operator()(ITERATOR const & i) const {
        return c_re(*i);
    }

        /// Read real part at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR const & i, DIFFERENCE d) const {
        return c_re(i[d]);
    }

        /// Write real part at iterator position.
    template <class ITERATOR>
    void set(value_type const & v, ITERATOR const & i) const {
        c_re(*i)= v;
    }

        /// Write real part at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    void set(value_type const & v, ITERATOR const & i, DIFFERENCE d) const {
        c_re(i[d])= v;
    }
};

/* documentation: see fftw3.hxx
*/
class FFTWImaginaryAccessor
{
  public:
        /// The accessor's value type.
    typedef fftw_real value_type;

        /// Read imaginary part at iterator position.
    template <class ITERATOR>
    value_type operator()(ITERATOR const & i) const {
        return c_im(*i);
    }

        /// Read imaginary part at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR const & i, DIFFERENCE d) const {
        return c_im(i[d]);
    }

        /// Write imaginary part at iterator position.
    template <class ITERATOR>
    void set(value_type const & v, ITERATOR const & i) const {
        c_im(*i)= v;
    }

        /// Write imaginary part at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    void set(value_type const & v, ITERATOR const & i, DIFFERENCE d) const {
        c_im(i[d])= v;
    }
};

/* documentation: see fftw3.hxx
*/
class FFTWWriteRealAccessor: public FFTWRealAccessor
{
  public:
        /// The accessor's value type.
    typedef fftw_real value_type;

        /** Write real number at iterator position. Set imaginary part
            to 0.
        */
    template <class ITERATOR>
    void set(value_type const & v, ITERATOR const & i) const {
        c_re(*i)= v;
        c_im(*i)= 0;
    }

        /** Write real number at offset from iterator position. Set imaginary part
            to 0.
        */
    template <class ITERATOR, class DIFFERENCE>
    void set(value_type const & v, ITERATOR const & i, DIFFERENCE d) const {
        c_re(i[d])= v;
        c_im(i[d])= 0;
    }
};

/* documentation: see fftw3.hxx
*/
class FFTWMagnitudeAccessor
{
  public:
        /// The accessor's value type.
    typedef fftw_real value_type;

        /// Read magnitude at iterator position.
    template <class ITERATOR>
    value_type operator()(ITERATOR const & i) const {
        return (*i).magnitude();
    }

        /// Read magnitude at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR const & i, DIFFERENCE d) const {
        return (i[d]).magnitude();
    }
};

/* documentation: see fftw3.hxx
*/
class FFTWPhaseAccessor
{
  public:
        /// The accessor's value type.
    typedef fftw_real value_type;

        /// Read phase at iterator position.
    template <class ITERATOR>
    value_type operator()(ITERATOR const & i) const {
        return (*i).phase();
    }

        /// Read phase at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR const & i, DIFFERENCE d) const {
        return (i[d]).phase();
    }
};

/********************************************************/
/*                                                      */
/*                    Fourier Transform                 */
/*                                                      */
/********************************************************/

/** \page FourierTransformFFTW2 Fast Fourier Transform
    
    This documentation describes the deprecated VIGRA interface to 
    FFTW 2. Use the \link FourierTransform interface to the newer
    version FFTW 3\endlink instead.
    
    VIGRA uses the <a href="http://www.fftw.org/">FFTW Fast Fourier
    Transform</a> package to perform Fourier transformations. VIGRA
    provides a wrapper for FFTW's complex number type (FFTWComplex),
    but FFTW's functions are used verbatim. If the image is stored as
    a FFTWComplexImage, a FFT is performed like this:

    \code
    vigra::FFTWComplexImage spatial(width,height), fourier(width,height);
    ... // fill image with data

    // create a plan for optimal performance
    fftwnd_plan forwardPlan=
        fftw2d_create_plan(height, width, FFTW_FORWARD, FFTW_ESTIMATE );

    // calculate FFT
    fftwnd_one(forwardPlan, spatial.begin(), fourier.begin());
    \endcode

    Note that in the creation of a plan, the height must be given
    first. Note also that <TT>spatial.begin()</TT> may only be passed
    to <TT>fftwnd_one</TT> if the transform shall be applied to the
    entire image. When you want to retrict operation to an ROI, you
    create a copy of the ROI in an image of appropriate size.

    More information on using FFTW can be found <a href="http://www.fftw.org/doc/">here</a>.

    FFTW produces fourier images that have the DC component (the
    origin of the Fourier space) in the upper left corner. Often, one
    wants the origin in the center of the image, so that frequencies
    always increase towards the border of the image. This can be
    achieved by calling \ref moveDCToCenter(). The inverse
    transformation is done by \ref moveDCToUpperLeft().

    <b>\#include</b> "<a href="fftw_8hxx-source.html">vigra/fftw.hxx</a>"<br>
    Namespace: vigra
*/

/********************************************************/
/*                                                      */
/*                     moveDCToCenter                   */
/*                                                      */
/********************************************************/

/* documentation: see fftw3.hxx
*/
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void moveDCToCenter(SrcImageIterator src_upperleft,
                               SrcImageIterator src_lowerright, SrcAccessor sa,
                               DestImageIterator dest_upperleft, DestAccessor da)
{
    int w= src_lowerright.x - src_upperleft.x;
    int h= src_lowerright.y - src_upperleft.y;
    int w1 = w/2;
    int h1 = h/2;
    int w2 = (w+1)/2;
    int h2 = (h+1)/2;

    // 2. Quadrant  zum 4.
    copyImage(srcIterRange(src_upperleft,
                           src_upperleft  + Diff2D(w2, h2), sa),
              destIter    (dest_upperleft + Diff2D(w1, h1), da));

    // 4. Quadrant zum 2.
    copyImage(srcIterRange(src_upperleft + Diff2D(w2, h2),
                           src_lowerright, sa),
              destIter    (dest_upperleft, da));

    // 1. Quadrant zum 3.
    copyImage(srcIterRange(src_upperleft  + Diff2D(w2, 0),
                           src_upperleft  + Diff2D(w,  h2), sa),
              destIter    (dest_upperleft + Diff2D(0,  h1), da));

    // 3. Quadrant zum 1.
    copyImage(srcIterRange(src_upperleft  + Diff2D(0,  h2),
                           src_upperleft  + Diff2D(w2, h), sa),
              destIter    (dest_upperleft + Diff2D(w1, 0), da));
}

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void moveDCToCenter(
    triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
    pair<DestImageIterator, DestAccessor> dest)
{
    moveDCToCenter(src.first, src.second, src.third,
                          dest.first, dest.second);
}

/********************************************************/
/*                                                      */
/*                   moveDCToUpperLeft                  */
/*                                                      */
/********************************************************/

/* documentation: see fftw3.hxx
*/
template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void moveDCToUpperLeft(SrcImageIterator src_upperleft,
                               SrcImageIterator src_lowerright, SrcAccessor sa,
                               DestImageIterator dest_upperleft, DestAccessor da)
{
    int w= src_lowerright.x - src_upperleft.x;
    int h= src_lowerright.y - src_upperleft.y;
    int w2 = w/2;
    int h2 = h/2;
    int w1 = (w+1)/2;
    int h1 = (h+1)/2;

    // 2. Quadrant  zum 4.
    copyImage(srcIterRange(src_upperleft,
                           src_upperleft  + Diff2D(w2, h2), sa),
              destIter    (dest_upperleft + Diff2D(w1, h1), da));

    // 4. Quadrant zum 2.
    copyImage(srcIterRange(src_upperleft + Diff2D(w2, h2),
                           src_lowerright, sa),
              destIter    (dest_upperleft, da));

    // 1. Quadrant zum 3.
    copyImage(srcIterRange(src_upperleft  + Diff2D(w2, 0),
                           src_upperleft  + Diff2D(w,  h2), sa),
              destIter    (dest_upperleft + Diff2D(0,  h1), da));

    // 3. Quadrant zum 1.
    copyImage(srcIterRange(src_upperleft  + Diff2D(0,  h2),
                           src_upperleft  + Diff2D(w2, h), sa),
              destIter    (dest_upperleft + Diff2D(w1, 0), da));
}

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
inline void moveDCToUpperLeft(
    triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
    pair<DestImageIterator, DestAccessor> dest)
{
    moveDCToUpperLeft(src.first, src.second, src.third,
                                          dest.first, dest.second);
}

/********************************************************/
/*                                                      */
/*                   applyFourierFilter                 */
/*                                                      */
/********************************************************/

/* documentation: see fftw3.hxx
*/

// applyFourierFilter versions without fftwnd_plans:
template <class SrcImageIterator, class SrcAccessor,
          class FilterImageIterator, class FilterAccessor,
          class DestImageIterator, class DestAccessor>
inline
void applyFourierFilter(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                        pair<FilterImageIterator, FilterAccessor> filter,
                        pair<DestImageIterator, DestAccessor> dest)
{
    applyFourierFilter(src.first, src.second, src.third,
                       filter.first, filter.second,
                       dest.first, dest.second);
}

template <class SrcImageIterator, class SrcAccessor,
          class FilterImageIterator, class FilterAccessor,
          class DestImageIterator, class DestAccessor>
void applyFourierFilter(SrcImageIterator srcUpperLeft,
                        SrcImageIterator srcLowerRight, SrcAccessor sa,
                        FilterImageIterator filterUpperLeft, FilterAccessor fa,
                        DestImageIterator destUpperLeft, DestAccessor da)
{
    // copy real input images into a complex one...
    int w= srcLowerRight.x - srcUpperLeft.x;
    int h= srcLowerRight.y - srcUpperLeft.y;

    FFTWComplexImage workImage(w, h);
    copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
              destImage(workImage, FFTWWriteRealAccessor()));

    // ...and call the impl
    FFTWComplexImage const & cworkImage = workImage;
    applyFourierFilterImpl(cworkImage.upperLeft(), cworkImage.lowerRight(), cworkImage.accessor(),
                           filterUpperLeft, fa,
                           destUpperLeft, da);
}

typedef FFTWComplexImage::const_traverser FFTWConstTraverser;

template <class FilterImageIterator, class FilterAccessor,
          class DestImageIterator, class DestAccessor>
inline
void applyFourierFilter(
    FFTWComplexImage::const_traverser srcUpperLeft,
    FFTWComplexImage::const_traverser srcLowerRight,
    FFTWComplexImage::ConstAccessor sa,
    FilterImageIterator filterUpperLeft, FilterAccessor fa,
    DestImageIterator destUpperLeft, DestAccessor da)
{
    int w= srcLowerRight.x - srcUpperLeft.x;
    int h= srcLowerRight.y - srcUpperLeft.y;

    // test for right memory layout (fftw expects a 2*width*height floats array)
    if (&(*(srcUpperLeft + Diff2D(w, 0))) == &(*(srcUpperLeft + Diff2D(0, 1))))
        applyFourierFilterImpl(srcUpperLeft, srcLowerRight, sa,
                               filterUpperLeft, fa,
                               destUpperLeft, da);
    else
    {
        FFTWComplexImage workImage(w, h);
        copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
                  destImage(workImage));

        FFTWComplexImage const & cworkImage = workImage;
        applyFourierFilterImpl(cworkImage.upperLeft(), cworkImage.lowerRight(), cworkImage.accessor(),
                               filterUpperLeft, fa,
                               destUpperLeft, da);
    }
}

template <class FilterImageIterator, class FilterAccessor,
          class DestImageIterator, class DestAccessor>
void applyFourierFilterImpl(
    FFTWComplexImage::const_traverser srcUpperLeft,
    FFTWComplexImage::const_traverser srcLowerRight,
    FFTWComplexImage::ConstAccessor sa,
    FilterImageIterator filterUpperLeft, FilterAccessor fa,
    DestImageIterator destUpperLeft, DestAccessor da)
{
    // create plans and call variant with plan parameters
    int w= srcLowerRight.x - srcUpperLeft.x;
    int h= srcLowerRight.y - srcUpperLeft.y;

    fftwnd_plan forwardPlan=
        fftw2d_create_plan(h, w, FFTW_FORWARD, FFTW_ESTIMATE );
    fftwnd_plan backwardPlan=
        fftw2d_create_plan(h, w, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);

    applyFourierFilterImpl(srcUpperLeft, srcLowerRight, sa,
                           filterUpperLeft, fa,
                           destUpperLeft, da,
                           forwardPlan, backwardPlan);

    fftwnd_destroy_plan(forwardPlan);
    fftwnd_destroy_plan(backwardPlan);
}

// applyFourierFilter versions with fftwnd_plans:
template <class SrcImageIterator, class SrcAccessor,
          class FilterImageIterator, class FilterAccessor,
          class DestImageIterator, class DestAccessor>
inline
void applyFourierFilter(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                        pair<FilterImageIterator, FilterAccessor> filter,
                        pair<DestImageIterator, DestAccessor> dest,
                        const fftwnd_plan &forwardPlan, const fftwnd_plan &backwardPlan)
{
    applyFourierFilter(src.first, src.second, src.third,
                       filter.first, filter.second,
                       dest.first, dest.second,
                       forwardPlan, backwardPlan);
}

template <class SrcImageIterator, class SrcAccessor,
          class FilterImageIterator, class FilterAccessor,
          class DestImageIterator, class DestAccessor>
void applyFourierFilter(SrcImageIterator srcUpperLeft,
                        SrcImageIterator srcLowerRight, SrcAccessor sa,
                        FilterImageIterator filterUpperLeft, FilterAccessor fa,
                        DestImageIterator destUpperLeft, DestAccessor da,
                        const fftwnd_plan &forwardPlan, const fftwnd_plan &backwardPlan)
{
    int w= srcLowerRight.x - srcUpperLeft.x;
    int h= srcLowerRight.y - srcUpperLeft.y;

    FFTWComplexImage workImage(w, h);
    copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
              destImage(workImage, FFTWWriteRealAccessor()));

    FFTWComplexImage const & cworkImage = workImage;
    applyFourierFilterImpl(cworkImage.upperLeft(), cworkImage.lowerRight(), cworkImage.accessor(),
                           filterUpperLeft, fa,
                           destUpperLeft, da,
                           forwardPlan, backwardPlan);
}

template <class FilterImageIterator, class FilterAccessor,
          class DestImageIterator, class DestAccessor>
inline
void applyFourierFilter(
    FFTWComplexImage::const_traverser srcUpperLeft,
    FFTWComplexImage::const_traverser srcLowerRight,
    FFTWComplexImage::ConstAccessor sa,
    FilterImageIterator filterUpperLeft, FilterAccessor fa,
    DestImageIterator destUpperLeft, DestAccessor da,
    const fftwnd_plan &forwardPlan, const fftwnd_plan &backwardPlan)
{
    applyFourierFilterImpl(srcUpperLeft, srcLowerRight, sa,
                           filterUpperLeft, fa,
                           destUpperLeft, da,
                           forwardPlan, backwardPlan);
}

template <class FilterImageIterator, class FilterAccessor,
          class DestImageIterator, class DestAccessor>
void applyFourierFilterImpl(
    FFTWComplexImage::const_traverser srcUpperLeft,
    FFTWComplexImage::const_traverser srcLowerRight,
    FFTWComplexImage::ConstAccessor sa,
    FilterImageIterator filterUpperLeft, FilterAccessor fa,
    DestImageIterator destUpperLeft, DestAccessor da,
    const fftwnd_plan &forwardPlan, const fftwnd_plan &backwardPlan)
{
    FFTWComplexImage complexResultImg(srcLowerRight - srcUpperLeft);

    // FFT from srcImage to complexResultImg
    fftwnd_one(forwardPlan, const_cast<FFTWComplex *>(&(*srcUpperLeft)),
               complexResultImg.begin());

    // convolve in freq. domain (in complexResultImg)
    combineTwoImages(srcImageRange(complexResultImg), srcIter(filterUpperLeft, fa),
                     destImage(complexResultImg), std::multiplies<FFTWComplex>());

    // FFT back into spatial domain (inplace in complexResultImg)
    fftwnd_one(backwardPlan, complexResultImg.begin(), 0);

    typedef typename
        NumericTraits<typename DestAccessor::value_type>::isScalar
        isScalarResult;

    // normalization (after FFTs), maybe stripping imaginary part
    applyFourierFilterImplNormalization(complexResultImg, destUpperLeft, da,
                                        isScalarResult());
}

template <class DestImageIterator, class DestAccessor>
void applyFourierFilterImplNormalization(FFTWComplexImage const &srcImage,
                                         DestImageIterator destUpperLeft,
                                         DestAccessor da,
                                         VigraFalseType)
{
    double normFactor= 1.0/(srcImage.width() * srcImage.height());

    for(int y=0; y<srcImage.height(); y++, destUpperLeft.y++)
    {
        DestImageIterator dIt= destUpperLeft;
        for(int x= 0; x< srcImage.width(); x++, dIt.x++)
        {
            da.setComponent(srcImage(x, y).re*normFactor, dIt, 0);
            da.setComponent(srcImage(x, y).im*normFactor, dIt, 1);
        }
    }
}

inline
void applyFourierFilterImplNormalization(FFTWComplexImage const & srcImage,
        FFTWComplexImage::traverser destUpperLeft,
        FFTWComplexImage::Accessor da,
        VigraFalseType)
{
    transformImage(srcImageRange(srcImage), destIter(destUpperLeft, da),
                   linearIntensityTransform<FFTWComplex>(1.0/(srcImage.width() * srcImage.height())));
}

template <class DestImageIterator, class DestAccessor>
void applyFourierFilterImplNormalization(FFTWComplexImage const & srcImage,
                                         DestImageIterator destUpperLeft,
                                         DestAccessor da,
                                         VigraTrueType)
{
    double normFactor= 1.0/(srcImage.width() * srcImage.height());

    for(int y=0; y<srcImage.height(); y++, destUpperLeft.y++)
    {
        DestImageIterator dIt= destUpperLeft;
        for(int x= 0; x< srcImage.width(); x++, dIt.x++)
            da.set(srcImage(x, y).re*normFactor, dIt);
    }
}

/**********************************************************/
/*                                                        */
/*                applyFourierFilterFamily                */
/*                                                        */
/**********************************************************/

/* documentation: see fftw3.hxx
*/

// applyFourierFilterFamily versions without fftwnd_plans:
template <class SrcImageIterator, class SrcAccessor,
          class FilterType, class DestImage>
inline
void applyFourierFilterFamily(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                              const ImageArray<FilterType> &filters,
                              ImageArray<DestImage> &results)
{
    applyFourierFilterFamily(src.first, src.second, src.third,
                             filters, results);
}

template <class SrcImageIterator, class SrcAccessor,
          class FilterType, class DestImage>
void applyFourierFilterFamily(SrcImageIterator srcUpperLeft,
                              SrcImageIterator srcLowerRight, SrcAccessor sa,
                              const ImageArray<FilterType> &filters,
                              ImageArray<DestImage> &results)
{
    int w= srcLowerRight.x - srcUpperLeft.x;
    int h= srcLowerRight.y - srcUpperLeft.y;

    FFTWComplexImage workImage(w, h);
    copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
              destImage(workImage, FFTWWriteRealAccessor()));

    FFTWComplexImage const & cworkImage = workImage;
    applyFourierFilterFamilyImpl(cworkImage.upperLeft(), cworkImage.lowerRight(), cworkImage.accessor(),
                                 filters, results);
}

template <class FilterType, class DestImage>
inline
void applyFourierFilterFamily(
    FFTWComplexImage::const_traverser srcUpperLeft,
    FFTWComplexImage::const_traverser srcLowerRight,
    FFTWComplexImage::ConstAccessor sa,
    const ImageArray<FilterType> &filters,
    ImageArray<DestImage> &results)
{
    applyFourierFilterFamilyImpl(srcUpperLeft, srcLowerRight, sa,
                                 filters, results);
}

template <class FilterType, class DestImage>
void applyFourierFilterFamilyImpl(
    FFTWComplexImage::const_traverser srcUpperLeft,
    FFTWComplexImage::const_traverser srcLowerRight,
    FFTWComplexImage::ConstAccessor sa,
    const ImageArray<FilterType> &filters,
    ImageArray<DestImage> &results)
{
    int w= srcLowerRight.x - srcUpperLeft.x;
    int h= srcLowerRight.y - srcUpperLeft.y;

    fftwnd_plan forwardPlan=
        fftw2d_create_plan(h, w, FFTW_FORWARD, FFTW_ESTIMATE );
    fftwnd_plan backwardPlan=
        fftw2d_create_plan(h, w, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);

    applyFourierFilterFamilyImpl(srcUpperLeft, srcLowerRight, sa,
                                 filters, results,
                                 forwardPlan, backwardPlan);

    fftwnd_destroy_plan(forwardPlan);
    fftwnd_destroy_plan(backwardPlan);
}

// applyFourierFilterFamily versions with fftwnd_plans:
template <class SrcImageIterator, class SrcAccessor,
          class FilterType, class DestImage>
inline
void applyFourierFilterFamily(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                              const ImageArray<FilterType> &filters,
                              ImageArray<DestImage> &results,
                              const fftwnd_plan &forwardPlan, const fftwnd_plan &backwardPlan)
{
    applyFourierFilterFamily(src.first, src.second, src.third,
                                 filters, results,
                                 forwardPlan, backwardPlan);
}

template <class SrcImageIterator, class SrcAccessor,
          class FilterType, class DestImage>
void applyFourierFilterFamily(SrcImageIterator srcUpperLeft,
                              SrcImageIterator srcLowerRight, SrcAccessor sa,
                              const ImageArray<FilterType> &filters,
                              ImageArray<DestImage> &results,
                              const fftwnd_plan &forwardPlan, const fftwnd_plan &backwardPlan)
{
    int w= srcLowerRight.x - srcUpperLeft.x;
    int h= srcLowerRight.y - srcUpperLeft.y;

    FFTWComplexImage workImage(w, h);
    copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
              destImage(workImage, FFTWWriteRealAccessor()));

    FFTWComplexImage const & cworkImage = workImage;
    applyFourierFilterFamilyImpl(cworkImage.upperLeft(), cworkImage.lowerRight(), cworkImage.accessor(),
                                 filters, results,
                                 forwardPlan, backwardPlan);
}

template <class FilterType, class DestImage>
inline
void applyFourierFilterFamily(
    FFTWComplexImage::const_traverser srcUpperLeft,
    FFTWComplexImage::const_traverser srcLowerRight,
    FFTWComplexImage::ConstAccessor sa,
    const ImageArray<FilterType> &filters,
    ImageArray<DestImage> &results,
    const fftwnd_plan &forwardPlan, const fftwnd_plan &backwardPlan)
{
    int w= srcLowerRight.x - srcUpperLeft.x;
    int h= srcLowerRight.y - srcUpperLeft.y;

    // test for right memory layout (fftw expects a 2*width*height floats array)
    if (&(*(srcUpperLeft + Diff2D(w, 0))) == &(*(srcUpperLeft + Diff2D(0, 1))))
        applyFourierFilterFamilyImpl(srcUpperLeft, srcLowerRight, sa,
                                     filters, results,
                                     forwardPlan, backwardPlan);
    else
    {
        FFTWComplexImage workImage(w, h);
        copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
                  destImage(workImage));

        FFTWComplexImage const & cworkImage = workImage;
        applyFourierFilterFamilyImpl(cworkImage.upperLeft(), cworkImage.lowerRight(), cworkImage.accessor(),
                                     filters, results,
                                     forwardPlan, backwardPlan);
    }
}

template <class FilterType, class DestImage>
void applyFourierFilterFamilyImpl(
    FFTWComplexImage::const_traverser srcUpperLeft,
    FFTWComplexImage::const_traverser srcLowerRight,
    FFTWComplexImage::ConstAccessor sa,
    const ImageArray<FilterType> &filters,
    ImageArray<DestImage> &results,
    const fftwnd_plan &forwardPlan, const fftwnd_plan &backwardPlan)
{
    // make sure the filter images have the right dimensions
    vigra_precondition((srcLowerRight - srcUpperLeft) == filters.imageSize(),
                       "applyFourierFilterFamily called with src image size != filters.imageSize()!");

    // make sure the result image array has the right dimensions
    results.resize(filters.size());
    results.resizeImages(filters.imageSize());

    // FFT from srcImage to freqImage
    int w= srcLowerRight.x - srcUpperLeft.x;
    int h= srcLowerRight.y - srcUpperLeft.y;

    FFTWComplexImage freqImage(w, h);
    FFTWComplexImage result(w, h);

    fftwnd_one(forwardPlan, const_cast<FFTWComplex *>(&(*srcUpperLeft)), freqImage.begin());

    typedef typename
        NumericTraits<typename DestImage::Accessor::value_type>::isScalar
        isScalarResult;

    // convolve with filters in freq. domain
    for (unsigned int i= 0;  i < filters.size(); i++)
    {
        combineTwoImages(srcImageRange(freqImage), srcImage(filters[i]),
                         destImage(result), std::multiplies<FFTWComplex>());

        // FFT back into spatial domain (inplace in destImage)
        fftwnd_one(backwardPlan, result.begin(), 0);

        // normalization (after FFTs), maybe stripping imaginary part
        applyFourierFilterImplNormalization(result,
                                            results[i].upperLeft(), results[i].accessor(),
                                            isScalarResult());
    }
}

} // namespace vigra

#endif // VIGRA_FFTW_HXX
