/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
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
 
#ifndef VIGRA_FFTW_HXX
#define VIGRA_FFTW_HXX

#include <cmath>

#include <fftw.h>
#include <vigra/stdimage.hxx>
#include <vigra/numerictraits.hxx>

namespace vigra {

/********************************************************/
/*                                                      */
/*                    FFTWComplex                       */
/*                                                      */
/********************************************************/

/** \brief Wrapper class for the FFTW type '<TT>fftw_complex</TT>'.

    This class provides constructors and other member functions 
    for the C struct '<TT>fftw_complex</TT>'. This struct is the basic 
    pixel type of the <a href="http://www.fftw.org/">FFTW Fast Fourier Transform</a>
    library. It inherits the data members '<TT>re</TT>' and '<TT>im</TT>'
    that denote the real and imaginary part of the number. In addition it
    defines transformations to polar coordinates, 
    as well as \ref FFTWComplexOperators "arithmetic operators" and 
    \ref FFTWComplexAccessors "accessors".
    
    <b>\#include</b> "<a href="fftw_8hxx-source.html">vigra/fftw.hxx</a>"<br>
    Namespace: vigra
*/
class FFTWComplex
: public fftw_complex
{
  public:
        /** The complex' component type, as defined in '<TT>fftw.h</TT>'
        */
    typedef fftw_real value_type;
    
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
    
        /** Assignment.
        */
    FFTWComplex& operator=(FFTWComplex const & o)
    {
        re = o.re;
        im = o.im;
        return *this;
    }
    
        /** Unary negation.
        */
    FFTWComplex operator-() const 
        { return FFTWComplex(-re, -im); }
    
        /** Squared magnitude x*conj(x)
        */
    value_type squaredMagnitude() const
        { return c_re(*this)*c_re(*this)+c_im(*this)*c_im(*this); }
        
        /** Magnitude (length of radius vector).
        */
    value_type magnitude() const
        { return VIGRA_CSTD::sqrt(squaredMagnitude()); }

        /** Phase angle.
        */
    value_type phase() const
        { return VIGRA_CSTD::atan2(c_im(*this),c_re(*this)); }

        /** Access components as if number were a vector.
        */
    value_type & operator[](int i)
        { return (&re)[i]; }

        /** Read components as if number were a vector.
        */
    value_type const & operator[](int i) const
        { return (&re)[i]; }

        /** Length of complex number (always 2).
        */
    int size() const
        { return 2; }
};

/********************************************************/
/*                                                      */
/*                    FFTWComplex Traits                */
/*                                                      */
/********************************************************/

/** \page FFTWComplexTraits Numeric and Promote Traits of FFTWComplex
    The numeric and promote traits for FFTWComplex follow 
    the general specifications for \ref NumericPromotionTraits. 
    They are implemented as follows:
    
    \code
    
    template <>
    struct NumericTraits<FFTWComplex >
    {
        typedef fftw_complex Type;
        typedef fftw_complex Promote;
        typedef fftw_complex RealPromote;
        typedef VigraFalseType isIntegral;
        typedef VigraFalseType isScalar;
        typedef VigraFalseType isOrdered;

        // etc.
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
    \endcode
    
    where the type '<TT>fftw_complex</TT>' is defined in '<TT>fftw.h</TT>'

    <b>\#include</b> "<a href="fftw_8hxx-source.html">vigra/fftw.hxx</a>"<br>
    Namespace: vigra
*/
template<>
struct NumericTraits<fftw_complex>
{
    typedef fftw_complex Type;
    typedef fftw_complex Promote;
    typedef fftw_complex RealPromote;
    typedef VigraFalseType isIntegral;
    typedef VigraFalseType isScalar;
    typedef VigraFalseType isOrdered;

    static FFTWComplex zero() { return FFTWComplex(0.0, 0.0); }
    static FFTWComplex one() { return FFTWComplex(1.0, 0.0); }
    static FFTWComplex nonZero() { return one(); }
};

template<>
struct NumericTraits<FFTWComplex>
: public NumericTraits<fftw_complex>
{
    typedef FFTWComplex Type;
    typedef FFTWComplex Promote;
    typedef FFTWComplex RealPromote;
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

/** \addtogroup FFTWComplexOperators Functions for FFTWComplex

    \brief <b>\#include</b> "<a href="fftw_8hxx-source.html">vigra/fftw.hxx</a>
    
    These functions fulfill the requirements of an Algebraic Field.
    Return types are determined according to \ref FFTWComplexTraits.

    Namespace: vigra
    <p>
    
 */
//@{
    /// equal
inline bool operator ==(FFTWComplex const &a, const FFTWComplex &b) {
    return c_re(a) == c_re(b) && c_im(a) == c_im(b);
}

    /// not equal
inline bool operator !=(FFTWComplex const &a, const FFTWComplex &b) {
    return c_re(a) != c_re(b) || c_im(a) != c_im(b);
}

    /// add-assignment
inline FFTWComplex &operator +=(FFTWComplex &a, const FFTWComplex &b) {
    c_re(a) += c_re(b);
    c_im(a) += c_im(b);
    return a;
}

    /// subtract-assignment
inline FFTWComplex &operator -=(FFTWComplex &a, const FFTWComplex &b) {
    c_re(a) -= c_re(b);
    c_im(a) -= c_im(b);
    return a;
}

    /// multiply-assignment
inline FFTWComplex &operator *=(FFTWComplex &a, const FFTWComplex &b) {
    FFTWComplex::value_type t = c_re(a)*c_re(b)-c_im(a)*c_im(b);
    c_im(a) = c_re(a)*c_im(b)+c_im(a)*c_re(b);
    c_re(a) = t;
    return a;
}

    /// divide-assignment
inline FFTWComplex &operator /=(FFTWComplex &a, const FFTWComplex &b) {
    FFTWComplex::value_type sm = b.squaredMagnitude();
    FFTWComplex::value_type t = (c_re(a)*c_re(b)+c_im(a)*c_im(b))/sm;
    c_im(a) = (c_re(b)*c_im(a)-c_re(a)*c_im(b))/sm;
    c_re(a) = t;
    return a;
}

    /// multiply-assignment with scalar double
inline FFTWComplex &operator *=(FFTWComplex &a, const double &b) {
    c_re(a) *= b;
    c_im(a) *= b;
    return a;
}

    /// divide-assignment with scalar double
inline FFTWComplex &operator /=(FFTWComplex &a, const double &b) {
    c_re(a) /= b;
    c_im(a) /= b;
    return a;
}

    /// addition
inline FFTWComplex operator +(FFTWComplex a, const FFTWComplex &b) {
    a += b;
    return a;
}

    /// subtraction
inline FFTWComplex operator -(FFTWComplex a, const FFTWComplex &b) {
    a -= b;
    return a;
}

    /// multiplication
inline FFTWComplex operator *(FFTWComplex a, const FFTWComplex &b) {
    a *= b;
    return a;
}

    /// right multiplication with scalar double
inline FFTWComplex operator *(FFTWComplex a, const double &b) {
    a *= b;
    return a;
}

    /// left multiplication with scalar double
inline FFTWComplex operator *(const double &a, FFTWComplex b) {
    b *= a;
    return b;
}

    /// division
inline FFTWComplex operator /(FFTWComplex a, const FFTWComplex &b) {
    a /= b;
    return a;
}

    /// right division with scalar double
inline FFTWComplex operator /(FFTWComplex a, const double &b) {
    a /= b;
    return a;
}

using VIGRA_CSTD::abs;

    /// absolute value (= magnitude)
inline FFTWComplex::value_type abs(const FFTWComplex &a)
{
    return a.magnitude();
}

    /// complex conjugate
inline FFTWComplex conj(const FFTWComplex &a)
{
    return FFTWComplex(a.re, -a.im);
}

//@}


/********************************************************/
/*                                                      */
/*                     FFTWComplexImage                 */
/*                                                      */
/********************************************************/

/** \addtogroup StandardImageTypes
*/
//@{

template<>
struct IteratorTraits<BasicImageIterator<FFTWComplex, FFTWComplex **> >
{
    typedef BasicImageIterator<FFTWComplex, FFTWComplex **>  Iterator;
    typedef Iterator                             iterator;
    typedef iterator::iterator_category          iterator_category;
    typedef iterator::value_type                 value_type;
    typedef iterator::reference                  reference;
    typedef iterator::index_reference            index_reference;
    typedef iterator::pointer                    pointer;
    typedef iterator::difference_type            difference_type;
    typedef iterator::row_iterator               row_iterator;
    typedef iterator::column_iterator            column_iterator;
    typedef StandardAccessor<FFTWComplex>        default_accessor; 
    typedef StandardAccessor<FFTWComplex>        DefaultAccessor; 
};  

template<>
struct IteratorTraits<ConstBasicImageIterator<FFTWComplex, FFTWComplex **> >
{
    typedef ConstBasicImageIterator<FFTWComplex, FFTWComplex **>    Iterator;
    typedef Iterator                             iterator;
    typedef iterator::iterator_category          iterator_category;
    typedef iterator::value_type                 value_type;
    typedef iterator::reference                  reference;
    typedef iterator::index_reference            index_reference;
    typedef iterator::pointer                    pointer;
    typedef iterator::difference_type            difference_type;
    typedef iterator::row_iterator               row_iterator;
    typedef iterator::column_iterator            column_iterator;
    typedef StandardConstAccessor<FFTWComplex>   default_accessor; 
    typedef StandardConstAccessor<FFTWComplex>   DefaultAccessor; 
};  

    /** Complex (FFTWComplex) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and 
        their const counterparts to access the data.
        
        <b>\#include</b> "<a href="fftw_8hxx-source.html">vigra/fftw.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<FFTWComplex> FFTWComplexImage;

//@}

/********************************************************/
/*                                                      */
/*                  FFTWComplex-Accessors               */
/*                                                      */
/********************************************************/

/** \defgroup FFTWComplexAccessors Accessors for FFTWComplex */
//@{
    /** Encapsulate access to the the real part of a complex number.

    <b>\#include</b> "<a href="fftw_8hxx-source.html">vigra/fftw.hxx</a>"<br>
    Namespace: vigra
    */
class FFTWRealAccessor
{
  public:
    
        /// The accessor's value type.
    typedef fftw_real value_type;

        /// Read real part at iterator position.
    template <class ITERATOR>
    value_type operator()(ITERATOR & i) const {
        return c_re(*i);
    }

        /// Read real part at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR & i, DIFFERENCE d) const {
        return c_re(i[d]);
    }

        /// Write real part at iterator position.
    template <class ITERATOR>
    void set(value_type const & v, ITERATOR & i) const {
        c_re(*i)= v;
    }

        /// Write real part at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    void set(value_type const & v, ITERATOR & i, DIFFERENCE d) const {
        c_re(i[d])= v;
    }
};

    /** Encapsulate access to the the imaginary part of a complex number.

    <b>\#include</b> "<a href="fftw_8hxx-source.html">vigra/fftw.hxx</a>"<br>
    Namespace: vigra
    */
class FFTWImaginaryAccessor
{
  public:
        /// The accessor's value type.
    typedef fftw_real value_type;

        /// Read imaginary part at iterator position.
    template <class ITERATOR>
    value_type operator()(ITERATOR & i) const {
        return c_im(*i);
    }

        /// Read imaginary part at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR & i, DIFFERENCE d) const {
        return c_im(i[d]);
    }

        /// Write imaginary part at iterator position.
    template <class ITERATOR>
    void set(value_type const & v, ITERATOR & i) const {
        c_im(*i)= v;
    }

        /// Write imaginary part at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    void set(value_type const & v, ITERATOR & i, DIFFERENCE d) const {
        c_im(i[d])= v;
    }
};

    /** Write a real number into a complex one. The imaginary part is set
        to 0.

    <b>\#include</b> "<a href="fftw_8hxx-source.html">vigra/fftw.hxx</a>"<br>
    Namespace: vigra
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
    void set(value_type const & v, ITERATOR & i) const {
        c_re(*i)= v;
        c_im(*i)= 0;
    }

        /** Write real number at offset from iterator position. Set imaginary part
            to 0.
        */
    template <class ITERATOR, class DIFFERENCE>
    void set(value_type const & v, ITERATOR & i, DIFFERENCE d) const {
        c_re(i[d])= v;
        c_im(i[d])= 0;
    }
};

    /** Calculate magnitude of complex number on the fly.

    <b>\#include</b> "<a href="fftw_8hxx-source.html">vigra/fftw.hxx</a>"<br>
    Namespace: vigra
    */
class FFTWMagnitudeAccessor
{
  public:
        /// The accessor's value type.
    typedef fftw_real value_type;

        /// Read magnitude at iterator position.
    template <class ITERATOR>
    value_type operator()(ITERATOR & i) const {
        return (*i).magnitude();
    }

        /// Read magnitude at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR & i, DIFFERENCE d) const {
        return (i[d]).magnitude();
    }
};

    /** Calculate phase of complex number on the fly.

    <b>\#include</b> "<a href="fftw_8hxx-source.html">vigra/fftw.hxx</a>"<br>
    Namespace: vigra
    */
class FFTWPhaseAccessor
{
  public:
        /// The accessor's value type.
    typedef fftw_real value_type;

        /// Read phase at iterator position.
    template <class ITERATOR>
    value_type operator()(ITERATOR & i) const {
        return (*i).phase();
    }

        /// Read phase at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR & i, DIFFERENCE d) const {
        return (i[d]).phase();
    }
};

//@}

/********************************************************/
/*                                                      */
/*                    Fourier Transform                 */
/*                                                      */
/********************************************************/

/** \page FourierTransform Fast Fourier Transform
    VIGRA uses the <a href="http://www.fftw.org/">FFTW Fast Fourier Transform</a>
    package to perform Fourier transformations. VIGRA provides a wrapper for 
    FFTW's complex number type (FFTWComplex), but FFTW's functions are used 
    verbatim. If the image is stored as a FFTWComplexImage, a FFT
    is performed like this:
    
    \code
    vigra::FFTWComplexImage spatial(width,height), fourier(width,height);
    ... // fill image with data
    
    // create a plan for optimal performance
    fftwnd_plan forwardPlan=
        fftw2d_create_plan(height, width, FFTW_FORWARD, FFTW_ESTIMATE );

    // calculate FFT
    fftwnd_one(forwardPlan, spatial.begin(), fourier.begin());    
    \endcode
    
    Note that in the creation of a plan, the height must be given first.
    Note also that <TT>spatial.begin()</TT> may only be passed to 
    <TT>fftwnd_one</TT> if the transform shall be applied to the entire
    image. When you want to retrict operation to an ROI, you create a 
    copy of the ROI in an image of appropriate size.
    
    More information on using FFTW can be found <a href="http://www.fftw.org/doc/">here</a>.

    FFTW produces fourier images that have the DC component (the origin
    of the Fourier space) in the upper left corner. Often, one wants the origin
    in the center of the image, so that frequencies always increase towards
    the border of the image. This can be achieved by calling 
    \ref moveDCToCenter(). The inverse transformation is done by 
    \ref moveDCToUpperLeft().

    <b>\#include</b> "<a href="fftw_8hxx-source.html">vigra/fftw.hxx</a>"<br>
    Namespace: vigra
*/
   
/** \addtogroup FourierHelpers Helper Functions for FFT
    Rearrange quadrants in a Fourier image.
*/
//@{

/********************************************************/
/*                                                      */
/*                     moveDCToCenter                   */
/*                                                      */
/********************************************************/

/** \brief Rearrange the quadrants of a Fourier image so that the origin is 
          in the image center.

    FFTW produces fourier images where the DC component (origin of fourier space)
    is located in the upper left corner of the image. The quadrants are placed
    like this (using a 4x4 image for example):
    
    \code
            DC 4 3 3
             4 4 3 3
             1 1 2 2
             1 1 2 2
    \endcode
    
    After applying the function, the quadrants are at their usual places:
             
    \code
            2 2  1 1
            2 2  1 1
            3 3 DC 4
            3 3  4 4
    \endcode
    
    This transformation can be reversed by \ref moveDCToUpperLeft(). Not that
    the transformation must not be executed in place - input and output images must
    be different.
    
    <b> Declarations:</b>
    
    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
        void moveDCToCenter(SrcImageIterator src_upperleft,
                               SrcImageIterator src_lowerright, SrcAccessor sa,
                               DestImageIterator dest_upperleft, DestAccessor da)
    }
    \endcode
    
    
    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        inline void moveDCToCenter(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest)
    }
    \endcode
    
    <b> Usage:</b>
    
        <b>\#include</b> "<a href="fftw_8hxx-source.html">vigra/fftw.hxx</a>"<br>
        Namespace: vigra
    
    \code
    vigra::FFTWComplexImage spatial(width,height), fourier(width,height);
    ... // fill image with data
    
    // create a plan for optimal performance
    fftwnd_plan forwardPlan=
        fftw2d_create_plan(height, width, FFTW_FORWARD, FFTW_ESTIMATE );

    // calculate FFT
    fftwnd_one(forwardPlan, spatial.begin(), fourier.begin());    
    
    vigra::FFTWComplexImage rearrangedFourier(width, height);
    moveDCToCenter(srcImageRange(fourier), destImage(rearrangedFourier));
    \endcode
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

/** \brief Rearrange the quadrants of a Fourier image so that the origin is 
          in the image's upper left.
          
     This function is the inversion of \ref moveDCToCenter(). See there
     for declarations and a usage example.

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

//@}

} // namespace vigra

#endif // VIGRA_FFTW_HXX
