/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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

#ifndef VIGRA_FFTW3_HXX
#define VIGRA_FFTW3_HXX

#include <cmath>
#include <functional>
#include "vigra/stdimage.hxx"
#include "vigra/copyimage.hxx"
#include "vigra/transformimage.hxx"
#include "vigra/combineimages.hxx"
#include "vigra/numerictraits.hxx"
#include "vigra/imagecontainer.hxx"
#include <fftw3.h>

namespace vigra {

typedef double fftw_real;

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

    FFTWComplex implements the concepts \ref AlgebraicField and
    \ref DivisionAlgebra.

    <b>\#include</b> "<a href="fftw3_8hxx-source.html">vigra/fftw3.hxx</a>"<br>
    Namespace: vigra
*/
class FFTWComplex
{
    fftw_complex data_;

  public:
        /** The complex' component type, as defined in '<TT>fftw3.h</TT>'
        */
    typedef fftw_real value_type;

        /** reference type (result of operator[])
        */
    typedef fftw_real & reference;

        /** const reference type (result of operator[] const)
        */
    typedef fftw_real const & const_reference;

        /** const reference type (result of operator[] const)
        */
    typedef fftw_real * iterator;

        /** const reference type (result of operator[] const)
        */
    typedef fftw_real const * const_iterator;

        /** Construct from real and imaginary part.
            Default: 0.
        */
    FFTWComplex(value_type const & re = 0.0, value_type const & im = 0.0)
    {
        data_[0] = re;
        data_[1] = im;
    }

        /** Copy constructor.
        */
    FFTWComplex(FFTWComplex const & o)
    {
        data_[0] = o.data_[0];
        data_[1] = o.data_[1];
    }

        /** Construct from plain <TT>fftw_complex</TT>.
        */
    FFTWComplex(fftw_complex const & o)
    {
        data_[0] = o[0];
        data_[1] = o[1];
    }

        /** Construct from TinyVector.
        */
    template <class T>
    FFTWComplex(TinyVector<T, 2> const & o)
    {
        data_[0] = o[0];
        data_[1] = o[1];
    }

        /** Assignment.
        */
    FFTWComplex& operator=(FFTWComplex const & o)
    {
        data_[0] = o.data_[0];
        data_[1] = o.data_[1];
        return *this;
    }

        /** Assignment.
        */
    FFTWComplex& operator=(fftw_complex const & o)
    {
        data_[0] = o[0];
        data_[1] = o[1];
        return *this;
    }

        /** Assignment.
        */
    FFTWComplex& operator=(fftw_real const & o)
    {
        data_[0] = o;
        data_[1] = 0.0;
        return *this;
    }

        /** Assignment.
        */
    template <class T>
    FFTWComplex& operator=(TinyVector<T, 2> const & o)
    {
        data_[0] = o[0];
        data_[1] = o[1];
        return *this;
    }

    reference re()
        { return data_[0]; }

    const_reference re() const
        { return data_[0]; }

    reference im()
        { return data_[1]; }

    const_reference im() const
        { return data_[1]; }

        /** Unary negation.
        */
    FFTWComplex operator-() const
        { return FFTWComplex(-data_[0], -data_[1]); }

        /** Squared magnitude x*conj(x)
        */
    value_type squaredMagnitude() const
        { return data_[0]*data_[0]+data_[1]*data_[1]; }

        /** Magnitude (length of radius vector).
        */
    value_type magnitude() const
        { return VIGRA_CSTD::sqrt(squaredMagnitude()); }

        /** Phase angle.
        */
    value_type phase() const
        { return VIGRA_CSTD::atan2(data_[1], data_[0]); }

        /** Access components as if number were a vector.
        */
    reference operator[](int i)
        { return data_[i]; }

        /** Read components as if number were a vector.
        */
    const_reference operator[](int i) const
        { return data_[i]; }

        /** Length of complex number (always 2).
        */
    int size() const
        { return 2; }

    iterator begin()
        { return data_; }

    iterator end()
        { return data_ + 2; }

    const_iterator begin() const
        { return data_; }

    const_iterator end() const
        { return data_ + 2; }
};

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
{
    typedef FFTWComplex Type;
    typedef FFTWComplex Promote;
    typedef FFTWComplex RealPromote;
    typedef FFTWComplex ComplexPromote;
    typedef fftw_real   ValueType;
    typedef VigraFalseType isIntegral;
    typedef VigraFalseType isScalar;
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

    \brief <b>\#include</b> "<a href="fftw3_8hxx-source.html">vigra/fftw3.hxx</a>

    These functions fulfill the requirements of an Algebraic Field.
    Return types are determined according to \ref FFTWComplexTraits.

    Namespace: vigra
    <p>

 */
//@{
    /// equal
inline bool operator ==(FFTWComplex const &a, const FFTWComplex &b) {
    return a.re() == b.re() && a.im() == b.im();
}

    /// not equal
inline bool operator !=(FFTWComplex const &a, const FFTWComplex &b) {
    return a.re() != b.re() || a.im() != b.im();
}

    /// add-assignment
inline FFTWComplex &operator +=(FFTWComplex &a, const FFTWComplex &b) {
    a.re() += b.re();
    a.im() += b.im();
    return a;
}

    /// subtract-assignment
inline FFTWComplex &operator -=(FFTWComplex &a, const FFTWComplex &b) {
    a.re() -= b.re();
    a.im() -= b.im();
    return a;
}

    /// multiply-assignment
inline FFTWComplex &operator *=(FFTWComplex &a, const FFTWComplex &b) {
    FFTWComplex::value_type t = a.re()*b.re()-a.im()*b.im();
    a.im() = a.re()*b.im()+a.im()*b.re();
    a.re() = t;
    return a;
}

    /// divide-assignment
inline FFTWComplex &operator /=(FFTWComplex &a, const FFTWComplex &b) {
    FFTWComplex::value_type sm = b.squaredMagnitude();
    FFTWComplex::value_type t = (a.re()*b.re()+a.im()*b.im())/sm;
    a.im() = (b.re()*a.im()-a.re()*b.im())/sm;
    a.re() = t;
    return a;
}

    /// multiply-assignment with scalar double
inline FFTWComplex &operator *=(FFTWComplex &a, const double &b) {
    a.re() *= b;
    a.im() *= b;
    return a;
}

    /// divide-assignment with scalar double
inline FFTWComplex &operator /=(FFTWComplex &a, const double &b) {
    a.re() /= b;
    a.im() /= b;
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
    return FFTWComplex(a.re(), -a.im());
}

//@}


/** \addtogroup StandardImageTypes
*/
//@{

/********************************************************/
/*                                                      */
/*                      FFTWRealImage                   */
/*                                                      */
/********************************************************/

    /** Float (<tt>fftw_real</tt>) image.
        The type <tt>fftw_real</tt> (either <tt>float</tt> or <tt>double</tt>)
        is defined during compilation of fftw and imported into VIGRA
        from <tt>fftw3.h</tt>. FFTWRealImage uses \ref vigra::BasicImageIterator
        and \ref vigra::StandardAccessor and
        their const counterparts to access the data.

        <b>\#include</b> "<a href="fftw3_8hxx-source.html">vigra/fftw3.hxx</a>"<br>
        Namespace: vigra
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
};

template<>
struct IteratorTraits<
        ConstBasicImageIterator<FFTWComplex, FFTWComplex **> >
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
    typedef VectorAccessor<FFTWComplex>          default_accessor;
    typedef VectorAccessor<FFTWComplex>          DefaultAccessor;
};

    /** Complex (FFTWComplex) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and
        their const counterparts to access the data.

        <b>\#include</b> "<a href="fftw3_8hxx-source.html">vigra/fftw3.hxx</a>"<br>
        Namespace: vigra
    */
typedef BasicImage<FFTWComplex> FFTWComplexImage;

//@}

/********************************************************/
/*                                                      */
/*                  FFTWComplex-Accessors               */
/*                                                      */
/********************************************************/

/** \addtogroup DataAccessors
*/
//@{
/** \defgroup FFTWComplexAccessors Accessors for FFTWComplex

    Encapsulate access to pixels of type FFTWComplex
*/
//@{
    /** Encapsulate access to the the real part of a complex number.

    <b>\#include</b> "<a href="fftw3_8hxx-source.html">vigra/fftw3.hxx</a>"<br>
    Namespace: vigra
    */
class FFTWRealAccessor
{
  public:

        /// The accessor's value type.
    typedef fftw_real value_type;

        /// Read real part at iterator position.
    template <class ITERATOR>
    value_type operator()(ITERATOR const & i) const {
        return (*i).re();
    }

        /// Read real part at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR const & i, DIFFERENCE d) const {
        return i[d].re();
    }

        /// Write real part at iterator position.
    template <class ITERATOR>
    void set(value_type const & v, ITERATOR const & i) const {
        (*i).re()= v;
    }

        /// Write real part at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    void set(value_type const & v, ITERATOR const & i, DIFFERENCE d) const {
        i[d].re()= v;
    }
};

    /** Encapsulate access to the the imaginary part of a complex number.

    <b>\#include</b> "<a href="fftw3_8hxx-source.html">vigra/fftw3.hxx</a>"<br>
    Namespace: vigra
    */
class FFTWImaginaryAccessor
{
  public:
        /// The accessor's value type.
    typedef fftw_real value_type;

        /// Read imaginary part at iterator position.
    template <class ITERATOR>
    value_type operator()(ITERATOR const & i) const {
        return (*i).im();
    }

        /// Read imaginary part at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR const & i, DIFFERENCE d) const {
        return i[d].im();
    }

        /// Write imaginary part at iterator position.
    template <class ITERATOR>
    void set(value_type const & v, ITERATOR const & i) const {
        (*i).im()= v;
    }

        /// Write imaginary part at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    void set(value_type const & v, ITERATOR const & i, DIFFERENCE d) const {
        i[d].im()= v;
    }
};

    /** Write a real number into a complex one. The imaginary part is set
        to 0.

    <b>\#include</b> "<a href="fftw3_8hxx-source.html">vigra/fftw3.hxx</a>"<br>
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
    void set(value_type const & v, ITERATOR const & i) const {
        (*i).re()= v;
        (*i).im()= 0;
    }

        /** Write real number at offset from iterator position. Set imaginary part
            to 0.
        */
    template <class ITERATOR, class DIFFERENCE>
    void set(value_type const & v, ITERATOR const & i, DIFFERENCE d) const {
        i[d].re()= v;
        i[d].im()= 0;
    }
};

    /** Calculate magnitude of complex number on the fly.

    <b>\#include</b> "<a href="fftw3_8hxx-source.html">vigra/fftw3.hxx</a>"<br>
    Namespace: vigra
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

    /** Calculate phase of complex number on the fly.

    <b>\#include</b> "<a href="fftw3_8hxx-source.html">vigra/fftw3.hxx</a>"<br>
    Namespace: vigra
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

//@}
//@}

/********************************************************/
/*                                                      */
/*                    Fourier Transform                 */
/*                                                      */
/********************************************************/

/** \page FourierTransform Fast Fourier Transform
    VIGRA uses the <a href="http://www.fftw.org/">FFTW Fast Fourier
    Transform</a> package to perform Fourier transformations. VIGRA
    provides a wrapper for FFTW's complex number type (FFTWComplex),
    but FFTW's functions are used verbatim. If the image is stored as
    a FFTWComplexImage, the simplest call to an FFT function is like this:

    \code
    vigra::FFTWComplexImage spatial(width,height), fourier(width,height);
    ... // fill image with data

    // create a plan with estimated performance optimization
    fftw_plan forwardPlan = fftw_plan_dft_2d(height, width, 
                                (fftw_complex *)spatial.begin(), (fftw_complex *)fourier.begin(), 
                                FFTW_FORWARD, FFTW_ESTIMATE );
    // calculate FFT (this can be repeated as often as needed, 
    //                with fresh data written into the source array)
    fftw_execute(forwardPlan);
    
    // release the plan memory
    fftw_destroy_plan(forwardPlan);
    
    // likewise for the inverse transform
    fftw_plan backwardPlan = fftw_plan_dft_2d(height, width, 
                                 (fftw_complex *)fourier.begin(), (fftw_complex *)spatial.begin(), 
                                 FFTW_BACKWARD, FFTW_ESTIMATE);        
    fftw_execute(backwardPlan);
    fftw_destroy_plan(backwardPlan);
    
    // do not forget to normalize the result according to the image size
    transformImage(srcImageRange(spatial), destImage(spatial),
                   std::bind1st(std::multiplies<FFTWComplex>(), 1.0 / width / height));
    \endcode

    Note that in the creation of a plan, the height must be given
    first. Note also that <TT>spatial.begin()</TT> may only be passed
    to <TT>fftw_plan_dft_2d</TT> if the transform shall be applied to the
    entire image. When you want to restrict operation to an ROI, you
    can create a copy of the ROI in an image of appropriate size, or
    you may use the Guru interface to FFTW. 

    More information on using FFTW can be found <a href="http://www.fftw.org/doc/">here</a>.

    FFTW produces fourier images that have the DC component (the
    origin of the Fourier space) in the upper left corner. Often, one
    wants the origin in the center of the image, so that frequencies
    always increase towards the border of the image. This can be
    achieved by calling \ref moveDCToCenter(). The inverse
    transformation is done by \ref moveDCToUpperLeft().

    <b>\#include</b> "<a href="fftw3_8hxx-source.html">vigra/fftw3.hxx</a>"<br>
    Namespace: vigra
*/

/** \addtogroup FourierHelpers Helper Functions for FFT
    Rearrange quadrants in a Fourier image or apply filters on images
    using the FFTW.
*/
//@{

/********************************************************/
/*                                                      */
/*                     moveDCToCenter                   */
/*                                                      */
/********************************************************/

/** \brief Rearrange the quadrants of a Fourier image so that the origin is
          in the image center.

    FFTW produces fourier images where the DC component (origin of
    fourier space) is located in the upper left corner of the
    image. The quadrants are placed like this (using a 4x4 image for
    example):

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

    This transformation can be reversed by \ref moveDCToUpperLeft().
    Note that the transformation must not be executed in place - input
    and output images must be different.

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

        <b>\#include</b> "<a href="fftw3_8hxx-source.html">vigra/fftw3.hxx</a>"<br>
        Namespace: vigra

    \code
    vigra::FFTWComplexImage spatial(width,height), fourier(width,height);
    ... // fill image with data

    // create a plan with estimated performance optimization
    fftw_plan forwardPlan = fftw_plan_dft_2d(height, width, 
                                (fftw_complex *)spatial.begin(), (fftw_complex *)fourier.begin(), 
                                FFTW_FORWARD, FFTW_ESTIMATE );
    // calculate FFT
    fftw_execute(forwardPlan);

    vigra::FFTWComplexImage rearrangedFourier(width, height);
    moveDCToCenter(srcImageRange(fourier), destImage(rearrangedFourier));

    // delete the plan
    fftw_destroy_plan(forwardPlan);
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

/********************************************************/
/*                                                      */
/*                   applyFourierFilter                 */
/*                                                      */
/********************************************************/

/** \brief Apply a filter (defined in the frequency domain) to an image.

    After transferring the image into the frequency domain, it is
    multiplied pixel-wise with the filter and transformed back. The
    result is always put into the given FFTWComplexImage
    <TT>destImg</TT> which must have the right size. Finally, the
    result will be normalized to compensate for the two FFTs.

    The input and filter images can be scalar images (more precisely,
    <TT>SrcAccessor::value_type</TT> must be scalar) or
    FFTWComplexImages (in this case, one conversion is saved).

    The DC entry of the filter must be in the upper left, which is the
    position where FFTW expects it (see \ref moveDCToUpperLeft()).

    You can optionally pass two fftwnd_plans as last parameters, the
    forward plan and the in-place backward plan. Both must have been
    created for the right image size (see sample code).

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
            class FilterImageIterator, class FilterAccessor>
        void applyFourierFilter(SrcImageIterator srcUpperLeft,
                                SrcImageIterator srcLowerRight, SrcAccessor sa,
                                FilterImageIterator filterUpperLeft, FilterAccessor fa,
                                FFTWComplexImage & destImg)

        template <class SrcImageIterator, class SrcAccessor,
            class FilterImageIterator, class FilterAccessor>
        void applyFourierFilter(SrcImageIterator srcUpperLeft,
                                SrcImageIterator srcLowerRight, SrcAccessor sa,
                                FilterImageIterator filterUpperLeft, FilterAccessor fa,
                                FFTWComplexImage & destImg,
                                const fftwnd_plan & forwardPlan, const fftwnd_plan & backwardPlan)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
            class FilterImageIterator, class FilterAccessor>
        inline
        void applyFourierFilter(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                                pair<FilterImageIterator, FilterAccessor> filter,
                                FFTWComplexImage &destImg)

        template <class SrcImageIterator, class SrcAccessor,
            class FilterImageIterator, class FilterAccessor>
        inline
        void applyFourierFilter(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                                pair<FilterImageIterator, FilterAccessor> filter,
                                FFTWComplexImage &destImg,
                                const fftwnd_plan &forwardPlan, const fftwnd_plan &backwardPlan)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="fftw3_8hxx-source.html">vigra/fftw3.hxx</a>"<br>
    Namespace: vigra

    \code
    // create a Gaussian filter in Fourier space
    vigra::FImage gaussFilter(w, h), filter(w, h);
    for(int y=0; y<h; ++y)
        for(int x=0; x<w; ++x)
        {
            xx = float(x - w / 2) / w;
            yy = float(y - h / 2) / h;

            gaussFilter(x,y) = std::exp(-(xx*xx + yy*yy) / 2.0 * scale);
        }

    // applyFourierFilter() expects the filter's DC in the upper left
    moveDCToUpperLeft(srcImageRange(gaussFilter), destImage(filter));

    vigra::FFTWComplexImage result(w, h);

    vigra::applyFourierFilter(srcImageRange(image), srcImage(filter), result);
    \endcode

    For inspection of the result, \ref FFTWMagnitudeAccessor might be
    useful. If you want to apply the same filter repeatedly, it may be more
    efficient to use the FFTW functions directly with FFTW plans optimized
    for good performance.
*/
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
    int w = srcLowerRight.x - srcUpperLeft.x;
    int h = srcLowerRight.y - srcUpperLeft.y;

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
    int w = srcLowerRight.x - srcUpperLeft.x;
    int h = srcLowerRight.y - srcUpperLeft.y;

    FFTWComplexImage complexResultImg(srcLowerRight - srcUpperLeft);

    // FFT from srcImage to complexResultImg
    fftw_plan forwardPlan=
        fftw_plan_dft_2d(h, w, (fftw_complex *)&(*srcUpperLeft),
                               (fftw_complex *)complexResultImg.begin(),
                               FFTW_FORWARD, FFTW_ESTIMATE );
    fftw_execute(forwardPlan);
    fftw_destroy_plan(forwardPlan);

    // convolve in freq. domain (in complexResultImg)
    combineTwoImages(srcImageRange(complexResultImg), srcIter(filterUpperLeft, fa),
                     destImage(complexResultImg), std::multiplies<FFTWComplex>());

    // FFT back into spatial domain (inplace in complexResultImg)
    fftw_plan backwardPlan=
        fftw_plan_dft_2d(h, w, (fftw_complex *)complexResultImg.begin(),
                               (fftw_complex *)complexResultImg.begin(),
                               FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(backwardPlan);
    fftw_destroy_plan(backwardPlan);

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
            da.setComponent(srcImage(x, y).re()*normFactor, dIt, 0);
            da.setComponent(srcImage(x, y).im()*normFactor, dIt, 1);
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
            da.set(srcImage(x, y).re()*normFactor, dIt);
    }
}

/**********************************************************/
/*                                                        */
/*                applyFourierFilterFamily                */
/*                                                        */
/**********************************************************/

/** \brief Apply an array of filters (defined in the frequency domain) to an image.

    This provides the same functionality as \ref applyFourierFilter(),
    but applying several filters at once allows to avoid
    repeated Fourier transforms of the source image.

    Filters and result images must be stored in \ref vigra::ImageArray data
    structures. In contrast to \ref applyFourierFilter(), this function adjusts
    the size of the result images and the the length of the array.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor, class FilterType>
        void applyFourierFilterFamily(SrcImageIterator srcUpperLeft,
                                      SrcImageIterator srcLowerRight, SrcAccessor sa,
                                      const ImageArray<FilterType> &filters,
                                      ImageArray<FFTWComplexImage> &results)
    }
    \endcode

    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor, class FilterType>
        inline
        void applyFourierFilterFamily(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                                      const ImageArray<FilterType> &filters,
                                      ImageArray<FFTWComplexImage> &results)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> "<a href="fftw3_8hxx-source.html">vigra/fftw3.hxx</a>"<br>
    Namespace: vigra

    \code
    // assuming the presence of a real-valued image named "image" to
    // be filtered in this example

    vigra::ImageArray<vigra::FImage> filters(16, image.size());

    for (int i=0; i<filters.size(); i++)
         // create some meaningful filters here
         createMyFilterOfScale(i, destImage(filters[i]));

    vigra::ImageArray<vigra::FFTWComplexImage> results();

    vigra::applyFourierFilterFamily(srcImageRange(image), filters, results);
    \endcode
*/
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
    int w= srcLowerRight.x - srcUpperLeft.x;

    // test for right memory layout (fftw expects a 2*width*height floats array)
    if (&(*(srcUpperLeft + Diff2D(w, 0))) == &(*(srcUpperLeft + Diff2D(0, 1))))
        applyFourierFilterFamilyImpl(srcUpperLeft, srcLowerRight, sa,
                                     filters, results);
    else
    {
        FFTWComplexImage workImage(w, h);
        copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
                  destImage(workImage));

        FFTWComplexImage const & cworkImage = workImage;
        applyFourierFilterFamilyImpl(cworkImage.upperLeft(), cworkImage.lowerRight(), cworkImage.accessor(),
                                     filters, results);
    }
}

template <class FilterType, class DestImage>
void applyFourierFilterFamilyImpl(
    FFTWComplexImage::const_traverser srcUpperLeft,
    FFTWComplexImage::const_traverser srcLowerRight,
    FFTWComplexImage::ConstAccessor sa,
    const ImageArray<FilterType> &filters,
    ImageArray<DestImage> &results)
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

    fftw_plan forwardPlan=
        fftw_plan_dft_2d(h, w, (fftw_complex *)&(*srcUpperLeft),
                               (fftw_complex *)freqImage.begin(),
                               FFTW_FORWARD, FFTW_ESTIMATE );
    fftw_execute(forwardPlan);
    fftw_destroy_plan(forwardPlan);

    fftw_plan backwardPlan=
        fftw_plan_dft_2d(h, w, (fftw_complex *)result.begin(),
                               (fftw_complex *)result.begin(),
                               FFTW_BACKWARD, FFTW_ESTIMATE );
    typedef typename
        NumericTraits<typename DestImage::Accessor::value_type>::isScalar
        isScalarResult;

    // convolve with filters in freq. domain
    for (unsigned int i= 0;  i < filters.size(); i++)
    {
        combineTwoImages(srcImageRange(freqImage), srcImage(filters[i]),
                         destImage(result), std::multiplies<FFTWComplex>());

        // FFT back into spatial domain (inplace in destImage)
        fftw_execute(backwardPlan);

        // normalization (after FFTs), maybe stripping imaginary part
        applyFourierFilterImplNormalization(result,
                                            results[i].upperLeft(), results[i].accessor(),
                                            isScalarResult());
    }
    fftw_destroy_plan(backwardPlan);
}

/********************************************************/
/*                                                      */
/*                fourierTransformReal                  */
/*                                                      */
/********************************************************/

/** \brief Real Fourier transforms for even and odd boundary conditions
           (aka. cosine and sine transforms).

    If the image is real and has even symmetry, its Fourier transform
    is also real and has even symmetry. The Fourier transform of a real image with odd
    symmetry is imaginary and has odd symmetry. In either case, only about a quarter
    of the pixels need to be stored because the rest can be calculated from the symmetry
    properties. This is especially useful, if the original image is implicitly assumed
    to have reflective or anti-reflective boundary conditions. Then the "negative"
    pixel locations are defined as

    \code
    even (reflective boundary conditions):      f[-x] = f[x]     (x = 1,...,N-1)
    odd (anti-reflective boundary conditions):  f[-1] = 0
                                                f[-x] = -f[x-2]  (x = 2,...,N-1)
    \endcode
    
    end similar at the other boundary (see the FFTW documentation for details). 
    This has the advantage that more efficient Fourier transforms that use only
    real numbers can be implemented. These are also known as cosine and sine transforms
    respectively. 
    
    If you use the odd transform it is important to note that in the Fourier domain,
    the DC component is always zero and is therefore dropped from the data structure.
    This means that index 0 in an odd symmetric Fourier domain image refers to
    the <i>first</i> harmonic. This is especially important if an image is first 
    cosine transformed (even symmetry), then in the Fourier domain multiplied 
    with an odd symmetric filter (e.g. a first derivative) and finally transformed
    back to the spatial domain with a sine transform (odd symmetric). For this to work
    properly the image must be shifted left or up by one pixel (depending on whether
    the x- or y-axis is odd symmetric) before the inverse transform can be applied.
    (see example below).
    
    The real Fourier transform functions are named <tt>fourierTransformReal??</tt>
    where the questions marks stand for either <tt>E</tt> or <tt>O</tt> indicating
    whether the x- and y-axis is to be transformed using even or odd symmetry.
    The same functions can be used for both the forward and inverse transforms,
    only the normalization changes. For signal processing, the following 
    normalization factors are most appropriate:
    
    \code
                          forward             inverse
    ------------------------------------------------------------
    X even, Y even           1.0         4.0 * (w-1) * (h-1)
    X even, Y odd           -1.0        -4.0 * (w-1) * (h+1)
    X odd,  Y even          -1.0        -4.0 * (w+1) * (h-1)
    X odd,  Y odd            1.0         4.0 * (w+1) * (h+1)
    \endcode

    where <tt>w</tt> and <tt>h</tt> denote the image width and height.

    <b> Declarations:</b>

    pass arguments explicitly:
    \code
    namespace vigra {
        template <class SrcTraverser, class SrcAccessor,
                  class DestTraverser, class DestAccessor>
        void
        fourierTransformRealEE(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                               DestTraverser dul, DestAccessor dest, fftw_real norm);
                               
        fourierTransformRealEO, fourierTransformRealOE, fourierTransformRealOO likewise
    }
    \endcode


    use argument objects in conjuction with \ref ArgumentObjectFactories:
    \code
    namespace vigra {
        template <class SrcTraverser, class SrcAccessor,
                  class DestTraverser, class DestAccessor>
        void
        fourierTransformRealEE(triple<SrcTraverser, SrcTraverser, SrcAccessor> src,
                               pair<DestTraverser, DestAccessor> dest, fftw_real norm);
                               
        fourierTransformRealEO, fourierTransformRealOE, fourierTransformRealOO likewise
    }
    \endcode

    <b> Usage:</b>

        <b>\#include</b> "<a href="fftw3_8hxx-source.html">vigra/fftw3.hxx</a>"<br>
        Namespace: vigra

    \code
    vigra::FImage spatial(width,height), fourier(width,height);
    ... // fill image with data

    // forward cosine transform == reflective boundary conditions
    fourierTransformRealEE(srcImageRange(spatial), destImage(fourier), (fftw_real)1.0);

    // multiply with a first derivative of Gaussian in x-direction
    for(int y = 0; y < height; ++y)
    {
        for(int x = 1; x < width; ++x)
        {
            double dx = x * M_PI / (width - 1);
            double dy = y * M_PI / (height - 1);
            fourier(x-1, y) = fourier(x, y) * dx * exp(-(dx*dx + dy*dy) * scale*scale / 2.0);
        }
        fourier(width-1, y) = 0.0;
    }
    
    // inverse transform -- odd symmetry in x-direction, even in y,
    //                      due to symmetry of the filter
    fourierTransformRealOE(srcImageRange(fourier), destImage(spatial), 
                           (fftw_real)-4.0 * (width+1) * (height-1));
    \endcode
*/
template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealEE(triple<SrcTraverser, SrcTraverser, SrcAccessor> src,
                               pair<DestTraverser, DestAccessor> dest, fftw_real norm)
{
    fourierTransformRealEE(src.first, src.second, src.third,
                                   dest.first, dest.second, norm);
}

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealEE(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                               DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                      norm, FFTW_REDFT00, FFTW_REDFT00);
}

template <class DestTraverser, class DestAccessor>
inline void
fourierTransformRealEE(
         FFTWRealImage::const_traverser sul,
         FFTWRealImage::const_traverser slr,
         FFTWRealImage::Accessor src,
         DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    int w = slr.x - sul.x;

    // test for right memory layout (fftw expects a width*height fftw_real array)
    if (&(*(sul + Diff2D(w, 0))) == &(*(sul + Diff2D(0, 1))))
        fourierTransformRealImpl(sul, slr, dul, dest,
                                 norm, FFTW_REDFT00, FFTW_REDFT00);
    else
        fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                 norm, FFTW_REDFT00, FFTW_REDFT00);
}

/********************************************************************/

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealOE(triple<SrcTraverser, SrcTraverser, SrcAccessor> src,
                               pair<DestTraverser, DestAccessor> dest, fftw_real norm)
{
    fourierTransformRealOE(src.first, src.second, src.third,
                                   dest.first, dest.second, norm);
}

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealOE(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                               DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                      norm, FFTW_RODFT00, FFTW_REDFT00);
}

template <class DestTraverser, class DestAccessor>
inline void
fourierTransformRealOE(
         FFTWRealImage::const_traverser sul,
         FFTWRealImage::const_traverser slr,
         FFTWRealImage::Accessor src,
         DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    int w = slr.x - sul.x;

    // test for right memory layout (fftw expects a width*height fftw_real array)
    if (&(*(sul + Diff2D(w, 0))) == &(*(sul + Diff2D(0, 1))))
        fourierTransformRealImpl(sul, slr, dul, dest,
                                 norm, FFTW_RODFT00, FFTW_REDFT00);
    else
        fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                 norm, FFTW_RODFT00, FFTW_REDFT00);
}

/********************************************************************/

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealEO(triple<SrcTraverser, SrcTraverser, SrcAccessor> src,
                               pair<DestTraverser, DestAccessor> dest, fftw_real norm)
{
    fourierTransformRealEO(src.first, src.second, src.third,
                                   dest.first, dest.second, norm);
}

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealEO(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                               DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                      norm, FFTW_REDFT00, FFTW_RODFT00);
}

template <class DestTraverser, class DestAccessor>
inline void
fourierTransformRealEO(
         FFTWRealImage::const_traverser sul,
         FFTWRealImage::const_traverser slr,
         FFTWRealImage::Accessor src,
         DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    int w = slr.x - sul.x;

    // test for right memory layout (fftw expects a width*height fftw_real array)
    if (&(*(sul + Diff2D(w, 0))) == &(*(sul + Diff2D(0, 1))))
        fourierTransformRealImpl(sul, slr, dul, dest,
                                 norm, FFTW_REDFT00, FFTW_RODFT00);
    else
        fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                 norm, FFTW_REDFT00, FFTW_RODFT00);
}

/********************************************************************/

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealOO(triple<SrcTraverser, SrcTraverser, SrcAccessor> src,
                               pair<DestTraverser, DestAccessor> dest, fftw_real norm)
{
    fourierTransformRealOO(src.first, src.second, src.third,
                                   dest.first, dest.second, norm);
}

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
inline void
fourierTransformRealOO(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                               DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                      norm, FFTW_RODFT00, FFTW_RODFT00);
}

template <class DestTraverser, class DestAccessor>
inline void
fourierTransformRealOO(
         FFTWRealImage::const_traverser sul,
         FFTWRealImage::const_traverser slr,
         FFTWRealImage::Accessor src,
         DestTraverser dul, DestAccessor dest, fftw_real norm)
{
    int w = slr.x - sul.x;

    // test for right memory layout (fftw expects a width*height fftw_real array)
    if (&(*(sul + Diff2D(w, 0))) == &(*(sul + Diff2D(0, 1))))
        fourierTransformRealImpl(sul, slr, dul, dest,
                                 norm, FFTW_RODFT00, FFTW_RODFT00);
    else
        fourierTransformRealWorkImageImpl(sul, slr, src, dul, dest,
                                 norm, FFTW_RODFT00, FFTW_RODFT00);
}

/*******************************************************************/

template <class SrcTraverser, class SrcAccessor,
          class DestTraverser, class DestAccessor>
void
fourierTransformRealWorkImageImpl(SrcTraverser sul, SrcTraverser slr, SrcAccessor src,
                                  DestTraverser dul, DestAccessor dest,
                                  fftw_real norm, fftw_r2r_kind kindx, fftw_r2r_kind kindy)
{
    FFTWRealImage workImage(slr - sul);
    copyImage(srcIterRange(sul, slr, src), destImage(workImage));
    FFTWRealImage const & cworkImage = workImage;
    fourierTransformRealImpl(cworkImage.upperLeft(), cworkImage.lowerRight(),
                             dul, dest, norm, kindx, kindy);
}


template <class DestTraverser, class DestAccessor>
void
fourierTransformRealImpl(
         FFTWRealImage::const_traverser sul,
         FFTWRealImage::const_traverser slr,
         DestTraverser dul, DestAccessor dest,
         fftw_real norm, fftw_r2r_kind kindx, fftw_r2r_kind kindy)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;
    BasicImage<fftw_real> res(w, h);
    
    fftw_plan plan = fftw_plan_r2r_2d(h, w, 
                         (fftw_real *)&(*sul), (fftw_real *)res.begin(),
                         kindy, kindx, FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    if(norm != 1.0)
        transformImage(srcImageRange(res), destIter(dul, dest),
                       std::bind1st(std::multiplies<fftw_real>(), 1.0 / norm));
    else
        copyImage(srcImageRange(res), destIter(dul, dest));
}


//@}

} // namespace vigra

#endif // VIGRA_FFTW3_HXX
