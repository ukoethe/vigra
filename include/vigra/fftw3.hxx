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

#ifndef VIGRA_FFTW3_HXX
#define VIGRA_FFTW3_HXX

#include <cmath>
#include <functional>
#include <complex>
#include "stdimage.hxx"
#include "copyimage.hxx"
#include "transformimage.hxx"
#include "combineimages.hxx"
#include "numerictraits.hxx"
#include "imagecontainer.hxx"
#include <fftw3.h>

namespace vigra {

typedef double fftw_real;

template <class T>
struct FFTWReal;

template <>
struct FFTWReal<fftw_complex>
{
    typedef double type;
};

template <>
struct FFTWReal<fftwf_complex>
{
    typedef float type;
};

template <>
struct FFTWReal<fftwl_complex>
{
    typedef long double type;
};

template <class T>
struct FFTWReal2Complex;

template <>
struct FFTWReal2Complex<double>
{
    typedef fftw_complex type;
    typedef fftw_plan plan_type;
};

template <>
struct FFTWReal2Complex<float>
{
    typedef fftwf_complex type;
    typedef fftwf_plan plan_type;
};

template <>
struct FFTWReal2Complex<long double>
{
    typedef fftwl_complex type;
    typedef fftwl_plan plan_type;
};

/********************************************************/
/*                                                      */
/*                    FFTWComplex                       */
/*                                                      */
/********************************************************/

/** \brief Wrapper class for the FFTW complex types '<TT>fftw_complex</TT>'.

    This class encapsulates the low-level complex number types provided by the 
    <a href="http://www.fftw.org/">FFTW Fast Fourier Transform</a> library (i.e. 
    '<TT>fftw_complex</TT>', '<TT>fftwf_complex</TT>', '<TT>fftwl_complex</TT>'). 
    In particular, it provides constructors, member functions and 
    \ref FFTWComplexOperators "arithmetic operators" that make FFTW complex numbers
    compatible with <tt>std::complex</tt>. In addition, the class defines 
    transformations to polar coordinates and \ref FFTWComplexAccessors "accessors".

    FFTWComplex implements the concepts \ref AlgebraicField and
    \ref DivisionAlgebra. The standard image types <tt>FFTWRealImage</tt>
    and <tt>FFTWComplexImage</tt> are defined.

    <b>See also:</b>
    <ul>
        <li> \ref FFTWComplexTraits
        <li> \ref FFTWComplexOperators
        <li> \ref FFTWComplexAccessors
    </ul>

    <b>\#include</b> \<vigra/fftw3.hxx\> (for FFTW 3) or<br>
    <b>\#include</b> \<vigra/fftw.hxx\> (for deprecated FFTW 2)<br>
    Namespace: vigra
*/
template <class Real = double>
class FFTWComplex
{
  public:
        /** The wrapped complex type
        */
      typedef typename FFTWReal2Complex<Real>::type complex_type;

        /** The complex' component type, as defined in '<TT>fftw3.h</TT>'
        */
    typedef Real value_type;

        /** reference type (result of operator[])
        */
    typedef value_type & reference;

        /** const reference type (result of operator[] const)
        */
    typedef value_type const & const_reference;

        /** iterator type (result of begin() )
        */
    typedef value_type * iterator;

        /** const iterator type (result of begin() const)
        */
    typedef value_type const * const_iterator;

        /** The norm type (result of magnitude())
        */
    typedef value_type NormType;

        /** The squared norm type (result of squaredMagnitde())
        */
    typedef value_type SquaredNormType;

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

        /** Copy constructor.
        */
    template <class U>
    FFTWComplex(FFTWComplex<U> const & o)
    {
        data_[0] = (Real)o.real();
        data_[1] = (Real)o.imag();
    }

        /** Construct from plain <TT>fftw_complex</TT>.
        */
    FFTWComplex(fftw_complex const & o)
    {
        data_[0] = (Real)o[0];
        data_[1] = (Real)o[1];
    }

        /** Construct from plain <TT>fftwf_complex</TT>.
        */
    FFTWComplex(fftwf_complex const & o)
    {
        data_[0] = (Real)o[0];
        data_[1] = (Real)o[1];
    }

        /** Construct from plain <TT>fftwl_complex</TT>.
        */
    FFTWComplex(fftwl_complex const & o)
    {
        data_[0] = (Real)o[0];
        data_[1] = (Real)o[1];
    }

        /** Construct from std::complex.
        */
    template <class T>
    FFTWComplex(std::complex<T> const & o)
    {
        data_[0] = (Real)o.real();
        data_[1] = (Real)o.imag();
    }

        /** Construct from TinyVector.
        */
    template <class T>
    FFTWComplex(TinyVector<T, 2> const & o)
    {
        data_[0] = (Real)o[0];
        data_[1] = (Real)o[1];
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
    template <class U>
    FFTWComplex& operator=(FFTWComplex<U> const & o)
    {
        data_[0] = (Real)o.real();
        data_[1] = (Real)o.imag();
        return *this;
    }

        /** Assignment.
        */
    FFTWComplex& operator=(fftw_complex const & o)
    {
        data_[0] = (Real)o[0];
        data_[1] = (Real)o[1];
        return *this;
    }

        /** Assignment.
        */
    FFTWComplex& operator=(fftwf_complex const & o)
    {
        data_[0] = (Real)o[0];
        data_[1] = (Real)o[1];
        return *this;
    }

        /** Assignment.
        */
    FFTWComplex& operator=(fftwl_complex const & o)
    {
        data_[0] = (Real)o[0];
        data_[1] = (Real)o[1];
        return *this;
    }

        /** Assignment.
        */
    FFTWComplex& operator=(double o)
    {
        data_[0] = (Real)o;
        data_[1] = 0.0;
        return *this;
    }

        /** Assignment.
        */
    FFTWComplex& operator=(float o)
    {
        data_[0] = (Real)o;
        data_[1] = 0.0;
        return *this;
    }

        /** Assignment.
        */
    FFTWComplex& operator=(long double o)
    {
        data_[0] = (Real)o;
        data_[1] = 0.0;
        return *this;
    }

        /** Assignment.
        */
    template <class T>
    FFTWComplex& operator=(TinyVector<T, 2> const & o)
    {
        data_[0] = (Real)o[0];
        data_[1] = (Real)o[1];
        return *this;
    }

        /** Assignment.
        */
    template <class T>
    FFTWComplex& operator=(std::complex<T> const & o)
    {
        data_[0] = (Real)o.real();
        data_[1] = (Real)o.imag();
        return *this;
    }

    reference re()
        { return data_[0]; }

    const_reference re() const
        { return data_[0]; }

    reference real()
        { return data_[0]; }

    const_reference real() const
        { return data_[0]; }

    reference im()
        { return data_[1]; }

    const_reference im() const
        { return data_[1]; }

    reference imag()
        { return data_[1]; }

    const_reference imag() const
        { return data_[1]; }

        /** Unary negation.
        */
    FFTWComplex operator-() const
        { return FFTWComplex(-data_[0], -data_[1]); }

        /** Squared magnitude x*conj(x)
        */
    SquaredNormType squaredMagnitude() const
        { return data_[0]*data_[0]+data_[1]*data_[1]; }

        /** Magnitude (length of radius vector).
        */
    NormType magnitude() const
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

  private:
    complex_type data_;
};

/********************************************************/
/*                                                      */
/*                     FFTWComplexTraits                */
/*                                                      */
/********************************************************/

/** \page FFTWComplexTraits Numeric and Promote Traits of FFTWComplex

    The numeric and promote traits for fftw_complex and FFTWComplex follow
    the general specifications for \ref NumericPromotionTraits and
    \ref AlgebraicField. They are explicitly specialized for the types
    involved:

    \code

    template<>
    struct NumericTraits<fftw_complex>
    {
        typedef fftw_complex Promote;
        typedef fftw_complex RealPromote;
        typedef fftw_complex ComplexPromote;
        typedef fftw_real    ValueType;

        typedef VigraFalseType isIntegral;
        typedef VigraFalseType isScalar;
        typedef VigraFalseType isOrdered;
        typedef VigraTrueType  isComplex;

        // etc.
    };

    template<class Real>
    struct NumericTraits<FFTWComplex<Real> >
    {
        typedef FFTWComplex<Real> Promote;
        typedef FFTWComplex<Real> RealPromote;
        typedef FFTWComplex<Real> ComplexPromote;
        typedef Real              ValueType;

        typedef VigraFalseType isIntegral;
        typedef VigraFalseType isScalar;
        typedef VigraFalseType isOrdered;
        typedef VigraTrueType  isComplex;

        // etc.
    };

    template<>
    struct NormTraits<fftw_complex>
    {
        typedef fftw_complex Type;
        typedef fftw_real    SquaredNormType;
        typedef fftw_real    NormType;
    };

    template<class Real>
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

    template <class Real>
    struct PromoteTraits<FFTWComplex<Real>, FFTWComplex<Real> >
    {
        typedef FFTWComplex<Real> Promote;
    };

    template <class Real>
    struct PromoteTraits<FFTWComplex<Real>, double>
    {
        typedef FFTWComplex<Real> Promote;
    };

    template <class Real>
    struct PromoteTraits<double, FFTWComplex<Real> >
    {
        typedef FFTWComplex<Real> Promote;
    };
    \endcode

    <b>\#include</b> \<vigra/fftw3.hxx\> (for FFTW 3) or<br>
    <b>\#include</b> \<vigra/fftw.hxx\> (for deprecated FFTW 2)<br>
    Namespace: vigra

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

    static FFTWComplex<> zero() { return FFTWComplex<>(0.0, 0.0); }
    static FFTWComplex<> one() { return FFTWComplex<>(1.0, 0.0); }
    static FFTWComplex<> nonZero() { return one(); }

    static const Promote & toPromote(const Type & v) { return v; }
    static const RealPromote & toRealPromote(const Type & v) { return v; }
    static const Type & fromPromote(const Promote & v) { return v; }
    static const Type & fromRealPromote(const RealPromote & v) { return v; }
};

template<class Real>
struct NumericTraits<FFTWComplex<Real> >
{
    typedef FFTWComplex<Real> Type;
    typedef FFTWComplex<Real> Promote;
    typedef FFTWComplex<Real> RealPromote;
    typedef FFTWComplex<Real> ComplexPromote;
    typedef typename Type::value_type ValueType;

    typedef VigraFalseType isIntegral;
    typedef VigraFalseType isScalar;
    typedef typename NumericTraits<ValueType>::isSigned isSigned;
    typedef VigraFalseType isOrdered;
    typedef VigraTrueType  isComplex;

    static Type zero() { return Type(0.0, 0.0); }
    static Type one() { return Type(1.0, 0.0); }
    static Type nonZero() { return one(); }

    static const Promote & toPromote(const Type & v) { return v; }
    static const RealPromote & toRealPromote(const Type & v) { return v; }
    static const Type & fromPromote(const Promote & v) { return v; }
    static const Type & fromRealPromote(const RealPromote & v) { return v; }
};

template<>
struct NormTraits<fftw_complex>
{
    typedef fftw_complex Type;
    typedef fftw_real    SquaredNormType;
    typedef fftw_real    NormType;
};

template<class Real>
struct NormTraits<FFTWComplex<Real> >
{
    typedef FFTWComplex<Real>  Type;
    typedef typename Type::SquaredNormType   SquaredNormType;
    typedef typename Type::NormType   NormType;
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

template <class Real>
struct PromoteTraits<FFTWComplex<Real>, FFTWComplex<Real> >
{
    typedef FFTWComplex<Real> Promote;
};

template <class Real>
struct PromoteTraits<FFTWComplex<Real>, double>
{
    typedef FFTWComplex<Real> Promote;
};

template <class Real>
struct PromoteTraits<double, FFTWComplex<Real> >
{
    typedef FFTWComplex<Real> Promote;
};

template<class T>
struct CanSkipInitialization<std::complex<T> >
{
    typedef typename CanSkipInitialization<T>::type type;
    static const bool value = type::asBool;
};

template<class Real>
struct CanSkipInitialization<FFTWComplex<Real> >
{
    typedef typename CanSkipInitialization<Real>::type type;
    static const bool value = type::asBool;
};

namespace multi_math {

template <class ARG>
struct MultiMathOperand;

template <class Real>
struct MultiMathOperand<FFTWComplex<Real> >
{
    typedef MultiMathOperand<FFTWComplex<Real> > AllowOverload;
    typedef FFTWComplex<Real> result_type;
    
    static const int ndim = 0;
    
    MultiMathOperand(FFTWComplex<Real> const & v)
    : v_(v)
    {}
    
    template <class SHAPE>
    bool checkShape(SHAPE const &) const
    {
        return true;
    }
    
    template <class SHAPE>
    FFTWComplex<Real> const & operator[](SHAPE const &) const
    {
        return v_;
    }

    void inc(unsigned int /*LEVEL*/) const
    {}

    void reset(unsigned int /*LEVEL*/) const
    {}
    
    FFTWComplex<Real> const & operator*() const
    {
        return v_;
    }
    
    FFTWComplex<Real> v_;
};

} // namespace multi_math

template<class Ty>
class FFTWAllocator
{
  public:
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef Ty *pointer;
    typedef const Ty *const_pointer;
    typedef Ty& reference;
    typedef const Ty& const_reference;
    typedef Ty value_type;
    
    pointer address(reference val) const
        { return &val; }
        
    const_pointer address(const_reference val) const
        { return &val; }
        
    template<class Other>
    struct rebind
    {
        typedef FFTWAllocator<Other> other;
    };
    
    FFTWAllocator() throw()
    {}
    
    template<class Other>
    FFTWAllocator(const FFTWAllocator<Other>& right) throw()
    {}
    
    template<class Other>
    FFTWAllocator& operator=(const FFTWAllocator<Other>& right)
    {
        return *this;
    }
    
    pointer allocate(size_type count, void * = 0)
    {
        return (pointer)fftw_malloc(count * sizeof(Ty));
    }
    
    void deallocate(pointer ptr, size_type count)
    {
        fftw_free(ptr);
    }
    
    void construct(pointer ptr, const Ty& val)
    {
        new(ptr) Ty(val);
        
    }
    
    void destroy(pointer ptr)
    {
        ptr->~Ty();
    }
    
    size_type max_size() const throw()
    {
        return NumericTraits<std::ptrdiff_t>::max() / sizeof(Ty);
    }
};

} // namespace vigra

namespace std {

template<class Real>
class allocator<vigra::FFTWComplex<Real> >
{
  public:
    typedef vigra::FFTWComplex<Real> value_type;
    typedef size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef value_type *pointer;
    typedef const value_type *const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    pointer address(reference val) const
        { return &val; }
        
    const_pointer address(const_reference val) const
        { return &val; }
        
    template<class Other>
    struct rebind
    {
        typedef allocator<Other> other;
    };
    
    allocator() throw()
    {}
    
    template<class Other>
    allocator(const allocator<Other>& right) throw()
    {}
    
    template<class Other>
    allocator& operator=(const allocator<Other>& right)
    {
        return *this;
    }
    
    pointer allocate(size_type count, void * = 0)
    {
        return (pointer)fftw_malloc(count * sizeof(value_type));
    }

    void deallocate(pointer ptr, size_type /*count*/)
    {
        fftw_free(ptr);
    }
    
    void construct(pointer ptr, const value_type& val)
    {
        new(ptr) value_type(val);
        
    }
    
    void destroy(pointer ptr)
    {
        ptr->~value_type();
    }
    
    size_type max_size() const throw()
    {
        return vigra::NumericTraits<std::ptrdiff_t>::max() / sizeof(value_type);
    }
};

} // namespace std

namespace vigra {

/********************************************************/
/*                                                      */
/*                    FFTWComplex Operations            */
/*                                                      */
/********************************************************/

/** \addtogroup FFTWComplexOperators Functions for FFTWComplex

    <b>\#include</b> \<vigra/fftw3.hxx\> (for FFTW 3) or<br>
    <b>\#include</b> \<vigra/fftw.hxx\> (for deprecated FFTW 2)<br>

    These functions fulfill the requirements of an Algebraic Field.
    Return types are determined according to \ref FFTWComplexTraits.

    Namespace: vigra
    <p>

 */
//@{
    /// equal
template <class R>
inline bool operator ==(FFTWComplex<R> const &a, const FFTWComplex<R> &b) {
    return a.re() == b.re() && a.im() == b.im();
}

template <class R>
inline bool operator ==(FFTWComplex<R> const &a,  double b) {
    return a.re() == b && a.im() == 0.0;
}

template <class R>
inline bool operator ==(double a, const FFTWComplex<R> &b) {
    return a == b.re() && 0.0 == b.im();
}

    /// not equal
template <class R>
inline bool operator !=(FFTWComplex<R> const &a, const FFTWComplex<R> &b) {
    return a.re() != b.re() || a.im() != b.im();
}

    /// not equal
template <class R>
inline bool operator !=(FFTWComplex<R> const &a,  double b) {
    return a.re() != b || a.im() != 0.0;
}

    /// not equal
template <class R>
inline bool operator !=(double a, const FFTWComplex<R> &b) {
    return a != b.re() || 0.0 != b.im();
}

    /// add-assignment
template <class R>
inline FFTWComplex<R> & operator +=(FFTWComplex<R> &a, const FFTWComplex<R> &b) {
    a.re() += b.re();
    a.im() += b.im();
    return a;
}

    /// subtract-assignment
template <class R>
inline FFTWComplex<R> & operator -=(FFTWComplex<R> &a, const FFTWComplex<R> &b) {
    a.re() -= b.re();
    a.im() -= b.im();
    return a;
}

    /// multiply-assignment
template <class R>
inline FFTWComplex<R> & operator *=(FFTWComplex<R> &a, const FFTWComplex<R> &b) {
    typename FFTWComplex<R>::value_type t = a.re()*b.re()-a.im()*b.im();
    a.im() = a.re()*b.im()+a.im()*b.re();
    a.re() = t;
    return a;
}

    /// divide-assignment
template <class R>
inline FFTWComplex<R> & operator /=(FFTWComplex<R> &a, const FFTWComplex<R> &b) {
    typename FFTWComplex<R>::value_type sm = b.squaredMagnitude();
    typename FFTWComplex<R>::value_type t = (a.re()*b.re()+a.im()*b.im())/sm;
    a.im() = (b.re()*a.im()-a.re()*b.im())/sm;
    a.re() = t;
    return a;
}

    /// add-assignment with scalar double
template <class R>
inline FFTWComplex<R> & operator +=(FFTWComplex<R> &a, double b) {
    a.re() += (R)b;
    return a;
}

    /// subtract-assignment with scalar double
template <class R>
inline FFTWComplex<R> & operator -=(FFTWComplex<R> &a, double b) {
    a.re() -= (R)b;
    return a;
}

    /// multiply-assignment with scalar double
template <class R>
inline FFTWComplex<R> & operator *=(FFTWComplex<R> &a, double b) {
    a.re() *= (R)b;
    a.im() *= (R)b;
    return a;
}

    /// divide-assignment with scalar double
template <class R>
inline FFTWComplex<R> & operator /=(FFTWComplex<R> &a, double b) {
    a.re() /= (R)b;
    a.im() /= (R)b;
    return a;
}

    /// addition
template <class R>
inline FFTWComplex<R> operator +(FFTWComplex<R> a, const FFTWComplex<R> &b) {
    a += b;
    return a;
}

    /// right addition with scalar double
template <class R>
inline FFTWComplex<R> operator +(FFTWComplex<R> a, double b) {
    a += b;
    return a;
}

    /// left addition with scalar double
template <class R>
inline FFTWComplex<R> operator +(double a, FFTWComplex<R> b) {
    b += a;
    return b;
}

    /// subtraction
template <class R>
inline FFTWComplex<R> operator -(FFTWComplex<R> a, const FFTWComplex<R> &b) {
    a -= b;
    return a;
}

    /// right subtraction with scalar double
template <class R>
inline FFTWComplex<R> operator -(FFTWComplex<R> a, double b) {
    a -= b;
    return a;
}

    /// left subtraction with scalar double
template <class R>
inline FFTWComplex<R> operator -(double a, FFTWComplex<R> const & b) {
    return (-b) += a;
}

    /// multiplication
template <class R>
inline FFTWComplex<R> operator *(FFTWComplex<R> a, const FFTWComplex<R> &b) {
    a *= b;
    return a;
}

    /// right multiplication with scalar double
template <class R>
inline FFTWComplex<R> operator *(FFTWComplex<R> a, double b) {
    a *= b;
    return a;
}

    /// left multiplication with scalar double
template <class R>
inline FFTWComplex<R> operator *(double a, FFTWComplex<R> b) {
    b *= a;
    return b;
}

    /// division
template <class R>
inline FFTWComplex<R> operator /(FFTWComplex<R> a, const FFTWComplex<R> &b) {
    a /= b;
    return a;
}

    /// right division with scalar double
template <class R>
inline FFTWComplex<R> operator /(FFTWComplex<R> a, double b) {
    a /= b;
    return a;
}

using VIGRA_CSTD::abs;

    /// absolute value (= magnitude)
template <class R>
inline typename FFTWComplex<R>::NormType abs(const FFTWComplex<R> &a)
{
    return a.magnitude();
}

    /// phase
template <class R>
inline R arg(const FFTWComplex<R> &a)
{
    return a.phase();
}

    /// real part
template <class R>
inline R real(const FFTWComplex<R> &a)
{
    return a.real();
}

    /// imaginary part
template <class R>
inline R imag(const FFTWComplex<R> &a)
{
    return a.imag();
}

    /// complex conjugate
template <class R>
inline FFTWComplex<R> conj(const FFTWComplex<R> &a)
{
    return FFTWComplex<R>(a.re(), -a.im());
}

    /// norm (= magnitude)
template <class R>
inline typename FFTWComplex<R>::NormType norm(const FFTWComplex<R> &a)
{
    return a.magnitude();
}

    /// squared norm (= squared magnitude)
template <class R>
inline typename FFTWComplex<R>::SquaredNormType squaredNorm(const FFTWComplex<R> &a)
{
    return a.squaredMagnitude();
}

#define VIGRA_DEFINE_FFTW_COMPLEX_FUNCTION(fct) \
template <class R> \
inline FFTWComplex<R> fct(const FFTWComplex<R> &a) \
{ \
    return std::fct(reinterpret_cast<std::complex<R> const &>(a)); \
}

VIGRA_DEFINE_FFTW_COMPLEX_FUNCTION(cos)
VIGRA_DEFINE_FFTW_COMPLEX_FUNCTION(cosh)
VIGRA_DEFINE_FFTW_COMPLEX_FUNCTION(exp)
VIGRA_DEFINE_FFTW_COMPLEX_FUNCTION(log)
VIGRA_DEFINE_FFTW_COMPLEX_FUNCTION(log10)
VIGRA_DEFINE_FFTW_COMPLEX_FUNCTION(sin)
VIGRA_DEFINE_FFTW_COMPLEX_FUNCTION(sinh)
VIGRA_DEFINE_FFTW_COMPLEX_FUNCTION(sqrt)
VIGRA_DEFINE_FFTW_COMPLEX_FUNCTION(tan)
VIGRA_DEFINE_FFTW_COMPLEX_FUNCTION(tanh)

#undef VIGRA_DEFINE_FFTW_COMPLEX_FUNCTION

template <class R>
inline FFTWComplex<R> pow(const FFTWComplex<R> &a, int e)
{
    return std::pow(reinterpret_cast<std::complex<R> const &>(a), e);
}

template <class R>
inline FFTWComplex<R> pow(const FFTWComplex<R> &a, R const & e)
{
    return std::pow(reinterpret_cast<std::complex<R> const &>(a), e);
}

template <class R>
inline FFTWComplex<R> pow(const FFTWComplex<R> &a, const FFTWComplex<R> & e)
{
    return std::pow(reinterpret_cast<std::complex<R> const &>(a), 
                     reinterpret_cast<std::complex<R> const &>(e));
}

template <class R>
inline FFTWComplex<R> pow(R const & a, const FFTWComplex<R> &e)
{
    return std::pow(a, reinterpret_cast<std::complex<R> const &>(e));
}

//@}

} // namespace vigra

namespace std {

template <class Real>
ostream & operator<<(ostream & s, vigra::FFTWComplex<Real> const & v)
{
    s << std::complex<Real>(v.re(), v.im());
    return s;
}

} // namespace std

namespace vigra {

/** \addtogroup StandardImageTypes
*/
//@{

/********************************************************/
/*                                                      */
/*                      FFTWRealImage                   */
/*                                                      */
/********************************************************/

    /** Float (<tt>fftw_real</tt>) image.

        The type <tt>fftw_real</tt> is defined as <tt>double</tt> (in FFTW 2 it used to be
        either <tt>float</tt> or <tt>double</tt>, as specified during compilation of FFTW).
        FFTWRealImage uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and
        their const counterparts to access the data.

        <b>\#include</b> \<vigra/fftw3.hxx\> (for FFTW 3) or<br>
        <b>\#include</b> \<vigra/fftw.hxx\> (for deprecated FFTW 2)<br>
        Namespace: vigra
    */
typedef BasicImage<fftw_real> FFTWRealImage;

/********************************************************/
/*                                                      */
/*                     FFTWComplexImage                 */
/*                                                      */
/********************************************************/

template<class R>
struct IteratorTraits<
        BasicImageIterator<FFTWComplex<R>, FFTWComplex<R> **> >
{
    typedef FFTWComplex<R> Type;
    typedef BasicImageIterator<Type, Type **>       Iterator;
    typedef Iterator                                iterator;
    typedef BasicImageIterator<Type, Type **>       mutable_iterator;
    typedef ConstBasicImageIterator<Type, Type **>  const_iterator;
    typedef typename iterator::iterator_category    iterator_category;
    typedef typename iterator::value_type           value_type;
    typedef typename iterator::reference            reference;
    typedef typename iterator::index_reference      index_reference;
    typedef typename iterator::pointer              pointer;
    typedef typename iterator::difference_type      difference_type;
    typedef typename iterator::row_iterator         row_iterator;
    typedef typename iterator::column_iterator      column_iterator;
    typedef VectorAccessor<Type>                    default_accessor;
    typedef VectorAccessor<Type>                    DefaultAccessor;
    typedef VigraTrueType                           hasConstantStrides;
};

template<class R>
struct IteratorTraits<
        ConstBasicImageIterator<FFTWComplex<R>, FFTWComplex<R> **> >
{
    typedef FFTWComplex<R> Type;
    typedef ConstBasicImageIterator<Type, Type **> Iterator;
    typedef Iterator                               iterator;
    typedef BasicImageIterator<Type, Type **>      mutable_iterator;
    typedef ConstBasicImageIterator<Type, Type **> const_iterator;
    typedef typename iterator::iterator_category   iterator_category;
    typedef typename iterator::value_type          value_type;
    typedef typename iterator::reference           reference;
    typedef typename iterator::index_reference     index_reference;
    typedef typename iterator::pointer             pointer;
    typedef typename iterator::difference_type     difference_type;
    typedef typename iterator::row_iterator        row_iterator;
    typedef typename iterator::column_iterator     column_iterator;
    typedef VectorAccessor<Type>                   default_accessor;
    typedef VectorAccessor<Type>                   DefaultAccessor;
    typedef VigraTrueType                          hasConstantStrides;
};

    /** Complex (FFTWComplex) image.
        It uses \ref vigra::BasicImageIterator and \ref vigra::StandardAccessor and
        their const counterparts to access the data.

        <b>\#include</b> \<vigra/fftw3.hxx\> (for FFTW 3) or<br>
        <b>\#include</b> \<vigra/fftw.hxx\> (for deprecated FFTW 2)<br>
        Namespace: vigra
    */
typedef BasicImage<FFTWComplex<> > FFTWComplexImage;

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

    <b>\#include</b> \<vigra/fftw3.hxx\> (for FFTW 3) or<br>
    <b>\#include</b> \<vigra/fftw.hxx\> (for deprecated FFTW 2)<br>
    Namespace: vigra
    */
template <class Real = double>
class FFTWRealAccessor
{
  public:

        /// The accessor's value type.
    typedef Real value_type;

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

        /// Write real part at iterator position from a scalar.
    template <class ITERATOR>
    void set(value_type const & v, ITERATOR const & i) const {
        (*i).re()= v;
    }

        /// Write real part at offset from iterator position from a scalar.
    template <class ITERATOR, class DIFFERENCE>
    void set(value_type const & v, ITERATOR const & i, DIFFERENCE d) const {
        i[d].re()= v;
    }

        /// Write real part at iterator position into a scalar.
    template <class R, class ITERATOR>
    void set(FFTWComplex<R> const & v, ITERATOR const & i) const {
        *i = v.re();
    }

        /// Write real part at offset from iterator position into a scalar.
    template <class R, class ITERATOR, class DIFFERENCE>
    void set(FFTWComplex<R> const & v, ITERATOR const & i, DIFFERENCE d) const {
        i[d] = v.re();
    }
};

    /** Encapsulate access to the the imaginary part of a complex number.

    <b>\#include</b> \<vigra/fftw3.hxx\> (for FFTW 3) or<br>
    <b>\#include</b> \<vigra/fftw.hxx\> (for deprecated FFTW 2)<br>
    Namespace: vigra
    */
template <class Real = double>
class FFTWImaginaryAccessor
{
  public:
        /// The accessor's value type.
    typedef Real value_type;

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

        /// Write imaginary part at iterator position from a scalar.
    template <class ITERATOR>
    void set(value_type const & v, ITERATOR const & i) const {
        (*i).im()= v;
    }

        /// Write imaginary part at offset from iterator position from a scalar.
    template <class ITERATOR, class DIFFERENCE>
    void set(value_type const & v, ITERATOR const & i, DIFFERENCE d) const {
        i[d].im()= v;
    }

        /// Write imaginary part at iterator position into a scalar.
    template <class R, class ITERATOR>
    void set(FFTWComplex<R> const & v, ITERATOR const & i) const {
        *i = v.im();
    }

        /// Write imaginary part at offset from iterator position into a scalar.
    template <class R, class ITERATOR, class DIFFERENCE>
    void set(FFTWComplex<R> const & v, ITERATOR const & i, DIFFERENCE d) const {
        i[d] = v.im();
    }
};

    /** Write a real number into a complex one. The imaginary part is set
        to 0.

    <b>\#include</b> \<vigra/fftw3.hxx\> (for FFTW 3) or<br>
    <b>\#include</b> \<vigra/fftw.hxx\> (for deprecated FFTW 2)<br>
    Namespace: vigra
    */
template <class Real = double>
class FFTWWriteRealAccessor
: public FFTWRealAccessor<Real>
{
  public:
        /// The accessor's value type.
    typedef Real value_type;

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

    /** Calculate squared magnitude of complex number on the fly.

    <b>\#include</b> \<vigra/fftw3.hxx\> (for FFTW 3) or<br>
    Namespace: vigra
    */
template <class Real = double>
class FFTWSquaredMagnitudeAccessor
{
  public:
        /// The accessor's value type.
    typedef Real value_type;

        /// Read squared magnitude at iterator position.
    template <class ITERATOR>
    value_type operator()(ITERATOR const & i) const {
        return (*i).squaredMagnitude();
    }

        /// Read squared magnitude at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR const & i, DIFFERENCE d) const {
        return (i[d]).squaredMagnitude();
    }
};

    /** Calculate magnitude of complex number on the fly.

    <b>\#include</b> \<vigra/fftw3.hxx\> (for FFTW 3) or<br>
    <b>\#include</b> \<vigra/fftw.hxx\> (for deprecated FFTW 2)<br>
    Namespace: vigra
    */
template <class Real = double>
class FFTWMagnitudeAccessor
{
  public:
        /// The accessor's value type.
    typedef Real value_type;

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

    /** Calculate natural logarithm of magnitude of complex number on the fly.

    <b>\#include</b> \<vigra/fftw3.hxx\> (for FFTW 3) or<br>
    <b>\#include</b> \<vigra/fftw.hxx\> (for deprecated FFTW 2)<br>
    Namespace: vigra
    */
template <class Real = double>
class FFTWLogMagnitudeAccessor
{
  public:
        /// The accessor's value type.
    typedef Real value_type;

        /// Read natural log of magnitude at iterator position.
    template <class ITERATOR>
    value_type operator()(ITERATOR const & i) const {
        return std::log((*i).magnitude() + 1);
    }

        /// Read natural log of magnitude at offset from iterator position.
    template <class ITERATOR, class DIFFERENCE>
    value_type operator()(ITERATOR const & i, DIFFERENCE d) const {
        return std::log((i[d]).magnitude() + 1);
    }
};

    /** Calculate phase of complex number on the fly.

    <b>\#include</b> \<vigra/fftw3.hxx\> (for FFTW 3) or<br>
    <b>\#include</b> \<vigra/fftw.hxx\> (for deprecated FFTW 2)<br>
    Namespace: vigra
    */
template <class Real = double>
class FFTWPhaseAccessor
{
  public:
        /// The accessor's value type.
    typedef Real value_type;

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

/* \addtogroup FourierTransform Fast Fourier Transform

    This documentation describes the VIGRA interface to FFTW version 3. The interface
    to the old FFTW version 2 (file "vigra/fftw.hxx") is deprecated.

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

    <b>\#include</b> \<vigra/fftw3.hxx\><br>
    Namespace: vigra
*/

/** \addtogroup FourierTransform Fast Fourier Transform

    VIGRA provides a powerful C++ API for the popular <a href="http://www.fftw.org/">FFTW library</a>
    for fast Fourier transforms. There are two versions of the API: an older one based on image 
    iterators (and therefore restricted to 2D) and a new one based on \ref MultiArrayView that
    works for arbitrary dimensions. In addition, the functions \ref convolveFFT() and 
    \ref applyFourierFilter() provide an easy-to-use interface for FFT-based convolution,
    a major application of Fourier transforms.
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
    Note that the 2D versions of this transformation must not be executed in place - input
    and output images must be different. In contrast, the nD version (with MultiArrayView
    argument) always works in-place.

    <b> Declarations:</b>

    use MultiArrayView (this works in-place, with arbitrary dimension N):
    \code
    namespace vigra {
        template <unsigned int N, class T, class Stride>
        void moveDCToCenter(MultiArrayView<N, T, Stride> a);
    }
    \endcode

    pass iterators explicitly (2D only, not in-place):
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void moveDCToCenter(SrcImageIterator src_upperleft,
                               SrcImageIterator src_lowerright, SrcAccessor sa,
                               DestImageIterator dest_upperleft, DestAccessor da);
    }
    \endcode


    use argument objects in conjunction with \ref ArgumentObjectFactories (2D only, not in-place):
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class DestImageIterator, class DestAccessor>
        void moveDCToCenter(
            triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
            pair<DestImageIterator, DestAccessor> dest);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/fftw3.hxx\> (for 2D variants) <br>
    <b>\#include</b> \<vigra/multi_fft.hxx\> (for nD variants) <br>
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
doxygen_overloaded_function(template <...> void moveDCToCenter)

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void moveDCToCenter(SrcImageIterator src_upperleft,
                               SrcImageIterator src_lowerright, SrcAccessor sa,
                               DestImageIterator dest_upperleft, DestAccessor da)
{
    int w = int(src_lowerright.x - src_upperleft.x);
    int h = int(src_lowerright.y - src_upperleft.y);
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
     for a detailed description and usage examples.

     <b> Declarations:</b>

    use MultiArrayView (this works in-place, with arbitrary dimension N):
    \code
    namespace vigra {
        template <unsigned int N, class T, class Stride>
        void moveDCToUpperLeft(MultiArrayView<N, T, Stride> a);
    }
    \endcode

    pass iterators explicitly (2D only, not in-place):
     \code
        namespace vigra {
            template <class SrcImageIterator, class SrcAccessor,
                      class DestImageIterator, class DestAccessor>
            void moveDCToUpperLeft(SrcImageIterator src_upperleft,
                                   SrcImageIterator src_lowerright, SrcAccessor sa,
                                   DestImageIterator dest_upperleft, DestAccessor da);
        }
     \endcode


     use argument objects in conjunction with \ref ArgumentObjectFactories (2D only, not in-place):
     \code
        namespace vigra {
            template <class SrcImageIterator, class SrcAccessor,
                      class DestImageIterator, class DestAccessor>
            void moveDCToUpperLeft(
                triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                pair<DestImageIterator, DestAccessor> dest);
        }
     \endcode
*/
doxygen_overloaded_function(template <...> void moveDCToUpperLeft)

template <class SrcImageIterator, class SrcAccessor,
          class DestImageIterator, class DestAccessor>
void moveDCToUpperLeft(SrcImageIterator src_upperleft,
                               SrcImageIterator src_lowerright, SrcAccessor sa,
                               DestImageIterator dest_upperleft, DestAccessor da)
{
    int w = int(src_lowerright.x - src_upperleft.x);
    int h = int(src_lowerright.y - src_upperleft.y);
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

namespace detail {

template <class T>
void
fourierTransformImpl(FFTWComplexImage::const_traverser sul,
                     FFTWComplexImage::const_traverser slr, FFTWComplexImage::ConstAccessor src,
                     FFTWComplexImage::traverser dul, FFTWComplexImage::Accessor dest, T sign)
{
    int w = int(slr.x - sul.x);
    int h = int(slr.y - sul.y);

    FFTWComplexImage sworkImage, dworkImage;

    fftw_complex * srcPtr = (fftw_complex *)(&*sul);
    fftw_complex * destPtr = (fftw_complex *)(&*dul);

    // test for right memory layout (fftw expects a 2*width*height floats array)
    if (h > 1 && &(*(sul + Diff2D(w, 0))) != &(*(sul + Diff2D(0, 1))))
    {
        sworkImage.resize(w, h);
        copyImage(srcIterRange(sul, slr, src), destImage(sworkImage));
        srcPtr = (fftw_complex *)(&(*sworkImage.upperLeft()));
    }
    if (h > 1 && &(*(dul + Diff2D(w, 0))) != &(*(dul + Diff2D(0, 1))))
    {
        dworkImage.resize(w, h);
        destPtr = (fftw_complex *)(&(*dworkImage.upperLeft()));
    }

    fftw_plan plan = fftw_plan_dft_2d(h, w, srcPtr, destPtr, sign, FFTW_ESTIMATE );
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    if (h > 1 && &(*(dul + Diff2D(w, 0))) != &(*(dul + Diff2D(0, 1))))
    {
        copyImage(srcImageRange(dworkImage), destIter(dul, dest));
    }
}

} // namespace detail

/********************************************************/
/*                                                      */
/*                   fourierTransform                   */
/*                                                      */
/********************************************************/

/** \brief Compute forward and inverse Fourier transforms.

    The array referring to the spatial domain (i.e. the input in a forward transform, 
    and the output in an inverse transform) may be scalar or complex. The array representing
    the frequency domain (i.e. output for forward transform, input for inverse transform) 
    must always be complex.
    
    The new implementations (those using MultiArrayView arguments) perform a normalized transform, 
    whereas the old ones (using 2D iterators or argument objects) perform an un-normalized 
    transform (i.e. the result of the inverse transform is scaled by the number of pixels).

    In general, input and output arrays must have the same shape, with the exception of the 
    special <a href="http://www.fftw.org/doc/Multi_002dDimensional-DFTs-of-Real-Data.html">R2C 
    and C2R modes</a> defined by FFTW.
    
    The R2C transform reduces the redundancy in the Fourier representation of a real-valued signal:
    Since the Fourier representation of a real signal is symmetric, about half of the Fourier coefficients 
    can simply be dropped. By convention, this reduction is applied to the first (innermost) dimension, 
    such that <tt>fourier.shape(0) == spatial.shape(0)/2 + 1</tt> holds. The correct frequency domain
    shape can be conveniently computed by means of the function \ref fftwCorrespondingShapeR2C().
    
    Note that your program must always link against <tt>libfftw3</tt>. If you want to compute Fourier 
    transforms for <tt>float</tt> or <tt>long double</tt> arrays, you must <i>additionally</i> link against <tt>libfftw3f</tt> and <tt>libfftw3l</tt> respectively. (Old-style functions only support <tt>double</tt>).
    
    The Fourier transform functions internally create <a href="http://www.fftw.org/doc/Using-Plans.html">FFTW plans</a>
    which control the algorithm details. The plans are creates with the flag <tt>FFTW_ESTIMATE</tt>, i.e.
    optimal settings are guessed or read from saved "wisdom" files. If you need more control over planning,
    you can use the class \ref FFTWPlan.
    
    <b> Declarations:</b>

    use complex-valued MultiArrayView arguments (this works for arbitrary dimension N):
    \code
    namespace vigra {
        template <unsigned int N, class Real, class C1, class C2>
        void 
        fourierTransform(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
                         MultiArrayView<N, FFTWComplex<Real>, C2> out);

        template <unsigned int N, class Real, class C1, class C2>
        void 
        fourierTransformInverse(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
                                MultiArrayView<N, FFTWComplex<Real>, C2> out);
    }
    \endcode

    use real-valued MultiArrayView in the spatial domain, complex-valued MultiArrayView 
    in the frequency domain (this works for arbitrary dimension N, and also supports
    the R2C and C2R transform, depending on the array shape in the frequency domain):
    \code
    namespace vigra {
        template <unsigned int N, class Real, class C1, class C2>
        void 
        fourierTransform(MultiArrayView<N, Real, C1> in, 
                         MultiArrayView<N, FFTWComplex<Real>, C2> out);

        template <unsigned int N, class Real, class C1, class C2>
        void 
        fourierTransformInverse(MultiArrayView<N, FFTWComplex<Real>, C1> in, 
                                MultiArrayView<N, Real, C2> out);
    }
    \endcode

    pass iterators explicitly (2D only, double only):
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor>
        void fourierTransform(SrcImageIterator srcUpperLeft,
                              SrcImageIterator srcLowerRight, SrcAccessor src,
                              FFTWComplexImage::traverser destUpperLeft, FFTWComplexImage::Accessor dest);

        void
        fourierTransformInverse(FFTWComplexImage::const_traverser sul,
                                FFTWComplexImage::const_traverser slr, FFTWComplexImage::ConstAccessor src,
                                FFTWComplexImage::traverser dul, FFTWComplexImage::Accessor dest)
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories (2D only, double only):
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor>
        void fourierTransform(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                              pair<FFTWComplexImage::traverser, FFTWComplexImage::Accessor> dest);

        void
        fourierTransformInverse(triple<FFTWComplexImage::const_traverser,
                                       FFTWComplexImage::const_traverser, FFTWComplexImage::ConstAccessor> src,
                                pair<FFTWComplexImage::traverser, FFTWComplexImage::Accessor> dest);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/fftw3.hxx\> (old-style 2D variants)<br>
    <b>\#include</b> \<vigra/multi_fft.hxx\> (new-style nD variants)<br>
    Namespace: vigra

    old-style example (using FFTWComplexImage):
    \code
    // compute complex Fourier transform of a real image, old-style implementation
    vigra::DImage src(w, h);
    vigra::FFTWComplexImage fourier(w, h);

    fourierTransform(srcImageRange(src), destImage(fourier));

    // compute inverse Fourier transform
    // note that both source and destination image must be of type vigra::FFTWComplexImage
    vigra::FFTWComplexImage inverseFourier(w, h);

    fourierTransformInverse(srcImageRange(fourier), destImage(inverseFourier));
    \endcode
    
    new-style examples (using MultiArray):
    \code
    // compute Fourier transform of a real array, using the R2C algorithm
    MultiArray<2, double> src(Shape2(w, h));
    MultiArray<2, FFTWComplex<double> > fourier(fftwCorrespondingShapeR2C(src.shape()));

    fourierTransform(src, fourier);

    // compute inverse Fourier transform, using the C2R algorithm
    MultiArray<2, double> dest(src.shape());
    fourierTransformInverse(fourier, dest);
    \endcode

    \code
    // compute Fourier transform of a real array with standard algorithm
    MultiArray<2, double> src(Shape2(w, h));
    MultiArray<2, FFTWComplex<double> > fourier(src.shape());

    fourierTransform(src, fourier);

    // compute inverse Fourier transform, using the C2R algorithm
    MultiArray<2, double> dest(src.shape());
    fourierTransformInverse(fourier, dest);
    \endcode
    Complex input arrays are handled in the same way. 
    
*/
doxygen_overloaded_function(template <...> void fourierTransform)

inline void
fourierTransform(FFTWComplexImage::const_traverser sul,
                 FFTWComplexImage::const_traverser slr, FFTWComplexImage::ConstAccessor src,
                 FFTWComplexImage::traverser dul, FFTWComplexImage::Accessor dest)
{
    detail::fourierTransformImpl(sul, slr, src, dul, dest, FFTW_FORWARD);
}

template <class SrcImageIterator, class SrcAccessor>
void fourierTransform(SrcImageIterator srcUpperLeft,
                      SrcImageIterator srcLowerRight, SrcAccessor sa,
                      FFTWComplexImage::traverser destUpperLeft, FFTWComplexImage::Accessor da)
{
    // copy real input images into a complex one...
    int w= srcLowerRight.x - srcUpperLeft.x;
    int h= srcLowerRight.y - srcUpperLeft.y;

    FFTWComplexImage workImage(w, h);
    copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
              destImage(workImage, FFTWWriteRealAccessor<>()));

    // ...and call the complex -> complex version of the algorithm
    FFTWComplexImage const & cworkImage = workImage;
    fourierTransform(cworkImage.upperLeft(), cworkImage.lowerRight(), cworkImage.accessor(),
                     destUpperLeft, da);
}

template <class SrcImageIterator, class SrcAccessor>
inline
void fourierTransform(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                      pair<FFTWComplexImage::traverser, FFTWComplexImage::Accessor> dest)
{
    fourierTransform(src.first, src.second, src.third, dest.first, dest.second);
}

/** \brief Compute inverse Fourier transforms.

    See \ref fourierTransform() for details.
*/
doxygen_overloaded_function(template <...> void fourierTransformInverse)

inline void
fourierTransformInverse(FFTWComplexImage::const_traverser sul,
                        FFTWComplexImage::const_traverser slr, FFTWComplexImage::ConstAccessor src,
                        FFTWComplexImage::traverser dul, FFTWComplexImage::Accessor dest)
{
    detail::fourierTransformImpl(sul, slr, src, dul, dest, FFTW_BACKWARD);
}

template <class DestImageIterator, class DestAccessor>
void fourierTransformInverse(FFTWComplexImage::const_traverser sul,
                             FFTWComplexImage::const_traverser slr, FFTWComplexImage::ConstAccessor src,
                             DestImageIterator dul, DestAccessor dest)
{
    int w = slr.x - sul.x;
    int h = slr.y - sul.y;

    FFTWComplexImage workImage(w, h);
    fourierTransformInverse(sul, slr, src, workImage.upperLeft(), workImage.accessor());
    copyImage(srcImageRange(workImage), destIter(dul, dest));
}


template <class DestImageIterator, class DestAccessor>
inline void
fourierTransformInverse(triple<FFTWComplexImage::const_traverser,
                               FFTWComplexImage::const_traverser, FFTWComplexImage::ConstAccessor> src,
                        pair<DestImageIterator, DestAccessor> dest)
{
    fourierTransformInverse(src.first, src.second, src.third, dest.first, dest.second);
}



/********************************************************/
/*                                                      */
/*                   applyFourierFilter                 */
/*                                                      */
/********************************************************/

/** \brief Apply a filter (defined in the frequency domain) to an image.

    After transferring the image into the frequency domain, it is
    multiplied pixel-wise with the filter and transformed back. The
    result is put into the given destination image which must have the right size.
    The result will be normalized to compensate for the two FFTs.

    If the destination image is scalar, only the real part of the result image is
    retained. In this case, you are responsible for choosing a filter image
    which ensures a zero imaginary part of the result (e.g. use a real, even symmetric
    filter image, or a purely imaginary, odd symmetric one).

    The DC entry of the filter must be in the upper left, which is the
    position where FFTW expects it (see \ref moveDCToUpperLeft()).
    
    See also \ref convolveFFT() for corresponding functionality on the basis of the
    \ref MultiArrayView interface.

    <b> Declarations:</b>

    pass 2D array views:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class FilterImageIterator, class FilterAccessor,
                  class DestImageIterator, class DestAccessor>
        void applyFourierFilter(SrcImageIterator srcUpperLeft,
                                SrcImageIterator srcLowerRight, SrcAccessor sa,
                                FilterImageIterator filterUpperLeft, FilterAccessor fa,
                                DestImageIterator destUpperLeft, DestAccessor da);
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class FilterImageIterator, class FilterAccessor,
                  class DestImageIterator, class DestAccessor>
        void applyFourierFilter(SrcImageIterator srcUpperLeft,
                                SrcImageIterator srcLowerRight, SrcAccessor sa,
                                FilterImageIterator filterUpperLeft, FilterAccessor fa,
                                DestImageIterator destUpperLeft, DestAccessor da);
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor,
                  class FilterImageIterator, class FilterAccessor,
                  class DestImageIterator, class DestAccessor>
        void applyFourierFilter(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                                pair<FilterImageIterator, FilterAccessor> filter,
                                pair<DestImageIterator, DestAccessor> dest);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/fftw3.hxx\><br>
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
doxygen_overloaded_function(template <...> void applyFourierFilter)

template <class SrcImageIterator, class SrcAccessor,
          class FilterImageIterator, class FilterAccessor,
          class DestImageIterator, class DestAccessor>
void applyFourierFilter(SrcImageIterator srcUpperLeft,
                        SrcImageIterator srcLowerRight, SrcAccessor sa,
                        FilterImageIterator filterUpperLeft, FilterAccessor fa,
                        DestImageIterator destUpperLeft, DestAccessor da)
{
    // copy real input images into a complex one...
    int w = int(srcLowerRight.x - srcUpperLeft.x);
    int h = int(srcLowerRight.y - srcUpperLeft.y);

    FFTWComplexImage workImage(w, h);
    copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
              destImage(workImage, FFTWWriteRealAccessor<>()));

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

template <class FilterImageIterator, class FilterAccessor,
          class DestImageIterator, class DestAccessor>
void applyFourierFilterImpl(
    FFTWComplexImage::const_traverser srcUpperLeft,
    FFTWComplexImage::const_traverser srcLowerRight,
    FFTWComplexImage::ConstAccessor,
    FilterImageIterator filterUpperLeft, FilterAccessor fa,
    DestImageIterator destUpperLeft, DestAccessor da)
{
    int w = int(srcLowerRight.x - srcUpperLeft.x);
    int h = int(srcLowerRight.y - srcUpperLeft.y);

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
                     destImage(complexResultImg), std::multiplies<FFTWComplex<> >());

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
                   linearIntensityTransform<FFTWComplex<> >(1.0/(srcImage.width() * srcImage.height())));
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

    pass 2D array views:
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor, class FilterType>
        void applyFourierFilterFamily(SrcImageIterator srcUpperLeft,
                                      SrcImageIterator srcLowerRight, SrcAccessor sa,
                                      const ImageArray<FilterType> &filters,
                                      ImageArray<FFTWComplexImage> &results)
    }
    \endcode

    pass \ref ImageIterators and \ref DataAccessors :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor, class FilterType>
        void applyFourierFilterFamily(SrcImageIterator srcUpperLeft,
                                      SrcImageIterator srcLowerRight, SrcAccessor sa,
                                      const ImageArray<FilterType> &filters,
                                      ImageArray<FFTWComplexImage> &results)
    }
    \endcode

    use argument objects in conjunction with \ref ArgumentObjectFactories :
    \code
    namespace vigra {
        template <class SrcImageIterator, class SrcAccessor, class FilterType>
        void applyFourierFilterFamily(triple<SrcImageIterator, SrcImageIterator, SrcAccessor> src,
                                      const ImageArray<FilterType> &filters,
                                      ImageArray<FFTWComplexImage> &results)
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/fftw3.hxx\><br>
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
doxygen_overloaded_function(template <...> void applyFourierFilterFamily)

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
    int w = int(srcLowerRight.x - srcUpperLeft.x);
    int h = int(srcLowerRight.y - srcUpperLeft.y);

    FFTWComplexImage workImage(w, h);
    copyImage(srcIterRange(srcUpperLeft, srcLowerRight, sa),
              destImage(workImage, FFTWWriteRealAccessor<>()));

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
        int h = srcLowerRight.y - srcUpperLeft.y;
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
    // FIXME: sa is not used
    // (maybe check if StandardAccessor, else copy?)    

    // make sure the filter images have the right dimensions
    vigra_precondition((srcLowerRight - srcUpperLeft) == filters.imageSize(),
                       "applyFourierFilterFamily called with src image size != filters.imageSize()!");

    // make sure the result image array has the right dimensions
    results.resize(filters.size());
    results.resizeImages(filters.imageSize());

    // FFT from srcImage to freqImage
    int w = int(srcLowerRight.x - srcUpperLeft.x);
    int h = int(srcLowerRight.y - srcUpperLeft.y);

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
                         destImage(result), std::multiplies<FFTWComplex<> >());

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

    pass 2D array views:
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

    pass \ref ImageIterators and \ref DataAccessors :
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


    use argument objects in conjunction with \ref ArgumentObjectFactories :
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

    <b>\#include</b> \<vigra/fftw3.hxx\><br>
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
            fourier(x-1, y) = fourier(x, y) * dx * std::exp(-(dx*dx + dy*dy) * scale*scale / 2.0);
        }
        fourier(width-1, y) = 0.0;
    }

    // inverse transform -- odd symmetry in x-direction, even in y,
    //                      due to symmetry of the filter
    fourierTransformRealOE(srcImageRange(fourier), destImage(spatial),
                           (fftw_real)-4.0 * (width+1) * (height-1));
    \endcode
*/
doxygen_overloaded_function(template <...> void fourierTransformReal)

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
