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
 
 
#ifndef VIGRA_COLORCONVERSIONS_HXX
#define VIGRA_COLORCONVERSIONS_HXX

#include <cmath>
#include "mathutil.hxx"
#include "rgbvalue.hxx"
#include "functortraits.hxx"

/** \page ColorConversions  Color Space Conversions

    Convert between RGB, sRGB, R'G'B', XYZ, L*a*b*, L*u*v*, Y'PbPr, Y'CbCr, Y'IQ, and Y'UV color spaces.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    <UL>
    <LI> <b>RGB/sRGB/R'G'B'</b><br>
        <em>linear and non-linear (gamma corrected) additive color</em>
        <p>
        <DL>
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::RGB2sRGBFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::sRGB2RGBFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::RGB2RGBPrimeFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::RGBPrime2RGBFunctor
        </DL><p>
    <LI> <b>XYZ</b><br>
        <em>device independent color representation 
               (according to Publication CIE  No  15.2 "Colorimetry"
                and ITU-R Recommendation BT.709)</em>
        <p>
        <DL>
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::RGB2XYZFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::RGBPrime2XYZFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::XYZ2RGBFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::XYZ2RGBPrimeFunctor
        </DL><p>
    <LI> <b>L*a*b* </b><br>
        <em>perceptually uniform color representation 
               (according to Publication CIE No 15.2 "Colorimetry" and
               ITU-R Recommendation BT.709)</em>
        <p>
        <DL>
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::RGB2LabFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::RGBPrime2LabFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::XYZ2LabFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::Lab2RGBFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::Lab2RGBPrimeFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::Lab2XYZFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             \ref polar2Lab "vigra::polar2Lab"()
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             \ref lab2Polar "vigra::lab2Polar"()
        </DL><p>
    <LI> <b>L*u*v* </b><br>
        <em>perceptually uniform color representation 
               (according to Publication CIE No 15.2 "Colorimetry" and
               ITU-R Recommendation BT.709)</em>
        <p>
        <DL>
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::RGB2LuvFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::RGBPrime2LuvFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::XYZ2LuvFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::Luv2RGBFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::Luv2RGBPrimeFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::Luv2XYZFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             \ref polar2Luv "vigra::polar2Luv"()
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             \ref luv2Polar "vigra::luv2Polar"()
        </DL><p>
    <LI> <b>Y'PbPr and Y'CbCr </b><br>
        <em>color difference coding
                (according to ITU-R Recommendation BT. 601)</em>
        <p>
        <DL>
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::RGBPrime2YPrimePbPrFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::YPrimePbPr2RGBPrimeFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             \ref polar2YPrimePbPr "vigra::polar2YPrimePbPr"()
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             \ref yPrimePbPr2Polar "vigra::yPrimePbPr2Polar"()
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::RGBPrime2YPrimeCbCrFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::YPrimeCbCr2RGBPrimeFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             \ref polar2YPrimeCbCr "vigra::polar2YPrimeCbCr"()
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             \ref yPrimeCbCr2Polar "vigra::yPrimeCbCr2Polar"()
        </DL><p>
    <LI> <b>Y'UV and Y'IQ </b><br>
        <em>analog video coding according to NTSC and PAL standards</em>
        <p>
        <DL>
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::RGBPrime2YPrimeUVFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::YPrimeUV2RGBPrimeFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             \ref polar2YPrimeUV "vigra::polar2YPrimeUV"()
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             \ref yPrimeUV2Polar "vigra::yPrimeUV2Polar"()
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::RGBPrime2YPrimeIQFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             vigra::YPrimeIQ2RGBPrimeFunctor
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             \ref polar2YPrimeIQ "vigra::polar2YPrimeIQ"()
        <DT> <IMG BORDER=0 ALT="-" SRC="documents/bullet.gif">
             \ref yPrimeIQ2Polar "vigra::yPrimeIQ2Polar"()
        </DL><p>
    </UL>
    
    \anchor _details
    This module provides conversion from RGB/R'G'B' into more perceptually uniform
    color spaces. In image analysis, colors are usually converted into another color space 
    in order to get good estimates of perceived color differences by just calculating 
    Euclidean distances between the transformed colors. The L*a*b* and L*u*v* were 
    designed with exactly this application in mind and thus give the best results. But these
    conversions are also the most computationally demanding. The Y'PbPr color difference
    space (designed for coding digital video) is computationally much cheaper, and 
    almost as good. Y'CbCr represents esentially the same transformation, but the color values 
    are scaled so that they can be stored with 8 bits per channel with minimal loss of 
    information. The other transformations are of lesser interest here: XYZ is a device independent
    (but not perceptually uniform) color representation, and Y'IQ and Y'UV are the color 
    spaces used by the PAL and NTSC analog video standards. Detailed information about
    these color spaces and their transformations can be found in 
    <a href="http://www.poynton.com/ColorFAQ.html">Charles Poynton's Color FAQ</a>
    
    When you want to perform a color conversion, you must first know in which
    color space the data are given. Although this sounds trivial, it is
    quite often done wrong, because the distinction between RGB and sRGB (still images) or R'G'B' 
    (digital video) is frequently overlooked: nowadays, most still images are stored in
    sRGB space, and treating them as RGB leads to wrong results (although the color primaries
    are named the same). RGB and R'G'B' are related by a so called <em>gamma correction</em>:
    
    \f[
        C' = C_{max} \left(\frac{C_{RGB}}{C_{max}} \right)^{0.45} \qquad
    \f]
    
    where C represents one of the color channels R, G, and B, and \f$ C_{max} \f$ usually equals 255. 
    The sRGB color space realizes a slight enhancement of this definition:
    
    \f[
        C_{sRGB} = \begin{cases}
        12.92 C_{RGB} & \textrm{if }\frac{C_{RGB}}{C_{max}} \le 0.00304 \\
        C_{max}\left( 1.055 \left(\frac{C_{RGB}}{C_{max}}\right)^{1/2.4}-0.055\right) & \textrm{otherwise}
        \end{cases}
    \f]
    
    sRGB has now become a widely accepted international standard (IEC 61966-2.1) which is used by most 
    consumer products (digital cameras, printers, and screens). In practice, you can 
    distinguish between linear and gamma-corrected red, green, and blue by displaying the images: if they look
    too dark, they are probably RGB, if they are OK, they are likely sRGB. (However, there are still a few older 
    graphics cards and display programs which silently apply an additional gamma correction to every image, 
    so that RGB appears correct and sRGB is too bright.) Whether or not the data are represented
    in the sRGB color space can also be seen in the color space tag of an image's EXIF data, if available.
    
    The distinction between RGB and R'G'B' is important because some conversions start at 
    RGB (XYZ, L*a*b*, L*u*v*), while others start at R'G'B' (Y'PbPr, Y'CbCr, Y'IQ, and Y'UV). 
    The names of VIGRA's color conversion functors always make clear to which color space 
    they must be applied.
   
    In addition VIGRA provides a <em>polar coordinate interface</em>
    to several color spaces (L*a*b*, L*u*v*, Y'PbPr, Y'CbCr, Y'IQ, and Y'UV). This
    interface makes use of the fact that these color spaces are conceptually similar:
    they represent colors by a "brightness" coordinate (L* or Y') and a pair of 
    "chromaticity" coordinates that span a plane of colors with equal brightness.
    The polar representation transforms chroma coordinates into a color "angle"
    (similar to hue in the HSV system) and a "saturation". The polar coordinates are 
    normalized so that a color angle of 0 degrees is always associated with red
    (green is at about 120 degrees, blue at about 240 degrees - exact values differ
    between color spaces). A saturation of 1 is the highest saturation that any RGB color 
    in the unit cube can have after transformation into the respective color space, 
    and saturation 0 corresponds to gray. Polar coordinates provide a more intuitive 
    interface to color specification by users and make different color spaces somewhat 
    comparable.
*/
namespace vigra {

namespace detail
{

inline double gammaCorrection(double value, double gamma)
{
    return (value < 0.0) ? 
            -VIGRA_CSTD::pow(-value, gamma) :
            VIGRA_CSTD::pow(value, gamma);
}

inline double gammaCorrection(double value, double gamma, double norm)
{
    return (value < 0.0) ? 
            -norm*VIGRA_CSTD::pow(-value/norm, gamma) :
            norm*VIGRA_CSTD::pow(value/norm, gamma);
}

inline double sRGBCorrection(double value, double norm)
{
    value /= norm;
    return (value <= 0.00304) 
               ? norm*12.92*value 
               : norm*(1.055*VIGRA_CSTD::pow(value, 0.41666666666666667) - 0.055);
}

inline double inverse_sRGBCorrection(double value, double norm)
{
    value /= norm;
    return (value <= 0.03928) 
               ? norm*value / 12.92
               : norm*VIGRA_CSTD::pow((value + 0.055)/1.055, 2.4);
}


} // namespace detail

/** \brief Convert linear (raw) RGB into non-linear (gamma corrected) R'G'B'.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the transformation
    
    \f[
        R' = R_{max} \left(\frac{R}{R_{max}} \right)^{0.45} \qquad
        G' = G_{max} \left(\frac{G}{G_{max}} \right)^{0.45} \qquad
        B' = B_{max} \left(\frac{B}{B_{max}} \right)^{0.45}
    \f]
    
    By default, \f$ R_{max} = G_{max} = B_{max} = 255 \f$. This default can be overridden
    in the constructor. If both source and target colors components are stored 
    as <tt>unsigned char</tt>, a look-up-table will be used to speed up the transformation.

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
    */
template <class From, class To = From>
class RGB2RGBPrimeFunctor
{
  public:
  
        /** the functor's argument type
        */
    typedef TinyVector<From, 3> argument_type;
  
        /** the functor's result type
        */
    typedef RGBValue<To> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef RGBValue<To> value_type;
  
        /** the result component's promote type
        */
    typedef typename NumericTraits<To>::RealPromote component_type;
    
        /** Default constructor.
            The maximum value for each RGB component defaults to 255
        */
    RGB2RGBPrimeFunctor()
    : max_(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    RGB2RGBPrimeFunctor(component_type max)
    : max_(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & rgb) const
    {
        return RGBValue<To>(
            NumericTraits<To>::fromRealPromote(detail::gammaCorrection(rgb[0], 0.45, max_)),
            NumericTraits<To>::fromRealPromote(detail::gammaCorrection(rgb[1], 0.45, max_)),
            NumericTraits<To>::fromRealPromote(detail::gammaCorrection(rgb[2], 0.45, max_)));
    }
    
  private:
    component_type max_;    
};

template <>
class RGB2RGBPrimeFunctor<unsigned char, unsigned char>
{
    unsigned char lut_[256];
        
  public:
  
    typedef RGBValue<unsigned char> value_type;
    
    RGB2RGBPrimeFunctor()
    {
        for(int i=0; i<256; ++i)
        {
            lut_[i] = NumericTraits<unsigned char>::fromRealPromote(detail::gammaCorrection(i, 0.45, 255.0));
        }
    }
    
    RGB2RGBPrimeFunctor(double max)
    {
        for(int i=0; i<256; ++i)
        {
            lut_[i] = NumericTraits<unsigned char>::fromRealPromote(detail::gammaCorrection(i, 0.45, max));
        }
    }
    
    template <class V>
    RGBValue<unsigned char> operator()(V const & rgb) const
    {
        return RGBValue<unsigned char>(lut_[rgb[0]], lut_[rgb[1]], lut_[rgb[2]]);
    }
};

template <class From, class To>
class FunctorTraits<RGB2RGBPrimeFunctor<From, To> >
: public FunctorTraitsBase<RGB2RGBPrimeFunctor<From, To> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert linear (raw) RGB into standardized sRGB.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The sRGB color space is a slight improvement over the R'G'B' space. Is is now a widely accepted 
    international standard (IEC 61966-2.1) which is used by most consumer products
    (digital cameras, printers, and screens). The functor realizes the transformation
    
    \f[
        C_{sRGB} = \left\{ \begin{array}{ll}
        12.92\,C_{RGB} & \textrm{ if }\frac{C_{RGB}}{C_{max}} \le 0.00304 \\
        C_{max}\left( 1.055 \left(\frac{C_{RGB}}{C_{max}}\right)^{1/2.4}-0.055\right) & \textrm{ otherwise}
        \end{array}  \right.
    \f]
    
    where C is any of the primaries R, G, and B. By default, \f$ C_{max} = 255 \f$ (this default can be
    overridden in the constructor). If both source and target color components are stored
    as <tt>unsigned char</tt>, a look-up-table will be used to speed up the transformation.

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
    */
template <class From, class To = From>
class RGB2sRGBFunctor
{
  public:
  
        /** the functor's argument type
        */
    typedef TinyVector<From, 3> argument_type;
  
        /** the functor's result type
        */
    typedef RGBValue<To> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef RGBValue<To> value_type;
  
        /** the result component's promote type
        */
    typedef typename NumericTraits<To>::RealPromote component_type;
    
        /** Default constructor.
            The maximum value for each RGB component defaults to 255
        */
    RGB2sRGBFunctor()
    : max_(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    RGB2sRGBFunctor(component_type max)
    : max_(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & rgb) const
    {
        return RGBValue<To>(
            NumericTraits<To>::fromRealPromote(detail::sRGBCorrection(rgb[0], max_)),
            NumericTraits<To>::fromRealPromote(detail::sRGBCorrection(rgb[1], max_)),
            NumericTraits<To>::fromRealPromote(detail::sRGBCorrection(rgb[2], max_)));
    }
    
  private:
    component_type max_;    
};

template <>
class RGB2sRGBFunctor<unsigned char, unsigned char>
{
    unsigned char lut_[256];
        
  public:
  
    typedef RGBValue<unsigned char> value_type;
    
    RGB2sRGBFunctor()
    {
        for(int i=0; i<256; ++i)
        {
            lut_[i] = NumericTraits<unsigned char>::fromRealPromote(detail::sRGBCorrection(i, 255.0));
        }
    }
    
    RGB2sRGBFunctor(double max)
    {
        for(int i=0; i<256; ++i)
        {
            lut_[i] = NumericTraits<unsigned char>::fromRealPromote(detail::sRGBCorrection(i, max));
        }
    }
    
    template <class V>
    RGBValue<unsigned char> operator()(V const & rgb) const
    {
        return RGBValue<unsigned char>(lut_[rgb[0]], lut_[rgb[1]], lut_[rgb[2]]);
    }
};

template <class From, class To>
class FunctorTraits<RGB2sRGBFunctor<From, To> >
: public FunctorTraitsBase<RGB2sRGBFunctor<From, To> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert non-linear (gamma corrected) R'G'B' into non-linear (raw) RGB.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the transformation
    
    \f[
        R = R_{max} \left(\frac{R'}{R_{max}} \right)^{1/0.45} \qquad
        G = G_{max} \left(\frac{G'}{G_{max}} \right)^{1/0.45} \qquad
        B = B_{max} \left(\frac{B'}{B_{max}} \right)^{1/0.45}
    \f]
    
    By default, \f$ R_{max} = G_{max} = B_{max} = 255 \f$. This default can be overridden
    in the constructor. If both source and target color components are stored 
    as <tt>unsigned char</tt>, a look-up-table will be used to speed up the transformation.

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class From, class To = From>
class RGBPrime2RGBFunctor
{
  public:
  
        /** the functor's argument type
        */
    typedef TinyVector<From, 3> argument_type;
  
        /** the functor's result type
        */
    typedef RGBValue<To> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef RGBValue<To> value_type;
  
        /** the result component's promote type
        */
    typedef typename NumericTraits<To>::RealPromote component_type;
    
        /** Default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    RGBPrime2RGBFunctor()
    : max_(255.0), gamma_(1.0/0.45)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    RGBPrime2RGBFunctor(component_type max)
    : max_(max), gamma_(1.0/0.45)
    {}
    
        /** apply the transformation
        */
    result_type operator()(argument_type const & rgb) const
    {
        return RGBValue<To>(
            NumericTraits<To>::fromRealPromote(detail::gammaCorrection(rgb[0], gamma_, max_)),
            NumericTraits<To>::fromRealPromote(detail::gammaCorrection(rgb[1], gamma_, max_)),
            NumericTraits<To>::fromRealPromote(detail::gammaCorrection(rgb[2], gamma_, max_)));
    }

  private:
    component_type max_;
    double gamma_;
};

template <>
class RGBPrime2RGBFunctor<unsigned char, unsigned char>
{    
    unsigned char lut_[256];
        
  public:
  
    typedef RGBValue<unsigned char> value_type;
    
    RGBPrime2RGBFunctor()
    {
        for(int i=0; i<256; ++i)
        {
            lut_[i] = NumericTraits<unsigned char>::fromRealPromote(detail::gammaCorrection(i, 1.0/0.45, 255.0));
        }
    }
    
    RGBPrime2RGBFunctor(double max)
    {
        for(int i=0; i<256; ++i)
        {
            lut_[i] = NumericTraits<unsigned char>::fromRealPromote(detail::gammaCorrection(i, 1.0/0.45, max));
        }
    }
    
    template <class V>
    RGBValue<unsigned char> operator()(V const & rgb) const
    {
        return RGBValue<unsigned char>(lut_[rgb[0]], lut_[rgb[1]], lut_[rgb[2]]);
    }
};

template <class From, class To>
class FunctorTraits<RGBPrime2RGBFunctor<From, To> >
: public FunctorTraitsBase<RGBPrime2RGBFunctor<From, To> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert standardized sRGB into non-linear (raw) RGB.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The sRGB color space is a slight improvement over the R'G'B' space. Is is now a widely accepted 
    international standard (IEC 61966-2.1) which is used by most consumer products
    (digital cameras, printers, and screens). The functor realizes the transformation
    
    \f[
        C_{RGB} = \begin{cases}
        C_{sRGB} / 12.92 & \textrm{if }\frac{C_{sRGB}}{C_{max}} \le 0.03928 \\
        C_{max}\left( \frac{C_{sRGB}/C_{max}+0.055}{1.055}\right)^{2.4} & \textrm{otherwise}
        \end{cases}
    \f]
    
    where C is one of the color channels R, G, or B, and \f$ C_{max}\f$ equals 255 by default (This default 
    can be overridden in the constructor). If both source and target color components are stored 
    as <tt>unsigned char</tt>, a look-up-table will be used to speed up the transformation.

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class From, class To = From>
class sRGB2RGBFunctor
{
  public:
  
        /** the functor's argument type
        */
    typedef TinyVector<From, 3> argument_type;
  
        /** the functor's result type
        */
    typedef RGBValue<To> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef RGBValue<To> value_type;
  
        /** the result component's promote type
        */
    typedef typename NumericTraits<To>::RealPromote component_type;
    
        /** Default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    sRGB2RGBFunctor()
    : max_(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    sRGB2RGBFunctor(component_type max)
    : max_(max)
    {}
    
        /** apply the transformation
        */
    result_type operator()(argument_type const & rgb) const
    {
        return RGBValue<To>(
            NumericTraits<To>::fromRealPromote(detail::inverse_sRGBCorrection(rgb[0], max_)),
            NumericTraits<To>::fromRealPromote(detail::inverse_sRGBCorrection(rgb[1], max_)),
            NumericTraits<To>::fromRealPromote(detail::inverse_sRGBCorrection(rgb[2], max_)));
    }

  private:
    component_type max_;
};

template <>
class sRGB2RGBFunctor<unsigned char, unsigned char>
{    
    unsigned char lut_[256];
        
  public:
  
    typedef RGBValue<unsigned char> value_type;
    
    sRGB2RGBFunctor()
    {
        for(int i=0; i<256; ++i)
        {
            lut_[i] = NumericTraits<unsigned char>::fromRealPromote(detail::inverse_sRGBCorrection(i, 255.0));
        }
    }
    
    sRGB2RGBFunctor(double max)
    {
        for(int i=0; i<256; ++i)
        {
            lut_[i] = NumericTraits<unsigned char>::fromRealPromote(detail::inverse_sRGBCorrection(i, max));
        }
    }
    
    template <class V>
    RGBValue<unsigned char> operator()(V const & rgb) const
    {
        return RGBValue<unsigned char>(lut_[rgb[0]], lut_[rgb[1]], lut_[rgb[2]]);
    }
};

template <class From, class To>
class FunctorTraits<sRGB2RGBFunctor<From, To> >
: public FunctorTraitsBase<sRGB2RGBFunctor<From, To> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert linear (raw) RGB into standardized tri-stimulus XYZ.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    According to ITU-R Recommendation BT.709, the functor realizes the transformation
    
    \f[
        \begin{array}{rcl}
        X & = & 0.412453\enspace R / R_{max} + 0.357580\enspace G / G_{max} + 0.180423\enspace B / B_{max}\\
        Y & = & 0.212671\enspace R / R_{max} + 0.715160\enspace G / G_{max} + 0.072169\enspace B / B_{max} \\
        Z & = & 0.019334\enspace R / R_{max} + 0.119193\enspace G / G_{max} + 0.950227\enspace B / B_{max}
        \end{array}
    \f]
    
    By default, \f$ R_{max} = G_{max} = B_{max} = 255 \f$. This default can be overridden
    in the constructor. X, Y, and Z are always positive and reach their maximum for white. 
    The white point is obtained by transforming RGB(255, 255, 255). It corresponds to the 
    D65 illuminant. Y represents the <em>luminance</em> ("brightness") of the color.

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class RGB2XYZFunctor
{
  public:
  
        /** the result's component type
        */
    typedef typename NumericTraits<T>::RealPromote component_type;

        /** the functor's argument type
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef TinyVector<component_type, 3> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef TinyVector<component_type, 3> value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    RGB2XYZFunctor()
    : max_(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    RGB2XYZFunctor(component_type max)
    : max_(max)
    {}
    
        /** apply the transformation
        */
    result_type operator()(argument_type const & rgb) const
    {
        component_type red = rgb[0] / max_;
        component_type green = rgb[1] / max_;
        component_type blue = rgb[2] / max_;
        result_type result;
        result[0] = 0.412453*red + 0.357580*green + 0.180423*blue;
        result[1] = 0.212671*red + 0.715160*green + 0.072169*blue;
        result[2] = 0.019334*red + 0.119193*green + 0.950227*blue;
        return result;
    }

  private:
    component_type max_;
};

template <class T>
class FunctorTraits<RGB2XYZFunctor<T> >
: public FunctorTraitsBase<RGB2XYZFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert non-linear (gamma corrected) R'G'B' into standardized tri-stimulus XYZ.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the transformation
    
    \f[
        R'G'B' \Rightarrow RGB \Rightarrow XYZ
    \f]
    
    See vigra::RGBPrime2RGBFunctor and vigra::RGB2XYZFunctor for a description of the two 
    steps.

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class RGBPrime2XYZFunctor
{
  public:
  
        /** the result's component type
        */
    typedef typename NumericTraits<T>::RealPromote component_type;

        /** the functor's argument type
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef TinyVector<component_type, 3> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef TinyVector<component_type, 3> value_type;
    
        /** default constructor
            The maximum value for each RGB component defaults to 255.
        */
    RGBPrime2XYZFunctor()
    : max_(255.0), gamma_(1.0/ 0.45)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    RGBPrime2XYZFunctor(component_type max)
    : max_(max), gamma_(1.0/ 0.45)
    {}
    
        /** apply the transformation
        */
    result_type operator()(argument_type const & rgb) const
    {
        component_type red = detail::gammaCorrection(rgb[0]/max_, gamma_);
        component_type green = detail::gammaCorrection(rgb[1]/max_, gamma_);
        component_type blue = detail::gammaCorrection(rgb[2]/max_, gamma_);
        result_type result;
        result[0] = 0.412453*red + 0.357580*green + 0.180423*blue;
        result[1] = 0.212671*red + 0.715160*green + 0.072169*blue;
        result[2] = 0.019334*red + 0.119193*green + 0.950227*blue;
        return result;
    }

  private:
    component_type max_, gamma_;
};

template <class T>
class FunctorTraits<RGBPrime2XYZFunctor<T> >
: public FunctorTraitsBase<RGBPrime2XYZFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert standardized tri-stimulus XYZ into linear (raw) RGB.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    According to ITU-R Recommendation BT.709, the functor realizes the transformation
    
    \f[
        \begin{array}{rcl}
        R & = & R_{max} (3.2404813432\enspace X - 1.5371515163\enspace Y - 0.4985363262\enspace Z) \\
        G & = & G_{max} (-0.9692549500\enspace X + 1.8759900015\enspace Y + 0.0415559266\enspace Z) \\
        B & = & B_{max} (0.0556466391\enspace X - 0.2040413384\enspace Y + 1.0573110696\enspace Z)
        \end{array}
    \f]
    
    By default, \f$ R_{max} = G_{max} = B_{max} = 255 \f$. This default can be overridden
    in the constructor. This is the inverse transform of vigra::RGB2XYZFunctor.

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class XYZ2RGBFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
        /** the functor's argument type. (Actually, the argument type
            is more general: <TT>V</TT> with arbitrary
            <TT>V</TT>. But this cannot be expressed in a typedef.)
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef RGBValue<T> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef RGBValue<T> value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    XYZ2RGBFunctor()
    : max_(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    XYZ2RGBFunctor(component_type max)
    : max_(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & xyz) const
    {
        component_type red =    3.2404813432*xyz[0] - 1.5371515163*xyz[1] - 0.4985363262*xyz[2];
        component_type green = -0.9692549500*xyz[0] + 1.8759900015*xyz[1] + 0.0415559266*xyz[2];
        component_type blue =   0.0556466391*xyz[0] - 0.2040413384*xyz[1] + 1.0573110696*xyz[2];
        return value_type(NumericTraits<T>::fromRealPromote(red * max_),
                          NumericTraits<T>::fromRealPromote(green * max_),
                          NumericTraits<T>::fromRealPromote(blue * max_));
    }
};

template <class T>
class FunctorTraits<XYZ2RGBFunctor<T> >
: public FunctorTraitsBase<XYZ2RGBFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert standardized tri-stimulus XYZ into non-linear (gamma corrected) R'G'B'.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the transformation
    
    \f[
        XYZ \Rightarrow RGB \Rightarrow R'G'B'
    \f]
    
    See vigra::XYZ2RGBFunctor and vigra::RGB2RGBPrimeFunctor for a description of the two 
    steps.

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class XYZ2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_, gamma_;
    
  public:
  
  public:
        /** the functor's argument type. (actually, the argument type
            can be any vector type with the same interface. 
            But this cannot be expressed in a typedef.)
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef RGBValue<T> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef RGBValue<T> value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    XYZ2RGBPrimeFunctor()
    : max_(255.0), gamma_(0.45)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    XYZ2RGBPrimeFunctor(component_type max)
    : max_(max), gamma_(0.45)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & xyz) const
    {
        component_type red =    3.2404813432*xyz[0] - 1.5371515163*xyz[1] - 0.4985363262*xyz[2];
        component_type green = -0.9692549500*xyz[0] + 1.8759900015*xyz[1] + 0.0415559266*xyz[2];
        component_type blue =   0.0556466391*xyz[0] - 0.2040413384*xyz[1] + 1.0573110696*xyz[2];
        return value_type(NumericTraits<T>::fromRealPromote(detail::gammaCorrection(red, gamma_) * max_),
                          NumericTraits<T>::fromRealPromote(detail::gammaCorrection(green, gamma_) * max_),
                          NumericTraits<T>::fromRealPromote(detail::gammaCorrection(blue, gamma_) * max_));
    }
};

template <class T>
class FunctorTraits<XYZ2RGBPrimeFunctor<T> >
: public FunctorTraitsBase<XYZ2RGBPrimeFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert standardized tri-stimulus XYZ into perceptual uniform CIE L*u*v*.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the transformation
    
    \f[
        \begin{array}{rcl}
        L^{*} & = & 116 \left( \frac{Y}{Y_n} \right)^\frac{1}{3}-16 \quad \mbox{if} \quad 0.008856 < \frac{Y}{Y_n}\\
        & & \\
        L^{*} & = & 903.3\enspace \frac{Y}{Y_n} \quad \mbox{otherwise} \\
        & & \\
        
        u' & = & \frac{4 X}{X+15 Y + 3 Z}, \quad 
             v' = \frac{9 Y}{X+15 Y + 3 Z}\\
        & & \\
        u^{*} & = & 13 L^{*} (u' - u_n'), \quad v^{*} = 13 L^{*} (v' - v_n')
        \end{array}
    \f]
    
    where \f$(X_n, Y_n, Z_n)\f$ is the reference white point, and 
    \f$u_n' = 0.197839, v_n'=0.468342\f$ are the quantities \f$u', v'\f$ calculated for this
    point. \f$L^{*}\f$ represents the
    <em>lighness</em> ("brightness") of the color, and \f$u^{*}, v^{*}\f$ code the 
    chromaticity.

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class XYZ2LuvFunctor
{
  public:
  
        /** the result's component type
        */
    typedef typename NumericTraits<T>::RealPromote component_type;

        /** the functor's argument type
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef TinyVector<component_type, 3> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef TinyVector<component_type, 3> value_type;
    
    XYZ2LuvFunctor()
    : gamma_(1.0/3.0)
    {}
    
    template <class V>
    result_type operator()(V const & xyz) const
    {
        result_type result;
        if(xyz[1] == NumericTraits<T>::zero())
        {
            result[0] = NumericTraits<component_type>::zero();
            result[1] = NumericTraits<component_type>::zero();
            result[2] = NumericTraits<component_type>::zero();
        }
        else
        {
            component_type L = xyz[1] < 0.008856 ?
                                  903.3 * xyz[1] :
                                  116.0 * VIGRA_CSTD::pow((double)xyz[1], gamma_) - 16.0;
            component_type denom = xyz[0] + 15.0*xyz[1] + 3.0*xyz[2];
            component_type uprime = 4.0 * xyz[0] / denom;
            component_type vprime = 9.0 * xyz[1] / denom;
            result[0] = L;
            result[1] = 13.0*L*(uprime - 0.197839);
            result[2] = 13.0*L*(vprime - 0.468342);
        }
        return result;
    }

  private:
    double gamma_;
};

template <class T>
class FunctorTraits<XYZ2LuvFunctor<T> >
: public FunctorTraitsBase<XYZ2LuvFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert perceptual uniform CIE L*u*v* into standardized tri-stimulus XYZ.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the inverse of the transformation described in vigra::XYZ2LuvFunctor

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class Luv2XYZFunctor
{
  public:
  
        /** the result's component type
        */
    typedef typename NumericTraits<T>::RealPromote component_type;

        /** the functor's argument type
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef TinyVector<component_type, 3> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef TinyVector<component_type, 3> value_type;
    
    Luv2XYZFunctor()
    : gamma_(3.0)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & luv) const
    {
        result_type result;
        if(luv[0] == NumericTraits<T>::zero())
        {
            result[0] = NumericTraits<component_type>::zero();
            result[1] = NumericTraits<component_type>::zero();
            result[2] = NumericTraits<component_type>::zero();
        }
        else
        {
            component_type uprime = luv[1] / 13.0 / luv[0] + 0.197839;
            component_type vprime = luv[2] / 13.0 / luv[0] + 0.468342;

            result[1] = luv[0] < 8.0 ?
                                  luv[0] / 903.3 :
                                  VIGRA_CSTD::pow((luv[0] + 16.0) / 116.0, gamma_);
            result[0] = 9.0*uprime*result[1] / 4.0 / vprime;
            result[2] = ((9.0 / vprime - 15.0)*result[1] - result[0])/ 3.0;
        }
        return result;
    }

  private:
    double gamma_;
};

template <class T>
class FunctorTraits<Luv2XYZFunctor<T> >
: public FunctorTraitsBase<Luv2XYZFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert standardized tri-stimulus XYZ into perceptual uniform CIE L*a*b*.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the transformation
    
    \f[
        \begin{array}{rcl}
        L^{*} & = & 116 \left( \frac{Y}{Y_n} \right)^\frac{1}{3}-16 \quad \mbox{if} \quad 0.008856 < \frac{Y}{Y_n}\\
        & & \\
        L^{*} & = & 903.3\enspace \frac{Y}{Y_n} \quad \mbox{otherwise} \\
        & & \\
        a^{*} & = & 500 \left[ \left( \frac{X}{X_n} \right)^\frac{1}{3} - \left( \frac{Y}{Y_n} \right)^\frac{1}{3} \right] \\
        & & \\
        b^{*} & = & 200 \left[ \left( \frac{Y}{Y_n} \right)^\frac{1}{3} - \left( \frac{Z}{Z_n} \right)^\frac{1}{3} \right] \\
        \end{array}
    \f]
    
    where \f$(X_n, Y_n, Z_n)\f$ is the reference white point. \f$L^{*}\f$ represents the
    <em>lighness</em> ("brightness") of the color, and \f$a^{*}, b^{*}\f$ code the 
    chromaticity.

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class XYZ2LabFunctor
{
  public:
  
        /** the result's component type
        */
    typedef typename NumericTraits<T>::RealPromote component_type;

        /** the functor's argument type
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef TinyVector<component_type, 3> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef TinyVector<component_type, 3> value_type;
    
    XYZ2LabFunctor()
    : gamma_(1.0/3.0)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & xyz) const
    {
        component_type xgamma = VIGRA_CSTD::pow(xyz[0] / 0.950456, gamma_);
        component_type ygamma = VIGRA_CSTD::pow((double)xyz[1], gamma_);
        component_type zgamma = VIGRA_CSTD::pow(xyz[2] / 1.088754, gamma_);
        component_type L = xyz[1] < 0.008856 ?
                              903.3 * xyz[1] :
                              116.0 * ygamma - 16.0;
        result_type result;
        result[0] = L;
        result[1] = 500.0*(xgamma - ygamma);
        result[2] = 200.0*(ygamma - zgamma);
        return result;
    }

  private:
    double gamma_;
};

template <class T>
class FunctorTraits<XYZ2LabFunctor<T> >
: public FunctorTraitsBase<XYZ2LabFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert perceptual uniform CIE L*a*b* into standardized tri-stimulus XYZ.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the inverse of the transformation described in vigra::XYZ2LabFunctor

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class Lab2XYZFunctor
{
  public:
  
        /** the result's component type
        */
    typedef typename NumericTraits<T>::RealPromote component_type;

        /** the functor's argument type
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef TinyVector<component_type, 3> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef TinyVector<component_type, 3> value_type;
    
        /** the functor's value type
        */
    Lab2XYZFunctor()
    : gamma_(3.0)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & lab) const
    {
        component_type Y = lab[0] < 8.0 ?
                              lab[0] / 903.3 :
                              VIGRA_CSTD::pow((lab[0] + 16.0) / 116.0, gamma_);
        component_type ygamma = VIGRA_CSTD::pow((double)Y, 1.0 / gamma_);
        component_type X = VIGRA_CSTD::pow(lab[1] / 500.0 + ygamma, gamma_) * 0.950456;
        component_type Z = VIGRA_CSTD::pow(-lab[2] / 200.0 + ygamma, gamma_) * 1.088754;
        result_type result;
        result[0] = X;
        result[1] = Y;
        result[2] = Z;
        return result;
    }

  private:
    double gamma_;
};

template <class T>
class FunctorTraits<Lab2XYZFunctor<T> >
: public FunctorTraitsBase<Lab2XYZFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert linear (raw) RGB into perceptual uniform CIE L*u*v*.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the transformation
    
    \f[
        RGB \Rightarrow XYZ \Rightarrow L^*u^*v^*
    \f]
    
    See vigra::RGB2XYZFunctor and vigra::XYZ2LuvFunctor for a description of the two 
    steps. The resulting color components will have the following bounds:
    
    \f[
        \begin{array}{rcl}
        0 \leq & L^* & \leq 100 \\
        -83.077 \leq & u^* & \leq 175.015 \\
        -134.101 \leq & v^* & \leq 107.393
        \end{array}
    \f]

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class RGB2LuvFunctor
{
    /*
    L in [0, 100]
    u in [-83.077, 175.015]
    v in [-134.101, 107.393]
    maximum saturation: 179.04 
    red = [53.2406, 175.015, 37.7522]
    */
  public:
  
        /** the result's component type
        */
    typedef typename NumericTraits<T>::RealPromote component_type;

        /** the functor's argument type
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef typename XYZ2LuvFunctor<component_type>::result_type result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef typename XYZ2LuvFunctor<component_type>::result_type value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    RGB2LuvFunctor()
    : rgb2xyz(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    RGB2LuvFunctor(component_type max)
    : rgb2xyz(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & rgb) const
    {
        return xyz2luv(rgb2xyz(rgb));
    }

  private:
    RGB2XYZFunctor<T> rgb2xyz;
    XYZ2LuvFunctor<component_type> xyz2luv;
};

template <class T>
class FunctorTraits<RGB2LuvFunctor<T> >
: public FunctorTraitsBase<RGB2LuvFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert linear (raw) RGB into perceptual uniform CIE L*a*b*.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the transformation
    
    \f[
        RGB \Rightarrow XYZ \Rightarrow L^*a^*b^*
    \f]
    
    See vigra::RGB2XYZFunctor and vigra::XYZ2LabFunctor for a description of the two 
    steps. The resulting color components will have the following bounds:
    
    \f[
        \begin{array}{rcl}
        0 \leq & L^* & \leq 100 \\
        -86.1813 \leq & u^* & \leq 98.2352 \\
        -107.862 \leq & v^* & \leq 94.4758
        \end{array}
    \f]

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class RGB2LabFunctor
{
    /*
    L in [0, 100]
    a in [-86.1813, 98.2352]
    b in [-107.862, 94.4758] 
    maximum saturation: 133.809
    red = [53.2406, 80.0942, 67.2015]
    */
  public:
  
        /** the result's component type
        */
    typedef typename NumericTraits<T>::RealPromote component_type;

        /** the functor's argument type
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef typename XYZ2LabFunctor<component_type>::result_type result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef typename XYZ2LabFunctor<component_type>::result_type value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    RGB2LabFunctor()
    : rgb2xyz(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    RGB2LabFunctor(component_type max)
    : rgb2xyz(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & rgb) const
    {
        return xyz2lab(rgb2xyz(rgb));
    }

  private:
    RGB2XYZFunctor<T> rgb2xyz;
    XYZ2LabFunctor<component_type> xyz2lab;
};

template <class T>
class FunctorTraits<RGB2LabFunctor<T> >
: public FunctorTraitsBase<RGB2LabFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert perceptual uniform CIE L*u*v* into linear (raw) RGB.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the inverse of the transformation described in vigra::RGB2LuvFunctor

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class Luv2RGBFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    XYZ2RGBFunctor<T> xyz2rgb;
    Luv2XYZFunctor<component_type> luv2xyz;
    
  public:
        /** the functor's argument type. (Actually, the argument type
            can be any vector type with the same interface. 
            But this cannot be expressed in a typedef.)
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef typename XYZ2RGBFunctor<T>::result_type result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef typename XYZ2RGBFunctor<T>::result_type value_type;
    
    Luv2RGBFunctor()
    : xyz2rgb(255.0)
    {}
    
    Luv2RGBFunctor(component_type max)
    : xyz2rgb(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & luv) const
    {
        return xyz2rgb(luv2xyz(luv));
    }
};

template <class T>
class FunctorTraits<Luv2RGBFunctor<T> >
: public FunctorTraitsBase<Luv2RGBFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert perceptual uniform CIE L*a*b* into linear (raw) RGB.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the inverse of the transformation described in vigra::RGB2LabFunctor

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class Lab2RGBFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    XYZ2RGBFunctor<T> xyz2rgb;
    Lab2XYZFunctor<component_type> lab2xyz;
    
  public:
  
        /** the functor's argument type. (Actually, the argument type
            can be any vector type with the same interface. 
            But this cannot be expressed in a typedef.)
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef typename XYZ2RGBFunctor<T>::result_type result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef typename XYZ2RGBFunctor<T>::result_type value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    Lab2RGBFunctor()
    : xyz2rgb(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    Lab2RGBFunctor(component_type max)
    : xyz2rgb(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & lab) const
    {
        return xyz2rgb(lab2xyz(lab));
    }
};

template <class T>
class FunctorTraits<Lab2RGBFunctor<T> >
: public FunctorTraitsBase<Lab2RGBFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert non-linear (gamma corrected) R'G'B' into perceptual uniform CIE L*u*v*.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the transformation
    
    \f[
        R'G'B' \Rightarrow RGB \Rightarrow XYZ \Rightarrow L^*u^*v^*
    \f]
    
    See vigra::RGBPrime2RGBFunctor, vigra::RGB2XYZFunctor and vigra::XYZ2LuvFunctor for a description of the three 
    steps. The resulting color components will have the following bounds:
    
    \f[
        \begin{array}{rcl}
        0 \leq & L^* & \leq 100 \\
        -83.077 \leq & u^* & \leq 175.015 \\
        -134.101 \leq & v^* & \leq 107.393
        \end{array}
    \f]

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class RGBPrime2LuvFunctor
{
  public:
  
        /** the result's component type
        */
    typedef typename NumericTraits<T>::RealPromote component_type;

        /** the functor's argument type
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef typename XYZ2LuvFunctor<component_type>::result_type result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef typename XYZ2LuvFunctor<component_type>::result_type value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    RGBPrime2LuvFunctor()
    : rgb2xyz(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    RGBPrime2LuvFunctor(component_type max)
    : rgb2xyz(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & rgb) const
    {
        return xyz2luv(rgb2xyz(rgb));
    }

  private:
    RGBPrime2XYZFunctor<T> rgb2xyz;
    XYZ2LuvFunctor<component_type> xyz2luv;
};

template <class T>
class FunctorTraits<RGBPrime2LuvFunctor<T> >
: public FunctorTraitsBase<RGBPrime2LuvFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert non-linear (gamma corrected) R'G'B' into perceptual uniform CIE L*a*b*.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the transformation
    
    \f[
        R'G'B' \Rightarrow RGB \Rightarrow XYZ \Rightarrow L^*a^*b^*
    \f]
    
    See vigra::RGBPrime2RGBFunctor, vigra::RGB2XYZFunctor and vigra::XYZ2LabFunctor for a description of the three 
    steps. The resulting color components will have the following bounds:
    
    \f[
        \begin{array}{rcl}
        0 \leq & L^* & \leq 100 \\
        -86.1813 \leq & u^* & \leq 98.2352 \\
        -107.862 \leq & v^* & \leq 94.4758
        \end{array}
    \f]

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class RGBPrime2LabFunctor
{
  public:
  
        /** the result's component type
        */
    typedef typename NumericTraits<T>::RealPromote component_type;

        /** the functor's argument type
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef typename XYZ2LabFunctor<component_type>::result_type result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef typename XYZ2LabFunctor<component_type>::result_type value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    RGBPrime2LabFunctor()
    : rgb2xyz(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    RGBPrime2LabFunctor(component_type max)
    : rgb2xyz(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & rgb) const
    {
        return xyz2lab(rgb2xyz(rgb));
    }

  private:
    RGBPrime2XYZFunctor<T> rgb2xyz;
    XYZ2LabFunctor<component_type> xyz2lab;
};

template <class T>
class FunctorTraits<RGBPrime2LabFunctor<T> >
: public FunctorTraitsBase<RGBPrime2LabFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert perceptual uniform CIE L*u*v* into non-linear (gamma corrected) R'G'B'.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the inverse of the transformation described in vigra::RGBPrime2LuvFunctor

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class Luv2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    XYZ2RGBPrimeFunctor<T> xyz2rgb;
    Luv2XYZFunctor<component_type> luv2xyz;
    
  public:
  
        /** the functor's argument type. (Actually, the argument type
            can be any vector type with the same interface. 
            But this cannot be expressed in a typedef.)
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef typename XYZ2RGBFunctor<T>::result_type result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef typename XYZ2RGBFunctor<T>::result_type value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    Luv2RGBPrimeFunctor()
    : xyz2rgb(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    Luv2RGBPrimeFunctor(component_type max)
    : xyz2rgb(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & luv) const
    {
        return xyz2rgb(luv2xyz(luv));
    }
};

template <class T>
class FunctorTraits<Luv2RGBPrimeFunctor<T> >
: public FunctorTraitsBase<Luv2RGBPrimeFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert perceptual uniform CIE L*a*b* into non-linear (gamma corrected) R'G'B'.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the inverse of the transformation described in vigra::RGBPrime2LabFunctor

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class Lab2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    XYZ2RGBPrimeFunctor<T> xyz2rgb;
    Lab2XYZFunctor<component_type> lab2xyz;
    
  public:
  
        /** the functor's argument type. (Actually, the argument type
            can be any vector type with the same interface. 
            But this cannot be expressed in a typedef.)
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef typename XYZ2RGBFunctor<T>::result_type result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef typename XYZ2RGBFunctor<T>::result_type value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    Lab2RGBPrimeFunctor()
    : xyz2rgb(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    Lab2RGBPrimeFunctor(component_type max)
    : xyz2rgb(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & lab) const
    {
        return xyz2rgb(lab2xyz(lab));
    }
};

template <class T>
class FunctorTraits<Lab2RGBPrimeFunctor<T> >
: public FunctorTraitsBase<Lab2RGBPrimeFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert non-linear (gamma corrected) R'G'B' into Y'PbPr color difference components.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    According to ITU-R Recommendation BT.601, the functor realizes the transformation
    
    \f[
        \begin{array}{rcl}
        Y' & = & 0.299\enspace R / R_{max} + 0.587\enspace G / G_{max} + 0.114\enspace B / B_{max}\\
        Pb & = & -0.1687358916\enspace R / R_{max} + 0.3312641084\enspace G / G_{max} + 0.5\enspace B / B_{max} \\
        Pr & = & 0.5\enspace R / R_{max} + 0.4186875892\enspace G / G_{max} + 0.0813124108\enspace B / B_{max}
        \end{array}
    \f]
    
    By default, \f$ R_{max} = G_{max} = B_{max} = 255 \f$. This default can be overridden
    in the constructor. Y' represents the <em>luminance</em> ("brightness") of the color, and
    Pb and Pr are the blue (B'-Y') and red (R'-Y') color difference components. 
    The transformation is scaled so that the following bounds apply:
    
    \f[
        \begin{array}{rcl}
        0 \leq & Y' & \leq 1 \\
        -0.5 \leq & Pb & \leq 0.5 \\
        -0.5 \leq & Pr & \leq 0.5
        \end{array}
    \f]

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class RGBPrime2YPrimePbPrFunctor
{
    /*
    Y in [0, 1]
    Pb in [-0.5, 0.5]
    Pr in [-0.5, 0.5]
    maximum saturation: 0.533887
    red = [0.299, -0.168736, 0.5]
    */
  public:
  
        /** the result's component type
        */
    typedef typename NumericTraits<T>::RealPromote component_type;

        /** the functor's argument type
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef TinyVector<component_type, 3> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef TinyVector<component_type, 3> value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    RGBPrime2YPrimePbPrFunctor()
    : max_(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    RGBPrime2YPrimePbPrFunctor(component_type max)
    : max_(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & rgb) const
    {
        component_type red = rgb[0] / max_;
        component_type green = rgb[1] / max_;
        component_type blue = rgb[2] / max_;
        
        result_type result;
        result[0] = 0.299*red + 0.587*green + 0.114*blue;
        result[1] = -0.1687358916*red - 0.3312641084*green + 0.5*blue;
        result[2] = 0.5*red - 0.4186875892*green - 0.0813124108*blue;
        return result;
    }

  private:
    component_type max_;
};

template <class T>
class FunctorTraits<RGBPrime2YPrimePbPrFunctor<T> >
: public FunctorTraitsBase<RGBPrime2YPrimePbPrFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert Y'PbPr color difference components into non-linear (gamma corrected) R'G'B'.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the inverse of the transformation described in vigra::RGBPrime2YPrimePbPrFunctor

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class YPrimePbPr2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
        /** the functor's argument type. (Actually, the argument type
            can be any vector type with the same interface. 
            But this cannot be expressed in a typedef.)
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef RGBValue<T> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef RGBValue<T> value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    YPrimePbPr2RGBPrimeFunctor()
    : max_(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    YPrimePbPr2RGBPrimeFunctor(component_type max)
    : max_(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & ypbpr) const
    {
        component_type nred =   ypbpr[0] + 1.402*ypbpr[2];
        component_type ngreen = ypbpr[0] - 0.3441362862*ypbpr[1] - 0.7141362862*ypbpr[2];
        component_type nblue =  ypbpr[0] + 1.772*ypbpr[1];
        return result_type(NumericTraits<T>::fromRealPromote(nred * max_),
                           NumericTraits<T>::fromRealPromote(ngreen * max_),
                           NumericTraits<T>::fromRealPromote(nblue * max_));
    }
};

template <class T>
class FunctorTraits<YPrimePbPr2RGBPrimeFunctor<T> >
: public FunctorTraitsBase<YPrimePbPr2RGBPrimeFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert non-linear (gamma corrected) R'G'B' into Y'IQ components.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    According to the PAL analog videa standard, the functor realizes the transformation
    
    \f[
        \begin{array}{rcl}
        Y' & = & 0.299\enspace R / R_{max} + 0.587\enspace G / G_{max} + 0.114\enspace B / B_{max}\\
        I & = & 0.596\enspace R / R_{max} - 0.274\enspace G / G_{max} - 0.322\enspace B / B_{max} \\
        Q & = & 0.212\enspace R / R_{max} - 0.523\enspace G / G_{max} + 0.311\enspace B / B_{max}
        \end{array}
    \f]
    
    By default, \f$ R_{max} = G_{max} = B_{max} = 255 \f$. This default can be overridden
    in the constructor. Y' represents the <em>luminance</em> ("brightness") of the color. 
    The transformation is scaled so that the following bounds apply:
    
    \f[
        \begin{array}{rcl}
        0 \leq & Y' & \leq 1 \\
        -0.596 \leq & I & \leq 0.596 \\
        -0.523 \leq & Q & \leq 0.523
        \end{array}
    \f]

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class RGBPrime2YPrimeIQFunctor
{
    /*
    Y in [0, 1]
    I in [-0.596, 0.596]
    Q in [-0.523, 0.523]
    maximum saturation: 0.632582
    red = [0.299, 0.596, 0.212]
    */
  public:
  
        /** the result's component type
        */
    typedef typename NumericTraits<T>::RealPromote component_type;

        /** the functor's argument type
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef TinyVector<component_type, 3> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef TinyVector<component_type, 3> value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    RGBPrime2YPrimeIQFunctor()
    : max_(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    RGBPrime2YPrimeIQFunctor(component_type max)
    : max_(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & rgb) const
    {
        component_type red = rgb[0] / max_;
        component_type green = rgb[1] / max_;
        component_type blue = rgb[2] / max_;
        
        result_type result;
        result[0] = 0.299*red + 0.587*green + 0.114*blue;
        result[1] = 0.596*red - 0.274*green - 0.322*blue;
        result[2] = 0.212*red - 0.523*green + 0.311*blue;
        return result;
    }

  private:
    component_type max_;
};

template <class T>
class FunctorTraits<RGBPrime2YPrimeIQFunctor<T> >
: public FunctorTraitsBase<RGBPrime2YPrimeIQFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert Y'IQ color components into non-linear (gamma corrected) R'G'B'.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the inverse of the transformation described in vigra::RGBPrime2YPrimeIQFunctor

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class YPrimeIQ2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
        /** the functor's argument type. (Actually, the argument type
            can be any vector type with the same interface. 
            But this cannot be expressed in a typedef.)
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef RGBValue<T> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef RGBValue<T> value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    YPrimeIQ2RGBPrimeFunctor()
    : max_(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    YPrimeIQ2RGBPrimeFunctor(component_type max)
    : max_(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & yiq) const
    {
        component_type nred =   yiq[0] + 0.9548892043*yiq[1] + 0.6221039350*yiq[2];
        component_type ngreen = yiq[0] - 0.2713547827*yiq[1] - 0.6475120259*yiq[2];
        component_type nblue =  yiq[0] - 1.1072510054*yiq[1] + 1.7024603738*yiq[2];
        return result_type(NumericTraits<T>::fromRealPromote(nred * max_),
                           NumericTraits<T>::fromRealPromote(ngreen * max_),
                           NumericTraits<T>::fromRealPromote(nblue * max_));
    }
};

template <class T>
class FunctorTraits<YPrimeIQ2RGBPrimeFunctor<T> >
: public FunctorTraitsBase<YPrimeIQ2RGBPrimeFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert non-linear (gamma corrected) R'G'B' into Y'UV components.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    According to the NTSC analog videa standard, the functor realizes the transformation
    
    \f[
        \begin{array}{rcl}
        Y' & = & 0.299\enspace R / R_{max} + 0.587\enspace G / G_{max} + 0.114\enspace B / B_{max}\\
        U & = & -0.147\enspace R / R_{max} - 0.289\enspace G / G_{max} + 0.436\enspace B / B_{max} \\
        V & = & 0.615\enspace R / R_{max} - 0.515\enspace G / G_{max} - 0.100\enspace B / B_{max}
        \end{array}
    \f]
    
    By default, \f$ R_{max} = G_{max} = B_{max} = 255 \f$. This default can be overridden
    in the constructor. Y' represents the <em>luminance</em> ("brightness") of the color. 
    The transformation is scaled so that the following bounds apply:
    
    \f[
        \begin{array}{rcl}
        0 \leq & Y' & \leq 1 \\
        -0.436 \leq & U & \leq 0.436 \\
        -0.615 \leq & V & \leq 0.615
        \end{array}
    \f]

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class RGBPrime2YPrimeUVFunctor
{
    /*
    Y in [0, 1]
    U in [-0.436, 0.436]
    V in [-0.615, 0.615]
    maximum saturation: 0.632324
    red = [0.299, -0.147, 0.615]
    */
  public:
  
        /** the result's component type
        */
    typedef typename NumericTraits<T>::RealPromote component_type;

        /** the functor's argument type
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef TinyVector<component_type, 3> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef TinyVector<component_type, 3> value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    RGBPrime2YPrimeUVFunctor()
    : max_(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    RGBPrime2YPrimeUVFunctor(component_type max)
    : max_(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & rgb) const
    {
        component_type red = rgb[0] / max_;
        component_type green = rgb[1] / max_;
        component_type blue = rgb[2] / max_;
        
        result_type result;
        result[0] = 0.299*red + 0.587*green + 0.114*blue;
        result[1] = -0.1471376975*red - 0.2888623025*green + 0.436*blue;
        result[2] = 0.6149122807*red - 0.5149122807*green - 0.100*blue;
        return result;
    }

  private:
    component_type max_;
};

template <class T>
class FunctorTraits<RGBPrime2YPrimeUVFunctor<T> >
: public FunctorTraitsBase<RGBPrime2YPrimeUVFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert Y'UV color components into non-linear (gamma corrected) R'G'B'.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the inverse of the transformation described in vigra::RGBPrime2YPrimeUVFunctor

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class YPrimeUV2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
        /** the functor's argument type. (Actually, the argument type
            can be any vector type with the same interface. 
            But this cannot be expressed in a typedef.)
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef RGBValue<T> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef RGBValue<T> value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    YPrimeUV2RGBPrimeFunctor()
    : max_(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    YPrimeUV2RGBPrimeFunctor(component_type max)
    : max_(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & yuv) const
    {
        component_type nred =   yuv[0] + 1.140*yuv[2];
        component_type ngreen = yuv[0] - 0.3946517044*yuv[1] - 0.580681431*yuv[2];
        component_type nblue =  yuv[0] + 2.0321100920*yuv[1];
        return result_type(NumericTraits<T>::fromRealPromote(nred * max_),
                           NumericTraits<T>::fromRealPromote(ngreen * max_),
                           NumericTraits<T>::fromRealPromote(nblue * max_));
    }
};

template <class T>
class FunctorTraits<YPrimeUV2RGBPrimeFunctor<T> >
: public FunctorTraitsBase<YPrimeUV2RGBPrimeFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert non-linear (gamma corrected) R'G'B' into Y'CbCr color difference components.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    This functor basically applies the same transformation as vigra::RGBPrime2YPrimePbPrFunctor
    but the color components are scaled so that they can be coded as 8 bit intergers with
    minimal loss of information:
    
    \f[
        \begin{array}{rcl}
        16\leq & Y' & \leq 235 \\
        16 \leq & Cb & \leq 240 \\
        16 \leq & Cr & \leq 240
        \end{array}
    \f]

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class RGBPrime2YPrimeCbCrFunctor
{
    /*
    Y in [16, 235]
    Cb in [16, 240]
    Cr in [16, 240]
    maximum saturation: 119.591
    red = [81.481, 90.203, 240]
    */
  public:
  
        /** the result's component type
        */
    typedef typename NumericTraits<T>::RealPromote component_type;

        /** the functor's argument type
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef TinyVector<component_type, 3> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef TinyVector<component_type, 3> value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    RGBPrime2YPrimeCbCrFunctor()
    : max_(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    RGBPrime2YPrimeCbCrFunctor(component_type max)
    : max_(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & rgb) const
    {
        component_type red = rgb[0] / max_;
        component_type green = rgb[1] / max_;
        component_type blue = rgb[2] / max_;
        
        result_type result;
        result[0] = 16.0 + 65.481*red + 128.553*green + 24.966*blue;
        result[1] = 128.0 - 37.79683972*red - 74.20316028*green + 112.0*blue;
        result[2] = 128.0 + 112.0*red - 93.78601998*green - 18.21398002*blue;
        return result;
    }

  private:
    component_type max_;
};

template <class T>
class FunctorTraits<RGBPrime2YPrimeCbCrFunctor<T> >
: public FunctorTraitsBase<RGBPrime2YPrimeCbCrFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/** \brief Convert Y'CbCr color difference components into non-linear (gamma corrected) R'G'B'.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    The functor realizes the inverse of the transformation described in vigra::RGBPrime2YPrimeCbCrFunctor

    <b> Traits defined:</b>
    
    <tt>FunctorTraits::isUnaryFunctor</tt> is true (<tt>VigraTrueType<tt>)
*/
template <class T>
class YPrimeCbCr2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
        /** the functor's argument type. (Actually, the argument type
            can be any vector type with the same interface. 
            But this cannot be expressed in a typedef.)
        */
    typedef TinyVector<T, 3> argument_type;
  
        /** the functor's result type
        */
    typedef RGBValue<T> result_type;
  
        /** \deprecated use argument_type and result_type
        */
    typedef RGBValue<T> value_type;
    
        /** default constructor.
            The maximum value for each RGB component defaults to 255.
        */
    YPrimeCbCr2RGBPrimeFunctor()
    : max_(255.0)
    {}
    
        /** constructor
            \arg max - the maximum value for each RGB component
        */
    YPrimeCbCr2RGBPrimeFunctor(component_type max)
    : max_(max)
    {}
    
        /** apply the transformation
        */
    template <class V>
    result_type operator()(V const & ycbcr) const
    {
        component_type y = ycbcr[0] - 16.0;
        component_type cb = ycbcr[1] - 128.0;
        component_type cr = ycbcr[2] - 128.0;
        
        component_type nred =   0.00456621*y + 0.006258928571*cr;
        component_type ngreen = 0.00456621*y - 0.001536322706*cb - 0.003188108420*cr;
        component_type nblue =  0.00456621*y + 0.007910714286*cb;
        return result_type(NumericTraits<T>::fromRealPromote(nred * max_),
                           NumericTraits<T>::fromRealPromote(ngreen * max_),
                           NumericTraits<T>::fromRealPromote(nblue * max_));
    }
};

template <class T>
class FunctorTraits<YPrimeCbCr2RGBPrimeFunctor<T> >
: public FunctorTraitsBase<YPrimeCbCr2RGBPrimeFunctor<T> >
{
  public:
    typedef VigraTrueType isUnaryFunctor;
};

/*
Polar coordinates of standard colors:
=====================================

Lab: black = [320.002, 0, 0]
Luv: black = [347.827, 0, 0]
YPbPr: black = [341.352, 0, 0]
YCbCr: black = [341.352, 0, 0]
YIQ: black = [19.5807, 0, 0]
YUV: black = [346.557, 0, 0]
Lab: red = [1.20391e-05, 0.532406, 0.781353]
Luv: red = [360, 0.532406, 1]
YPbPr: red = [360, 0.299, 0.988419]
YCbCr: red = [360, 0.299, 0.988417]
YIQ: red = [360, 0.299, 1]
YUV: red = [360, 0.299, 1]
Lab: green = [96.0184, 0.877351, 0.895108]
Luv: green = [115.552, 0.877351, 0.758352]
YPbPr: green = [123.001, 0.587, 1]
YCbCr: green = [123.001, 0.587, 0.999996]
YIQ: green = [137.231, 0.587, 0.933362]
YUV: green = [137.257, 0.587, 0.933931]
Lab: blue = [266.287, 0.322957, 0.999997]
Luv: blue = [253.7, 0.322957, 0.729883]
YPbPr: blue = [242.115, 0.114, 0.948831]
YCbCr: blue = [242.115, 0.114, 0.948829]
YIQ: blue = [243.585, 0.114, 0.707681]
YUV: blue = [243.639, 0.114, 0.707424]
Lab: yellow = [62.8531, 0.971395, 0.724189]
Luv: yellow = [73.7, 0.971395, 0.597953]
YPbPr: yellow = [62.1151, 0.886, 0.948831]
YCbCr: yellow = [62.1149, 0.886, 0.948829]
YIQ: yellow = [63.5851, 0.886, 0.707681]
YUV: yellow = [63.6393, 0.886, 0.707424]
Lab: magenta = [288.237, 0.603235, 0.863482]
Luv: magenta = [295.553, 0.603235, 0.767457]
YPbPr: magenta = [303.001, 0.413, 1]
YCbCr: magenta = [303.001, 0.413, 0.999996]
YIQ: magenta = [317.231, 0.413, 0.933362]
YUV: magenta = [317.257, 0.413, 0.933931]
Lab: cyan = [156.378, 0.911133, 0.374577]
Luv: cyan = [180, 0.911133, 0.402694]
YPbPr: cyan = [180, 0.701, 0.988419]
YCbCr: cyan = [180, 0.701, 0.988417]
YIQ: cyan = [180, 0.701, 1]
YUV: cyan = [180, 0.701, 1]
Lab: white = [320.002, 1, 0]
Luv: white = [14.3606, 1, 3.26357e-06]
YPbPr: white = [341.352, 1, 0]
YCbCr: white = [341.352, 1, 0]
YIQ: white = [154.581, 1, 1.24102e-16]
YUV: white = [229.992, 1, 9.81512e-17]

*/

/** \addtogroup PolarColors Polar Color Coordinates
    
    Transform colors from/to a polar representation (hue, brighness, saturation).
    In many situations, this is more inituitive than direct initialization in a 
    particular color space. The polar coordinates are 
    normalized so that a color angle of 0 degrees is always associated with red
    (green is at about 120 degrees, blue at about 240 degrees - exact values differ
    between color spaces). A saturation of 1 is the highest saturation that any RGB color 
    gets after transformation into the respective color space, and saturation 0 corresponds to
    gray. Thus, different color spaces become somewhat comparable.
*/
//@{
/** \brief Init L*a*b* color triple from polar representation.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    <b> Declarations:</b>
    
    \code
    TinyVector<float, 3>
    polar2Lab(double color, double brightness, double saturation);
    
    TinyVector<float, 3>
    polar2Lab(TinyVector<float, 3> const & polar);
    \endcode
    
    \arg color - the color angle in degrees
    \arg brightness - between 0 and 1
    \arg saturation - between 0 and 1
    
    L*a*b* polar coordinates of some important colors:
    
    \code
    black   = [*, 0, 0]    * - arbitrary
    white   = [*, 1, 0]    * - arbitrary
    
    red     = [      0, 0.532406, 0.781353]
    yellow  = [62.8531, 0.971395, 0.724189]
    green   = [96.0184, 0.877351, 0.895108]
    cyan    = [156.378, 0.911133, 0.374577]
    blue    = [266.287, 0.322957, 0.999997]
    magenta = [288.237, 0.603235, 0.863482]
    \endcode
*/
inline TinyVector<float, 3>
polar2Lab(double color, double brightness, double saturation)
{
    double angle = (color+39.9977)/180.0*M_PI;
    double normsat = saturation*133.809;
    
    TinyVector<float, 3> result;
    result[0] = 100.0*brightness;
    result[1] = normsat*VIGRA_CSTD::cos(angle);
    result[2] = normsat*VIGRA_CSTD::sin(angle);
    return result;
}


template <class V>
TinyVector<float, 3>
polar2Lab(V const & polar)
{
    return polar2Lab(polar[0], polar[1], polar[2]);
}

/** \brief Create polar representation form L*a*b*

    <b> Declaration:</b>
    
    \code
    namespace vigra {
        TinyVector<float, 3> lab2Polar(TinyVector<float, 3> const & lab);
    }
    \endcode
    
    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    This realizes the inverse of the transformation described in 
    \link PolarColors#polar2Lab polar2Lab\endlink().
*/
template <class V>
TinyVector<float, 3>
lab2Polar(V const & lab)
{
    TinyVector<float, 3> result;
    result[1] = lab[0]/100.0;
    double angle = (lab[1] == 0.0 && lab[2] == 0.0)
        ? 0.0
        : VIGRA_CSTD::atan2(lab[2], lab[1])/M_PI*180.0-39.9977;
    result[0] = angle < 0.0 ?
                    angle + 360.0 :
                    angle;
    result[2] = VIGRA_CSTD::sqrt(lab[1]*lab[1] + lab[2]*lab[2])/133.809;
    return result;
}

/** \brief Init L*u*v* color triple from polar representation.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    <b> Declarations:</b>
    
    \code
    TinyVector<float, 3>
    polar2Luv(double color, double brightness, double saturation);
    
    TinyVector<float, 3>
    polar2Luv(TinyVector<float, 3> const & polar);
    \endcode
    
    \arg color - the color angle in degrees
    \arg brightness - between 0 and 1
    \arg saturation - between 0 and 1
    
    L*u*v* polar coordinates of some important colors:
    
    \code
    black   = [*, 0, 0]    * - arbitrary
    white   = [*, 1, 0]    * - arbitrary
    
    red     = [      0, 0.532406,        1]
    yellow  = [   73.7, 0.971395, 0.597953]
    green   = [115.552, 0.877351, 0.758352]
    cyan    = [  180.0, 0.911133, 0.402694]
    blue    = [  253.7, 0.322957, 0.729883]
    magenta = [295.553, 0.603235, 0.767457]
    \endcode
*/
inline TinyVector<float, 3>
polar2Luv(double color, double brightness, double saturation)
{
    double angle = (color+12.1727)/180.0*M_PI;
    double normsat = saturation*179.04;
    
    TinyVector<float, 3> result;
    result[0] = 100.0*brightness;
    result[1] = normsat*VIGRA_CSTD::cos(angle);
    result[2] = normsat*VIGRA_CSTD::sin(angle);
    return result;
}

template <class V>
TinyVector<float, 3>
polar2Luv(V const & polar)
{
    return polar2Luv(polar[0], polar[1], polar[2]);
}

/** \brief Create polar representation form L*u*v*

    <b> Declaration:</b>
    
    \code
    namespace vigra {
        TinyVector<float, 3> luv2Polar(TinyVector<float, 3> const & luv);
    }
    \endcode
    
    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    This realizes the inverse of the transformation described in 
    \link PolarColors#polar2Luv polar2Luv\endlink().
*/
template <class V>
TinyVector<float, 3>
luv2Polar(V const & luv)
{
    TinyVector<float, 3> result;
    result[1] = luv[0]/100.0;
    double angle = (luv[1] == 0.0 && luv[2] == 0.0)
        ? 0.0
        : VIGRA_CSTD::atan2(luv[2], luv[1])/M_PI*180.0-12.1727;
    result[0] = angle < 0.0 ?
                    angle + 360.0 :
                    angle;
    result[2] = VIGRA_CSTD::sqrt(luv[1]*luv[1] + luv[2]*luv[2])/179.04;
    return result;
}

/** \brief Init Y'PbPr color triple from polar representation.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    <b> Declarations:</b>
    
    \code
    TinyVector<float, 3>
    polar2YPrimePbPr(double color, double brightness, double saturation);
    
    TinyVector<float, 3>
    polar2YPrimePbPr(TinyVector<float, 3> const & polar);
    \endcode
    
    \arg color - the color angle in degrees
    \arg brightness - between 0 and 1
    \arg saturation - between 0 and 1
    
    Y'PbPr polar coordinates of some important colors:
    
    \code
    black   = [*, 0, 0]    * - arbitrary
    white   = [*, 1, 0]    * - arbitrary
    
    red     = [      0,  0.299, 0.988419]
    yellow  = [62.1151,  0.886, 0.948831]
    green   = [123.001,  0.587,        1]
    cyan    = [  180.0,  0.701, 0.988419]
    blue    = [242.115,  0.114, 0.948831]
    magenta = [303.001,  0.413,        1]
    \endcode
*/
inline TinyVector<float, 3>
polar2YPrimePbPr(double color, double brightness, double saturation)
{
    double angle = (color+18.6481)/180.0*M_PI;
    double normsat = saturation*0.533887;
    
    TinyVector<float, 3> result;
    result[0] = brightness;
    result[1] = -normsat*VIGRA_CSTD::sin(angle);
    result[2] = normsat*VIGRA_CSTD::cos(angle);
    return result;
}

template <class V>
TinyVector<float, 3>
polar2YPrimePbPr(V const & polar)
{
    return polar2YPrimePbPr(polar[0], polar[1], polar[2]);
}

/** \brief Create polar representation form Y'PbPr

    <b> Declaration:</b>
    
    \code
    namespace vigra {
        TinyVector<float, 3> yPrimePbPr2Polar(TinyVector<float, 3> const & ypbpr);
    }
    \endcode
    
    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    This realizes the inverse of the transformation described in 
    \link PolarColors#polar2YPrimePbPr polar2YPrimePbPr\endlink().
*/
template <class V>
TinyVector<float, 3>
yPrimePbPr2Polar(V const & ypbpr)
{
    TinyVector<float, 3> result;
    result[1] = ypbpr[0];
    double angle = (ypbpr[1] == 0.0 && ypbpr[2] == 0.0)
        ? 0.0
        : VIGRA_CSTD::atan2(-ypbpr[1], ypbpr[2])/M_PI*180.0-18.6481;
    result[0] = angle < 0.0 ?
                    angle + 360.0 :
                    angle;
    result[2] = VIGRA_CSTD::sqrt(ypbpr[1]*ypbpr[1] + ypbpr[2]*ypbpr[2])/0.533887;
    return result;
}

/** \brief Init Y'CbCr color triple from polar representation.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    <b> Declarations:</b>
    
    \code
    TinyVector<float, 3>
    polar2YPrimeCbCr(double color, double brightness, double saturation);
    
    TinyVector<float, 3>
    polar2YPrimeCbCr(TinyVector<float, 3> const & polar);
    \endcode
    
    \arg color - the color angle in degrees
    \arg brightness - between 0 and 1
    \arg saturation - between 0 and 1
    
    Y'CbCr polar coordinates of some important colors:
    
    \code
    black   = [*, 0, 0]    * - arbitrary
    white   = [*, 1, 0]    * - arbitrary
    
    red     = [      0,  0.299, 0.988419]
    yellow  = [62.1151,  0.886, 0.948831]
    green   = [123.001,  0.587,        1]
    cyan    = [  180.0,  0.701, 0.988419]
    blue    = [242.115,  0.114, 0.948831]
    magenta = [303.001,  0.413,        1]
    \endcode
*/
inline TinyVector<float, 3>
polar2YPrimeCbCr(double color, double brightness, double saturation)
{
    double angle = (color+18.6482)/180.0*M_PI;
    double normsat = saturation*119.591;
    
    TinyVector<float, 3> result;
    result[0] = brightness*219.0 + 16.0;
    result[1] = -normsat*VIGRA_CSTD::sin(angle)+128.0;
    result[2] = normsat*VIGRA_CSTD::cos(angle)+128.0;
    return result;
}

template <class V>
TinyVector<float, 3>
polar2YPrimeCbCr(V const & polar)
{
    return polar2YPrimeCbCr(polar[0], polar[1], polar[2]);
}

/** \brief Create polar representation form Y'CbCr

    <b> Declaration:</b>
    
    \code
    namespace vigra {
        TinyVector<float, 3> yPrimeCbCr2Polar(TinyVector<float, 3> const & ycbcr);
    }
    \endcode
    
    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    This realizes the inverse of the transformation described in 
    \link PolarColors#polar2YPrimeCbCr polar2YPrimeCbCr\endlink().
*/
template <class V>
TinyVector<float, 3>
yPrimeCbCr2Polar(V const & ycbcr)
{
    TinyVector<float, 3> result;
    result[1] = (ycbcr[0]-16.0)/219.0;
    double cb = ycbcr[1]-128.0;
    double cr = ycbcr[2]-128.0;
    double angle = (cb == 0.0 && cr == 0.0)
        ? 0.0
        : VIGRA_CSTD::atan2(-cb, cr)/M_PI*180.0-18.6482;
    result[0] = angle < 0.0 ?
                    angle + 360.0 :
                    angle;
    result[2] = VIGRA_CSTD::sqrt(cb*cb + cr*cr)/119.591;
    return result;
}

/** \brief Init Y'IQ color triple from polar representation.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    <b> Declarations:</b>
    
    \code
    TinyVector<float, 3>
    polar2YPrimeIQ(double color, double brightness, double saturation);
    
    TinyVector<float, 3>
    polar2YPrimeIQ(TinyVector<float, 3> const & polar);
    \endcode
    
    \arg color - the color angle in degrees
    \arg brightness - between 0 and 1
    \arg saturation - between 0 and 1
    
    Y'IQ polar coordinates of some important colors:
    
    \code
    black   = [*, 0, 0]    * - arbitrary
    white   = [*, 1, 0]    * - arbitrary
    
    red     = [      0, 0.299,        1]
    yellow  = [63.5851, 0.886, 0.707681]
    green   = [137.231, 0.587, 0.933362]
    cyan    = [  180.0, 0.701,        1]
    blue    = [243.585, 0.114, 0.707681]
    magenta = [317.231, 0.413, 0.933362]
    \endcode
*/
inline TinyVector<float, 3>
polar2YPrimeIQ(double color, double brightness, double saturation)
{
    double angle = (color-19.5807)/180.0*M_PI;
    double normsat = saturation*0.632582;
    
    TinyVector<float, 3> result;
    result[0] = brightness;
    result[1] = normsat*VIGRA_CSTD::cos(angle);
    result[2] = -normsat*VIGRA_CSTD::sin(angle);
    return result;
}

template <class V>
TinyVector<float, 3>
polar2YPrimeIQ(V const & polar)
{
    return polar2YPrimeIQ(polar[0], polar[1], polar[2]);
}

/** \brief Create polar representation form Y'IQ

    <b> Declaration:</b>
    
    \code
    namespace vigra {
        TinyVector<float, 3> yPrimeIQ2Polar(TinyVector<float, 3> const & yiq);
    }
    \endcode
    
    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    This realizes the inverse of the transformation described in 
    \link PolarColors#polar2YPrimeIQ polar2YPrimeIQ\endlink().
*/
template <class V>
TinyVector<float, 3>
yPrimeIQ2Polar(V const & yiq)
{
    TinyVector<float, 3> result;
    result[1] = yiq[0];
    double angle = (yiq[1] == 0.0 && yiq[2] == 0.0)
        ? 0.0
        : VIGRA_CSTD::atan2(-yiq[2], yiq[1])/M_PI*180.0+19.5807;
    result[0] = angle < 0.0 ?
                    angle + 360.0 :
                    angle;
    result[2] = VIGRA_CSTD::sqrt(yiq[1]*yiq[1] + yiq[2]*yiq[2])/0.632582;
    return result;
}

/** \brief Init Y'UV color triple from polar representation.

    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    <b> Declarations:</b>
    
    \code
    TinyVector<float, 3>
    polar2YPrimeUV(double color, double brightness, double saturation);
    
    TinyVector<float, 3>
    polar2YPrimeUV(TinyVector<float, 3> const & polar);
    \endcode
    
    \arg color - the color angle in degrees
    \arg brightness - between 0 and 1
    \arg saturation - between 0 and 1
    
    Y'UV polar coordinates of some important colors:
    
    \code
    black   = [*, 0, 0]    * - arbitrary
    white   = [*, 1, 0]    * - arbitrary
    
    red     = [      0, 0.299,        1]
    yellow  = [63.5851, 0.886, 0.707681]
    green   = [137.231, 0.587, 0.933362]
    cyan    = [  180.0, 0.701,        1]
    blue    = [243.585, 0.114, 0.707681]
    magenta = [317.231, 0.413, 0.933362]
    \endcode
*/
inline TinyVector<float, 3>
polar2YPrimeUV(double color, double brightness, double saturation)
{
    double angle = (color+13.4569)/180.0*M_PI;
    double normsat = saturation*0.632324;
    
    TinyVector<float, 3> result;
    result[0] = brightness;
    result[1] = -normsat*VIGRA_CSTD::sin(angle);
    result[2] = normsat*VIGRA_CSTD::cos(angle);
    return result;
}

template <class V>
TinyVector<float, 3>
polar2YPrimeUV(V const & polar)
{
    return polar2YPrimeUV(polar[0], polar[1], polar[2]);
}

/** \brief Create polar representation form Y'UV

    <b> Declaration:</b>
    
    \code
    namespace vigra {
        TinyVector<float, 3> yPrimeUV2Polar(TinyVector<float, 3> const & yuv);
    }
    \endcode
    
    <b>\#include</b> "<a href="colorconversions_8hxx-source.html">vigra/colorconversions.hxx</a>"<br>
    Namespace: vigra
    
    This realizes the inverse of the transformation described in 
    \link PolarColors#polar2YPrimeUV polar2YPrimeUV\endlink().
*/
template <class V>
TinyVector<float, 3>
yPrimeUV2Polar(V const & yuv)
{
    TinyVector<float, 3> result;
    result[1] = yuv[0];
    double angle = (yuv[1] == 0.0 && yuv[2] == 0.0)
        ? 0.0
        : VIGRA_CSTD::atan2(-yuv[1], yuv[2])/M_PI*180.0-13.4569;
    result[0] = angle < 0.0 ?
                    angle + 360.0 :
                    angle;
    result[2] = VIGRA_CSTD::sqrt(yuv[1]*yuv[1] + yuv[2]*yuv[2])/0.632324;
    return result;
}

//@}

} // namespace vigra 

#endif /* VIGRA_COLORCONVERSIONS_HXX */
