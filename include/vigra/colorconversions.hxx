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
 
 
#ifndef VIGRA_COLORCONVERSIONS_HXX
#define VIGRA_COLORCONVERSIONS_HXX

#include <cmath>
#include "vigra/rgbvalue.hxx"

namespace vigra {

namespace detail
{

double gammaCorrection(double value, double gamma)
{
    return (value < 0.0) ? 
            -std::pow(-value, gamma) :
            std::pow(value, gamma);
}

double gammaCorrection(double value, double gamma, double norm)
{
    return (value < 0.0) ? 
            -norm*std::pow(-value/norm, gamma) :
            norm*std::pow(value/norm, gamma);
}

} // namespace detail

template <class From, class To>
class RGB2RGBPrimeFunctor
{
    typedef typename NumericTraits<To>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef RGBValue<To> value_type;
    
    RGB2RGBPrimeFunctor(component_type max = 255.0)
    : max_(max)
    {}
    
    value_type operator()(TinyVector<From, 3> const & rgb)
    {
        return RGBValue<To>(
            NumericTraits<To>::fromRealPromote(detail::gammaCorrection(rgb[0], 0.45, max_)),
            NumericTraits<To>::fromRealPromote(detail::gammaCorrection(rgb[1], 0.45, max_)),
            NumericTraits<To>::fromRealPromote(detail::gammaCorrection(rgb[2], 0.45, max_)));
    }
};

template <>
class RGB2RGBPrimeFunctor<unsigned char, unsigned char>
{
    unsigned char lut_[256];
        
  public:
  
    typedef RGBValue<unsigned char> value_type;
    
    RGB2RGBPrimeFunctor(double max = 255.0)
    {
        for(int i=0; i<256; ++i)
        {
            lut_[i] = NumericTraits<unsigned char>::fromRealPromote(detail::gammaCorrection(i, 0.45, max));
        }
    }
    
    RGBValue<unsigned char> operator()(TinyVector<unsigned char, 3> const & rgb)
    {
        return RGBValue<unsigned char>(lut_[rgb[0]], lut_[rgb[1]], lut_[rgb[2]]);
    }
};

template <class From, class To>
class RGBPrime2RGBFunctor
{
    typedef typename NumericTraits<To>::RealPromote component_type;
    
    component_type max_;
    double gamma_;
    
  public:
  
    typedef RGBValue<To> value_type;
    
    RGBPrime2RGBFunctor(component_type max = 255.0)
    : max_(max), gamma_(1.0/0.45)
    {}
    
    value_type operator()(TinyVector<From, 3> const & rgb)
    {
        return RGBValue<To>(
            NumericTraits<To>::fromRealPromote(detail::gammaCorrection(rgb[0], gamma_, max_)),
            NumericTraits<To>::fromRealPromote(detail::gammaCorrection(rgb[1], gamma_, max_)),
            NumericTraits<To>::fromRealPromote(detail::gammaCorrection(rgb[2], gamma_, max_)));
    }
};

template <>
class RGBPrime2RGBFunctor<unsigned char, unsigned char>
{    
    unsigned char lut_[256];
        
  public:
  
    typedef RGBValue<unsigned char> value_type;
    
    RGBPrime2RGBFunctor(double max = 255.0)
    {
        for(int i=0; i<256; ++i)
        {
            lut_[i] = NumericTraits<unsigned char>::fromRealPromote(detail::gammaCorrection(i, 1.0/0.45, max));
        }
    }
    
    RGBValue<unsigned char> operator()(TinyVector<unsigned char, 3> const & rgb)
    {
        return RGBValue<unsigned char>(lut_[rgb[0]], lut_[rgb[1]], lut_[rgb[2]]);
    }
};

template <class T>
class RGB2XYZFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef TinyVector<component_type, 3> value_type;
    
    RGB2XYZFunctor(component_type max = 255.0)
    : max_(max)
    {}
    
    value_type operator()(TinyVector<T, 3> const & rgb)
    {
        component_type red = rgb[0] / max_;
        component_type green = rgb[1] / max_;
        component_type blue = rgb[2] / max_;
        value_type result;
        result[0] = 0.412453*red + 0.357580*green + 0.180423*blue;
        result[1] = 0.212671*red + 0.715160*green + 0.072169*blue;
        result[2] = 0.019334*red + 0.119193*green + 0.950227*blue;
        return result;
    }
};

template <class T>
class RGBPrime2XYZFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_, gamma_;
    
  public:
  
    typedef TinyVector<component_type, 3> value_type;
    
    RGBPrime2XYZFunctor(component_type max = 255.0)
    : max_(max), gamma_(1.0/ 0.45)
    {}
    
    value_type operator()(TinyVector<T, 3> const & rgb)
    {
        component_type red = detail::gammaCorrection(rgb[0]/max_, gamma_);
        component_type green = detail::gammaCorrection(rgb[1]/max_, gamma_);
        component_type blue = detail::gammaCorrection(rgb[2]/max_, gamma_);
        value_type result;
        result[0] = 0.412453*red + 0.357580*green + 0.180423*blue;
        result[1] = 0.212671*red + 0.715160*green + 0.072169*blue;
        result[2] = 0.019334*red + 0.119193*green + 0.950227*blue;
        return result;
    }
};

template <class T>
class XYZ2RGBFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef RGBValue<T> value_type;
    
    XYZ2RGBFunctor(component_type max = 255.0)
    : max_(max)
    {}
    
    template <class V>
    value_type operator()(TinyVector<V, 3> const & xyz)
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
class XYZ2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_, gamma_;
    
  public:
  
    typedef RGBValue<T> value_type;
    
    XYZ2RGBPrimeFunctor(component_type max = 255.0)
    : max_(max), gamma_(0.45)
    {}
    
    template <class V>
    value_type operator()(TinyVector<V, 3> const & xyz)
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
class XYZ2LuvFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type gamma_;
    
  public:
  
    typedef TinyVector<component_type, 3> value_type;
    
    XYZ2LuvFunctor()
    : gamma_(1.0/3.0)
    {}
    
    value_type operator()(TinyVector<T, 3> const & xyz)
    {
        value_type result;
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
                                  116.0 * std::pow(xyz[1], gamma_) - 16.0;
            component_type denom = xyz[0] + 15.0*xyz[1] + 3.0*xyz[2];
            component_type uprime = 4.0 * xyz[0] / denom;
            component_type vprime = 9.0 * xyz[1] / denom;
            result[0] = L;
            result[1] = 13.0*L*(uprime - 0.197839);
            result[2] = 13.0*L*(vprime - 0.468342);
        }
        return result;
    }
};

template <class T>
class Luv2XYZFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type gamma_;
    
  public:
  
    typedef TinyVector<component_type, 3> value_type;
    
    Luv2XYZFunctor()
    : gamma_(3.0)
    {}
    
    value_type operator()(TinyVector<T, 3> const & luv)
    {
        value_type result;
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
                                  std::pow((luv[0] + 16.0) / 116.0, gamma_);
            result[0] = 9.0*uprime*result[1] / 4.0 / vprime;
            result[2] = ((9.0 / vprime - 15.0)*result[1] - result[0])/ 3.0;
        }
        return result;
    }
};

template <class T>
class XYZ2LabFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type gamma_;
    
  public:
  
    typedef TinyVector<component_type, 3> value_type;
    
    XYZ2LabFunctor()
    : gamma_(1.0/3.0)
    {}
    
    value_type operator()(TinyVector<T, 3> const & xyz)
    {
        component_type xgamma = std::pow(xyz[0] / 0.950456, gamma_);
        component_type ygamma = std::pow(xyz[1], gamma_);
        component_type zgamma = std::pow(xyz[2] / 1.088754, gamma_);
        component_type L = xyz[1] < 0.008856 ?
                              903.3 * xyz[1] :
                              116.0 * ygamma - 16.0;
        value_type result;
        result[0] = L;
        result[1] = 500.0*(xgamma - ygamma);
        result[2] = 200.0*(ygamma - zgamma);
        return result;
    }
};

template <class T>
class Lab2XYZFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type gamma_;
    
  public:
  
    typedef TinyVector<component_type, 3> value_type;
    
    Lab2XYZFunctor()
    : gamma_(3.0)
    {}
    
    value_type operator()(TinyVector<T, 3> const & lab)
    {
        component_type Y = lab[0] < 8.0 ?
                              lab[0] / 903.3 :
                              std::pow((lab[0] + 16.0) / 116.0, gamma_);
        component_type ygamma = std::pow(Y, 1.0 / gamma_);
        component_type X = std::pow(lab[1] / 500.0 + ygamma, gamma_) * 0.950456;
        component_type Z = std::pow(-lab[2] / 200.0 + ygamma, gamma_) * 1.088754;
        value_type result;
        result[0] = X;
        result[1] = Y;
        result[2] = Z;
        return result;
    }
};


/*
L in [0, 100]
u in [-83.077, 175.015]
v in [-134.101, 107.393]
maximum saturation: 179.04 
red = [53.2406, 175.015, 37.7522]
*/
template <class T>
class RGB2LuvFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    RGB2XYZFunctor<T> rgb2xyz;
    XYZ2LuvFunctor<component_type> xyz2luv;
    
  public:
  
    typedef typename XYZ2LuvFunctor<component_type>::value_type value_type;
    
    RGB2LuvFunctor(component_type max = 255.0)
    : rgb2xyz(max)
    {}
    
    value_type operator()(TinyVector<T, 3> const & rgb)
    {
        return xyz2luv(rgb2xyz(rgb));
    }
};

/*
L in [0, 100]
a in [-86.1813, 98.2352]
b in [-107.862, 94.4758] 
maximum saturation: 133.809
red = [53.2406, 80.0942, 67.2015]
*/
template <class T>
class RGB2LabFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    RGB2XYZFunctor<T> rgb2xyz;
    XYZ2LabFunctor<component_type> xyz2lab;
    
  public:
  
    typedef typename XYZ2LabFunctor<component_type>::value_type value_type;
    
    RGB2LabFunctor(component_type max = 255.0)
    : rgb2xyz(max)
    {}
    
    value_type operator()(TinyVector<T, 3> const & rgb)
    {
        return xyz2lab(rgb2xyz(rgb));
    }
};

template <class T>
class Luv2RGBFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    XYZ2RGBFunctor<T> xyz2rgb;
    Luv2XYZFunctor<component_type> luv2xyz;
    
  public:
  
    typedef typename XYZ2RGBFunctor<T>::value_type value_type;
    
    Luv2RGBFunctor(component_type max = 255.0)
    : xyz2rgb(max)
    {}
    
    template <class V>
    value_type operator()(TinyVector<V, 3> const & luv)
    {
        return xyz2rgb(luv2xyz(luv));
    }
};

template <class T>
class Lab2RGBFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    XYZ2RGBFunctor<T> xyz2rgb;
    Lab2XYZFunctor<component_type> lab2xyz;
    
  public:
  
    typedef typename XYZ2RGBFunctor<T>::value_type value_type;
    
    Lab2RGBFunctor(component_type max = 255.0)
    : xyz2rgb(max)
    {}
    
    template <class V>
    value_type operator()(TinyVector<V, 3> const & lab)
    {
        return xyz2rgb(lab2xyz(lab));
    }
};

template <class T>
class RGBPrime2LuvFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    RGBPrime2XYZFunctor<T> rgb2xyz;
    XYZ2LuvFunctor<component_type> xyz2luv;
    
  public:
  
    typedef typename XYZ2LuvFunctor<component_type>::value_type value_type;
    
    RGBPrime2LuvFunctor(component_type max = 255.0)
    : rgb2xyz(max)
    {}
    
    value_type operator()(TinyVector<T, 3> const & rgb)
    {
        return xyz2luv(rgb2xyz(rgb));
    }
};

template <class T>
class RGBPrime2LabFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    RGBPrime2XYZFunctor<T> rgb2xyz;
    XYZ2LabFunctor<component_type> xyz2lab;
    
  public:
  
    typedef typename XYZ2LabFunctor<component_type>::value_type value_type;
    
    RGBPrime2LabFunctor(component_type max = 255.0)
    : rgb2xyz(max)
    {}
    
    value_type operator()(TinyVector<T, 3> const & rgb)
    {
        return xyz2lab(rgb2xyz(rgb));
    }
};

template <class T>
class Luv2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    XYZ2RGBPrimeFunctor<T> xyz2rgb;
    Luv2XYZFunctor<component_type> luv2xyz;
    
  public:
  
    typedef typename XYZ2RGBFunctor<T>::value_type value_type;
    
    Luv2RGBPrimeFunctor(component_type max = 255.0)
    : xyz2rgb(max)
    {}
    
    template <class V>
    value_type operator()(TinyVector<V, 3> const & luv)
    {
        return xyz2rgb(luv2xyz(luv));
    }
};

template <class T>
class Lab2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    XYZ2RGBPrimeFunctor<T> xyz2rgb;
    Lab2XYZFunctor<component_type> lab2xyz;
    
  public:
  
    typedef typename XYZ2RGBFunctor<T>::value_type value_type;
    
    Lab2RGBPrimeFunctor(component_type max = 255.0)
    : xyz2rgb(max)
    {}
    
    template <class V>
    value_type operator()(TinyVector<V, 3> const & lab)
    {
        return xyz2rgb(lab2xyz(lab));
    }
};

/*
Y in [0, 1]
Pb in [-0.5, 0.5]
Pr in [-0.5, 0.5]
maximum saturation: 0.533887
red = [0.299, -0.168736, 0.5]
*/
template <class T>
class RGBPrime2YPrimePbPrFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef TinyVector<component_type, 3> value_type;
    
    RGBPrime2YPrimePbPrFunctor(component_type max = 255.0)
    : max_(max)
    {}
    
    value_type operator()(TinyVector<T, 3> const & rgb)
    {
        component_type red = rgb[0] / max_;
        component_type green = rgb[1] / max_;
        component_type blue = rgb[2] / max_;
        
        value_type result;
        result[0] = 0.299*red + 0.587*green + 0.114*blue;
        result[1] = -0.1687358916*red - 0.3312641084*green + 0.5*blue;
        result[2] = 0.5*red - 0.4186875892*green - 0.0813124108*blue;
        return result;
    }
};

template <class T>
class YPrimePbPr2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef RGBValue<T> value_type;
    
    YPrimePbPr2RGBPrimeFunctor(component_type max = 255.0)
    : max_(max)
    {}
    
    template <class V>
    value_type operator()(TinyVector<V, 3> const & ypbpr)
    {
        component_type nred =   ypbpr[0] + 1.402*ypbpr[2];
        component_type ngreen = ypbpr[0] - 0.3441362862*ypbpr[1] - 0.7141362862*ypbpr[2];
        component_type nblue =  ypbpr[0] + 1.772*ypbpr[1];
        return value_type(NumericTraits<T>::fromRealPromote(nred * max_),
                          NumericTraits<T>::fromRealPromote(ngreen * max_),
                          NumericTraits<T>::fromRealPromote(nblue * max_));
    }
};

/*
Y in [0, 1]
I in [-0.596, 0.596]
Q in [-0.523, 0.523]
maximum saturation: 0.632582
red = [0.299, 0.596, 0.212]
*/
template <class T>
class RGBPrime2YPrimeIQFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef TinyVector<component_type, 3> value_type;
    
    RGBPrime2YPrimeIQFunctor(component_type max = 255.0)
    : max_(max)
    {}
    
    value_type operator()(TinyVector<T, 3> const & rgb)
    {
        component_type red = rgb[0] / max_;
        component_type green = rgb[1] / max_;
        component_type blue = rgb[2] / max_;
        
        value_type result;
        result[0] = 0.299*red + 0.587*green + 0.114*blue;
        result[1] = 0.596*red - 0.274*green - 0.322*blue;
        result[2] = 0.212*red - 0.523*green + 0.311*blue;
        return result;
    }
};

template <class T>
class YPrimeIQ2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef RGBValue<T> value_type;
    
    YPrimeIQ2RGBPrimeFunctor(component_type max = 255.0)
    : max_(max)
    {}
    
    template <class V>
    value_type operator()(TinyVector<V, 3> const & yiq)
    {
        component_type nred =   yiq[0] + 0.9548892043*yiq[1] + 0.6221039350*yiq[2];
        component_type ngreen = yiq[0] - 0.2713547827*yiq[1] - 0.6475120259*yiq[2];
        component_type nblue =  yiq[0] - 1.1072510054*yiq[1] + 1.7024603738*yiq[2];
        return value_type(NumericTraits<T>::fromRealPromote(nred * max_),
                          NumericTraits<T>::fromRealPromote(ngreen * max_),
                          NumericTraits<T>::fromRealPromote(nblue * max_));
    }
};

/*
Y in [0, 1]
U in [-0.436, 0.436]
V in [-0.615, 0.615]
maximum saturation: 0.632324
red = [0.299, -0.147, 0.615]
*/
template <class T>
class RGBPrime2YPrimeUVFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef TinyVector<component_type, 3> value_type;
    
    RGBPrime2YPrimeUVFunctor(component_type max = 255.0)
    : max_(max)
    {}
    
    value_type operator()(TinyVector<T, 3> const & rgb)
    {
        component_type red = rgb[0] / max_;
        component_type green = rgb[1] / max_;
        component_type blue = rgb[2] / max_;
        
        value_type result;
        result[0] = 0.299*red + 0.587*green + 0.114*blue;
        result[1] = -0.1471376975*red - 0.2888623025*green + 0.436*blue;
        result[2] = 0.6149122807*red - 0.5149122807*green - 0.100*blue;
        return result;
    }
};

template <class T>
class YPrimeUV2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef RGBValue<T> value_type;
    
    YPrimeUV2RGBPrimeFunctor(component_type max = 255.0)
    : max_(max)
    {}
    
    template <class V>
    value_type operator()(TinyVector<V, 3> const & yuv)
    {
        component_type nred =   yuv[0] + 1.140*yuv[2];
        component_type ngreen = yuv[0] - 0.3946517044*yuv[1] - 0.580681431*yuv[2];
        component_type nblue =  yuv[0] + 2.0321100920*yuv[1];
        return value_type(NumericTraits<T>::fromRealPromote(nred * max_),
                          NumericTraits<T>::fromRealPromote(ngreen * max_),
                          NumericTraits<T>::fromRealPromote(nblue * max_));
    }
};

/*
Y in [16, 235]
Cb in [16, 240]
Cr in [16, 240]
maximum saturation: 119.591
red = [81.481, 90.203, 240]
*/
template <class T>
class RGBPrime2YPrimeCbCrFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef TinyVector<component_type, 3> value_type;
    
    RGBPrime2YPrimeCbCrFunctor(component_type max = 255.0)
    : max_(max)
    {}
    
    value_type operator()(TinyVector<T, 3> const & rgb)
    {
        component_type red = rgb[0] / max_;
        component_type green = rgb[1] / max_;
        component_type blue = rgb[2] / max_;
        
        value_type result;
        result[0] = 16.0 + 65.481*red + 128.553*green + 24.966*blue;
        result[1] = 128.0 - 37.79683972*red - 74.20316028*green + 112.0*blue;
        result[2] = 128.0 + 112.0*red - 93.78601998*green - 18.21398002*blue;
        return result;
    }
};

template <class T>
class YPrimeCbCr2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef RGBValue<T> value_type;
    
    YPrimeCbCr2RGBPrimeFunctor(component_type max = 255.0)
    : max_(max)
    {}
    
    template <class V>
    value_type operator()(TinyVector<V, 3> const & ycbcr)
    {
        component_type y = ycbcr[0] - 16.0;
        component_type cb = ycbcr[1] - 128.0;
        component_type cr = ycbcr[2] - 128.0;
        
        component_type nred =   0.00456621*y + 0.006258928571*cr;
        component_type ngreen = 0.00456621*y - 0.001536322706*cb - 0.003188108420*cr;
        component_type nblue =  0.00456621*y + 0.007910714286*cb;
        return value_type(NumericTraits<T>::fromRealPromote(nred * max_),
                          NumericTraits<T>::fromRealPromote(ngreen * max_),
                          NumericTraits<T>::fromRealPromote(nblue * max_));
    }
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

inline TinyVector<float, 3>
polar2Lab(double color, double brightness, double saturation)
{
    double angle = (color+39.9977)/180.0*M_PI;
    double normsat = saturation*133.809;
    
    TinyVector<float, 3> result;
    result[0] = 100.0*brightness;
    result[1] = normsat*std::cos(angle);
    result[2] = normsat*std::sin(angle);
    return result;
}

inline TinyVector<float, 3>
polar2Lab(TinyVector<float, 3> const & polar)
{
    return polar2Lab(polar[0], polar[1], polar[2]);
}

inline TinyVector<float, 3>
lab2Polar(TinyVector<float, 3> const & lab)
{
    TinyVector<float, 3> result;
    result[1] = lab[0]/100.0;
    double angle = std::atan2(lab[2], lab[1])/M_PI*180.0-39.9977;
    result[0] = angle < 0.0 ?
                    angle + 360.0 :
                    angle;
    result[2] = std::sqrt(lab[1]*lab[1] + lab[2]*lab[2])/133.809;
    return result;
}

inline TinyVector<float, 3>
polar2Luv(double color, double brightness, double saturation)
{
    double angle = (color+12.1727)/180.0*M_PI;
    double normsat = saturation*179.04;
    
    TinyVector<float, 3> result;
    result[0] = 100.0*brightness;
    result[1] = normsat*std::cos(angle);
    result[2] = normsat*std::sin(angle);
    return result;
}

inline TinyVector<float, 3>
polar2Luv(TinyVector<float, 3> const & polar)
{
    return polar2Luv(polar[0], polar[1], polar[2]);
}

inline TinyVector<float, 3>
luv2Polar(TinyVector<float, 3> const & luv)
{
    TinyVector<float, 3> result;
    result[1] = luv[0]/100.0;
    double angle = std::atan2(luv[2], luv[1])/M_PI*180.0-12.1727;
    result[0] = angle < 0.0 ?
                    angle + 360.0 :
                    angle;
    result[2] = std::sqrt(luv[1]*luv[1] + luv[2]*luv[2])/179.04;
    return result;
}

inline TinyVector<float, 3>
polar2YPrimePbPr(double color, double brightness, double saturation)
{
    double angle = (color+18.6481)/180.0*M_PI;
    double normsat = saturation*0.533887;
    
    TinyVector<float, 3> result;
    result[0] = brightness;
    result[1] = -normsat*std::sin(angle);
    result[2] = normsat*std::cos(angle);
    return result;
}

inline TinyVector<float, 3>
polar2YPrimePbPr(TinyVector<float, 3> const & polar)
{
    return polar2YPrimePbPr(polar[0], polar[1], polar[2]);
}

inline TinyVector<float, 3>
yPrimePbPr2Polar(TinyVector<float, 3> const & ypbpr)
{
    TinyVector<float, 3> result;
    result[1] = ypbpr[0];
    double angle = std::atan2(-ypbpr[1], ypbpr[2])/M_PI*180.0-18.6481;
    result[0] = angle < 0.0 ?
                    angle + 360.0 :
                    angle;
    result[2] = std::sqrt(ypbpr[1]*ypbpr[1] + ypbpr[2]*ypbpr[2])/0.533887;
    return result;
}

inline TinyVector<float, 3>
polar2YPrimeCbCr(double color, double brightness, double saturation)
{
    double angle = (color+18.6482)/180.0*M_PI;
    double normsat = saturation*119.591;
    
    TinyVector<float, 3> result;
    result[0] = brightness*219.0 + 16.0;
    result[1] = -normsat*std::sin(angle)+128.0;
    result[2] = normsat*std::cos(angle)+128.0;
    return result;
}

inline TinyVector<float, 3>
polar2YPrimeCbCr(TinyVector<float, 3> const & polar)
{
    return polar2YPrimeCbCr(polar[0], polar[1], polar[2]);
}

inline TinyVector<float, 3>
yPrimeCbCr2Polar(TinyVector<float, 3> const & ycbcr)
{
    TinyVector<float, 3> result;
    result[1] = (ycbcr[0]-16.0)/219.0;
    double cb = ycbcr[1]-128.0;
    double cr = ycbcr[2]-128.0;
    double angle = std::atan2(-cb, cr)/M_PI*180.0-18.6482;
    result[0] = angle < 0.0 ?
                    angle + 360.0 :
                    angle;
    result[2] = std::sqrt(cb*cb + cr*cr)/119.591;
    return result;
}

inline TinyVector<float, 3>
polar2YPrimeIQ(double color, double brightness, double saturation)
{
    double angle = (color-19.5807)/180.0*M_PI;
    double normsat = saturation*0.632582;
    
    TinyVector<float, 3> result;
    result[0] = brightness;
    result[1] = normsat*std::cos(angle);
    result[2] = -normsat*std::sin(angle);
    return result;
}

inline TinyVector<float, 3>
polar2YPrimeIQ(TinyVector<float, 3> const & polar)
{
    return polar2YPrimeIQ(polar[0], polar[1], polar[2]);
}

inline TinyVector<float, 3>
yPrimeIQ2Polar(TinyVector<float, 3> const & yiq)
{
    TinyVector<float, 3> result;
    result[1] = yiq[0];
    double angle = std::atan2(-yiq[2], yiq[1])/M_PI*180.0+19.5807;
    result[0] = angle < 0.0 ?
                    angle + 360.0 :
                    angle;
    result[2] = std::sqrt(yiq[1]*yiq[1] + yiq[2]*yiq[2])/0.632582;
    return result;
}

inline TinyVector<float, 3>
polar2YPrimeUV(double color, double brightness, double saturation)
{
    double angle = (color+13.4429)/180.0*M_PI;
    double normsat = saturation*0.632324;
    
    TinyVector<float, 3> result;
    result[0] = brightness;
    result[1] = -normsat*std::sin(angle);
    result[2] = normsat*std::cos(angle);
    return result;
}

inline TinyVector<float, 3>
polar2YPrimeUV(TinyVector<float, 3> const & polar)
{
    return polar2YPrimeUV(polar[0], polar[1], polar[2]);
}

inline TinyVector<float, 3>
yPrimeUV2Polar(TinyVector<float, 3> const & yuv)
{
    TinyVector<float, 3> result;
    result[1] = yuv[0];
    double angle = std::atan2(-yuv[1], yuv[2])/M_PI*180.0-13.4569;
    result[0] = angle < 0.0 ?
                    angle + 360.0 :
                    angle;
    result[2] = std::sqrt(yuv[1]*yuv[1] + yuv[2]*yuv[2])/0.632324;
    return result;
}

} // namespace vigra 

#endif /* VIGRA_COLORCONVERSIONS_HXX */
