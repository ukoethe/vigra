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

double gammaCorrection(double value, double gamma, double norm)
{
    return norm*std::pow(value/norm, gamma);
}

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
    
    value_type operator()(RGBValue<From> const & rgb)
    {
        return RGBValue<T>(
            NumericTraits<T>::fromRealPromote(gammaCorrection(rgb.red(), 0.45, max_)),
            NumericTraits<T>::fromRealPromote(gammaCorrection(rgb.green(), 0.45, max_)),
            NumericTraits<T>::fromRealPromote(gammaCorrection(rgb.blue(), 0.45, max_)));
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
            lut_[i] = NumericTraits<unsigned char>::fromRealPromote(gammaCorrection(i, 0.45, max)
        }
    }
    
    RGBValue<unsigned char> operator()(RGBValue<unsigned char> const & rgb)
    {
        return RGBValue<unsigned char>(lut_[rgb.red()], lut_[rgb.green()], lut_[rgb.blue()]);
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
    
    value_type operator()(RGBValue<From> const & rgb)
    {
        return RGBValue<T>(
            NumericTraits<T>::fromRealPromote(gammaCorrection(rgb.red(), gamma_, max_)),
            NumericTraits<T>::fromRealPromote(gammaCorrection(rgb.green(), gamma_, max_)),
            NumericTraits<T>::fromRealPromote(gammaCorrection(rgb.blue(), gamma_, max_)));
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
            lut_[i] = NumericTraits<unsigned char>::fromRealPromote(gammaCorrection(i, 1.0/0.45, max)
        }
    }
    
    RGBValue<unsigned char> operator()(RGBValue<unsigned char> const & rgb)
    {
        return RGBValue<unsigned char>(lut_[rgb.red()], lut_[rgb.green()], lut_[rgb.blue()]);
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
    
    value_type operator()(RGBValue<T> const & rgb)
    {
        component_type red = rgb.red() / max_;
        component_type green = rgb.green() / max_;
        component_type blue = rgb.blue() / max_;
        return value_type(0.412453*red + 0.357580*green + 0.180423*blue,
                          0.212671*red + 0.715160*green + 0.072169*blue,
                          0.019334*red + 0.119193*green + 0.950227*blue);
    }
};

template <class T>
class RGBPrime2XYZFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef TinyVector<component_type, 3> value_type;
    
    RGBPrime2XYZFunctor(component_type max = 255.0)
    : max_(max), gamma_(1.0/ 0.45)
    {}
    
    value_type operator()(RGBValue<T> const & rgb)
    {
        component_type red = std::pow(rgb.red() / max_, gamma_);
        component_type green = std::pow(rgb.green() / max_, gamma_);
        component_type blue = std::pow(rgb.blue() / max_, gamma_);
        return value_type(0.412453*red + 0.357580*green + 0.180423*blue,
                          0.212671*red + 0.715160*green + 0.072169*blue,
                          0.019334*red + 0.119193*green + 0.950227*blue);
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
        component_type red =    3.240479*xyz[0] - 1.537150*xyz[1] - 0.498535*xyz[2];
        component_type green = -0.969256*xyz[0] + 1.875992*xyz[1] + 0.041556*xyz[2];
        component_type blue =   0.055648*xyz[0] - 0.204043*xyz[1] + 1.057311*xyz[2];
        return value_type(NumericTraits<T>::fromRealPromote(red * max_),
                          NumericTraits<T>::fromRealPromote(green * max_),
                          NumericTraits<T>::fromRealPromote(blue * max_));
    }
};

template <class T>
class XYZ2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef RGBValue<T> value_type;
    
    XYZ2RGBPrimeFunctor(component_type max = 255.0)
    : max_(max), gamma_(0.45)
    {}
    
    template <class V>
    value_type operator()(TinyVector<V, 3> const & xyz)
    {
        component_type red =    3.240479*xyz[0] - 1.537150*xyz[1] - 0.498535*xyz[2];
        component_type green = -0.969256*xyz[0] + 1.875992*xyz[1] + 0.041556*xyz[2];
        component_type blue =   0.055648*xyz[0] - 0.204043*xyz[1] + 1.057311*xyz[2];
        return value_type(NumericTraits<T>::fromRealPromote(std::pow(red, gamma_) * max_),
                          NumericTraits<T>::fromRealPromote(std::pow(green, gamma_) * max_),
                          NumericTraits<T>::fromRealPromote(std::pow(blue, gamma_) * max_));
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
        component_type L = xyz[1] < 0.008856 ?
                              903.3 * xyz[1] :
                              116.0 * std:pow(xyz[1], gamma_) - 16.0;
        component_type denom = xyz[0] + 15.0*xyz[1] + 3.0*xyz[2];
        component_type uprime = 4.0 * xyz[0] / denom;
        component_type vprime = 9.0 * xyz[1] / denom;
        return value_type(L, 13.0*L*(uprime - 0.197839), 13.0*L*(vprime - 0.468342));
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
        component_type uprime = luv[1] / 13.0 / luv[0] + 0.197839;
        component_type vprime = luv[2] / 13.0 / luv[0] + 0.468342;
        
        component_type Y = luv[0] < 8.0 ?
                              luv[0] / 903.3 :
                              std::pow((luv[0] + 16.0) / 116.0, gamma_);
        component_type X = 9.0*uprime*Y / 4.0 / vprime;
        component_type Z = ((9.0 / vprime - 15.0)*Y - X)/ 3.0; 
        return value_type(X, Y, Z);
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
        component_type xgamma = std:pow(xyz[0] / 0.950456, gamma_);
        component_type ygamma = std:pow(xyz[1], gamma_);
        component_type zgamma = std:pow(xyz[2] / 1.088754, gamma_);
        component_type L = xyz[1] < 0.008856 ?
                              903.3 * xyz[1] :
                              116.0 * ygamma - 16.0;
        return value_type(L, 500.0*(xgamma - ygamma), 200.0*(ygamma - zgamma));
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
        component_type ygamma = std:pow(Y, 1.0 / gamma_);
        component_type X = std::pow(lab[1] / 500.0 + ygamma, gamma_) * 0.950456;
        component_type Z = std::pow(-lab[2] / 200.0 + ygamma, gamma_) * 1.088754;
        return value_type(X, Y, Z);
    }
};

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
    
    value_type operator()(RGBValue<T> const & rgb)
    {
        return xyz2luv(rgb2xyz(rgb));
    }
};

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
    
    value_type operator()(RGBValue<T> const & rgb)
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
    Labv2XYZFunctor<component_type> lab2xyz;
    
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
    
    value_type operator()(RGBValue<T> const & rgb)
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
    
    value_type operator()(RGBValue<T> const & rgb)
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
    Labv2XYZFunctor<component_type> lab2xyz;
    
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

template <class T>
class RGBPrime2YprimePbPrFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef TinyVector<component_type, 3> value_type;
    
    RGBPrime2YprimePbPrFunctor(component_type max = 255.0)
    : max_(max)
    {}
    
    value_type operator()(RGBValue<T> const & rgb)
    {
        component_type red = rgb.red() / max_;
        component_type green = rgb.green() / max_;
        component_type blue = rgb.blue() / max_;
        
        return value_type(0.299*red + 0.587*green + 0.114*blue,
                          -0.168736*red - 0.331264*green + 0.5*blue,
                          0.5*red - 0.418688*green - 0.081312*blue);
    }
};

template <class T>
class YprimePbPr2RGBPrimeFunctor
{
    typedef typename NumericTraits<T>::RealPromote component_type;
    
    component_type max_;
    
  public:
  
    typedef RGBValue<T> value_type;
    
    YprimePbPr2RGBPrimeFunctor(component_type max = 255.0)
    : max_(max)
    {}
    
    template <class V>
    value_type operator()(TinyVector<V, 3> const & ypbpr)
    {
        component_type nred =   ypbpr[0] + 1.402*ypbpr[2];
        component_type ngreen = ypbpr[0] - 0.344136*ypbpr[1] - 0.714136*ypbpr[2];
        component_type nblue =  ypbpr[0] + 1.772*ypbpr[1];
        return value_type(NumericTraits<T>::fromRealPromote(nred * max_),
                          NumericTraits<T>::fromRealPromote(ngree * max_),
                          NumericTraits<T>::fromRealPromote(nblue * max_));
    }
};



} // namespace vigra 

#endif /* VIGRA_COLORCONVERSIONS_HXX */
