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
/*    The VIGRA Website is http://www.egd.igd.fhg.de/~ulli/vigra/       */
/*    Please direct questions, bug reports, and contributions to        */
/*                        ulli@egd.igd.fhg.de                           */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#include <iostream>
#include <algorithm>
#include "vigra/stdimage.hxx"
#include "vigra/stdimagefunctions.hxx"
#include "vigra/impex.hxx"
#include "vigra/colorconversions.hxx"

using namespace vigra;

template<class Polar2Cartesian, class Cartesian2RGB, class RGB2RGBPrime>
void createColorVsSaturation(BRGBImage & result, double brightness, 
                 Polar2Cartesian polar2Cartesian, Cartesian2RGB cartesian2RGB, RGB2RGBPrime rgb2RGBPrime)
{
    int w = result.width(); 
    int h = result.height();
    
    for(int y=0; y<h; ++y)
    {
        for(int x=0; x<w; ++x)
        {
            double saturation = (float)x / (w-1);
            double color = (float)y / (h-1) * 360.0;
            
            RGBValue<float> rgb = cartesian2RGB(
               polar2Cartesian(color, brightness, saturation));
                        
            if(saturation > 1.0 ||
               rgb.red() < 0.0 || rgb.red() > 255.0 ||
               rgb.green() < 0.0 || rgb.green() > 255.0 ||
               rgb.blue() < 0.0 || rgb.blue() > 255.0)
            {
                result(x,y) = RGBValue<unsigned char>(170.0);
            }
            else
            {
                result(x,y) = rgb2RGBPrime(rgb);
            }
        }
    }
}

template<class Polar2Cartesian, class Cartesian2RGB, class RGB2RGBPrime>
void createColorVsBrightness(BRGBImage & result, double saturation, 
                 Polar2Cartesian polar2Cartesian, Cartesian2RGB cartesian2RGB, RGB2RGBPrime rgb2RGBPrime)
{
    int w = result.width(); 
    int h = result.height();
    
    for(int y=0; y<h; ++y)
    {
        for(int x=0; x<w; ++x)
        {
            double brightness = (float)x / (w-1);
            double color = (float)y / (h-1) * 360.0;
            
            RGBValue<float> rgb = cartesian2RGB(
               polar2Cartesian(color, brightness, saturation));
            
            if(saturation > 1.0 ||
               rgb.red() < 0.0 || rgb.red() > 255.0 ||
               rgb.green() < 0.0 || rgb.green() > 255.0 ||
               rgb.blue() < 0.0 || rgb.blue() > 255.0)
            {
                result(x,y) = RGBValue<unsigned char>(170.0);
            }
            else
            {
                result(x,y) = rgb2RGBPrime(rgb);
            }
        }
    }
}

template<class Polar2Cartesian, class Cartesian2RGB, class RGB2RGBPrime>
void createSaturationVsBrightness(BRGBImage & result, double color, 
                 Polar2Cartesian polar2Cartesian, Cartesian2RGB cartesian2RGB, RGB2RGBPrime rgb2RGBPrime)
{
    int w = result.width(); 
    int h = result.height();
    
    for(int y=0; y<h; ++y)
    {
        for(int x=0; x<w; ++x)
        {
            double brightness = (float)x / (w-1);
            double saturation = (float)y / (h-1);
            
            RGBValue<float> rgb = cartesian2RGB(
               polar2Cartesian(color, brightness, saturation));
            
            if(saturation > 1.0 ||
               rgb.red() < 0.0 || rgb.red() > 255.0 ||
               rgb.green() < 0.0 || rgb.green() > 255.0 ||
               rgb.blue() < 0.0 || rgb.blue() > 255.0)
            {
                result(x,y) = RGBValue<unsigned char>(170.0);
            }
            else
            {
                result(x,y) = rgb2RGBPrime(rgb);
            }
        }
    }
}

template<class Polar2Cartesian, class Cartesian2RGB, class RGB2RGBPrime>
void createColorCircle(BRGBImage & result, double brightness, 
                 Polar2Cartesian polar2Cartesian, Cartesian2RGB cartesian2RGB, RGB2RGBPrime rgb2RGBPrime)
{
    int w = result.width(); 
    int h = result.height();
    
    for(int y=0; y<h; ++y)
    {
        for(int x=0; x<w; ++x)
        {
            double dx = x/128.0 - 1.0;
            double dy = -y/128.0 + 1.0;
            double color = 180.0/M_PI*std::atan2(dy,dx);
            double saturation = std::sqrt(dx*dx+dy*dy);
            
            RGBValue<float> rgb = cartesian2RGB(
               polar2Cartesian(color, brightness, saturation));
            
            if(saturation > 1.0 ||
               rgb.red() < 0.0 || rgb.red() > 255.0 ||
               rgb.green() < 0.0 || rgb.green() > 255.0 ||
               rgb.blue() < 0.0 || rgb.blue() > 255.0)
            {
                result(x,y) = RGBValue<unsigned char>(170.0);
            }
            else
            {
                result(x,y) = rgb2RGBPrime(rgb);
            }
        }
    }
}

void write(char const * colorspace, char const * diagram, int i, BRGBImage const & img)
{
    char buf[1000];
    if(i < 10)
        sprintf(buf, "%s_%s_0%d.gif", colorspace, diagram, i);
    else
        sprintf(buf, "%s_%s_%d.gif", colorspace, diagram, i);
    exportImage(srcImageRange(img), ImageExportInfo(buf));
    std::cout << "Wrote " << buf << std::endl;
}

template<class Polar2Cartesian, class Cartesian2RGB, class RGB2RGBPrime>
void createColorSpaceSlices(char const * colorspace,
    Polar2Cartesian polar2Cartesian, Cartesian2RGB cartesian2RGB, RGB2RGBPrime rgb2RGBPrime)
{
    int w = 257; 
    int h = 257;
    int Ymax = 10;
    
    for(int i=0; i<=Ymax; ++i)
    {

        BRGBImage result(w, h);

        createColorVsSaturation(result, (float)i/Ymax, 
                        polar2Cartesian, cartesian2RGB, rgb2RGBPrime);
        write(colorspace, "ColorVsSaturation", i, result);

        createColorVsBrightness(result, (float)i/Ymax, 
                        polar2Cartesian, cartesian2RGB, rgb2RGBPrime);
        write(colorspace, "ColorVsBrightness", i, result);

        createSaturationVsBrightness(result, (float)i/Ymax*360.0, 
                        polar2Cartesian, cartesian2RGB, rgb2RGBPrime);
        write(colorspace, "SaturationVsBrightness", i, result);

        createColorCircle(result, (float)i/Ymax, 
                        polar2Cartesian, cartesian2RGB, rgb2RGBPrime);
        write(colorspace, "ColorCircle", i, result);
    }
}

void usage(char const * prog)
{
    std::cerr << "Usage: " << prog << " colorspace\n"
                 "with colorspace in [lab luv ypbpr ycbcr yiq yuv]\n\n";
    std::cerr << "This programm calculates slices through the given color space\n"
                 "Images are named 'lab_SaturationVsBrightness_01.gif' etc.\n"
                 "where the first part of the name designates the colorspace used,\n"
                 "the second part says what is varied on the image\n"
                 "and the number codes the value of the quantity that is kept\n"
                 "constant in the image - 01 in the example means that the color\n"
                 "angle is 36 degrees = 1 * 360 degrees / 10\n";
}

int main(int argc, char ** argv)
{
    if(argc <2)
    {
        usage(argv[0]);
        return 1;
    }
    
    try
    {
        typedef TinyVector<float, 3> (*PolarFct)(double, double, double);
        std::string colorspace(argv[1]);
        
        if(colorspace == "lab")
        {
            createColorSpaceSlices("lab", 
                         (PolarFct)&polar2Lab, Lab2RGBFunctor<float>(),
                         RGB2RGBPrimeFunctor<float, unsigned char>());
        }
        else if(colorspace == "luv")
        {
            createColorSpaceSlices("luv", 
                         (PolarFct)&polar2Luv, Luv2RGBFunctor<float>(),
                         RGB2RGBPrimeFunctor<float, unsigned char>());
        }
        else if(colorspace == "ypbpr")
        {
            createColorSpaceSlices("ypbpr", 
                         (PolarFct)&polar2YPrimePbPr, YPrimePbPr2RGBPrimeFunctor<float>(),
                         &NumericTraits<RGBValue<unsigned char> >::fromRealPromote);
        }
        else if(colorspace == "ycbcr")
        {
            createColorSpaceSlices("ycbcr", 
                         (PolarFct)&polar2YPrimeCbCr, YPrimeCbCr2RGBPrimeFunctor<float>(),
                         &NumericTraits<RGBValue<unsigned char> >::fromRealPromote);
        }
        else if(colorspace == "yiq")
        {
            createColorSpaceSlices("yiq", 
                         (PolarFct)&polar2YPrimeIQ, YPrimeIQ2RGBPrimeFunctor<float>(),
                         &NumericTraits<RGBValue<unsigned char> >::fromRealPromote);
        }
        else if(colorspace == "yuv")
        {
            createColorSpaceSlices("yuv", 
                         (PolarFct)&polar2YPrimeUV, YPrimeUV2RGBPrimeFunctor<float>(),
                         &NumericTraits<RGBValue<unsigned char> >::fromRealPromote);
        }
        else
        {
            std::cerr << "Unknown colorspace: " << colorspace << std::endl;
            usage(argv[0]);
            return 1;
        }
    }
    catch (StdException & e)
    {
        std::cout << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
