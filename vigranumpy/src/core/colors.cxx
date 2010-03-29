
/************************************************************************/
/*                                                                      */
/*                 Copyright 2009 by Ullrich Koethe                     */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpycolors_PyArray_API
//#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/multi_pointoperators.hxx>
#include <vigra/colorconversions.hxx>
#include <vigra/mathutil.hxx>

namespace python = boost::python;

namespace vigra
{

template<class ValueType>
struct BrightnessFunctor
{
    double factor, min, max, diff;
    
    typedef ValueType argument_type;
    typedef ValueType result_type;
    
    BrightnessFunctor(double f, double mi, double ma)
    : factor(0.0), min(mi), max(ma), diff(ma-mi)
    {
        vigra_precondition(f > 0.0,
            "brightness(): Factor must be positive.");
        vigra_precondition(diff > 0.0,
            "brightness(): Range upper bound must be greater than lower bound.");
        factor = 0.25 * diff * std::log(f);
    }
    
    result_type operator()(argument_type v) const
    {
        double r = v + factor;
        
        return detail::RequiresExplicitCast<result_type>::cast(
                   r < min
                       ? min
                       : r > max
                            ? max
                            : r);
    }
};

template<class ValueType>
struct ContrastFunctor
{
    double factor, min, max, diff, offset;
    
    typedef ValueType argument_type;
    typedef ValueType result_type;
    
    ContrastFunctor(double f, double mi, double ma)
    : factor(f), min(mi), max(ma), diff(0.5*(ma-mi)), offset(diff*(1.0-f))
    {
        vigra_precondition(f > 0.0,
            "contrast(): Factor must be positive.");
        vigra_precondition(diff > 0.0,
            "contrast(): Range upper bound must be greater than lower bound.");
    }
   
    result_type operator()(argument_type v) const
    {
        double r = v * factor + offset;
        
        return detail::RequiresExplicitCast<result_type>::cast(
                   r < min
                       ? min
                       : r > max
                            ? max
                            : r);
    }
};

inline
bool parseRange(python::object range, double * min, double * max,
                const char * message)
{
    if(!range)
    {
        return false;
    }
    python::extract<std::string> isString(range);
    if(isString.check())
    {
        std::string text(isString());
        if(text == "" || text == "auto")
            return false;
        vigra_precondition(false, message);
    }
    python::extract<python::tuple> isTuple(range);
    if(isTuple.check())
    {
        python::extract<double> arg1(isTuple()[0]),
                                arg2(isTuple()[1]);
        if(arg1.check() && arg2.check())
        {
            *min = arg1();
            *max = arg2();
            return true;
        }
    }
    vigra_precondition(false, message);
    return false;
}

template < class PixelType, unsigned int N>
NumpyAnyArray 
pythonBrightnessTransform(NumpyArray<N, Multiband<PixelType> > image,
                          double factor,
                          python::object range,
                          NumpyArray<N, Multiband<PixelType> > res)
{
    double min = 0.0, max = 0.0;
    if(!parseRange(range, &min, &max, "brightness(): Invalid range argument."))
    {
        FindMinMax<PixelType> minmax;
        inspectMultiArray(srcMultiArrayRange(image), minmax);
        min = minmax.min;
        max = minmax.max;
    }
    
    vigra_precondition(min < max,
          "brightness(): Range upper bound must be greater than lower bound.");
    
    res.reshapeIfEmpty(image.shape(),"brightness(): Output images has wrong dimensions");

    transformMultiArray(srcMultiArrayRange(image), destMultiArray(res), 
                        BrightnessFunctor<PixelType>(factor, min, max));
    return res;
}

template < class PixelType, unsigned int N>
NumpyAnyArray 
pythonContrastTransform(NumpyArray<N, Multiband<PixelType> > image,
                        double factor,
                        python::object range,
                        NumpyArray<N, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(image.shape(),"contrast(): Output images has wrong dimensions");

    double min = 0.0, max = 0.0;
    if(!parseRange(range, &min, &max, "contrast(): Invalid range argument."))
    {
        FindMinMax<PixelType> minmax;
        inspectMultiArray(srcMultiArrayRange(image), minmax);
        min = minmax.min;
        max = minmax.max;
    }
    
    vigra_precondition(min < max,
          "contrast(): Range upper bound must be greater than lower bound.");
    
    transformMultiArray(srcMultiArrayRange(image), destMultiArray(res), 
                        ContrastFunctor<PixelType>(factor, min, max));
    return res;
}

template < class PixelType, unsigned int N>
NumpyAnyArray 
pythonGammaTransform(NumpyArray<N, Multiband<PixelType> > image,
                     double gamma,
                     python::object range,
                     NumpyArray<N, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(image.shape(),"gamma_correction(): Output images has wrong dimensions");

    double min = 0.0, max = 0.0;
    if(!parseRange(range, &min, &max, "gamma_correction(): Invalid range argument."))
    {
        FindMinMax<PixelType> minmax;
        inspectMultiArray(srcMultiArrayRange(image), minmax);
        min = minmax.min;
        max = minmax.max;
    }
    
    vigra_precondition(min < max,
          "gamma_correction(): Range upper bound must be greater than lower bound.");
    
    transformMultiArray(srcMultiArrayRange(image), destMultiArray(res), 
                        GammaFunctor<PixelType>(1.0 / gamma, PixelType(min), PixelType(max)));
    return res;
}

template < class SrcPixelType, class DestPixelType, unsigned int N>
NumpyAnyArray 
pythonLinearRangeMapping(NumpyArray<N, Multiband<SrcPixelType> > image,
                         python::object oldRange,
                         python::object newRange,
                         NumpyArray<N, Multiband<DestPixelType> > res)
{
    res.reshapeIfEmpty(image.shape(),"linearRangeMapping(): Output images has wrong dimensions");

    double oldMin = 0.0, oldMax = 0.0,
           newMin = 0.0, newMax = 0.0;
    if(!parseRange(oldRange, &oldMin, &oldMax, "linearRangeMapping(): Argument 'oldRange' is invalid."))
    {
        FindMinMax<SrcPixelType> minmax;
        inspectMultiArray(srcMultiArrayRange(image), minmax);
        oldMin = minmax.min;
        oldMax = minmax.max;
    }
    if(!parseRange(newRange, &newMin, &newMax, "linearRangeMapping(): Argument 'newRange' is invalid."))
    {
        newMin = 0.0;
        newMax = 255.0;
    }
    
    vigra_precondition(oldMin < oldMax && newMin < newMax,
          "linearRangeMapping(): Range upper bound must be greater than lower bound.");
    
    transformMultiArray(srcMultiArrayRange(image), destMultiArray(res), 
                        linearRangeMapping(oldMin, oldMax, newMin, newMax));
    return res;
}


template < class PixelType, unsigned int N, class Functor >
NumpyAnyArray 
pythonColorTransform(NumpyArray<N, TinyVector<PixelType, 3> > image,
                     NumpyArray<N, TinyVector<PixelType, 3> > res)
{
    res.reshapeIfEmpty(image.shape(),"colorTransform(): Output images has wrong dimensions");

    transformMultiArray(srcMultiArrayRange(image), destMultiArray(res), Functor());
    return res;
}

#define exportColorTransform(name) \
def("transform_" #name, registerConverters(&pythonColorTransform<float, 2, name##Functor<float> >), \
    (arg("image"), arg("out")=object()), \
    "Convert the colors of the given 'image' using " #name "Functor.\n" \
    "\n" \
    "For details see " #name "Functor_ in the C++ documentation.\n")



void defineColors()
{
    using namespace python;
    
    def("brightness", 
         registerConverters(&pythonBrightnessTransform<float, 3>),
         (arg("image"), arg("factor"), arg("range")=make_tuple(0.0, 255.0), arg("out")=object()),
        "Adjust the brightness of an image or volume. The function applies the formula::\n"
        "\n"
        "   out = image + 0.25 * log(factor) * (range[1] - range[0])\n"
        "\n"
        "to each element of the array. 'factor' and 'range[1] - range[0]' must be "
        "positive. Elements outside the given range are clipped at the range borders. "
        "If 'range' is None or \"\" or \"auto\", the range is set to the actual range "
        "of 'image'::\n"
        "\n"
        "   range = image.min(), image.max()\n\n");

    def("brightness", 
         registerConverters(&pythonBrightnessTransform<float, 4>),
         (arg("volume"), arg("factor"), arg("range")=make_tuple(0.0, 255.0), arg("out")=object()));
         
    def("contrast", 
         registerConverters(&pythonContrastTransform<float, 3>),
         (arg("image"), arg("factor"), arg("range")=make_tuple(0.0, 255.0), arg("out")=object()),
        "Adjust the contrast of an image or volume. The function applies the formula::\n"
        "\n"
        "    out = factor * image + (1.0 - factor) * (range[1] - range[0]) / 2.0\n"
        "\n"
        "to each element of the array. 'factor' and 'range[1] - range[0]' must be "
        "positive. Elements outside the given range are clipped at the range borders. "
        "If 'range' is None or \"\" or \"auto\", the range is set to the actual range "
        "of 'image'::\n"
        "\n"
        "   range = image.min(), image.max()\n\n");

    def("contrast", 
         registerConverters(&pythonContrastTransform<float, 4>),
         (arg("volume"), arg("factor"), arg("range")=make_tuple(0.0, 255.0), arg("out")=object()));
          
    def("gammaCorrection", 
         registerConverters(&pythonGammaTransform<float, 3>),
         (arg("image"), arg("factor"), arg("range")=make_tuple(0.0, 255.0), arg("out")=object()),
        "Adjust gamma correction to an image or volume. The function applies the formula::\n"
        "\n"
        "    diff = range[1] - range[0]\n"
        "    out = pow((image - range[0]) / diff, 1.0 / gamma) * diff + range[0]\n"
        "\n"
        "to each element of the array. 'factor' and 'range[1] - range[0]' must be "
        "positive. Elements outside the given range are clipped at the range borders. "
        "If 'range' is None or \"\" or \"auto\", the range is set to the actual range "
        "of 'image'::\n"
        "\n"
        "   range = image.min(), image.max()\n\n");

    def("gammaCorrection", 
         registerConverters(&pythonGammaTransform<float, 4>),
         (arg("volume"), arg("factor"), arg("range")=make_tuple(0.0, 255.0), arg("out")=object()));
          
    def("linearRangeMapping", 
         registerConverters(&pythonLinearRangeMapping<float, UInt8, 3>),
         (arg("image"), arg("oldRange")="auto", arg("newRange")=make_tuple(0.0, 255.0), arg("out")=object()),
        "Convert the intensity range of an image or volume. The function applies a linear transformation "
        "to the intensities such that the value oldRange[0] is mapped onto newRange[0], "
        "and oldRange[1] is mapped onto newRange[1]. That is, the algorithm applies the formula::\n"
        "\n"
        "    oldDiff = oldRange[1] - oldRange[0]\n"
        "    newDiff = newRange[1] - newRange[0]\n"
        "    out = (image - oldRange[0]) / oldDiff * newDiff + newRange[0]\n"
        "\n"
        "to each element of the array. 'oldDiff' and 'newDiff' must be "
        "positive. If 'oldRange' is None or \"\" or \"auto\" (the default), the range is set to " 
        "the actual range of 'image'::\n"
        "\n"
        "   range = image.min(), image.max()\n\n"
        "If 'newRange' is None or \"\" or \"auto\", it is set to (0, 255.0).\n");

    def("linearRangeMapping", 
         registerConverters(&pythonLinearRangeMapping<float, float, 3>),
         (arg("image"), arg("oldRange")="auto", arg("newRange")=make_tuple(0.0, 255.0), arg("out")=object()));
          
    def("linearRangeMapping", 
         registerConverters(&pythonLinearRangeMapping<float, UInt8, 4>),
         (arg("volume"), arg("oldRange")="auto", arg("newRange")=make_tuple(0.0, 255.0), arg("out")=object()));
          
    def("linearRangeMapping", 
         registerConverters(&pythonLinearRangeMapping<float, float, 4>),
         (arg("volume"), arg("oldRange")="auto", arg("newRange")=make_tuple(0.0, 255.0), arg("out")=object()));
          

    exportColorTransform(RGB2sRGB);
    exportColorTransform(sRGB2RGB);
    exportColorTransform(RGB2RGBPrime);
    exportColorTransform(RGBPrime2RGB);

    exportColorTransform(RGB2XYZ);
    exportColorTransform(RGBPrime2XYZ);
    exportColorTransform(XYZ2RGB);
    exportColorTransform(XYZ2RGBPrime);

    exportColorTransform(RGB2Lab);
    exportColorTransform(RGBPrime2Lab);
    exportColorTransform(XYZ2Lab);
    exportColorTransform(Lab2RGB);
    exportColorTransform(Lab2RGBPrime);
    exportColorTransform(Lab2XYZ);

    exportColorTransform(RGB2Luv);
    exportColorTransform(RGBPrime2Luv);
    exportColorTransform(XYZ2Luv);
    exportColorTransform(Luv2RGB);
    exportColorTransform(Luv2RGBPrime);
    exportColorTransform(Luv2XYZ);

    exportColorTransform(RGBPrime2YPrimePbPr);
    exportColorTransform(YPrimePbPr2RGBPrime);
    exportColorTransform(RGBPrime2YPrimeCbCr);
    exportColorTransform(YPrimeCbCr2RGBPrime);

    exportColorTransform(RGBPrime2YPrimeUV);
    exportColorTransform(YPrimeUV2RGBPrime);
    exportColorTransform(RGBPrime2YPrimeIQ);
    exportColorTransform(YPrimeIQ2RGBPrime);
}

} // namespace vigra

using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(colors)
{
    import_vigranumpy();
    defineColors();
}
