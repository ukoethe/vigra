
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
    res.reshapeIfEmpty(image.taggedShape(),
            "brightness(): Output images has wrong dimensions");

    double min = 0.0, max = 0.0;
    bool computeRange = !parseRange(range, &min, &max, "brightness(): Invalid range argument.");
    
    {
        PyAllowThreads _pythread;
        if(computeRange)
        {
            FindMinMax<PixelType> minmax;
            inspectMultiArray(srcMultiArrayRange(image), minmax);
            min = minmax.min;
            max = minmax.max;
        }

        vigra_precondition(min < max,
              "brightness(): Range upper bound must be greater than lower bound.");

        transformMultiArray(srcMultiArrayRange(image), destMultiArray(res),
                            BrightnessFunctor<PixelType>(factor, min, max));
    }
    return res;
}

template < class PixelType, unsigned int N>
NumpyAnyArray
pythonContrastTransform(NumpyArray<N, Multiband<PixelType> > image,
                        double factor,
                        python::object range,
                        NumpyArray<N, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(image.taggedShape(),
            "contrast(): Output images has wrong dimensions");

    double min = 0.0, max = 0.0;
    bool computeRange = !parseRange(range, &min, &max, "contrast(): Invalid range argument.");
    {
        PyAllowThreads _pythread;
        if(computeRange)
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
    }
    return res;
}

template < class PixelType, unsigned int N>
NumpyAnyArray
pythonGammaTransform(NumpyArray<N, Multiband<PixelType> > image,
                     double gamma,
                     python::object range,
                     NumpyArray<N, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(image.taggedShape(),
            "gamma_correction(): Output images has wrong dimensions");

    double min = 0.0, max = 0.0;
    bool computeRange = !parseRange(range, &min, &max, "gamma_correction(): Invalid range argument.");    
    {
        PyAllowThreads _pythread;
        if(computeRange)
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
    }
    return res;
}

template < class SrcPixelType, class DestPixelType, unsigned int N>
NumpyAnyArray
pythonLinearRangeMapping(NumpyArray<N, Multiband<SrcPixelType> > image,
                         python::object oldRange,
                         python::object newRange,
                         NumpyArray<N, Multiband<DestPixelType> > res)
{
    res.reshapeIfEmpty(image.taggedShape(),
        "linearRangeMapping(): Output images has wrong dimensions");

    double oldMin = 0.0, oldMax = 0.0,
           newMin = 0.0, newMax = 0.0;
    bool computeRange = !parseRange(oldRange, &oldMin, &oldMax, 
                                    "linearRangeMapping(): Argument 'oldRange' is invalid.");
    if(!parseRange(newRange, &newMin, &newMax, "linearRangeMapping(): Argument 'newRange' is invalid."))
    {
        newMin = 0.0;
        newMax = 255.0;
    }
    
    {
        PyAllowThreads _pythread;
        if(computeRange)
        {
            FindMinMax<SrcPixelType> minmax;
            inspectMultiArray(srcMultiArrayRange(image), minmax);
            oldMin = minmax.min;
            oldMax = minmax.max;
        }
        
        vigra_precondition(oldMin < oldMax && newMin < newMax,
              "linearRangeMapping(): Range upper bound must be greater than lower bound.");

        transformMultiArray(srcMultiArrayRange(image), destMultiArray(res),
                            linearRangeMapping(oldMin, oldMax, newMin, newMax));
    }
    
    return res;
}

template < class SrcPixelType>
inline NumpyAnyArray
pythonLinearRangeMapping2D(NumpyArray<3, Multiband<SrcPixelType> > image,
                           python::object oldRange,
                           python::object newRange,
                           NumpyArray<3, Multiband<UInt8> > res)
{
    return pythonLinearRangeMapping(image, oldRange, newRange, res);
}

VIGRA_PYTHON_MULTITYPE_FUNCTOR(pyLinearRangeMapping2D, pythonLinearRangeMapping2D)

template < class PixelType, unsigned int N, class Functor >
NumpyAnyArray
pythonColorTransform(NumpyArray<N, TinyVector<PixelType, 3> > image,
                     NumpyArray<N, TinyVector<PixelType, 3> > res)
{
    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(Functor::targetColorSpace()),
        "colorTransform(): Output images has wrong dimensions");

    {
        PyAllowThreads _pythread;
        transformMultiArray(srcMultiArrayRange(image), destMultiArray(res), Functor());
    }
    return res;
}

#define exportColorTransform(name) \
def("transform_" #name, registerConverters(&pythonColorTransform<float, 2, name##Functor<float> >), \
    (arg("image"), arg("out")=object()), \
    "Convert the colors of the given 'image' using " #name "Functor.\n" \
    "\n" \
    "For details see " #name "Functor_ in the C++ documentation.\n")

template<class T>
NumpyAnyArray pythonApplyColortable(const NumpyArray<2, Singleband<T> >& valueImage, 
                                    const NumpyArray<2, UInt8>& colortable,
                                    NumpyArray<3, Multiband<npy_uint8> > res =  NumpyArray<3, Multiband<npy_uint8> >())
{
    vigra_precondition(!colortable.axistags(),
                       "applyColortable(): colortable must not have axistags\n"
                       "(use 'array.view(numpy.ndarray)' to remove them).");
    
    // Singleband: there is only a singleton channel axis (which is removed when converted from python numpy array to C++
    // Multiband: channel axis is allowed to be singleband, but does not have to be,
    //            will be last when converted Python -> C++ and channel axis is counted in the dimension ('3')
    typedef NumpyArray<2, Singleband<T> > InputType;
    
    res.reshapeIfEmpty(valueImage.taggedShape().setChannelCount(colortable.shape(1)),
                       "pythonApplyColortable: shape of res is wrong");
    
    const unsigned int N = colortable.shape(0);

    bool startsWithTransparent = (colortable(0,3) == 0);

    for(MultiArrayIndex c=0; c<colortable.shape(1); ++c)
    {
        MultiArrayView<2, UInt8>::iterator channelIter = res.bind<2>(c).begin();
        
        //make an unstrided copy of the current column of the colortable
        ArrayVector<UInt8> ctable(colortable.bind<1>(c).begin(), colortable.bind<1>(c).end());
        
        for(typename InputType::const_iterator v = valueImage.begin(); v != valueImage.end(); ++v, ++channelIter)
        {
            if (*v == 0)
            {
                *channelIter = ctable[0];
            }
            else if (startsWithTransparent)
            {
                // Special behavior: If the colortable has too few values for the image,
                // we simply repeat the table for the higher indexes (see below).
                // BUT:
                // It's common for the colortable to start with a transparent value for "background".
                // In that case, we only repeat the remaining values in the colortable (don't repeat the transparent color)
                *channelIter = ctable[((*v-1) % (N-1))+1];
            }
            else
            {
                // Use % to repeat the colortable if there aren't enough values in the colortable.
                *channelIter = ctable[*v % N];
            }
        }
    }
    
    return res;
}
VIGRA_PYTHON_MULTITYPE_FUNCTOR(pyApplyColortable, pythonApplyColortable)

template<class T>
void pythonGray2QImage_ARGB32Premultiplied(
    const NumpyArray<2, Singleband<T> >& image, 
    NumpyArray<3, Multiband<npy_uint8> > qimageView,
    NumpyArray<1, T> normalize = boost::python::object()
) 
{
    vigra_precondition(image.isUnstrided() || image.transpose().isUnstrided(),
        "gray2qimage_ARGB32Premultiplied(): Can only handle arrays with contiguous memory.");
    typedef typename NumericTraits<T>::RealPromote TmpType;
    
    T* data = image.data(); 
    const T* dataEnd = data+image.size();
    UInt8* imgData = qimageView.data();
    UInt8 pixel = 0;
    TmpType pixelF;
    
    TmpType normalizeLow, normalizeHigh; 
    if(normalize.pyObject() != Py_None)
    {
        vigra_precondition(normalize.shape(0) == 2,
            "gray2qimage_ARGB32Premultiplied(): normalize.shape[0] == 2 required.");
            
        //normalize = None
        normalizeLow = normalize[0];
        normalizeHigh = normalize[1];
        
        vigra_precondition(normalizeHigh > normalizeLow,
            "gray2qimage_ARGB32Premultiplied(): normalize[0] < normalize[1] is required.");
            
        const TmpType f = TmpType(255) / static_cast<TmpType>(normalizeHigh-normalizeLow);
        
        while(data < dataEnd) 
        {
            pixelF = detail::RequiresExplicitCast<TmpType>::cast(*data);
            
            if(pixelF < normalizeLow) 
            {
                pixel = 0;
            }
            else if(pixelF > normalizeHigh) 
            {
                pixel = 255;
            }
            else
            {
                pixel = detail::RequiresExplicitCast<UInt8>::cast((pixelF-normalizeLow)*f);
            }
            *imgData = pixel; ++imgData; //B
            *imgData = pixel; ++imgData; //G
            *imgData = pixel; ++imgData; //R
            *imgData = 255;   ++imgData; //A
            ++data;
        }
    }
    else 
    {
        while(data < dataEnd) 
        {
            pixel = detail::RequiresExplicitCast<UInt8>::cast(*data);
            *imgData = pixel; ++imgData; //B
            *imgData = pixel; ++imgData; //G
            *imgData = pixel; ++imgData; //R
            *imgData = 255;   ++imgData; //A
            ++data;
        }
    }
}
VIGRA_PYTHON_MULTITYPE_FUNCTOR(pyGray2QImage_ARGB32Premultiplied, pythonGray2QImage_ARGB32Premultiplied)

template<class T>
void pythonAlphaModulated2QImage_ARGB32Premultiplied(
    const NumpyArray<2, Singleband<T> >& image, 
    NumpyArray<3, Multiband<npy_uint8> > qimageView,
    NumpyArray<1, float> tintColor,
    NumpyArray<1, T> normalize
) 
{
    vigra_precondition(image.isUnstrided() || image.transpose().isUnstrided(),
        "alphamodulated2qimage_ARGB32Premultiplied(): Can only handle arrays with contiguous memory.");
    typedef typename NumericTraits<T>::RealPromote TmpType;
    
    vigra_precondition(normalize.shape(0) == 2,
        "alphamodulated2qimage_ARGB32Premultiplied(): normalize.shape[0] == 2 required.");
    vigra_precondition(tintColor.shape(0) == 3,
        "alphamodulated2qimage_ARGB32Premultiplied(): tintColor.shape[0] == 3 required.");
            
    const TmpType l = normalize[0];
    const TmpType h = normalize[1];
    
    vigra_precondition(h > l,
        "alphamodulated2qimage_ARGB32Premultiplied(): normalize[0] < normalize[1] is required.");
    
    const TmpType r = tintColor[0];
    const TmpType g = tintColor[1];
    const TmpType b = tintColor[2];
    
    T* data = image.data();
    const T* dataEnd = image.data()+image.size();
    unsigned char* imgData = qimageView.data();
    TmpType pixelF;
    const TmpType f = TmpType(255) / static_cast<TmpType>(h-l);
    while(data < dataEnd) 
    {
        pixelF = detail::RequiresExplicitCast<TmpType>::cast(*data);
        if(pixelF < l)
        {
            pixelF = TmpType();
        }
        else if(pixelF > h)
        {
            pixelF = TmpType(255);
        }
        else
        {
            pixelF = (pixelF-l)*f;
        }
        *imgData = detail::RequiresExplicitCast<UInt8>::cast(pixelF*b); ++imgData; //B
        *imgData = detail::RequiresExplicitCast<UInt8>::cast(pixelF*g); ++imgData; //G
        *imgData = detail::RequiresExplicitCast<UInt8>::cast(pixelF*r); ++imgData; //R
        *imgData = detail::RequiresExplicitCast<UInt8>::cast(pixelF);   ++imgData; //A
        ++data;
    }
}
VIGRA_PYTHON_MULTITYPE_FUNCTOR(pyAlphaModulated2QImage_ARGB32Premultiplied, pythonAlphaModulated2QImage_ARGB32Premultiplied)

void defineColors()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    multidef("applyColortable", pyApplyColortable<vigra::Int8, vigra::UInt8, vigra::Int16, vigra::UInt16, vigra::Int32, vigra::UInt32>(),
        (arg("valueImage"), 
        arg("colortable"),
        arg("out")=python::object()), 
        "Applies a colortable to the given 2D valueImage.\n\n"
        "Colortable must have 4 columns, each row represents a color (for example, RGBA). \n"
        "Values in valueImage are first taken modulo the length of the colortable. \n"
        "In the special case where the first color in the table is transparent, that value "
        "is NOT repeated for values outside the colortable length.\n\n"
        "Returns: uint8 image with 4 channels\n");
    
    multidef("gray2qimage_ARGB32Premultiplied", pyGray2QImage_ARGB32Premultiplied<vigra::Int8, vigra::UInt8, vigra::Int16, vigra::UInt16, vigra::Int32, vigra::UInt32, float, double>(),
        (arg("image"), 
        arg("qimage"),
        arg("normalize")=python::object()), 
        "Convert the image (single-band) into a QImage of format Format_ARGB32_Premultiplied.\n"
        "\n"
        "import qimage2ndarray\n"
        "qimg = QImage(a.shape[0], a.shape[1], QImage.Format_ARGB32_Premultiplied)\n"
        "normalize = numpy.asarray([10, 217], dtype=image.dtype)\n"
        "vigra.colors.gray2qimage_ARGB32Premultiplied(a, qimage2ndarray.byte_view(qimg), normalize)\n");
    
    multidef("alphamodulated2qimage_ARGB32Premultiplied", pyAlphaModulated2QImage_ARGB32Premultiplied<vigra::Int8, vigra::UInt8, vigra::Int16, vigra::UInt16, vigra::Int32, vigra::UInt32, float, double>(),
        (arg("image"), 
        arg("qimage"),
        arg("tintColor"),
        arg("normalize")), 
        "Convert the image (single-band) into a QImage of format Format_ARGB32_Premultiplied.\n"
        "\n"
        "import qimage2ndarray\n"
        "qimg = QImage(a.shape[0], a.shape[1], QImage.Format_ARGB32_Premultiplied)\n"
        "normalize = numpy.asarray([10, 217], dtype=image.dtype)\n"
        "tintColor = numpy.asarray([1.0, 0.0, 0.0], dtype=numpy.float32) #RGB\n"
        "vigra.colors.alphamodulated2qimage_ARGB32Premultiplied(a, qimage2ndarray.byte_view(qimg), tintColor, normalize)\n");
    
    def("brightness",
         registerConverters(&pythonBrightnessTransform<float, 3>),
         (arg("image"), arg("factor"), arg("range")=make_tuple(0.0, 255.0), arg("out")=object()),
        "Adjust the brightness of a 2D scalar or multiband image. The function applies the formula::\n"
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
         (arg("volume"), arg("factor"), arg("range")=make_tuple(0.0, 255.0), arg("out")=object()),
         "Likewise for a 3D scalar or multiband volume.\n");

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
         (arg("volume"), arg("factor"), arg("range")=make_tuple(0.0, 255.0), arg("out")=object()),
         "Likewise for a 3D scalar or multiband volume.\n");

    def("gammaCorrection",
         registerConverters(&pythonGammaTransform<float, 3>),
         (arg("image"), arg("gamma"), arg("range")=make_tuple(0.0, 255.0), arg("out")=object()),
        "Adjust gamma correction to an image or volume. The function applies the formula::\n"
        "\n"
        "    diff = range[1] - range[0]\n"
        "    out = pow((image - range[0]) / diff, 1.0 / gamma) * diff + range[0]\n"
        "\n"
        "to each element of the array. 'gamma' and 'range[1] - range[0]' must be "
        "positive. Elements outside the given range are clipped at the range borders. "
        "If 'range' is None or \"\" or \"auto\", the range is set to the actual range "
        "of 'image'::\n"
        "\n"
        "   range = image.min(), image.max()\n\n");

    def("gammaCorrection",
         registerConverters(&pythonGammaTransform<float, 4>),
         (arg("volume"), arg("gamma"), arg("range")=make_tuple(0.0, 255.0), arg("out")=object()),
         "Likewise for a 3D scalar or multiband volume.\n");

    multidef("linearRangeMapping", pyLinearRangeMapping2D<vigra::Int8, vigra::UInt8, vigra::Int16, vigra::UInt16, vigra::Int32, vigra::UInt32, float, double>(),
         (arg("image"), arg("oldRange")="auto", arg("newRange")=make_tuple(0.0, 255.0), arg("out")=object()),
        "Convert the intensity range of a 2D scalar or multiband image. The function applies a linear transformation "
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
        "If 'newRange' is None or \"\" or \"auto\", it is set to (0, 255.0). "
        "If 'out' is explicitly passed, it must be a uint8 image.\n");

    def("linearRangeMapping",
         registerConverters(&pythonLinearRangeMapping<float, float, 3>),
         (arg("image"), arg("oldRange")="auto", arg("newRange")=make_tuple(0.0, 255.0), arg("out")=object()),
         "Likewise, but #in' and 'out' are float32 images.\n");

    def("linearRangeMapping",
         registerConverters(&pythonLinearRangeMapping<float, UInt8, 4>),
         (arg("volume"), arg("oldRange")="auto", arg("newRange")=make_tuple(0.0, 255.0), arg("out")=object()),
         "Likewise for a 3D scalar or multiband volume, when 'in' is a float32 and 'out' a unit8 volume.\n");

    def("linearRangeMapping",
         registerConverters(&pythonLinearRangeMapping<float, float, 4>),
         (arg("volume"), arg("oldRange")="auto", arg("newRange")=make_tuple(0.0, 255.0), arg("out")=object()),
         "Likewise, but 'in' and 'out' are float32 volumes.\n");


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
