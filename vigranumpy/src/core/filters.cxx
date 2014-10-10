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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyfilters_PyArray_API
//#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/nonlineardiffusion.hxx>
#include <vigra/rank_filter.hxx>
#include <vigra/symmetry.hxx>
#include <vigra/tv_filter.hxx>

#include <utility>

namespace python = boost::python;

namespace vigra
{

template <class InValue, class OutValue>
NumpyAnyArray
pythonNonlinearDiffusion2D(NumpyArray<3, Multiband<InValue> > image,
                           double edgeThreshold, double scale,
                           NumpyArray<3, Multiband<OutValue> > res=NumpyArray<3, Multiband<float> >())
{
    res.reshapeIfEmpty(image.taggedShape(),
        "nonlinearDiffusion2D(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        for(int k=0; k<image.shape(2); ++k)
        {
            MultiArrayView<2, OutValue, StridedArrayTag> bres = res.bindOuter(k);
            nonlinearDiffusion(srcImageRange(image.bindOuter(k)),
                               destImage(bres),
                               DiffusivityFunctor< double >(edgeThreshold), scale);
        }
    }
    return res;
}


template <class InValue, class OutValue>
NumpyAnyArray
pythonTotalVariationFilter2D(NumpyArray<2, Singleband<InValue> > image,
                           double alpha, int steps, double eps = 0,
                           NumpyArray<2, Singleband<OutValue> > res = python::object())
{
    std::string description("totalVariationFilter, alpha, steps, eps=");
    description += asString(eps);

    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description),
            "totalVariationFilter(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        totalVariationFilter(image, res, alpha, steps, eps);
    }
    return res;
}

template <class InValue, class InValue2, class OutValue>
NumpyAnyArray
pythonTotalVariationFilter2D(NumpyArray<2, Singleband<InValue> > image,
                             NumpyArray<2, Singleband<InValue2> > weight,
                           double alpha, int steps, double eps = 0,
                           NumpyArray<2, Singleband<OutValue> > res = python::object())
{
    std::string description("totalVariationFilter, weight, alpha, steps, eps=");
    description += asString(eps);

    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description),
            "totalVariationFilter(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        totalVariationFilter(image, weight, res, alpha, steps, eps);
    }
    return res;
}


template < class SrcPixelType >
NumpyAnyArray
pythonRadialSymmetryTransform2D(NumpyArray<2, Singleband<SrcPixelType> > image,
                                double scale = 1.0,
                                NumpyArray<2, Singleband<SrcPixelType> > res = python::object())
{
    std::string description("radial symmetry transform, scale=");
    description += asString(scale);

    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description),
            "radialSymmetryTransform2D(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;
        radialSymmetryTransform(srcImageRange(image), destImage(res), scale);
    }
    return res;
}


template < unsigned int N, class SrcPixelType >
NumpyAnyArray pythonLineRankOrderFilter(const NumpyArray< N, Singleband<SrcPixelType> > & image,
        unsigned long half_length, double rank, unsigned int axis = N - 1,
        NumpyArray< N, Singleband<SrcPixelType> > res = python::object())
{
    std::string description("rank order filter over 1-Dimension, axis=");
    description += asString(axis);

    vigra_precondition((2*half_length + 1) <= image.shape(axis),
            "lineRankOrderFilter(): Window must be no bigger than the image.");

    res.reshapeIfEmpty(image.taggedShape().setChannelDescription(description),
            "lineRankOrderFilter(): Output array has wrong shape.");

    {
        PyAllowThreads _pythread;

        MultiArrayView<N, SrcPixelType, StridedArrayTag> image_view = image;
        MultiArrayView<N, SrcPixelType, StridedArrayTag> res_view = res;
        lineRankOrderFilter(image_view, res_view, half_length, rank, axis);
    }
    return res;
}



void defineFilters2D()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    def("nonlinearDiffusion",
        registerConverters(&pythonNonlinearDiffusion2D<float, float>),
        (arg("image"), arg("edgeThreshold"), arg("scale"), arg("out")=python::object()),
        "Perform edge-preserving smoothing at the given scale."
        "\n\n"
        "For details see nonlinearDiffusion_ in the vigra C++ documentation.\n");

    def("totalVariationFilter",
        registerConverters(&pythonTotalVariationFilter2D<double,double>),
        (arg("image"), arg("alpha"), arg("steps"), arg("eps"), arg("out")=python::object()),
        "Perform total variation filter on 2D single band images."
        "\n\n"
        "For details see totalVariationFilter in the vigra C++ documentation.\n");

    def("totalVariationFilter",
        registerConverters(&pythonTotalVariationFilter2D<double,double,double>),
        (arg("image"), arg("weight"), arg("alpha"), arg("steps"), arg("eps"), arg("out")=python::object()),
        "Perform weighted total variation filter on 2D single band images."
        "\n\n"
        "For details see totalVariationFilter in the vigra C++ documentation.\n");

    def("radialSymmetryTransform2D",
        registerConverters(&pythonRadialSymmetryTransform2D<float>),
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Find centers of radial symmetry in an 2D image.\n\n"
        "This algorithm implements the Fast Radial Symmetry Transform according to "
        "[G. Loy, A. Zelinsky: \"A Fast Radial Symmetry Transform for Detecting Points of Interest\", "
        "in: A. Heyden et al. (Eds.): Proc. of 7th European Conf. on Computer Vision, Part 1, pp. 358-368, Springer LNCS 2350, 2002]\n\n"
        "For details see radialSymmetryTransform_ in the vigra C++ documentation.\n");
}


void defineRankFilter()
{
    using namespace python;

    docstring_options doc_options(true, true, false);

    def("lineRankOrderFilter",
            registerConverters(&pythonLineRankOrderFilter<1, float>),
            (arg("image"), arg("half_length"), arg("rank"), arg("axis")=(1-1), arg("out")=python::object()),
            "Convolution along a single dimension of a 2D scalar or multiband image. "
                    "'kernel' must be an instance of Kernel1D.\n"
                    "\n"
                    "For details see convolveMultiArrayOneDimension_ in the vigra C++ documentation.\n");

    def("lineRankOrderFilter",
            registerConverters(&pythonLineRankOrderFilter<2, float>),
            (arg("image"), arg("half_length"), arg("rank"), arg("axis")=(2-1), arg("out")=python::object()),
            "Convolution along a single dimension of a 2D scalar or multiband image. "
                    "'kernel' must be an instance of Kernel1D.\n"
                    "\n"
                    "For details see convolveMultiArrayOneDimension_ in the vigra C++ documentation.\n");

    def("lineRankOrderFilter",
            registerConverters(&pythonLineRankOrderFilter<3, float>),
            (arg("image"), arg("half_length"), arg("rank"), arg("axis")=(3-1), arg("out")=python::object()),
            "Convolution along a single dimension of a 2D scalar or multiband image. "
                    "'kernel' must be an instance of Kernel1D.\n"
                    "\n"
                    "For details see convolveMultiArrayOneDimension_ in the vigra C++ documentation.\n");

    def("lineRankOrderFilter",
            registerConverters(&pythonLineRankOrderFilter<4, float>),
            (arg("image"), arg("half_length"), arg("rank"), arg("axis")=(4-1), arg("out")=python::object()),
            "Convolution along a single dimension of a 2D scalar or multiband image. "
                    "'kernel' must be an instance of Kernel1D.\n"
                    "\n"
                    "For details see convolveMultiArrayOneDimension_ in the vigra C++ documentation.\n");

    def("lineRankOrderFilter",
            registerConverters(&pythonLineRankOrderFilter<5, float>),
            (arg("image"), arg("half_length"), arg("rank"), arg("axis")=(5-1), arg("out")=python::object()),
            "Convolution along a single dimension of a 2D scalar or multiband image. "
                    "'kernel' must be an instance of Kernel1D.\n"
                    "\n"
                    "For details see convolveMultiArrayOneDimension_ in the vigra C++ documentation.\n");

    def("lineRankOrderFilter",
            registerConverters(&pythonLineRankOrderFilter<1, double>),
            (arg("image"), arg("half_length"), arg("rank"), arg("axis")=(1-1), arg("out")=python::object()),
            "Convolution along a single dimension of a 2D scalar or multiband image. "
                    "'kernel' must be an instance of Kernel1D.\n"
                    "\n"
                    "For details see convolveMultiArrayOneDimension_ in the vigra C++ documentation.\n");

    def("lineRankOrderFilter",
            registerConverters(&pythonLineRankOrderFilter<2, double>),
            (arg("image"), arg("half_length"), arg("rank"), arg("axis")=(2-1), arg("out")=python::object()),
            "Convolution along a single dimension of a 2D scalar or multiband image. "
                    "'kernel' must be an instance of Kernel1D.\n"
                    "\n"
                    "For details see convolveMultiArrayOneDimension_ in the vigra C++ documentation.\n");

    def("lineRankOrderFilter",
            registerConverters(&pythonLineRankOrderFilter<3, double>),
            (arg("image"), arg("half_length"), arg("rank"), arg("axis")=(3-1), arg("out")=python::object()),
            "Convolution along a single dimension of a 2D scalar or multiband image. "
                    "'kernel' must be an instance of Kernel1D.\n"
                    "\n"
                    "For details see convolveMultiArrayOneDimension_ in the vigra C++ documentation.\n");

    def("lineRankOrderFilter",
            registerConverters(&pythonLineRankOrderFilter<4, double>),
            (arg("image"), arg("half_length"), arg("rank"), arg("axis")=(4-1), arg("out")=python::object()),
            "Convolution along a single dimension of a 2D scalar or multiband image. "
                    "'kernel' must be an instance of Kernel1D.\n"
                    "\n"
                    "For details see convolveMultiArrayOneDimension_ in the vigra C++ documentation.\n");

    def("lineRankOrderFilter",
            registerConverters(&pythonLineRankOrderFilter<5, double>),
            (arg("image"), arg("half_length"), arg("rank"), arg("axis")=(5-1), arg("out")=python::object()),
            "Convolution along a single dimension of a 2D scalar or multiband image. "
                    "'kernel' must be an instance of Kernel1D.\n"
                    "\n"
                    "For details see convolveMultiArrayOneDimension_ in the vigra C++ documentation.\n");
}

void defineKernels();
void defineConvolutionFunctions();
void defineMorphology();
void defineTensor();
void defineNonLocalMean();

} // namespace vigra

using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(filters)
{
    import_vigranumpy();
    defineFilters2D();
    defineKernels();
    defineConvolutionFunctions();
    defineRankFilter();
    defineMorphology();
    defineTensor();
    defineNonLocalMean();
}
