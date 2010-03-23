
/************************************************************************/
/*                                                                      */
/*                 Copyright 2009 by Ullrich Koethe                     */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpynoise_PyArray_API
//#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/noise_normalization.hxx>

#include <cmath>

#include "tensors.hxx"
namespace python = boost::python;

namespace vigra
{

//TODO
NumpyAnyArray vectorToArray(std::vector< TinyVector< double, 2 > > & result)
{
    NumpyArray<2,double> res(MultiArrayShape<2>::type(result.size(),2));

    for(size_t ii=0;ii<result.size();++ii)
    {
        res(ii,0)=result[ii][0];
        res(ii,1)=result[ii][1];
    }
    return res;
}

template < class PixelType >
NumpyAnyArray pythonNoiseVarianceEstimation(NumpyArray<2, Singleband<PixelType> > image,
                                            bool useGradient=true, unsigned int windowRadius=6,
                                            unsigned int clusterCount=10, double averagingQuantile=0.8,
                                            double noiseEstimationQuantile=1.5, double noiseVarianceInitialGuess=10.0,
                                            NumpyArray<3, Multiband<PixelType> > res=python::object())
{
    NoiseNormalizationOptions noiseNormalizationOptions;
    noiseNormalizationOptions.useGradient(useGradient).windowRadius(windowRadius).clusterCount(clusterCount).averagingQuantile(averagingQuantile).noiseEstimationQuantile(noiseEstimationQuantile).noiseVarianceInitialGuess(noiseVarianceInitialGuess);
    std::vector< TinyVector< double, 2 > > result;

    noiseVarianceEstimation(srcImageRange(image), result, noiseNormalizationOptions);

    return vectorToArray(result);
}

template < class PixelType >
NumpyAnyArray pythonNoiseVarianceClustering(NumpyArray<2, Singleband<PixelType> > image,
                                            bool useGradient=true, unsigned int windowRadius=6,
                                            unsigned int clusterCount=10, double averagingQuantile=0.8,
                                            double noiseEstimationQuantile=1.5, double noiseVarianceInitialGuess=10.0,
                                            NumpyArray<3, Multiband<PixelType> > res=python::object())
{
    NoiseNormalizationOptions noiseNormalizationOptions;
    noiseNormalizationOptions.useGradient(useGradient).windowRadius(windowRadius).clusterCount(clusterCount).averagingQuantile(averagingQuantile).noiseEstimationQuantile(noiseEstimationQuantile).noiseVarianceInitialGuess(noiseVarianceInitialGuess);
    std::vector< TinyVector< double, 2 > > result;
    noiseVarianceClustering(srcImageRange(image), result,
        noiseNormalizationOptions);
    return vectorToArray(result);
}

template < class PixelType >
NumpyAnyArray pythonNonparametricNoiseNormalization(NumpyArray<3, Multiband<PixelType> > image,
                                                    bool useGradient=true, unsigned int windowRadius=6,
                                                    unsigned int clusterCount=10, double averagingQuantile=0.8,
                                                    double noiseEstimationQuantile=1.5, double noiseVarianceInitialGuess=10.0,
                                                    NumpyArray<3, Multiband<PixelType> > res=python::object())
{
    NoiseNormalizationOptions noiseNormalizationOptions;
    noiseNormalizationOptions.useGradient(useGradient).windowRadius(windowRadius).clusterCount(clusterCount).averagingQuantile(averagingQuantile).noiseEstimationQuantile(noiseEstimationQuantile).noiseVarianceInitialGuess(noiseVarianceInitialGuess);
    
    res.reshapeIfEmpty(image.shape(),"nonparametricNoiseNormalization(): Output images has wrong dimensions");
    
    for(int k=0;k<image.shape(2);++k)
    {
        nonparametricNoiseNormalization(srcImageRange(image),
                                        destImage(res), noiseNormalizationOptions);
    }
    return res;
}

template < class PixelType >
NumpyAnyArray pythonQuadraticNoiseNormalizationEstimated( NumpyArray<3, Multiband<PixelType> > image,
                                                          bool useGradient=true, unsigned int windowRadius=6,
                                                          unsigned int clusterCount=10, double averagingQuantile=0.8,
                                                          double noiseEstimationQuantile=1.5, double noiseVarianceInitialGuess=10.0,
                                                          NumpyArray<3, Multiband<PixelType> > res=python::object())
{
    NoiseNormalizationOptions noiseNormalizationOptions;
    noiseNormalizationOptions.useGradient(useGradient).windowRadius(windowRadius).clusterCount(clusterCount).averagingQuantile(averagingQuantile).noiseEstimationQuantile(noiseEstimationQuantile).noiseVarianceInitialGuess(noiseVarianceInitialGuess);

    res.reshapeIfEmpty(image.shape(),"quadraticNoiseNormalizationEstimated(): Output images has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        quadraticNoiseNormalization(srcImageRange(bimage),
                                    destImage(bres), noiseNormalizationOptions);
    }
    return res;
}

template < class PixelType >
NumpyAnyArray pythonLinearNoiseNormalizationEstimated(NumpyArray<3, Multiband<PixelType> > image,
                                                      bool useGradient=true, unsigned int windowRadius=6,
                                                      unsigned int clusterCount=10, double averagingQuantile=0.8,
                                                      double noiseEstimationQuantile=1.5, double noiseVarianceInitialGuess=10.0,
                                                      NumpyArray<3, Multiband<PixelType> > res=python::object())
{
    NoiseNormalizationOptions noiseNormalizationOptions;
    noiseNormalizationOptions.useGradient(useGradient).windowRadius(windowRadius).clusterCount(clusterCount).averagingQuantile(averagingQuantile).noiseEstimationQuantile(noiseEstimationQuantile).noiseVarianceInitialGuess(noiseVarianceInitialGuess);
    
    res.reshapeIfEmpty(image.shape(),"linearNoiseNormalizationEstimated(): Output images has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        linearNoiseNormalization(srcImageRange(bimage),
                                 destImage(bres), noiseNormalizationOptions);
    }
    return res;
}


template < class PixelType >
NumpyAnyArray pythonQuadraticNoiseNormalization(NumpyArray<3, Multiband<PixelType> > image,
                                                double a0, double a1, double a2,
                                                NumpyArray<3, Multiband<PixelType> > res=python::object())
{
    res.reshapeIfEmpty(image.shape(),"quadraticNoiseNormalization(): Output images has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        quadraticNoiseNormalization(srcImageRange(bimage), destImage(bres),a0, a1, a2);
    }
    
    return res;
}

template < class PixelType >
NumpyAnyArray pythonLinearNoiseNormalization(NumpyArray<3, Multiband<PixelType> > image,
                                             double a0, double a1, NumpyArray<3, Multiband<PixelType> > res)
{
    res.reshapeIfEmpty(image.shape(),"linearNoiseNormalization(): Output images has wrong dimensions");

    for(int k=0;k<image.shape(2);++k)
    {
        MultiArrayView<2, PixelType, StridedArrayTag> bimage = image.bindOuter(k);
        MultiArrayView<2, PixelType, StridedArrayTag> bres = res.bindOuter(k);
        linearNoiseNormalization(srcImageRange(bimage), destImage(bres),a0, a1);
    } 
    return res;
}


void defineNoise()
{
    using namespace python;

    def("noiseVarianceEstimation",
        registerConverters(&pythonNoiseVarianceEstimation<float>),
        (arg("image"), arg("useGradient")=true, arg("windowRadius")=6,
        arg("clusterCount")=10, arg("averagingQuantile")=0.8,
        arg("noiseEstimationQuantile")=1.5,
        arg("noiseVarianceInitialGuess")=10.0, arg("out")=object()),
        "Determine the noise variance as a function of the image intensity.\n"
        "\n"
        "Returns an array with the means in the first column and the variances in the second column.\n"
        "Since the length of the resulting array is not known beforhand, it can not be written into an preallocated array\n"
        "(the \"out\" argument in most other vigra python functions.\n\n"
        "For details see the vigra documentation noiseVarianceEstimation_.\n"
        );

    def("noiseVarianceClustering",
        registerConverters(&pythonNoiseVarianceClustering<float>),
        (arg("image"), arg("useGradient")=true, arg("windowRadius")=6,
        arg("clusterCount")=10, arg("averagingQuantile")=0.8,
        arg("noiseEstimationQuantile")=1.5,
        arg("noiseVarianceInitialGuess")=10.0, arg("out")=object()),
        "Determine the noise variance as a function of the image intensity and cluster the results.\n"
        "This operator first calls noiseVarianceEstimation() to obtain a sequence of intensity/variance pairs,\n"
        "which are then clustered using the median cut algorithm. Then the cluster centers (i.e. average variance vs. average intensity)\n"
        "are determined and returned in the result sequence.\n"
        "\n"
        "Since the length of the resulting array is not known beforhand, it cannot be written into an preallocated array\n"
        "(the \"out\" argument in most other vigra python functions)\n.\n"
        "For details see the vigra documentation noiseVarianceClustering_.\n"
        );

    def("nonparametricNoiseNormalization",
        registerConverters(&pythonNonparametricNoiseNormalization<float>),    // also multiband
        (arg("image"), arg("useGradient")=true, arg("windowRadius")=6,
        arg("clusterCount")=10, arg("averagingQuantile")=0.8,
        arg("noiseEstimationQuantile")=1.5,
        arg("noiseVarianceInitialGuess")=10.0, arg("out")=object()),
        "Noise normalization by means of an estimated non-parametric noise model.\n\n"
        "For details see nonparametricNoiseNormalization_ in the vigra C++ documentation.\n");

    def("quadraticNoiseNormalizationEstimated",
        registerConverters(&pythonQuadraticNoiseNormalizationEstimated<float>),    // also multiband
        (arg("image"), arg("useGradient")=true, arg("windowRadius")=6,
        arg("clusterCount")=10, arg("averagingQuantile")=0.8,
        arg("noiseEstimationQuantile")=1.5,
        arg("noiseVarianceInitialGuess")=10.0, arg("out")=object()));

    def("linearNoiseNormalizationEstimated",
        registerConverters(&pythonLinearNoiseNormalizationEstimated<float>),    // also multiband
        (arg("image"), arg("useGradient")=true, arg("windowRadius")=6,
        arg("clusterCount")=10, arg("averagingQuantile")=0.8,
        arg("noiseEstimationQuantile")=1.5,
        arg("noiseVarianceInitialGuess")=10.0, arg("out")=object()));

    def("quadraticNoiseNormalization",
        registerConverters(&pythonQuadraticNoiseNormalization<float>),    // also multiband
        (arg("image"), arg("a0"), arg("a1"), arg("a2"), arg("out")=object()),
        "Noise normalization by means of an estimated quadratic noise model.\n\n"
        "For details see quadraticNoiseNormalization_ in the vigra C++ documentation.\n"
        );

    def("linearNoiseNormalization",
        registerConverters(&pythonLinearNoiseNormalization<float>),    // also multiband
        (arg("image"), arg("a0"), arg("a1"), arg("out")=object()),
        "Noise normalization by means of an estimated linear noise model.\n\n"
        "For details see linearNoiseNormalization_ in the vigra C++ documentation.\n"
       );
}

} // namespace vigra

using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(noise)
{
    import_vigranumpy();
    defineNoise();
}
