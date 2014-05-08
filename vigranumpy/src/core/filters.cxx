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
#include <vigra/symmetry.hxx>
#include <vigra/non_local_mean.hxx>

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

    def("radialSymmetryTransform2D",
        registerConverters(&pythonRadialSymmetryTransform2D<float>),
        (arg("image"), arg("scale"),arg("out")=python::object()),
        "Find centers of radial symmetry in an 2D image.\n\n"
        "This algorithm implements the Fast Radial Symmetry Transform according to "
        "[G. Loy, A. Zelinsky: \"A Fast Radial Symmetry Transform for Detecting Points of Interest\", "
        "in: A. Heyden et al. (Eds.): Proc. of 7th European Conf. on Computer Vision, Part 1, pp. 358-368, Springer LNCS 2350, 2002]\n\n"
        "For details see radialSymmetryTransform_ in the vigra C++ documentation.\n");
}

template<int DIM,class PIXEL_TYPE>
NumpyAnyArray  pyNonLocalMean(
    NumpyArray<DIM,PIXEL_TYPE> image,
    const float sigma,
    const int searchRadius,
    const int patchRadius,
    const bool gaussNoise,
    const double sigmaMean,
    const int nThreads,
    const double epsilon,
    const double mu1,
    const double var1,
    const int stepSize,
    const bool verbose,
    NumpyArray<DIM,PIXEL_TYPE> out = NumpyArray<DIM,PIXEL_TYPE>()
){
    NonLocalMeanParameter param;

    param.sigma_=sigma;
    param.searchRadius_=searchRadius;
    param.patchRadius_=patchRadius;
    param.gaussNoise_=gaussNoise;
    param.sigmaMean_=sigmaMean;
    param.nThreads_ = nThreads;
    param.epsilon_=epsilon;
    param.mu1_=mu1;
    param.var1_=var1;
    param.stepSize_=stepSize;
    param.verbose_=verbose;
    out.reshapeIfEmpty(image.shape());
    nonLocalMean<DIM,PIXEL_TYPE>(image,param,out);
    return out;
}

template<int DIM,class PIXEL_TYPE>
void exportNonLocalMean(const std::string name){

    // export the function to python
    python::def(name.c_str(), registerConverters(&pyNonLocalMean<DIM,PIXEL_TYPE>) ,
        (
            python::arg("image"),
            python::arg("sigma")=1.0,
            python::arg("searchRadius")=3,
            python::arg("patchRadius")=1,
            python::arg("gaussNoise")=true,
            python::arg("sigmaMean")=1.0,
            python::arg("nThreads")=8,
            python::arg("epsilon")=0.00001,
            python::arg("mu1")=0.95,
            python::arg("var1")=0.5,
            python::arg("stepSize")=2,
            python::arg("verbose")=true,
            python::arg("out") = boost::python::object()
        ),
        "loop over an image and do something with each pixels\n\n"
        "Args:\n\n"
        "   image : input image\n\n"
        "returns an an image with the same shape as the input image"
    );
}

void defineNonLocalMean(){
    exportNonLocalMean<2,vigra::TinyVector<float,3> >("nonLocalMean2d");
    exportNonLocalMean<2,float>("nonLocalMean2d");
    exportNonLocalMean<3,float>("nonLocalMean3d");
    exportNonLocalMean<4,float>("nonLocalMean4d");
}

void defineKernels();
void defineConvolutionFunctions();
void defineMorphology();
void defineTensor();

} // namespace vigra

using namespace vigra;
using namespace boost::python;

BOOST_PYTHON_MODULE_INIT(filters)
{
    import_vigranumpy();
    defineFilters2D();
    defineKernels();
    defineConvolutionFunctions();
    defineMorphology();
    defineTensor();
    defineNonLocalMean();
}
