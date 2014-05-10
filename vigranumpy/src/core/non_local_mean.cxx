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
#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/non_local_mean.hxx>

namespace python = boost::python;

namespace vigra
{






template<int DIM,class PIXEL_TYPE,class SMOOTH_POLICY>
NumpyAnyArray  pyNonLocalMean(
    NumpyArray<DIM,PIXEL_TYPE> image,
    const typename SMOOTH_POLICY::ParameterType & policyParam,
    const double sigmaSpatial,
    const int searchRadius,
    const int patchRadius,
    const double sigmaMean,
    const int stepSize,
    const int iterations,
    const int nThreads,
    const bool verbose,
    NumpyArray<DIM,PIXEL_TYPE> out = NumpyArray<DIM,PIXEL_TYPE>()
){

    SMOOTH_POLICY smoothPolicy(policyParam);
    NonLocalMeanParameter param;
    param.sigmaSpatial_=sigmaSpatial;
    param.searchRadius_=searchRadius;
    param.patchRadius_=patchRadius;
    param.sigmaMean_=sigmaMean;
    param.stepSize_=stepSize;
    param.iterations_=iterations;
    param.nThreads_ = nThreads;
    param.verbose_=verbose;
    out.reshapeIfEmpty(image.shape());
    nonLocalMean<DIM,PIXEL_TYPE>(image,smoothPolicy,param,out);
    return out;
}


void exportNonLocalMeanPolicyParameterObjects(){

    {
        typedef RatioPolicyParameter ParamType;

        python::class_<ParamType>(
            "RatioPolicy",
            python::init<const double,const double,const double,const double>(
                (
                    python::arg("sigma"),
                    python::arg("meanRatio")=0.95,
                    python::arg("varRatio")=0.5,
                    python::arg("epsilon")=0.00001
                )
            )
        )
        .def_readwrite("sigma", &ParamType::sigma_)
        .def_readwrite("meanRatio", &ParamType::meanRatio_)
        .def_readwrite("varRatio", &ParamType::varRatio_)
        .def_readwrite("epsilon", &ParamType::epsilon_)
        ;
            
    }

    {
        typedef NormPolicyParameter ParamType;

        python::class_<ParamType>(
            "NormPolicy",
            python::init<const double,const double,const double>(
                (
                    python::arg("sigma"),
                    python::arg("meanDist"),
                    python::arg("varRatio")
                )
            )
        )
        .def_readwrite("sigma", &ParamType::sigma_)
        .def_readwrite("meanDist", &ParamType::meanDist_)
        .def_readwrite("varRatio", &ParamType::varRatio_)
        ;
            
    }


}




template<int DIM,class PIXEL_TYPE, class POLICY>
void exportNonLocalMean(const std::string name){

    typedef POLICY SmoothPolicyType;
    typedef typename SmoothPolicyType::ParameterType SmoothPolicyParameterType;
    // export the function to python
    python::def(name.c_str(), registerConverters(&pyNonLocalMean<DIM,PIXEL_TYPE,SmoothPolicyType>) ,
        (
            python::arg("image"),
            python::arg("policy"),
            python::arg("sigmaSpatial")=2.0,
            python::arg("searchRadius")=3,
            python::arg("patchRadius")=1,
            python::arg("sigmaMean")=1.0,
            python::arg("stepSize")=2,
            python::arg("iterations")=1,
            python::arg("nThreads")=8,
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
    using namespace python;
    docstring_options doc_options(true, true, false);

    // export different parameter objects
    exportNonLocalMeanPolicyParameterObjects();

    {
        exportNonLocalMean<2,TinyVector<float,3>, RatioPolicy<TinyVector<float,3> > >("nonLocalMean2d");
        exportNonLocalMean<2,float, RatioPolicy<float> >("nonLocalMean2d");
        exportNonLocalMean<3,float, RatioPolicy<float> >("nonLocalMean3d");
        exportNonLocalMean<4,float, RatioPolicy<float> >("nonLocalMean4d");
    }
    {
        exportNonLocalMean<2,TinyVector<float,3>, NormPolicy<TinyVector<float,3> > >("nonLocalMean2d");
        exportNonLocalMean<2,float, NormPolicy<float> >("nonLocalMean2d");
        exportNonLocalMean<3,float, NormPolicy<float> >("nonLocalMean3d");
        exportNonLocalMean<4,float, NormPolicy<float> >("nonLocalMean4d");
    }
}


} // namespace vigra
