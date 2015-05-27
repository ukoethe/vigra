/************************************************************************/
/*                                                                      */
/*                 Copyright 2011 by Ullrich Koethe                     */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyhistogram_PyArray_API
//#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/multi_histogram.hxx>
namespace python = boost::python;

namespace vigra{


    template<unsigned int DIM,unsigned int CHANNELS>
    NumpyAnyArray pyMultiGaussianHistogram(
        NumpyArray<DIM,TinyVector<float,CHANNELS> > image,
        const TinyVector<float,CHANNELS> minVals,
        const TinyVector<float,CHANNELS> maxVals,
        const size_t bins,
        const float  sigma,
        const float sigmaBin,
        NumpyArray<DIM+2,float> histogram = NumpyArray<DIM+2,float>()
    ){
        typename NumpyArray<DIM+2,float>::difference_type  outShape;
        for(size_t d=0;d<DIM;++d){
            outShape[d]=image.shape(d);
        }
        outShape[DIM]=bins;
        outShape[DIM+1]=CHANNELS;
        histogram.reshapeIfEmpty(outShape);
        {
             PyAllowThreads _pythread;
            multi_gaussian_histogram<DIM,float,CHANNELS,float>(image,minVals,maxVals,bins,
                sigma,sigmaBin,histogram);
        }
        return histogram;
    }


    template<unsigned int DIM>
    NumpyAnyArray pyMultiGaussianCoHistogram(
        NumpyArray<DIM,float > imageA,
        NumpyArray<DIM,float > imageB,
        const TinyVector<float,2> minVals,
        const TinyVector<float,2> maxVals,
        const TinyVector<int,2> bins,
        const TinyVector<float,3> sigma,
        NumpyArray<DIM+2,float> histogram = NumpyArray<DIM+2,float>()
    ){
        typename NumpyArray<DIM+2,float>::difference_type  outShape;
        for(size_t d=0;d<DIM;++d){
           outShape[d]=imageA.shape(d);
        }
        outShape[DIM]=bins[0];
        outShape[DIM+1]=bins[1];
        histogram.reshapeIfEmpty(outShape);
        {
            PyAllowThreads _pythread;
            multi_gaussian_co_histogram<DIM,float,float>(imageA,imageB,minVals,maxVals,bins,
            sigma,histogram);
        }
        return histogram;
    }




    template<unsigned int DIM>
    NumpyAnyArray pyMultiGaussianRankOrder(
        const NumpyArray<DIM, float>      & image,
        const float                         minVal,
        const float                         maxVal,
        const size_t                        bins,
        const NumpyArray<1, float>        & sigmas,  
        const NumpyArray<1, float>        & ranks,
        NumpyArray<DIM+1, float>            out
    ){
        // reshape 
        typename NumpyArray<DIM+1,float>::difference_type  outShape;
        for(size_t d=0;d<DIM;++d){
           outShape[d]=image.shape(d);
        }
        outShape[DIM]=ranks.size();
        out.reshapeIfEmpty(outShape);

        {
             PyAllowThreads _pythread;
            TinyVector<double, DIM+1> sigmaVec;
            std::copy(sigmas.begin(),sigmas.end(),sigmaVec.begin());

            multi_gaussian_rank_order(image, minVal, maxVal,
                                bins, sigmaVec, ranks, out);
        }
        return out;
    }



    template<unsigned int DIM,unsigned int CHANNELS>
    void defineMultiGaussianHistogram(){

        python::def("gaussianHistogram_",registerConverters(&pyMultiGaussianHistogram<DIM,CHANNELS>),
            (
                python::arg("image"),
                python::arg("minVals"),
                python::arg("maxVals"),
                python::arg("bins")=30,
                python::arg("sigma")=3.0,
                python::arg("sigmaBin")=2.0,
                python::arg("out")=python::object()
            )
        );
    }


    template<unsigned int DIM>
    void defineMultiGaussianCoHistogram(){

        python::def("gaussianCoHistogram",registerConverters(&pyMultiGaussianCoHistogram<DIM>),
            (
                python::arg("imageA"),
                python::arg("imageB"),
                python::arg("minVals"),
                python::arg("maxVals"),
                python::arg("bins"),
                python::arg("sigma"),
                python::arg("out")=python::object()
            )
        );
    }


    template<unsigned int DIM>
    void defineMultiGaussianRank(){

        python::def("_gaussianRankOrder",registerConverters(&pyMultiGaussianRankOrder<DIM>),
            (
                python::arg("image"),
                python::arg("minVal"),
                python::arg("maxVal"),
                python::arg("bins"),
                python::arg("sigmas"),
                python::arg("ranks"),
                python::arg("out")=python::object()
            )
        );
    }
    


} // namespace vigra

using namespace vigra;
using namespace boost::python;




BOOST_PYTHON_MODULE_INIT(histogram)
{
    import_vigranumpy();

    // all exporters needed for graph exporters (like lemon::INVALID)
    defineMultiGaussianHistogram<2,3>();
    defineMultiGaussianHistogram<3,1>();
    defineMultiGaussianHistogram<3,3>();
    defineMultiGaussianHistogram<3,10>();

    defineMultiGaussianCoHistogram<2>();
    defineMultiGaussianCoHistogram<3>();

    defineMultiGaussianRank<2>();
    defineMultiGaussianRank<3>();
}
