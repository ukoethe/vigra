/************************************************************************/
/*                                                                      */
/*      Copyright 2011 by Ullrich Koethe and Kevin Kiefer               */
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

#define PY_ARRAY_UNIQUE_SYMBOL vigranumpyblockwise_PyArray_API
//#define NO_IMPORT_ARRAY

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/multi_blocking.hxx>
#include <vigra/multi_blockwise.hxx>
#include <vigra/blockwise_labeling.hxx>
#include <vigra/blockwise_watersheds.hxx>


namespace python = boost::python;


namespace vigra{


#define VIGRA_NUMPY_FILTERS(FUNCTOR, FUNCTION)      \
template<unsigned int DIM, class T_IN, class T_OUT> \
NumpyAnyArray FUNCTOR(                              \
    const NumpyArray<DIM, T_IN> &  source,          \
    const BlockwiseConvolutionOptions<DIM>  & opt,  \
    NumpyArray<DIM, T_OUT> dest                     \
){                                                  \
    dest.reshapeIfEmpty(source.shape());            \
    FUNCTION(source, dest, opt);                    \
    return dest;                                    \
}

VIGRA_NUMPY_FILTERS(pyGaussianSmooth,                     gaussianSmoothMultiArray)
VIGRA_NUMPY_FILTERS(pyGaussianGradient,                   gaussianGradientMultiArray)
VIGRA_NUMPY_FILTERS(pyGaussianGradientMagnitude,          gaussianGradientMagnitudeMultiArray)
VIGRA_NUMPY_FILTERS(pyGaussianDivergence,                 gaussianDivergenceMultiArray)
VIGRA_NUMPY_FILTERS(pyHessianOfGaussian,                  hessianOfGaussianMultiArray)
VIGRA_NUMPY_FILTERS(pyHessianOfGaussianEigenvalues,       hessianOfGaussianEigenvaluesMultiArray)
VIGRA_NUMPY_FILTERS(pyHessianOfGaussianFirstEigenvalue,   hessianOfGaussianFirstEigenvalueMultiArray)
VIGRA_NUMPY_FILTERS(pyHessianOfGaussianLastEigenvalue,    hessianOfGaussianLastEigenvalueMultiArray)
VIGRA_NUMPY_FILTERS(pyLaplacianOfGaussian,                laplacianOfGaussianMultiArray)
VIGRA_NUMPY_FILTERS(pySymmetricGradient,                  symmetricGradientMultiArray)
VIGRA_NUMPY_FILTERS(pyStructureTensor,                    structureTensorMultiArray)

#undef VIGRA_NUMPY_FILTERS



#define VIGRA_NUMPY_FILTER_BINDINGS(NAME_STR, FUNCTOR, T_IN, T_OUT) \
python::def(NAME_STR, registerConverters(&FUNCTOR<N,T_IN,T_OUT>),   \
    (                                                               \
        python::arg("source"),                                      \
        python::arg("options"),                                     \
        python::arg("out") = python::object()                       \
    )                                                               \
);

template<unsigned int N, class T_IN, class T_OUT>
void defineBlockwiseFilters(){
    typedef TinyVector<T_IN,N> in_type;
    typedef TinyVector<T_OUT,N> out_type;
    typedef TinyVector<T_OUT,N*(N+1)/2> out_type_2;

    VIGRA_NUMPY_FILTER_BINDINGS("_gaussianSmooth",                  pyGaussianSmooth,                   T_IN,    T_OUT)
    VIGRA_NUMPY_FILTER_BINDINGS("_gaussianGradient",                pyGaussianGradient,                 T_IN,    out_type)
    VIGRA_NUMPY_FILTER_BINDINGS("_gaussianGradientMagnitude",       pyGaussianGradientMagnitude,        T_IN,    T_OUT)
    VIGRA_NUMPY_FILTER_BINDINGS("_gaussianDivergence",              pyGaussianDivergence,               in_type, T_OUT)
    VIGRA_NUMPY_FILTER_BINDINGS("_hessianOfGaussian",               pyHessianOfGaussian,                T_IN,    out_type_2)
    VIGRA_NUMPY_FILTER_BINDINGS("_hessianOfGaussianEigenvalues",    pyHessianOfGaussianEigenvalues,     T_IN,    out_type)
    VIGRA_NUMPY_FILTER_BINDINGS("_hessianOfGaussianFirstEigenvalue",pyHessianOfGaussianFirstEigenvalue, T_IN,    T_OUT)
    VIGRA_NUMPY_FILTER_BINDINGS("_hessianOfGaussianLastEigenvalue", pyHessianOfGaussianLastEigenvalue,  T_IN,    T_OUT)
    VIGRA_NUMPY_FILTER_BINDINGS("_laplacianOfGaussian",             pyLaplacianOfGaussian,              T_IN,    T_OUT)
    VIGRA_NUMPY_FILTER_BINDINGS("_symmetricGradient",               pySymmetricGradient,                T_IN,    out_type)
    VIGRA_NUMPY_FILTER_BINDINGS("_structureTensor",                 pyStructureTensor,                  T_IN,    out_type_2)
}

#undef VIGRA_NUMPY_FILTER_BINDINGS



template<unsigned int N, class T_IN, class T_OUT>
python::tuple pyUnionFindWatersheds(
    const NumpyArray<N, T_IN> & data,
    const BlockwiseLabelOptions & opt,
    NumpyArray<N, T_OUT> labels
){
    labels.reshapeIfEmpty(data.shape());
    auto res = unionFindWatershedsBlockwise(data, labels, opt);
    return python::make_tuple(labels, res);
}

template<unsigned int N, class Data, class Label>
void defineUnionFindWatershedsImpl()
{
    python::def("_unionFindWatersheds", registerConverters(&pyUnionFindWatersheds<N,Data,Label>),
        (
            python::arg("data"),
            python::arg("options"),
            python::arg("labels") = python::object()
        )
    );
}



template<unsigned int N, class T_IN, class T_OUT>
python::tuple pyLabelArray(
    const NumpyArray<N, T_IN> & data,
    const BlockwiseLabelOptions & opt,
    NumpyArray<N, T_OUT> labels
){
    labels.reshapeIfEmpty(data.shape());
    auto res = labelMultiArrayBlockwise(data, labels, opt);
    return python::make_tuple(labels, res);
}

template<unsigned int N, class Data, class Label>
void defineLabelArrayImpl()
{
    python::def("_labelArray", registerConverters(&pyLabelArray<N,Data,Label>),
        (
            python::arg("data"),
            python::arg("options"),
            python::arg("labels") = python::object()
        )
    );
}



template<class  MB>
NumpyAnyArray intersectingBlocks(
    const MB & mb,
    const typename MB::Shape begin,
    const typename MB::Shape end,
    NumpyArray<1, UInt32> out
){
    std::vector<UInt32> outVec = mb.intersectingBlocks(begin,end);
    out.reshapeIfEmpty(typename NumpyArray<1,UInt32>::difference_type(outVec.size()));
    std::copy(outVec.begin(),outVec.end(), out.begin());
    return out;
}

template<class  MB>
python::tuple getBlock(
    const MB & mb,
    const UInt32 blockIndex
){
    const auto iter = mb.blockBegin();
    const auto & block = iter[blockIndex];
    auto tl = block.begin();
    auto br = block.end();
    return python::make_tuple(tl,br);
}


template<class  MB>
python::tuple getBlock2(
    const MB & mb,
    const typename  MB::BlockDesc desc
){
    const auto block = mb.blockDescToBlock(desc);
    auto tl = block.begin();
    auto br = block.end();
    return python::make_tuple(tl,br);
}

template<class BLOCK>
typename BLOCK::Vector
blockBegin(const BLOCK & b){
    return b.begin();
}
template<class BLOCK>
typename BLOCK::Vector
blockEnd(const BLOCK & b){
    return b.end();
}

template<class BLOCK>
typename BLOCK::Vector
blockShape(const BLOCK & b){
    return b.size();
}


template<unsigned int DIM>
void defineMultiBlocking(const std::string & clsName){

    typedef MultiBlocking<DIM> Blocking;
    typedef typename Blocking::Shape Shape;
    typedef typename Blocking::Block Block;

    python::class_<Blocking>(clsName.c_str(), python::init<const Shape &, const Shape &>())
        .def("intersectingBlocks",registerConverters(&intersectingBlocks<Blocking>),
            (
                python::arg("begin"),
                python::arg("end"),
                python::arg("out") = python::object()
            )
        )
        .def("__len__", &Blocking::numBlocks)
        .def("__getitem__", &getBlock<Blocking>)
        .def("__getitem__", &getBlock2<Blocking>)
    ;

    const std::string blockName = clsName + std::string("Block");

    python::class_<Block>(blockName.c_str())
        .add_property("begin",&blockBegin<Block>)
        .add_property("end",  &blockEnd<Block>)
        .add_property("shape",&blockShape<Block>)
    ;
}



template<unsigned int DIM>
void defineBlockwiseConvolutionOptions(const std::string & clsName){

    typedef BlockwiseConvolutionOptions<DIM> Opt;
    python::class_<Opt>(clsName.c_str(), python::init<>())
    .add_property("stdDev", &Opt::getStdDev, &Opt::setStdDev)
    //.add_property("scale", &Opt::getScale, &Opt::setScale)
    .add_property("innerScale", &Opt::getInnerScale,  &Opt::setInnerScale)
    .add_property("outerScale", &Opt::getOuterScale,  &Opt::setOuterScale)
    .add_property("blockShape", &Opt::readBlockShape, &Opt::setBlockShape)
    .add_property("numThreads", &Opt::getNumThreads,  &Opt::setNumThreads)
    ;
}


void defineNeighborhoodType()
{
    python::enum_<NeighborhoodType>("NeighborhoodType")
    .value("DirectNeighborhood", DirectNeighborhood)
    .value("IndirectNeighborhood", IndirectNeighborhood)
    .export_values()
    ;
}


void defineBlockwiseLabelOptions()
{
    typedef BlockwiseLabelOptions Opt;

    python::class_<Opt>("BlockwiseLabelOptions", python::init<>())
    .add_property("blockShape", &Opt::readBlockShape, &Opt::setBlockShape)
    .add_property("numThreads", &Opt::getNumThreads, &Opt::setNumThreads)
    .add_property("backgroundValue", &Opt::getBackgroundValue<double>,
            python::make_function(&Opt::ignoreBackgroundValue<double>, python::return_internal_reference<>()))
    .add_property("neighbourhood", &Opt::getNeighborhood,
            python::make_function(&Opt::neighborhood, python::return_internal_reference<>()))
    .def("hasBackgroundValue", &Opt::hasBackgroundValue)
    ;
}

}
using namespace vigra;
using namespace boost::python;







BOOST_PYTHON_MODULE_INIT(blockwise)
{
    import_vigranumpy();

    python::docstring_options doc_options(true, true, false);

    defineMultiBlocking<2>("Blocking2D");
    defineMultiBlocking<3>("Blocking3D");

    defineBlockwiseConvolutionOptions<2>("BlockwiseConvolutionOptions2D");
    defineBlockwiseConvolutionOptions<3>("BlockwiseConvolutionOptions3D");
    defineBlockwiseConvolutionOptions<4>("BlockwiseConvolutionOptions4D");
    defineBlockwiseConvolutionOptions<5>("BlockwiseConvolutionOptions5D");

    defineNeighborhoodType();
    defineBlockwiseLabelOptions();

    defineBlockwiseFilters<2, npy_float32, npy_float32>();
    defineBlockwiseFilters<3, npy_float32, npy_float32>();

    defineUnionFindWatershedsImpl<2, npy_uint8, npy_uint32>();
    defineUnionFindWatershedsImpl<3, npy_uint8, npy_uint32>();
    defineUnionFindWatershedsImpl<2, npy_uint32, npy_uint32>();
    defineUnionFindWatershedsImpl<3, npy_uint32, npy_uint32>();

    defineLabelArrayImpl<2, npy_uint32, npy_uint32>();
    defineLabelArrayImpl<3, npy_uint32, npy_uint32>();
}
