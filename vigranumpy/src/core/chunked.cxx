#include <boost/python.hpp>
#include <vigra/python_utility.hxx>
#include <vigra/numpy_array.hxx>
#include <vigra/multi_array_chunked.hxx>
#include <vigra/multi_blockwise.hxx>
#include <vigra/multi_chunked.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/blockwise_watersheds.hxx>
#include <vigra/blockwise_labeling.hxx>
#include <vigra/error.hxx>
#include <string>


namespace python = boost::python;


namespace vigra {


namespace detail {

    template<class T>
    python::object makeOwningHolder(T* p)
    {
        typedef python::detail::make_owning_holder owning;
        auto nonzero = python_ptr::new_nonzero_reference;

        python_ptr ptr(python::to_python_indirect<T*,owning>()(p), nonzero);
        python::handle<> h(ptr.release());
        return python::object(h);
    }

    template<unsigned int N, class T1, class T2>
    python::object makeChunkedArray(const ChunkedArray<N,T1> & source)
    {
        std::string backend = source.backend();
        ChunkedArrayOptions opt;
        opt.cacheMax(source.cacheMaxSize());
        ChunkedArray<N,T2> * out = nullptr;

        if (backend == "ChunkedArrayFull") {
            out = new ChunkedArrayFull<N,T2>(source.shape(), opt);
        }
        else if (backend == "ChunkedArrayLazy") {
            out = new ChunkedArrayLazy<N,T2>(source.shape(), source.chunkShape(), opt);
        }
        else if (backend == "ChunkedArrayTmpFile") {
            out = new ChunkedArrayTmpFile<N,T2>(source.shape(), source.chunkShape(), opt);
        }
        else if (backend.find("ChunkedArrayCompressed") != std::string::npos) {
            const ChunkedArrayCompressed<N,T1> & sourceComp = dynamic_cast<const ChunkedArrayCompressed<N,T1> &>(source);
            opt.compression(sourceComp.compression_method_);
            out = new ChunkedArrayCompressed<N,T2>(sourceComp.shape(), sourceComp.chunkShape(), opt);
        }
        else if (backend.find("ChunkedArrayHDF5") != std::string::npos) {
            vigra_precondition(false, "The 'out' parameter is mandatory for ChunkedArrayHDF5");
        }
        else {
            vigra_fail("Unable to derive backend from 'input' parameter");
        }

        return makeOwningHolder(out);
    }

} // END NAMESPACE detail



#define PY_VIGRA_FILTERS(FUNCTOR, FUNCTION)                                                 \
template<unsigned int N, class T1, class T2>                                                \
python::object FUNCTOR(                                                                     \
        const ChunkedArray<N,T1> & source,                                                  \
        const BlockwiseConvolutionOptions<N> & opt,                                         \
        python::object out                                                                  \
){                                                                                          \
    if (out == python::object()) out = detail::makeChunkedArray<N,T1,T2>(source);           \
    ChunkedArray<N,T2> & dest = python::extract<ChunkedArray<N,T2>&>(out)();                \
    FUNCTION(source, dest, opt);                                                            \
\
    return out;                                                                             \
}

PY_VIGRA_FILTERS(pyGaussianSmooth,                  gaussianSmoothMultiArray)
PY_VIGRA_FILTERS(pyGaussianGradient,                gaussianGradientMultiArray)
PY_VIGRA_FILTERS(pyGaussianGradientMagnitude,       gaussianGradientMagnitudeMultiArray)
PY_VIGRA_FILTERS(pyGaussianDivergence,              gaussianDivergenceMultiArray)
PY_VIGRA_FILTERS(pyHessianOfGaussian,               hessianOfGaussianMultiArray)
PY_VIGRA_FILTERS(pyHessianOfGaussianEigenvalues,    hessianOfGaussianEigenvaluesMultiArray)
PY_VIGRA_FILTERS(pyHessianOfGaussianFirstEigenvalue,hessianOfGaussianFirstEigenvalueMultiArray)
PY_VIGRA_FILTERS(pyHessianOfGaussianLastEigenvalue, hessianOfGaussianLastEigenvalueMultiArray)
PY_VIGRA_FILTERS(pyLaplacianOfGaussian,             laplacianOfGaussianMultiArray)
PY_VIGRA_FILTERS(pySymmetricGradient,               symmetricGradientMultiArray)
PY_VIGRA_FILTERS(pyStructureTensor,                 structureTensorMultiArray)

#undef PY_VIGRA_CONVOLUTION


#define PY_VIGRA_BINDINGS(PY_NAME, FUNCTOR, T_IN, T_OUT) \
python::def(PY_NAME, &FUNCTOR<N,T_IN,T_OUT>,             \
    (                                                    \
        python::arg("source"),                           \
        python::arg("options"),                          \
        python::arg("out") = python::object()            \
    )                                                    \
);

template<unsigned int N, class T_IN, class T_OUT>
void defineChunkedFiltersImpl()
{
    typedef TinyVector<T_IN,N> in_type;
    typedef TinyVector<T_OUT,N> out_type;
    typedef TinyVector<T_OUT,N*(N+1)/2> out_type_2;

    PY_VIGRA_BINDINGS("_gaussianSmooth",                    pyGaussianSmooth, T_IN, T_OUT)
    PY_VIGRA_BINDINGS("_gaussianGradient",                  pyGaussianGradient, T_IN, out_type)
    PY_VIGRA_BINDINGS("_gaussianGradientMagnitude",         pyGaussianGradientMagnitude, T_IN, T_OUT)
    PY_VIGRA_BINDINGS("_gaussianDivergence",                pyGaussianDivergence, in_type, T_OUT)
    PY_VIGRA_BINDINGS("_hessianOfGaussian",                 pyHessianOfGaussian, T_IN, out_type_2)
    PY_VIGRA_BINDINGS("_hessianOfGaussianEigenvalues",      pyHessianOfGaussianEigenvalues, T_IN, out_type)
    PY_VIGRA_BINDINGS("_hessianOfGaussianFirstEigenvalue",  pyHessianOfGaussianFirstEigenvalue, T_IN, T_OUT)
    PY_VIGRA_BINDINGS("_hessianOfGaussianLastEigenvalue",   pyHessianOfGaussianLastEigenvalue, T_IN, T_OUT)
    PY_VIGRA_BINDINGS("_laplacianOfGaussian",               pyLaplacianOfGaussian, T_IN, T_OUT)
    PY_VIGRA_BINDINGS("_symmetricGradient",                 pySymmetricGradient, T_IN, out_type)
    PY_VIGRA_BINDINGS("_structureTensor",                   pyStructureTensor, T_IN, out_type_2)
}

#undef PY_VIGRA_BINDINGS



template<unsigned int N, class T1, class T2>
python::tuple pyUnionFindWatersheds(
        const ChunkedArray<N,T1> & source,
        const BlockwiseLabelOptions & opt,
        python::object out
){
    if (out == python::object()) out = detail::makeChunkedArray<N,T1,T2>(source);
    ChunkedArray<N,T2> & dest = python::extract<ChunkedArray<N,T2>&>(out)();
    auto res = unionFindWatershedsBlockwise(source, dest, opt);

    return python::make_tuple(out, res);
}

template<unsigned int N, class Data, class Label>
void defineChunkedWatershedsImpl()
{
    python::def("_unionFindWatersheds", &pyUnionFindWatersheds<N,Data,Label>,
        (
            python::arg("data"),
            python::arg("options"),
            python::arg("labels") = python::object()
        )
    );
}


template<unsigned int N, class T1, class T2>
python::tuple pyLabelArray(
        const ChunkedArray<N,T1> & source,
        const BlockwiseLabelOptions & opt,
        python::object out
){
    if (out == python::object()) out = detail::makeChunkedArray<N,T1,T2>(source);
    ChunkedArray<N,T2> & dest = python::extract<ChunkedArray<N,T2>&>(out)();
    auto res = labelMultiArrayBlockwise(source, dest, opt);

    return python::make_tuple(out, res);
}

template<unsigned int N, class Data, class Label>
void defineChunkedLabelImpl()
{
    python::def("_labelArray", &pyLabelArray<N,Data,Label>,
        (
            python::arg("data"),
            python::arg("options"),
            python::arg("labels") = python::object()
        )
    );
}



void defineChunkedFunctions()
{
    defineChunkedFiltersImpl<2, npy_float32, npy_float32>();
    defineChunkedFiltersImpl<3, npy_float32, npy_float32>();

    defineChunkedWatershedsImpl<2, npy_uint32, npy_uint32>();
    defineChunkedWatershedsImpl<3, npy_uint32, npy_uint32>();

    defineChunkedLabelImpl<2, npy_uint32, npy_uint32>();
    defineChunkedLabelImpl<3, npy_uint32, npy_uint32>();
}


} // END NAMESPACE VIGRA
