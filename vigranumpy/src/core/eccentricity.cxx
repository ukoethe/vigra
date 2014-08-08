#define PY_ARRAY_UNIQUE_SYMBOL vigranumpygraphs_PyArray_API
#define NO_IMPORT_ARRAY

/*boost python before anything else*/
#include <boost/python.hpp>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>
#include <vigra/eccentricitytransform.hxx>

namespace python = boost::python;

namespace vigra
{

    template < int N, class PixelType1, class PixelType2 >
    NumpyAnyArray
    pythonEccentricityTransformOnLabels(
            const NumpyArray<N, PixelType1 > & src,
            NumpyArray<N, PixelType2 > dest
    ){
        dest.reshapeIfEmpty(src.taggedShape(),
                           "eccentricityTransformOnLabels(): Output image has wrong dimensions");
        eccentricityTransformOnLabels(src, dest);
        return dest;
    }

    template< int N, class PixelType1, class PixelType2 >
    void defineEcc()
    {
        python::def("_eccentricityTransform",registerConverters(&pythonEccentricityTransformOnLabels< N, PixelType1, PixelType2 >),
                    (
                        python::arg("labels"),
                        python::arg("out")=python::object()
                    )
        );
    }

    void defineEccentricity()
    {
        defineEcc<2, UInt32, float>();
    }

}
