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

    template< int N, class T, class S >
    NumpyAnyArray
    pythonEccentricityTransformOnLabels(
            const NumpyArray< N, T > & src,
            NumpyArray< N, S > dest
    ){
        dest.reshapeIfEmpty(src.taggedShape(),
                           "eccentricityTransformOnLabels(): Output image has wrong dimensions");
        eccentricityTransformOnLabels(src, dest);
        return dest;
    }

    template< int N, class T >
    NumpyAnyArray
    pythonEccentricityCenters(
            const NumpyArray< N, T > & src,
            NumpyArray< 2, MultiArrayIndex > centers
    ){
        T maxLabel = *std::max_element(src.begin(), src.end());
        centers.reshapeIfEmpty(Shape2(maxLabel, N),
                               "eccentricityCenters(): Output has wrong dimensions");
        findEccentricityCenters(src, centers, maxLabel);
        return centers;
    }

    template< int N, class T, class S >
    python::tuple
    pythonEccentricityTransformWithCenters(
            const NumpyArray< N, T > & src,
            NumpyArray< N, S > dest,
            NumpyArray< 2, MultiArrayIndex > centers
    ){
        dest.reshapeIfEmpty(src.taggedShape(),
                            "eccentricityTransformWithCenters(): Output (dest) has wrong dimensions");
        T maxLabel = *std::max_element(src.begin(), src.end());
        centers.reshapeIfEmpty(Shape2(maxLabel, N),
                               "eccentricityTransformWithCenters(): Output (centers) has wrong dimensions");
        eccentricityTransformOnLabels(src, dest, centers);
        return python::make_tuple(dest, centers);
    }

    template< int N, class T, class S >
    void defineEcc()
    {
        python::def("_eccentricityTransform",registerConverters(&pythonEccentricityTransformOnLabels< N, T, S >),
                    (
                        python::arg("labels"),
                        python::arg("out")=python::object()
                    )
        );
    }

    template< int N, class T >
    void defineCenters()
    {
        python::def("_eccentricityCenters",registerConverters(&pythonEccentricityCenters< N, T >),
                    (
                        python::arg("labels"),
                        python::arg("out")=python::object()
                    )
        );
    }

    template< int N, class T, class S >
    void defineEccWithCenters()
    {
        python::def("_eccentricityTransformWithCenters",registerConverters(&pythonEccentricityTransformWithCenters<N, T, S >),
                    (
                        python::arg("labels"),
                        python::arg("ecc")=python::object(),
                        python::arg("centers")=python::object()
                    )
        );
    }

    void defineEccentricity()
    {
        defineEcc< 2, UInt32, float >();
        defineEcc< 3, UInt32, float >();
        defineCenters< 2, UInt32 >();
        defineCenters< 3, UInt32 >();
        defineEccWithCenters< 2, UInt32, float >();
        defineEccWithCenters< 3, UInt32, float >();
    }

}
