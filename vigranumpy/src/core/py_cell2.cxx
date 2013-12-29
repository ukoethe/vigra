#define PY_ARRAY_UNIQUE_SYMBOL vigranumpycgp_PyArray_API
#define NO_IMPORT_ARRAY

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

#include <vigra/cgp/cgp.hxx>
#include <vigra/cgp/cgp_python.hxx>
#include "py_cell_visitor.hxx"

namespace vigra {

namespace python = boost::python;

void export_cell2()
{
    using namespace python;
    
    docstring_options doc_options(true, true, false);


    ////////////////////////////////////////
    // Region Graph
    ////////////////////////////////////////
    // basic types
    // tgrid and input image type
    typedef Cgp<CgpCoordinateType,CgpLabelType> CgpType;
    typedef CgpType::TopologicalGridType TopologicalGridType;

    typedef  vigra::NumpyArray<2 ,vigra::Singleband < CgpLabelType > > InputLabelImageType;
    // cgp type and cell types
    typedef CgpType::PointType PointType;
    // bound vector
    typedef std::vector<float> FloatVectorType;
    typedef std::vector<CgpLabelType> LabelVectorType;
    // point vector
    typedef std::vector<PointType> PointVectorType;
    // geo cells 
    typedef CgpType::Cell2 Cell2Type;
    typedef CgpType::Cells2 Cell2VectorType;

    // cells
    python::class_<Cell2Type>("Cell2",python::init<>())
        .def(CellTypeSuite<Cell2Type>())
    ;
}

} // namespace vigra

