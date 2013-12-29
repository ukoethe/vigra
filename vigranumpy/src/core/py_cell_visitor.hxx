#ifndef CGP2D_PY_CELL_VISITOR_HXX
#define CGP2D_PY_CELL_VISITOR_HXX

#include <boost/python/detail/wrap_python.hpp>
#include <boost/python.hpp>
#include <string>
#include <cmath>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

#include <vigra/cgp/cgp.hxx>
#include <vigra/cgp/cgp_python.hxx>

#include <boost/python/detail/wrap_python.hpp>
#include <boost/python.hpp>

#include <sstream>
//#include <Python.h>
#include <numpy/arrayobject.h>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

/*
#include <numpy/noprefix.h>
#ifdef UInt32
#undef UInt32
#endif 
*/


namespace vigra {

namespace python = boost::python;




template<class CELL>
python::tuple pointNumpyTupel(const CELL & cell){
    const size_t numPoints=cell.points().size();
    typedef vigra::NumpyArray<1, vigra::UInt32 >  SingleCoordArrayType;
    typedef typename SingleCoordArrayType::difference_type ShapeType;
    const ShapeType shape(numPoints);

    SingleCoordArrayType cx(shape),cy(shape);
    for(size_t i=0;i<numPoints;++i){
        cx(i)=cell.points()[i][0];
        cy(i)=cell.points()[i][1];
    }
    vigra::NumpyAnyArray ax=cx,ay=cy;
    return python::make_tuple(ax,ay);
}


template<class CELL>
std::vector<float> * getAngles(const CELL & cell,const typename CELL::TopologicalGridType & grid,const int radius){
    std::vector<float> * angles=new std::vector<float>();
    cell.getAngles(grid,radius,*angles);
    return angles;
}



template<class CELL_TYPE>
class CellTypeSuite : public boost::python::def_visitor<CellTypeSuite<CELL_TYPE> >
{
public:
    friend class boost::python::def_visitor_access;
    typedef CELL_TYPE CellType;
    typedef typename  CellType::LabelType LabelType;
    typedef typename  CellType::PointType PointType;
    typedef typename  CellType::FloatPointType FloatPointType;
    typedef std::vector<LabelType> LabelVector;
    typedef std::vector<PointType> PointVector;

    typedef typename  CellType::CgpType  CgpType;

    template <class classT>
    void visit(classT& c) const
    {
        c
            //////////////////////
            //      properties
            //////////////////////
            .add_property("cgp", python::make_function(
                (const CgpType & (CellType::*)())  &CellType::cgp, python::return_internal_reference<>() )
            )
            .add_property("bounds", python::make_function(  
                (const LabelVector & (CellType::*)()) &CellType::bounds, python::return_internal_reference<>() 
            ))
            .add_property("boundedBy", python::make_function(
                (const LabelVector & (CellType::*)()) &CellType::boundedBy, python::return_internal_reference<>() 
            ))
            .add_property("points", python::make_function(
                (const PointVector & (CellType::*)()) &CellType::points, python::return_internal_reference<>() 
            ))
            .add_property("label", python::make_function(
                (LabelType (CellType::*)()) &CellType::label, python::return_value_policy<python::return_by_value>() 
            ))
            .add_property("cellType", python::make_function(
                (size_t (CellType::*)()) &CellType::cellType, python::return_value_policy<python::return_by_value>() 
            ))
            //////////////////////
            //      functions
            //////////////////////
            .def("centerPoint",&CellType::centerCoordinate,python::return_value_policy<python::return_by_value>() )
            .def("boundingBox",&boundingBox)
            .def("pointArray",vigra::registerConverters(&pointNumpyTupel<CellType>))
            .def("boundView",&boundView
                ,  
                python::with_custodian_and_ward_postcall<0,1>() 
            )
        ;
    }

    static python::tuple boundingBox(const CellType & cell){
        std::pair<PointType,PointType> point=cell.boundingBox();
        return python::make_tuple(point.first,point.second);
    }
    

    static boost::python::object   boundView(CellType & cell){
        LabelType  * ptr = & cell.__bounds__()[0];
        return numpyView1d(ptr,cell.bounds().size(),false);
    }

};





} // namespace vigra




#endif // CGP2D_PY_CELL_VISITOR_HXX