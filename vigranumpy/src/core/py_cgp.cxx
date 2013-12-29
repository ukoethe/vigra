#define PY_ARRAY_UNIQUE_SYMBOL vigranumpycgp_PyArray_API
#define NO_IMPORT_ARRAY

#include <string>
#include <cmath>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


#include <boost/array.hpp>

#include <boost/accumulators/accumulators.hpp>

#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/extended_p_square_quantile.hpp>
#include <boost/accumulators/statistics/tail_quantile.hpp>

#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

#include <vigra/cgp/cgp.hxx>
#include <vigra/cgp/cgp_python.hxx>

namespace python = boost::python;

namespace vigra
{





template<class CGP>
vigra::NumpyAnyArray cell1BoundsArray(
    const CGP & cgp
){
    typedef vigra::NumpyArray<2,int>  ResultArray;
    typedef typename ResultArray::difference_type ShapeType;

    typedef typename CGP::Cells0 Cells0;
    typedef typename CGP::Cells1 Cells1;
    typedef typename CGP::Cells2 Cells2;

    typedef typename CGP::Cell0 Cell0;
    typedef typename CGP::Cell1 Cell1;
    typedef typename CGP::Cell2 Cell2;

    const size_t numBoundaries = cgp.numCells(1);
    const size_t numRegion     = cgp.numCells(2);
    const Cells1 & cells1=cgp.geometry1();
    CGP_ASSERT_OP(numBoundaries,==,cells1.size());

    ResultArray resultArray = ResultArray(ShapeType(numBoundaries,2));

    for(size_t b=0;b<numBoundaries;++b){
        const Cell1 & cell1 = cells1[b];
        const int la = cell1.bounds()[0];
        const int lb = cell1.bounds()[1];
        CGP_ASSERT_OP(la,!=,0);
        CGP_ASSERT_OP(lb,!=,0);
        resultArray(b,0)=la;
        resultArray(b,1)=lb;
    }
    return resultArray;
}


template<class CELL>
python::tuple pointNumpyTupel(const CELL & cell){
    const size_t numPoints=cell.points().size();
    typedef vigra::NumpyArray<1,vigra::UInt32>  SingleCoordArrayType;
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


template<class TGRID>
vigra::NumpyAnyArray getCellLabelGrid(
                                const TGRID & tgrid,
                                int cellType,
                                bool useTopologicalShape,
                                vigra::NumpyArray<2, vigra::Singleband<npy_uint32> > res = vigra::NumpyArray<2,vigra::Singleband<npy_uint32> >()){

    if(useTopologicalShape){
        res.reshapeIfEmpty(tgrid.tgrid().shape());
        std::fill(res.begin(),res.end(),0);

        if(cellType==0){
            for(size_t y=1;y<tgrid.shape(1);y+=2)
            for(size_t x=1;x<tgrid.shape(0);x+=2){
                if( tgrid(x,y)!=0 )
                    res(x,y)=tgrid(x,y);
            }
        }
        else if(cellType==1){
            for(size_t y=0;y<tgrid.shape(1);++y)
            for(size_t x=0;x<tgrid.shape(0);++x){
                if(  (  (x%2==0 && y%2!=0) || (x%2!=0 && y%2==0) ) && tgrid(x,y)!=0 )
                    res(x,y)=tgrid(x,y);
            }
        }
        else if(cellType==2){
            for(size_t y=0;y<tgrid.shape(1);y+=2)
            for(size_t x=0;x<tgrid.shape(0);x+=2){
                res(x,y)=tgrid(x,y);
            }
        }
    }
    else{
        typedef typename TGRID::LabelImageType::difference_type ShapeType;
        const ShapeType shape( (tgrid.shape(0)+1)/2,(tgrid.shape(1)+1)/2   );
        res.reshapeIfEmpty(shape);
        std::fill(res.begin(),res.end(),0);

        if(cellType==0){
            for(size_t y=1;y<tgrid.shape(1);y+=2)
            for(size_t x=1;x<tgrid.shape(0);x+=2){
                if( tgrid(x,y)!=0 ){
                    // for junction-pixel:
                    // 1 pixel in t-grid => 4 pixels in grid 
                    res( (x+1)/2,(y+1)/2 )=tgrid(x,y);
                    res( (x+1)/2,(y-1)/2 )=tgrid(x,y);
                    res( (x-1)/2,(y+1)/2 )=tgrid(x,y);
                    res( (x-1)/2,(y+1)/2 )=tgrid(x,y);
                }
            }
        }

        else if(cellType==1){
            for(size_t y=0;y<tgrid.shape(1);++y)
            for(size_t x=0;x<tgrid.shape(0);++x){

                if(tgrid(x,y)!=0 ){

                    // - boundary 
                    if( x%2==0 && y%2!=0 ){
                        // for boundary-pixel:
                        // 1 pixel in t-grid => 2 pixels in grid 
                        res( x/2,(y+1)/2 )=tgrid(x,y);
                        res( x/2,(y-1)/2 )=tgrid(x,y);
                    }
                    //  |  boundary
                    else if(x%2!=0 && y%2==0){
                        // for boundary-pixel:
                        // 1 pixel in t-grid => 2 pixels in grid 
                        res( (x-1)/2,y/2 )=tgrid(x,y);
                        res( (x+1)/2,y/2 )=tgrid(x,y);
                    }
                }
            }
        }

        else if(cellType==2){
            for(size_t y=0;y<tgrid.shape(1);y+=2)
            for(size_t x=0;x<tgrid.shape(0);x+=2){
                // for region-pixel:
                // 1 pixel in t-grid => 4 pixels in grid 
                res(x/2,y/2)=tgrid(x,y);
            }
        }   
    }
    return res;
}


template<class CGP>
vigra::NumpyArray<1, unsigned int> pyCgpSerialize(
    const CGP& cgp
) {
    vigra::NumpyArray<1, unsigned int> result;
    std::vector<unsigned int> res = cgp.serialize();
    result.reshape(vigra::Shape1(res.size()));
    std::copy(res.begin(), res.end(), result.begin());
    return result;
}

template<class CGP>
const typename CGP::TopologicalGridType * merge2Cells(
    const CGP & cgp,
    vigra::NumpyArray<1,npy_uint32> cell1States
){
    typename CGP::TopologicalGridType * tgrid = new typename CGP::TopologicalGridType();
    cgp.merge2Cells(cell1States.begin(),cell1States.end(),*tgrid);
    return tgrid;
}





template<class CGP>
vigra::NumpyAnyArray featuresToFeatureImage(
    const CGP & cgp,
    int cellType,
    vigra::NumpyArray<1,float> features,
    const bool ignoreInactive,
    const float inactiveValue,
    const bool useTopologicalShape=true,
    vigra::NumpyArray<2, vigra::Singleband<float> > res = vigra::NumpyArray<2,vigra::Singleband<float> >()
){
    typedef typename  CGP::LabelType LabelType;
    if(useTopologicalShape)
        res.reshapeIfEmpty(cgp.tgrid().tgrid().shape());
    else
        res.reshapeIfEmpty(cgp.tgrid().shapeLabeling());
        

    if(ignoreInactive==false)
        std::fill(res.begin(),res.end(),inactiveValue);

    if(useTopologicalShape){
        if (cellType==0){
            for(size_t y=0;y<cgp.shape(1);++y)
            for(size_t x=0;x<cgp.shape(0);++x){
                const LabelType cellLabel=cgp(x,y);
                if(cellType==0)
                if( (x%2!=0 && y%2!=0) && cellLabel!=0)
                    res(x,y)=features[cellLabel-1];
            }
        }
        if(cellType==1){
            for(size_t y=0;y<cgp.shape(1);++y)
            for(size_t x=0;x<cgp.shape(0);++x){
                const LabelType cellLabel=cgp(x,y);
                if( ( (x%2==0 && y%2!=0)  || (x%2!=0 && y%2==0 ) ) && cellLabel!=0)
                    res(x,y)=features[cellLabel-1];
            }
        }
        else if(cellType==2){
            for(size_t y=0;y<cgp.shape(1);++y)
            for(size_t x=0;x<cgp.shape(0);++x){
                const LabelType cellLabel=cgp(x,y);
                if( (x%2==0 && y%2==0))
                    res(x,y)=features[cellLabel-1];
            }
        }
    }


    else{

        if(cellType==0){
            for(size_t y=1;y<cgp.shape(1);y+=2)
            for(size_t x=1;x<cgp.shape(0);x+=2){
                if( cgp(x,y)!=0 ){
                    // for junction-pixel:
                    // 1 pixel in t-grid => 4 pixels in grid 
                    res( (x+1)/2,(y+1)/2 )=features[cgp(x,y)-1];
                    res( (x+1)/2,(y-1)/2 )=features[cgp(x,y)-1];
                    res( (x-1)/2,(y+1)/2 )=features[cgp(x,y)-1];
                    res( (x-1)/2,(y+1)/2 )=features[cgp(x,y)-1];
                }
            }
        }

        else if(cellType==1){
            for(size_t y=0;y<cgp.shape(1);++y)
            for(size_t x=0;x<cgp.shape(0);++x){

                if(cgp(x,y)!=0 ){

                    // - boundary 
                    if( x%2==0 && y%2!=0 ){
                        // for boundary-pixel:
                        // 1 pixel in t-grid => 2 pixels in grid 
                        res( x/2,(y+1)/2 )=features[cgp(x,y)-1];
                        res( x/2,(y-1)/2 )=features[cgp(x,y)-1];
                    }
                    //  |  boundary
                    else if(x%2!=0 && y%2==0){
                        // for boundary-pixel:
                        // 1 pixel in t-grid => 2 pixels in grid 
                        res( (x-1)/2,y/2 )=features[cgp(x,y)-1];
                        res( (x+1)/2,y/2 )=features[cgp(x,y)-1];
                    }
                }
            }
        }

        else if(cellType==2){
            for(size_t y=0;y<cgp.shape(1);y+=2)
            for(size_t x=0;x<cgp.shape(0);x+=2){
                // for region-pixel:
                // 1 pixel in t-grid => 4 pixels in grid 
                res(x/2,y/2)=features[cgp(x,y)-1];
            }
        }   
    }
    return res;
}



template<class CGP>
vigra::NumpyAnyArray cellSizes(
    const CGP & cgp,
    const size_t cellType,
    vigra::NumpyArray<1, vigra::UInt32 > res = vigra::NumpyArray<1,vigra::UInt32> ()
){
    const size_t nCells = cgp.numCells(cellType);
    res.reshapeIfEmpty(typename vigra::NumpyArray<1,vigra::UInt32>::difference_type(nCells));
    for(size_t  c=0;c<nCells;++c){
        res(c)=cgp.cellSize(cellType,c);
    }
    return res;
}

void export_cgp()
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
    typedef CgpType::Cell0 Cell0Type;
    typedef CgpType::Cell1 Cell1Type;
    typedef CgpType::Cell2 Cell2Type;

    typedef CgpType::Cells0 Cell0VectorType;
    typedef CgpType::Cells1 Cell1VectorType;
    typedef CgpType::Cells2 Cell2VectorType;

    // cell vectors
    python::class_<TopologicalGridType>("TopologicalGrid",python::init<const InputLabelImageType & >())
    .add_property("shape", python::make_function(&TopologicalGridType::shapeTopologicalGrid, python::return_value_policy<return_by_value>()) )
    .add_property("shapeLabeling", python::make_function(&TopologicalGridType::shapeLabeling, python::return_value_policy<return_by_value>()) )
    .def("numCells",&TopologicalGridType::numCells)
    .def("labelGrid",vigra::registerConverters(&getCellLabelGrid<TopologicalGridType> ) ,
        (
            arg("cellType"),
            arg("useTopologicalShape")=true,
            arg("out")=python::object() 
        )  
    )
    ;

    // float vector
    python::class_<FloatVectorType>("FloatVector",init<>())
        .def(vector_indexing_suite<FloatVectorType ,true >())
    ;
    // bound / bounded by vector
    python::class_<LabelVectorType> exporter = python::class_<LabelVectorType>("LabelVector",init<>())
        .def(vector_indexing_suite<LabelVectorType >())
    ;
    // point   vector
    python::class_<PointVectorType>("PointVector",init<>())
        .def(vector_indexing_suite<PointVectorType ,true>())
    ;



    /************************************************************************/
    /* C e l l B a s e                                                      */
    /************************************************************************/

    python::class_<CgpType>("Cgp",python::init<const TopologicalGridType & >()[with_custodian_and_ward<1 /*custodian == self*/, 2 /*ward == const TopologicalGridType& */>()] )

        .add_property("shape", python::make_function(&CgpType::shapeTopologicalGrid, python::return_value_policy<return_by_value>()) )
        .add_property("shapeLabeling", python::make_function(&CgpType::shapeLabeling, python::return_value_policy<return_by_value>()))

        
        //.def("shape",&CgpType::shape)
        
        .add_property("tgrid", python::make_function(&CgpType::tgrid, return_internal_reference<>() ))
        .add_property("cells0", python::make_function(&CgpType::geometry0, return_internal_reference<>() ))
        .add_property("cells1", python::make_function(&CgpType::geometry1, return_internal_reference<>() ))
        .add_property("cells2", python::make_function(&CgpType::geometry2, return_internal_reference<>() ))
        .def("cell1BoundsArray",vigra::registerConverters(&cell1BoundsArray<CgpType>))



        .def("cellSizes",vigra::registerConverters(&cellSizes<CgpType>),
            (
                arg("cellType"),
                arg("out")=python::object()
            )
        )
        .def("serialize", &pyCgpSerialize<CgpType>)
        .def("numCells",&CgpType::numCells)
        .def("featureToImage",vigra::registerConverters(&featuresToFeatureImage<CgpType>),
            (
                arg("cellType"),
                arg("features"),
                arg("ignoreInactive")=false,
                arg("inactiveValue")=0.0f,
                arg("useTopologicalShape")=true,
                arg("out")=python::object()
            )
        )
        .def("merge2Cells",vigra::registerConverters( &merge2Cells<CgpType> ) ,python::return_value_policy<python::manage_new_object>(),
            (
                arg("cell1States")
            )
        )
    ;
}

} // namespace vigra

