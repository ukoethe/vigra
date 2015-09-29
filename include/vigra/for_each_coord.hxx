#ifndef VIGRA_FOR_EACH_COORD_HXX
#define  VIGRA_FOR_EACH_COORD_HXX

#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/accumulator.hxx>





namespace vigra{


template<unsigned int DIM>
struct ForEachCoord{
    typedef TinyVector<MultiArrayIndex, DIM> Coord;
    template<class F>
    static void forEachCoord(const Coord & shape, F &&f){
        MultiCoordinateIterator<DIM> i(shape), end = i.getEndIterator();
        for(; i != end; ++i){
            f(*i);
        }
    }
};

}


#endif /* VIGRA_FOR_EACH_COORD_HXX */
