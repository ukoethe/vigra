#ifndef VIGRA_GRID_RAG_VISUALIZATION_HXX
#define  VIGRA_GRID_RAG_VISUALIZATION_HXX

#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/accumulator.hxx>
#include <vigra/for_each_coord.hxx>

namespace vigra{

template<class LABEL_TYPE>
void findSlicesEdges(const MultiArrayView<2, LABEL_TYPE> & labels){

    const auto & shape = labels.shape();
    typedef TinyVector<MultiArrayIndex, 2> Coord;

    Coord coord;

    typedef std::pair<UInt64, UInt64> LPair;
    std::map<LPair, std::vector<Coord> > Edges; 
    Edges edges;

    // fill the coord map
    ForEachCoord<2>::forEachCoord(shape, [&](const Coord & coord){
        for(size_t d=0; d<2; ++d){
            const auto lu = labels[coord];
            coordV = coord;
            coordV[d] += 1;
            if(coordV<shape[d]){
                const auto lv = labels[coordV];
                if(lu!=lv){
                    const auto lMin = std::min(lu, lv);
                    const auto lMax = std::max(lu, lv);
                    edges[LPair(lMin,Lmax)].push_back(coord+coordV);
                }
            }
        }
    });
}


#endif /*VIGRA_GRID_RAG_VISUALIZATION_HXX*/
