#ifndef VIGRA_GRID_RAG_VISUALIZATION_HXX
#define  VIGRA_GRID_RAG_VISUALIZATION_HXX

#include <list>
#include <map>
#include <functional>

#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/accumulator.hxx>
#include <vigra/multi_blocking.hxx>
#include <vigra/for_each_coord.hxx>

namespace vigra{




class SliceEdges{

public:
    typedef AdjacencyListGraph Graph;
    typedef typename Graph::Edge Edge;
    typedef TinyVector<MultiArrayIndex, 2> Coord;
    typedef std::list<Coord> CoordList;

    SliceEdges(const AdjacencyListGraph & graph)
    : graph_(graph){

    }

    template<class LABEL_TYPE>
    void findSlicesEdges(const MultiArrayView<2, LABEL_TYPE> & labels){

        lines_.clear();
        edgeToLines_.clear();

        const auto & shape = labels.shape();


        Coord coord;


        typedef std::map<Edge, std::set<Coord> > Edges; 
        
        Edges edges;

        // fill the coord map
        ForEachCoord<2>::forEachCoord(shape, [&](const Coord & coord){
            for(size_t d=0; d<2; ++d){
                const auto lu = labels[coord];
                Coord coordV = coord;
                coordV[d] += 1;
                if(coordV[d]<shape[d]){
                    const auto lv = labels[coordV];
                    if(lu!=lv){
                        edges[graph_.findEdge(lu,lv)].insert(coord+coordV);
                    }
                }
            }
        });


        for(auto & kv : edges){

            const auto edge = kv.first;
            auto & edgeSet = kv.second;

            

            CoordList cl;

            // starting point
            // auto startEdgeIter = edgeSet.begin();
            // auto startEdge = *startEdgeIter;
            // edgeSet.erase(startEdgeIter);

            

            auto appendCoord = [&](const Coord & start, CoordList & coordList, const bool front)->bool{
                // horizontal edge
                if(start[0]%2==0){
                    const Coord offsets[6] ={
                        Coord(-1,-1),
                        Coord( 1,-1),
                        Coord(-2, 0),
                        Coord( 2, 0),
                        Coord(-1, 1),
                        Coord( 1, 1)
                    };
                    for(size_t i=0; i<6; ++i){
                        auto fRes = edgeSet.find(start+offsets[i]);
                        if(fRes!=edgeSet.end()){
                            Coord ecoord = *fRes;
                            if(front)
                                coordList.push_front(ecoord);
                            else
                                coordList.push_back(ecoord);
                            edgeSet.erase(fRes);
                            return true;
                        }
                    }
                }
                else{
                    const Coord offsets[6] ={
                        Coord(-1,-1),
                        Coord( 0,-2),
                        Coord( 1,-1),
                        Coord(-1, 1),
                        Coord( 1, 1),
                        Coord( 0, 2)
                    };
                    for(size_t i=0; i<6; ++i){
                        auto fRes = edgeSet.find(start+offsets[i]);
                        if(fRes!=edgeSet.end()){
                            Coord ecoord = *fRes;
                            if(front)
                                coordList.push_front(ecoord);
                            else
                                coordList.push_back(ecoord);
                            edgeSet.erase(fRes);
                            return true;
                        }
                    }
                }
                return false;
            };

            while(!edgeSet.empty()){
                // starting point
                auto startEdgeIter = edgeSet.begin();
                auto startEdge = *startEdgeIter;
                edgeSet.erase(startEdgeIter);
                CoordList cl;
                cl.push_back(startEdge);

                while(true){
                    auto f = appendCoord(cl.front(), cl, true);
                    auto b = appendCoord(cl.back(), cl, false);
                    if(!(f||b)){
                        break;
                    }
                }
                //coordListVec.push_back(cl);

                edgeToLines_[edge].push_back(lines_.size());
                lines_.push_back(cl);
            }

            //sortedCoords[key] = coordListVec;
        }
    }


    UInt64 nVisibleEdges()const{
        return edgeToLines_.size();
    }


    UInt64 nLines()const{
        return lines_.size();
    }

    UInt64 nLines(const Edge edge)const{
        const auto fIter = edgeToLines_.find(edge);
        return fIter==edgeToLines_.end() ? 0 :  fIter->second.size();
    }

    UInt64 nLines(const UInt64 edge)const{
        return nLines(graph_.edgeFromId(edge));
    }

    void visibleEdges(MultiArrayView<1,Int32> edges){
        auto  c = 0;
        for(const auto kv: edgeToLines_){
            edges[c++] = graph_.id(kv.first);
        }
    }

    UInt64 lineSize(const UInt64 lineIndex)const{
        return lines_[lineIndex].size();
    }

    void line(
        const UInt64 lineIndex, 
        MultiArrayView<1, Coord > & out
    )const{
        const auto line = lines_[lineIndex];
        auto c = 0 ;
        for(const auto lc : line){
            out[c++] = lc;
        }
    }

    void lineIds(const Edge & edge, MultiArrayView<1, UInt64> & out)const{
        const auto lineIds = edgeToLines_.find(edge)->second;
        auto c=0;
        for(const auto li : lineIds){
            out[c++] = li;
        }
    }
    void lineIds(const UInt64 & edge, MultiArrayView<1, UInt64> & out)const{
        return lineIds(graph_.edgeFromId(edge), out);
    }

private:

    std::vector< CoordList > lines_; 
    std::map<Edge, std::vector<UInt64> > edgeToLines_;
    const AdjacencyListGraph & graph_;
};




struct SortedLine{
};

class EdgeTileManager2D{

public:

    typedef typename AdjacencyListGraph::Edge Edge;
    typedef MultiBlocking<2> Blocking;
    typedef typename Blocking::Point Point;
    typedef typename Blocking::Block Block;
    typedef Point Shape;

    typedef MultiArray<2, UInt32>     LabelSliceArray;
    typedef MultiArrayView<2, UInt32> LabelSliceView;

    typedef std::set<UInt32> BlockIndexSet;

    typedef std::function< void (const Point &, const Point &, LabelSliceView & ) > GetLabelsSliceCallback;



    EdgeTileManager2D(
        const AdjacencyListGraph & graph,   
        const Shape & shape, 
        const Shape & blockShape
    )
    :   graph_(graph),
        blocking_(shape, blockShape),
        viewRect_(Point(0),Point(0)),
        isVisible_(blocking_.numBlocks(), false),
        visibleBlocks_(),
        appearedBlocks_(),
        disappearedBlocks_(),
        getLabelSliceCallback_()
    {
        const auto nb = blocking_.numBlocks();
        visibleBlocks_.reserve(nb);
        appearedBlocks_.reserve(nb);
        disappearedBlocks_.reserve(nb);
    }


    void setViewRect(const Point & begin, const Point & end){

        viewRect_ = Block(begin, end);
        visibleBlocks_.resize(0);
        appearedBlocks_.resize(0);
        disappearedBlocks_.resize(0);

        auto  c=0;
        for(auto iter=blocking_.blockBegin(); iter!=blocking_.blockEnd(); ++iter){
            if(viewRect_.intersects(*(iter))){
                visibleBlocks_.push_back(c);
                if(isVisible_[c] == false){
                    appearedBlocks_.push_back(c);
                    blockAppeared(c);
                }
                isVisible_[c] = true;
            }
            else{
                if(isVisible_[c] == true){
                    disappearedBlocks_.push_back(c);
                    blockDisappeared(c);
                }
                isVisible_[c] = false;
            }
        }
    }

    /**
     * @brief this will notify that the viewer was scrolled.
     *  As a result, all edges for visible blocks will 
     *  be recomputed
     */
    void scrolled(){
        // todo
    }


private:

    void blockAppeared(const UInt32 bi){

    }

    void blockDisappeared(const UInt32 bi){

    }


    const AdjacencyListGraph & graph_;
    Blocking blocking_;
    Block viewRect_;

    std::vector<bool>   isVisible_;
    std::vector<UInt32> visibleBlocks_;
    std::vector<UInt32> appearedBlocks_;
    std::vector<UInt32> disappearedBlocks_;
    GetLabelsSliceCallback getLabelSliceCallback_;


    // which edges can i find in a certain block
    std::vector<  std::vector<Edge>  >   blockToEdges_;

    // for each edge in each block there is a line
    std::vector<  std::vector<SortedLine>  >  blocksLines_;

    // in which blocks can i find a certain edge
    std::map< Edge, BlockIndexSet  > edgeToBlocks_;
};






}


#endif /*VIGRA_GRID_RAG_VISUALIZATION_HXX*/
