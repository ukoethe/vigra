#ifndef VIGRA_GRID_RAG_2_HXX
#define  VIGRA_GRID_RAG_2_HXX

#include <list>
#include <map>
#include <functional>
#include <iostream>
#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/accumulator.hxx>
#include <vigra/multi_blocking.hxx>
#include <vigra/for_each_coord.hxx>
#include "vigra/threadpool.hxx"

namespace vigra{




class SortedLines{

public:
    typedef AdjacencyListGraph Graph;
    typedef typename Graph::Edge Edge;
    typedef TinyVector<MultiArrayIndex, 2> Coord;
    typedef std::list<Coord> CoordList;

    template<class LABEL_TYPE>
    SortedLines(const AdjacencyListGraph & graph, const MultiArrayView<2, LABEL_TYPE> & labels)
    :   graph_(graph),
        edgeToLines_(graph_.maxEdgeId()+1){
        this->compute(labels);
    }



    UInt64 nLines()const{
        return lines_.size();
    }

    UInt64 nLines(const Edge edge)const{
        return edgeToLines_[graph_.id(edge)].size();
    }

    UInt64 nLines(const UInt64 edge)const{
        return nLines(graph_.edgeFromId(edge));
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

    const  CoordList  & line(const UInt64 lineIndex)const{
        return lines_[lineIndex];
    }

    void lineIds(const Edge & edge, MultiArrayView<1, UInt64> & out)const{
        const auto lineIds = edgeToLines_[graph_.id(edge)];
        auto c=0;
        for(const auto li : lineIds){
            out[c++] = li;
        }
    }
    void lineIds(const UInt64 & edge, MultiArrayView<1, UInt64> & out)const{
        return lineIds(graph_.edgeFromId(edge), out);
    }
    const AdjacencyListGraph & graph() const{
        return graph_;
    }
    const std::vector<UInt64> & edgeToLines(const UInt64 id)const{
        return edgeToLines_[id];
    }
private:

    template<class LABEL_TYPE>
    void compute(const MultiArrayView<2, LABEL_TYPE> & labels){

        lines_.clear();


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

                edgeToLines_[graph_.id(edge)].push_back(lines_.size());
                lines_.push_back(cl);
            }

            //sortedCoords[key] = coordListVec;
        }
    }

    const AdjacencyListGraph & graph_;
    std::vector< CoordList > lines_; 
    std::vector< std::vector<UInt64> > edgeToLines_;
    
};




inline void curvature(
    const SortedLines & sl, 
    const int nThreads,
    MultiArrayView<2, float> out
){
    typedef SortedLines::Coord Coord;
    typedef TinyVector<long double, 2> FCoord;

    using namespace vigra::acc;

    typedef StandardQuantiles<AutoRangeHistogram<0> > Quantiles;
    typedef AccumulatorChain<double,
    Select<Mean, Sum, Minimum, Maximum, Variance, Skewness, Kurtosis, Quantiles> > Acc;

    
    const auto & graph = sl.graph();

    parallel_foreach( nThreads, graph.edgeNum(),
    [&](size_t /*thread_id*/, int id){

        const long double stencilDX[7] ={ 
           -1.0/60.0, 
            3.0/20.0, 
           -3.0/4.0, 
            0.0, 
            3.0/4.0,
           -3.0/20.0, 
            1.0/60.0
        }; 

        const long double stencilDXX[7] ={ 
             1.0/90.0,    
            -3.0/20.0,   
             3.0/2.0,    
            -49/18.0,  
             3/2.0,     
            -3.0/20.0,   
             1.0/90.0
        };   

        const size_t paddingSize = 3;
        for(const auto lineIndex : sl.edgeToLines(id)){
            const auto & line  = sl.line(lineIndex);

            std::vector<Coord> paddedLineVec;
            paddedLineVec.reserve(line.size()+2*paddingSize);

            std::vector<FCoord> dxLine(line.size());
            std::vector<FCoord> dxxLine(line.size());
            std::vector<long double> rLine(line.size());

            // padded line
            for(size_t i=0; i<paddingSize; ++i)
                paddedLineVec.push_back(line.front());
            for(const auto & coord : line)
                paddedLineVec.push_back(coord);
            for(size_t i=0; i<paddingSize; ++i)
                paddedLineVec.push_back(line.back());

            for(size_t i=0; i<line.size(); ++i){
                auto vecI = i + paddingSize;
                const Coord & c = paddedLineVec[i+paddingSize];
                
                FCoord valDX(0.0), valDXX(0.0);
                for(size_t s=0; s<7; ++s){
                    FCoord otherVal  = paddedLineVec[ vecI + s-3];
                    //otherVal *= stencilDXX[s];
                    valDXX += otherVal*FCoord(stencilDXX[s]);
                    valDX  += otherVal*FCoord(stencilDX[s]);
                }
                dxLine[i] = valDX;
                dxxLine[i] = valDXX;

                auto va = std::pow(vigra::norm(valDX),3);
                auto vb = vigra::squaredNorm(valDX);
                auto vc = vigra::squaredNorm(valDXX);
                auto vd = vigra::sum(valDX*valDXX);
                auto r  = va/std::sqrt(std::max(vb*vc-vd*vd, std::numeric_limits<long double>::epsilon() * 40.0 ));
                rLine[i] = r;
            }
            //if(id==10){
            //    std::cout<<"\n";
            //    for(size_t i=0; i<line.size(); ++i){
            //        std::cout<<"l  "<<paddedLineVec[i+3]<<"\n";
            //        std::cout<<"dxx"<<dxLine[i]<<"\n";
            //        std::cout<<"dxx"<<dxxLine[i]<<"\n";
            //        std::cout<<"r  "<<1.0/rLine[i]<<"\n";
            //    }
            //}

            Acc accR;
            Acc accDxx;

            size_t n_bins = std::pow( line.size(), 1. / 2.5); 
            n_bins = std::max( 2, std::min(int(n_bins), 64) );
            accDxx.setHistogramOptions(HistogramOptions().setBinCount(n_bins));
            accR.setHistogramOptions(HistogramOptions().setBinCount(n_bins));

            for(unsigned int k=1; k <= accR.passesRequired(); ++k){
                for(size_t i=0;i<line.size();++i){
                    accR.updatePassN( 1.0/rLine[i], k);
                    accR.updatePassN( norm(rLine[i]), k);
                }
            }
        }

        

    });

}




}


#endif /*VIGRA_GRID_RAG_2_HXX*/
