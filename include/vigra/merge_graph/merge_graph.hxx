    
#ifndef VIGRA_MERGE_GRAPH_HXX
#define VIGRA_MERGE_GRAPH_HXX

/* std library */
#include <set>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <deque>
#include <map>
#include <stdexcept>
#include <sstream>

/* boost */
#include <boost/function.hpp>
#include <boost/signals2.hpp>


/* vigra */
#include <vigra/multi_array.hxx>
#include <vigra/tinyvector.hxx>
#include <vigra/multi_array.hxx>


/* this project*/
#include "partition.hxx"
#include "macros.hxx"
#include "merge_graph_node.hxx"
#include "merge_graph_edge.hxx"
#include "merge_graph_callbacks.hxx"

/* this project (TO BE REFACTORED) */
#include "partition.hxx"
#include "macros.hxx"

namespace vigra {


template<class LABEL_TYPE>
class MergeGraph : public MergeGraphCallbacks<LABEL_TYPE> {
    public:
        typedef LABEL_TYPE                  LabelType;
        typedef MergeGraphEdge<LabelType>   EdgeType;
        typedef MergeGraphNode<LabelType>   NodeType;
        typedef std::map<vigra::UInt64 , std::vector<LabelType>  > DoubleMap;
        // setup
        
    private:
        MergeGraph();                               // non empty-construction
        MergeGraph( const MergeGraph& other );      // non construction-copyable
        MergeGraph& operator=( const MergeGraph& ); // non copyable
    public:


        
        typedef detail_merge_graph::Partition<LabelType> UfdType;
        typedef typename UfdType::const_iterator ConstUdfIter;

        typedef ConstUdfIter EdgeIterator;
        typedef ConstUdfIter NodeIterator;

        EdgeIterator edgesBegin()const{
            return edgeUfd_.begin();
        }
        EdgeIterator edgesEnd()const{
            return edgeUfd_.end();
        }

        NodeIterator nodesBegin()const{
            return nodeUfd_.begin();
        }
        NodeIterator nodesEnd()const{
            return nodeUfd_.end();
        }
        

        MergeGraph(const size_t nNodes,const size_t nEdges);
        void   setInitalEdge(const size_t initEdge,const size_t initNode0,const size_t initNode1);

        // query
        size_t numberOfNodes()const;
        size_t numberOfEdges()const;
        size_t initNumberOfNodes()const;
        size_t initNumberOfEdges()const;


        const EdgeType & getInitalEdge(const LabelType index){
            return initEdges_[index];
        }

        //const EdgeType & getEdge_OLD(const LabelType index){
        //    return dynamicEdges_[index];
        //}

        // _NEW means better implementation
        EdgeType getEdge(const LabelType index)const{
            EdgeType edge = initEdges_[index];
            edge[0]=reprNode(edge[0]);
            edge[1]=reprNode(edge[1]);
            return edge;
        }





        const NodeType & getNode(const LabelType index){
            return dynamicNodes_[index];
        }

        //bool hasEdge_OLD(const LabelType edgeIndex)const{
        //    const bool hasEdge  = dynamicEdges_.find(edgeIndex)!=dynamicEdges_.end();
        //    return hasEdge;
        //}

        bool hasEdge(const LabelType edgeIndex)const{
            const LabelType rep = reprEdge(edgeIndex);
            if(rep!=edgeIndex){
                return false;
            }
            else{
                const EdgeType edge=getEdge(rep);
                return( edge[0]!=edge[1] );
            }
        }

        bool hasNode(const LabelType nodeIndex)const{
            const bool hasNode  = dynamicNodes_.find(nodeIndex)!=dynamicNodes_.end();
            return hasNode;
        }

        LabelType reprEdge(const LabelType edgeIndex)const{
            return edgeUfd_.find(edgeIndex);
        }
        LabelType reprNode(const LabelType nodeIndex)const{
            return nodeUfd_.find(nodeIndex);
        }


        bool stateOfInitalEdge(const LabelType initalEdge)const{
            const EdgeType edge=getEdge(initalEdge);
            return (edge[0]!=edge[1]);
        }


        template<class OUT_ITER>
        void stateOfInitalEdges(OUT_ITER begin,OUT_ITER end)const{
            const size_t d = std::distance(begin,end);
            for(size_t ie=0;ie<initNumberOfEdges();++ie){
                const EdgeType edge=getEdge(ie);
                if(edge[0]!=edge[1]){
                    begin[ie]=1;
                }
                else{
                    begin[ie]=0;
                }
            }
        }


        // modification
        void mergeParallelEdges();

        void mergeRegions(const LabelType edgeIndex);



  
    private:

        //typedef std::map<LabelType, EdgeType > EdgeMap;
        typedef std::map<LabelType, NodeType > NodeMap;


        //typedef typename EdgeMap::const_iterator ConstEdgeMapIterator;
        typedef typename NodeMap::const_iterator ConstNodeMapIterator;


        void combineDoubleEdges(const std::vector<LabelType> & ,const LabelType ,const LabelType );

        


        void searchLocalDoubleEdges(const NodeType & node , DoubleMap & doubleMap,const LabelType relabelFrom,const LabelType relabelTo);

        size_t nInitNodes_;
        size_t nInitEdges_;

        UfdType nodeUfd_;
        UfdType edgeUfd_;

        std::vector< EdgeType >     initEdges_;
        

        //EdgeMap dynamicEdges_;
        NodeMap dynamicNodes_;
};



template<class LABEL_TYPE>
MergeGraph<LABEL_TYPE>::MergeGraph(const size_t nNodes,const size_t nEdges)
:   MergeGraphCallbacks<LABEL_TYPE>(),
    nInitNodes_(nNodes),
    nInitEdges_(nEdges),
    nodeUfd_(nNodes),
    edgeUfd_(nEdges),
    initEdges_(nEdges)
{
     // fill nodes
    for(size_t n=0;n<nNodes;++n){
        dynamicNodes_[n]=NodeType();
        CGP_ASSERT_OP(nodeUfd_.find(n),==,n);
    }
}

template<class LABEL_TYPE>
void MergeGraph<LABEL_TYPE>::setInitalEdge(const size_t initEdge,const size_t initNode0,const size_t initNode1){

    // set up inital and dynamic edges
    initEdges_[initEdge].first =initNode0;
    initEdges_[initEdge].second=initNode1;
    //dynamicEdges_[initEdge] = initEdges_[initEdge];

    // set up the edges of a given region mapping
    dynamicNodes_[initNode0].edges_.insert(initEdge);
    dynamicNodes_[initNode1].edges_.insert(initEdge);
}

template<class LABEL_TYPE>
inline size_t MergeGraph<LABEL_TYPE>::numberOfNodes()const{
    return dynamicNodes_.size();
}

template<class LABEL_TYPE>
inline size_t MergeGraph<LABEL_TYPE>::numberOfEdges()const{
    //return dynamicEdges_.size();
    return edgeUfd_.numberOfSets();
}

template<class LABEL_TYPE>
inline size_t MergeGraph<LABEL_TYPE>::initNumberOfNodes()const{
    return nInitNodes_;
}

template<class LABEL_TYPE>
inline size_t MergeGraph<LABEL_TYPE>::initNumberOfEdges()const{
    return nInitEdges_;
}


template<class LABEL_TYPE>
void MergeGraph<LABEL_TYPE>::mergeParallelEdges(){
    typedef typename DoubleMap::const_iterator MapIter;
    DoubleMap pEdgeFinder;
    for(size_t e=0;e<nInitEdges_;++e){
        LabelType n0=initEdges_[e].first;
        LabelType n1=initEdges_[e].second;
        if(n0<n1){
            std::swap(n0,n1);
        }
        const size_t key = n0 + nInitNodes_*n1;
        pEdgeFinder[key].push_back(e);
    }

    for(MapIter iter=pEdgeFinder.begin();iter!=pEdgeFinder.end();++iter){
        const std::vector<LabelType> & dEdges = iter->second;
        CGP_ASSERT_OP(dEdges.size(),!=,0);

        if(dEdges.size()>1){
            //std::cout<<"found double edges "<<dEdges.size()<<"\n";
            const size_t key = iter->first;
            const size_t r1  = key/nInitNodes_;
            const size_t r0  = key - nInitNodes_*r1;
            this->combineDoubleEdges(dEdges,r0,r1);
        }
    }
}

template<class LABEL_TYPE>
void MergeGraph<LABEL_TYPE>::combineDoubleEdges(const std::vector<LABEL_TYPE> & toCombine,const LABEL_TYPE r0,const LABEL_TYPE r1){
    //std::set<LabelType> toCombineSet(toCombine.begin(),toCombine.end());
    //CGP_ASSERT_OP(toCombine.size(),==,toCombineSet.size());

    const LabelType newIndex = edgeUfd_.multiMerge(toCombine.front(),toCombine.begin()+1,toCombine.end());
   

    // update the two region between the double edge 
    const LabelType regions[2]={r0,r1};
    for(size_t r=0;r<2;++r){
        const size_t ri=regions[r];
        std::set<LabelType> & nodesEdges = dynamicNodes_[ri].edges_;
        for(size_t i=0;i<toCombine.size();++i){
            if(toCombine[i]!=newIndex){
                const bool found = static_cast<bool>(nodesEdges.find(toCombine[i])!=nodesEdges.end());
                CGP_ASSERT_OP(found,==,true);
                const size_t nErased = nodesEdges.erase(toCombine[i]);
                CGP_ASSERT_OP(nErased,==,1);
            }
        }
    }

    // call the registerd callbacks to merge the edges
    for(size_t i=0;i<toCombine.size();++i){
        if(toCombine[i]!=newIndex){
            this->callMergeEdgeCallbacks(newIndex,toCombine[i]);
        }
    }


    //CGP_ASSERT_OP(dynamicEdges_.size(),==,edgeUfd_.numberOfSets());
}

template<class LABEL_TYPE>
void MergeGraph<LABEL_TYPE>::searchLocalDoubleEdges(
    const MergeGraph<LABEL_TYPE>::NodeType & node , 
    MergeGraph<LABEL_TYPE>::DoubleMap & doubleMap,
    const LABEL_TYPE relabelFrom,
    const LABEL_TYPE relabelTo
){
    // loop over all edges of the new formed region
    for(
        typename std::set<LabelType>::const_iterator  edgeIter = node.edges_.begin();
        edgeIter!=node.edges_.end();
        ++edgeIter
    ){
        const LabelType outEdgeIndex = *edgeIter;
        //CGP_ASSERT_OP(outEdgeIndex,!=,edgeIndex);

        const EdgeType  oldEdge    = this->getEdge(outEdgeIndex);
        //const LabelType oldNodes[2]= {dynamicEdges_[outEdgeIndex].first,dynamicEdges_[outEdgeIndex].second };
        // do the relabling 
        LabelType newNodes[2]={
            oldEdge[0]==relabelFrom ? relabelTo : oldEdge[0] , 
            oldEdge[1]==relabelFrom ? relabelTo : oldEdge[1]
        };
        if(newNodes[1]<newNodes[0]){
            std::swap(newNodes[1],newNodes[0]);
        }
        const size_t  key = newNodes[0] + newNodes[1]*this->initNumberOfNodes();
        doubleMap[key].push_back(outEdgeIndex);
    }
}

template<class LABEL_TYPE>
void MergeGraph<LABEL_TYPE>::mergeRegions(const LABEL_TYPE toDeleteEdgeIndex){
    //std::cout<<"merge edge "<<toDeleteEdgeIndex<<"\n";
    const size_t preNumNodes = this->numberOfNodes();

    // assertions that edge is active and
    // its own repr.
    CGP_ASSERT_OP(reprEdge(toDeleteEdgeIndex),==,toDeleteEdgeIndex);
    CGP_ASSERT_OP(hasEdge(toDeleteEdgeIndex),==,true);

    const EdgeType toDeleteEdge = getEdge(toDeleteEdgeIndex);
    //const size_t nodes[2]= {dynamicEdges_[toDeleteEdgeIndex].first,dynamicEdges_[toDeleteEdgeIndex].second };
    std::vector<size_t> nodes(2);
    nodes[0]=toDeleteEdge[0];
    nodes[1]=toDeleteEdge[1];
    CGP_ASSERT_OP(nodes[0],!=,nodes[1]);

    for(size_t n=0;n<2;++n){
        // assertions that node is active and
        // its own repr.
        const size_t  ni=nodes[n];
        CGP_ASSERT_OP(reprNode(ni),==,ni);
        CGP_ASSERT_OP(hasNode(ni),==,true);
    }



    // merge the two nodes
    nodeUfd_.merge(nodes[0],nodes[1]);
    const LabelType newNodeRep    = reprNode(nodes[0]);
    const LabelType notNewNodeRep =  (newNodeRep == nodes[0] ? nodes[1] : nodes[0] );

    const size_t  edgeSizeRep    = dynamicNodes_[newNodeRep].numberOfEdges();
    const size_t  edgeSizeNotRep = dynamicNodes_[notNewNodeRep].numberOfEdges();

    // the new region wich is the result of the merge
    NodeType & newFormedNode = dynamicNodes_[newNodeRep];

    // merge the edges of the nodes
    newFormedNode.mergeEdges(dynamicNodes_[notNewNodeRep]);
    CGP_ASSERT_OP(newFormedNode.numberOfEdges(),==,edgeSizeRep+edgeSizeNotRep-1);

    // delete the old region
    dynamicNodes_.erase(notNewNodeRep);

    // delete the edge which has been between those two regions
    // which we merge (since this edge is the one getting deleted)
    newFormedNode.eraseEdge(toDeleteEdgeIndex);
    //dynamicEdges_.erase(toDeleteEdgeIndex);
    CGP_ASSERT_OP(newFormedNode.numberOfEdges(),==,edgeSizeRep+edgeSizeNotRep-2);


    // bevore processing with merging the edges we call the "merge" of the node maps
    // - we need to do this bevore any "merge" within the nodeMaps such that
    //   we can guarantee that the nodes maps are tidy when the edge-maps mergers
    //   are called
    this->callMergeNodeCallbacks(newNodeRep,notNewNodeRep);


    edgeUfd_.eraseElement(toDeleteEdgeIndex);

    // construct the "DoubleMap"
    // - if an vector in the map has a size >=2 
    //   this means that there are multiple edges
    //   between a pair of regions which needs to be merged
    DoubleMap doubleEdgeMap;
    this->searchLocalDoubleEdges(newFormedNode,doubleEdgeMap,notNewNodeRep,newNodeRep);

    // loop over the double map
    // if an vector in the map has a size >=2 
    // this means that there are multiple edges
    // between a pair of regions which needs to be merged
    for( typename DoubleMap::const_iterator dIter = doubleEdgeMap.begin();dIter!=doubleEdgeMap.end();++dIter){

        // if this vector has a size >=2 this means we have multiple
        // edges between 2 regions
        // the 2 regions are encoded in the key (dIter->first)
        // but we do not need them here
        const std::vector<LabelType> & edgeVec = dIter->second;
        if(edgeVec.size()>=2){

            // merge all these edges in the ufd and get the new representative
            //CGP_ASSERT_OP(hasEdge(toMergeEdgeIndex),==,true);
            const LabelType newEdgeRep = edgeUfd_.multiMerge(edgeVec.front(),edgeVec.begin()+1,edgeVec.end());
            //CGP_ASSERT_OP(hasEdge(toMergeEdgeIndex),==,false);
            // delte all edges which are not needed any more
            //  - edgeVec.size() -1 edges will be deleted 
            //  - (all edges except the new representative "newEdgeRep")
            // furthermore  the edge-sets all nodes adjacent to the "newFormedNode"
            // must be visited since they might refere to nodes which are deleted /merged
            for(size_t td=0;td<edgeVec.size();++td){

                // index of the edge which is considered for deletion
                const LabelType toMergeEdgeIndex = edgeVec[td];
                // delte this edge only if it is NOT the new representative edge
                if(toMergeEdgeIndex!=newEdgeRep){

                    // delete the edge from the new formed region
                    newFormedNode.edges_.erase(toMergeEdgeIndex);

                    //  not true any more
                    //CGP_ASSERT_OP(hasEdge(toMergeEdgeIndex),==,true);
                    

                    // at least one of the nodes of the edge "toMergeEdgeIndex" must be the "newFormedNode"
                    //  - we want to get the nodes adjacent to the "newFormedNode"
                    //CGP_ASSERT_OP(dynamicEdges_[toMergeEdgeIndex].hasNode(newNodeRep),==,true);
                    const size_t adjacentNodeIndex = getEdge(toMergeEdgeIndex).otherNode(newNodeRep);

                    dynamicNodes_[adjacentNodeIndex].eraseAndInsert(toMergeEdgeIndex,newEdgeRep);  
                    
                    // finaly delete the unneeded edge
                    //dynamicEdges_.erase(toMergeEdgeIndex);
                    //CGP_ASSERT_OP(hasEdge_OLD(toMergeEdgeIndex),==,false);
                }
            }
            CGP_ASSERT_OP(edgeVec.size(),==,2)
            CGP_ASSERT_OP(edgeVec[0],!=,toDeleteEdgeIndex);
            CGP_ASSERT_OP(edgeVec[1],!=,toDeleteEdgeIndex);
            CGP_ASSERT_OP(edgeVec[0],!=,edgeVec[1]);
            CGP_ASSERT_OP(hasEdge(newEdgeRep),==,true);
            // CALL CALLBACKS TO MERGE EDGES
            this->callMergeEdgeCallbacks(newEdgeRep, (newEdgeRep==edgeVec[0] ? edgeVec[1] : edgeVec[0]));
        } 
    }

    // CALL CALLBACKS TO ERASE EDGE
    this->callEraseEdgeCallbacks(toDeleteEdgeIndex);

    CGP_ASSERT_OP(dynamicNodes_.size(),==,preNumNodes-1);
    CGP_ASSERT_OP(nodeUfd_.numberOfSets(),==,preNumNodes-1);
    CGP_ASSERT_OP(this->numberOfNodes(),==,preNumNodes-1);
}




}



#endif //VIGRA_MERGE_GRAPH_HXX