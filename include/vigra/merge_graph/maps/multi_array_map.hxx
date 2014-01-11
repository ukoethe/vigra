#ifndef VIGRA_MULTI_ARRAY_MAP_HXX
#define VIGRA_MULTI_ARRAY_MAP_HXX

/* boost */
#include <boost/function.hpp>


/* vigra */
#include <vigra/multi_array.hxx>
#include <vigra/accumulator.hxx>

/* vigra - merge graph */
#include <vigra/merge_graph/min_indexed_pq.hxx>


namespace vigra{
namespace view_maps {


    template<unsigned int DIM,class T>
    class SumMap;


    template<class T>
    class SumMap<0,T> : public MultiArrayView<1,T> {
    private:
        SumMap();                               // non empty-construction
        SumMap( const SumMap& other );          // non construction-copyable
        SumMap& operator=( const SumMap& );     // non copyable
    public:
        typedef MultiArrayView<1,T> ArrayViewType;
        typedef Int32 IdType;
        typedef T value_type;

        //const value_type & operator[](const IdType label)const{
        //  return this->operator(label);
        //}

        void merge(const IdType a,const IdType b){
            this->operator()(a)+=this->operator()(b);
        } 

        template<class CB>
        CB mergeCallback(){
            CB cb;
            cb = boost::bind(boost::mem_fn(&SumMap<0,T>::merge), this , _1,_2);
            return cb;
        }

        SumMap(ArrayViewType & array) 
        :   ArrayViewType(array){
        }
    };

    template<unsigned int DIM, class T>
    class SumMap : public MultiArrayView<DIM+1,T> {
    private:
        SumMap();                               // non empty-construction
        SumMap( const SumMap& other );          // non construction-copyable
        SumMap& operator=( const SumMap& );     // non copyable
    public:
        typedef MultiArrayView<DIM+1,T> ArrayViewType;
        typedef Int32 IdType;
        typedef T value_type;

        //const value_type & operator[](const IdType label)const{
        //  return this->operator(label);
        //}

        void merge(const IdType a,const IdType b){
            this->bindInner(a)+=this->bindInner(b);
        } 

        template<class CB>
        CB mergeCallback(){
            CB cb;
            cb = boost::bind(boost::mem_fn(&SumMap<DIM,T>::merge), this , _1,_2);
            return cb;
        }

        SumMap(ArrayViewType & array) 
        :   ArrayViewType(array){
        }
    };



    template<unsigned int DIM,class T,class WEIGHT_MAP>
    class WeightedMeanMap;


    template<class T,class WEIGHT_MAP>
    class WeightedMeanMap<0,T,WEIGHT_MAP> : public MultiArrayView<1,T> {
    private:
        WeightedMeanMap();                                      // non empty-construction
        WeightedMeanMap( const WeightedMeanMap& other );        // non construction-copyable
        WeightedMeanMap& operator=( const WeightedMeanMap& );   // non copyable
    public:
        typedef WEIGHT_MAP WeightMapType;
        typedef MultiArrayView<1,T> ArrayViewType;
        typedef Int32 IdType;
        typedef T value_type;

        void merge(const IdType a,const IdType b){
            const T va=this->operator()(a);
            const T vb=this->operator()(b);
            const T wa=weightMap_(a);
            const T wb=weightMap_(b);
            this->operator()(a)=(va*wa + vb*wb)/(wa+wb);
        } 


        template<class CB>
        CB mergeCallback(){
            return  boost::bind(boost::mem_fn(&WeightedMeanMap<0,T,WEIGHT_MAP>::merge), this , _1,_2);
        }

        WeightedMeanMap(ArrayViewType & array,const WeightMapType & weightMap) 
        :   ArrayViewType(array),
            weightMap_(weightMap){
        }
    private:
        const WeightMapType & weightMap_;
    };


    template<unsigned int DIM,class T,class WEIGHT_MAP>
    class WeightedMeanMap : public MultiArrayView<DIM+1,T> {
    private:
        WeightedMeanMap();                                      // non empty-construction
        WeightedMeanMap( const WeightedMeanMap& other );        // non construction-copyable
        WeightedMeanMap& operator=( const WeightedMeanMap& );   // non copyable
    public:
        typedef WEIGHT_MAP WeightMapType;
        typedef MultiArrayView<DIM+1,T> ArrayViewType;
        typedef Int32 IdType;
        typedef T value_type;

        void merge(const IdType a,const IdType b){
            const T wa=weightMap_(a);
            const T wb=weightMap_(b);
            vigra::MultiArrayView<DIM,T> va = this->bindInner(a);
            vigra::MultiArrayView<DIM,T>     vb = this->bindInner(b);
            //a=(a*wa+b*wb)
            va*=wa;
            vb*=wb;
            va+=vb;
            va/=(wa+wb);
            vb/=wb;
        } 


        template<class CB>
        CB mergeCallback(){
            return  boost::bind(boost::mem_fn(&WeightedMeanMap<DIM,T,WEIGHT_MAP>::merge), this , _1,_2);
        }

        WeightedMeanMap(ArrayViewType & array,const WeightMapType & weightMap) 
        :   ArrayViewType(array),
            weightMap_(weightMap){
        }
    private:
        const WeightMapType & weightMap_;
    };





    template<class T,class MERGE_GRAPH,class EDGE_MAP,class NODE_MAP>
    class MinWeightEdgeMap : public MultiArrayView<1,T> {
    private:
        MinWeightEdgeMap();                                         // non empty-construction
        MinWeightEdgeMap( const MinWeightEdgeMap& other );          // non construction-copyable
        MinWeightEdgeMap& operator=( const MinWeightEdgeMap& );     // non copyable
    public:
        typedef MERGE_GRAPH MergeGraphType;
        typedef EDGE_MAP    EdgeMapType;
        typedef NODE_MAP    NodeMapType;
        typedef MultiArrayView<1,T> ArrayViewType;
        typedef Int32 IdType;
        typedef T value_type;

        void mergeEdges(const IdType a,const IdType b){
            this->operator()(a)=getMixedWeight(a);
            pq_.deleteValue(b);
        } 
        void eraseEdge(const IdType label){
            pq_.deleteValue(label);
        }

        template<class CB>
        CB eraseEdgeCallback(){
            return  boost::bind(boost::mem_fn(&MinWeightEdgeMap<T,MERGE_GRAPH,EDGE_MAP,NODE_MAP>::eraseEdge), this , _1);
        }

        template<class CB>
        CB mergeEdgeCallback(){
            return  boost::bind(boost::mem_fn(&MinWeightEdgeMap<T,MERGE_GRAPH,EDGE_MAP,NODE_MAP>::mergeEdges), this , _1,_2);
        }

        MinWeightEdgeMap(const MergeGraphType & mergeGraph,ArrayViewType & array,const EdgeMapType & edgeMap,const NodeMapType & nodeMap) 
        :   ArrayViewType(array),
            mergeGraph_(mergeGraph),
            edgeMap_(edgeMap),
            nodeMap_(nodeMap),
            pq_(mergeGraph.initNumberOfEdges()){

            // initalize mixed weights  and initalize pq
            for(IdType l=0;l<mergeGraph_.initNumberOfEdges();++l){
                const T mixedWeight = getMixedWeight(l);
                this->operator()(l)=mixedWeight;
                pq_.insert(l,mixedWeight);
            }
        }

        IdType minWeightEdgeLabel(){
            IdType minLabel = pq_.minIndex();
            while(mergeGraph_.hasEdgeId(minLabel)==false){
                pq_.deleteValue(minLabel);
                IdType minLabel = pq_.minIndex();
            }
            return minLabel;
        }

    private:
        T getMixedWeight(const IdType label)const{
            return edgeMap_(label);
        }

        void changePqWeight(const IdType l,const T newWeigt){
            pq_.changeValue(l,newWeigt);
        }


        const MergeGraphType & mergeGraph_;
        const EdgeMapType    & edgeMap_;
        const NodeMapType    & nodeMap_;
        vigra::MinIndexedPQ<T>     pq_;
    };






    template<class T,class MERGE_GRAPH,class EDGE_MAP>
    class MinWeightEdgeMapSimple : public MultiArrayView<1,T> {
    private:
        MinWeightEdgeMapSimple();                                           // non empty-construction
        MinWeightEdgeMapSimple( const MinWeightEdgeMapSimple& other );          // non construction-copyable
        MinWeightEdgeMapSimple& operator=( const MinWeightEdgeMapSimple& );     // non copyable
    public:
        typedef MERGE_GRAPH MergeGraphType;
        typedef typename MergeGraphType::EdgeIt EdgeIt;
        typedef EDGE_MAP    EdgeMapType;
        typedef MultiArrayView<1,T> ArrayViewType;
        typedef Int32 IdType;
        typedef T value_type;

        void mergeEdges(const IdType a,const IdType b){
            this->operator()(a)=getMixedWeight(a);
            pq_.deleteValue(b);
        } 
        void eraseEdge(const IdType label){
            pq_.deleteValue(label);
        }

        template<class CB>
        CB eraseEdgeCallback(){
            return  boost::bind(boost::mem_fn(&MinWeightEdgeMapSimple<T,MERGE_GRAPH,EDGE_MAP>::eraseEdge), this , _1);
        }

        template<class CB>
        CB mergeEdgeCallback(){
            return  boost::bind(boost::mem_fn(&MinWeightEdgeMapSimple<T,MERGE_GRAPH,EDGE_MAP>::mergeEdges), this , _1,_2);
        }

        MinWeightEdgeMapSimple(const MergeGraphType & mergeGraph,ArrayViewType & array,const EdgeMapType & edgeMap) 
        :   ArrayViewType(array),
            mergeGraph_(mergeGraph),
            edgeMap_(edgeMap),
            pq_(mergeGraph.maxEdgeId()+1){



            // initalize mixed weights  and initalize pq
            for(IdType l=0;l<mergeGraph_.maxEdgeId()+1;++l){

                if(mergeGraph.edgeFromId(l)!=lemon::INVALID){

                    const T mixedWeight = getMixedWeight(l);
                    this->operator()(l)=mixedWeight;
                    pq_.insert(l,mixedWeight);
                }


            }
        }

        IdType minWeightEdgeLabel(){
            IdType minLabel = pq_.minIndex();
            while(mergeGraph_.hasEdgeId(minLabel)==false){
                pq_.deleteValue(minLabel);
                IdType minLabel = pq_.minIndex();
            }
            return minLabel;
        }

    private:
        T getMixedWeight(const IdType label)const{
            return edgeMap_(label);
        }

        void changePqWeight(const IdType l,const T newWeigt){
            pq_.changeValue(l,newWeigt);
        }


        const MergeGraphType & mergeGraph_;
        const EdgeMapType    & edgeMap_;
        vigra::MinIndexedPQ<T>     pq_;
    };

} // end namespace view_maps
} // end namespace vigra

#endif // VIGRA_MULTI_ARRAY_MAP_HXX