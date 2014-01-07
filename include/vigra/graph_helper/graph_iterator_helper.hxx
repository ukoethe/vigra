#ifndef VIGRA_GRAPH_ITERATOR_HELPER_HXX
#define VIGRA_GRAPH_ITERATOR_HELPER_HXX

namespace vigra{
namespace detail{   



    template<class GRAPH,class INC_EDGE_ID_ITER>
    class GenericIncEdgeIt
    :  public boost::iterator_facade<
            RagIncEdgeIt
            T const,
            boost::bidirectional_traversal_tag
       >

    {
    private:
        friend class boost::iterator_core_access;
        typedef GRAPH Graph;
        typedef typename  Graph::Edge Edge;
         typedef typename Graph::index_type index_type;
    public:



        // default constructor
        GenericIncEdgeIt()
        :   graph_(NULL),
            idIter_(),
            currentEdge_(lemon::Invalid)  
        {
        }
        GenericIncEdgeIt(const GenericIncEdgeIt & other)
        :   graph_(other.graph_),
            idIter_(other.idIter_),
            currentEdge_(other.currentEdge_)  
        {
        }

        // Invalid constructor & conversion. 
        GenericIncEdgeIt(const lemon::Invalid & invalid)
        :   graph_(NULL),
            idIter_(),
            currentEdge_(lemon::Invalid)  
        {
        }
        GenericIncEdgeIt(const Graph & g)
        :   BaseIterType(0, GraphItemHelper<GRAPH,ITEM>::itemNum(g) ){
        }
        GenericIncEdgeIt(const Graph & g,const ITEM & item)
        :   BaseIterType(g.id(item), GraphItemHelper<GRAPH,ITEM>::itemNum(g) ){
        }


        bool isBegin()const{
            return idIter_.isBegin();
        }
        bool isEnd()const{
            return idIter_.isEnd();
        }

    private:
        bool equal(const ConstRepIter & other)const{
            return idIter_ == other.idIter_;
        }

        void increment(){
            ++idIter_;
        }

        void decrement(){
            ++idIter_;
        }

        const T & dereference()const{
            currentEdge_=graph_.edgeFromId(*idIter_);
            return currentEdge_;
        }

    private:
        GRAPH graph_;
        INC_EDGE_ID_ITER  idIter_;
        Edge currentEdge_;

    };


}
}


#endif