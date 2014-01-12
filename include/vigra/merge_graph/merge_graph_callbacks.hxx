#ifndef VIGRA_MERGE_GRAPH_CALLBACKS_HXX
#define VIGRA_MERGE_GRAPH_CALLBACKS_HXX

/* std library */
#include <vector>

/* boost */
#include <boost/function.hpp>
#include <boost/signals2.hpp>

namespace vigra {

/*
template<class ID_TYPE>
class MergeGraphCallbacks{
    public:
        typedef ID_TYPE  IdType;
        //callbacks typedefs
        typedef boost::function<void (const IdType,const IdType)>  MergeItemsCallBackType;
        typedef MergeItemsCallBackType                                   MergeNodeCallBackType;
        typedef MergeItemsCallBackType                                   MergeEdgeCallBackType;
        typedef boost::function<void (const IdType)>                  EraseEdgeCallBackType;

        MergeGraphCallbacks(){}


        template<class OBJ,class F>
        void registerMergeNodeCallBack(OBJ & obj,F  f){
            MergeNodeCallBackType internalF ;
            internalF = boost::bind(boost::mem_fn(f), &obj , _1,_2);
            mergeNodeCallbacks_.push_back(internalF);
        }
        template<class OBJ,class F>
        void registerMergeEdgeCallBack(OBJ & obj,F  f){
            MergeEdgeCallBackType internalF ;
            internalF = boost::bind(boost::mem_fn(f), &obj , _1,_2);
            mergeEdgeCallbacks_.push_back(internalF);
        }
        template<class OBJ,class F>
        void registerEraseEdgeCallBack(OBJ & obj,F  f){
            EraseEdgeCallBackType internalF ;
            internalF = boost::bind(boost::mem_fn(f), &obj , _1);
            eraseEdgeCallbacks_.push_back(internalF);
        }


        void registerMergeNodeCallBack(MergeNodeCallBackType  f){
            mergeNodeCallbacks_.push_back(f);
        }
        void registerMergeEdgeCallBack(MergeEdgeCallBackType  f){
            mergeEdgeCallbacks_.push_back(f);
        }
        void registerEraseEdgeCallBack(EraseEdgeCallBackType  f){
            eraseEdgeCallbacks_.push_back(f);
        }

    protected:
        void callMergeNodeCallbacks(const IdType a,const IdType b){
            for(size_t i=0;i<mergeNodeCallbacks_.size();++i)
                mergeNodeCallbacks_[i](a,b);
        }
        void callMergeEdgeCallbacks(const IdType a,const IdType b){
            for(size_t i=0;i<mergeEdgeCallbacks_.size();++i)
                mergeEdgeCallbacks_[i](a,b);
        }
        void callEraseEdgeCallbacks(const IdType a){
            for(size_t i=0;i<eraseEdgeCallbacks_.size();++i)
                eraseEdgeCallbacks_[i](a);
        }
    private:

        // callback vectors
        std::vector<MergeNodeCallBackType> mergeNodeCallbacks_;
        std::vector<MergeEdgeCallBackType> mergeEdgeCallbacks_;
        std::vector<EraseEdgeCallBackType> eraseEdgeCallbacks_;
};
*/


template<class NODE,class EDGE>
class NewMergeGraphCallbacks{
    public:
        //callbacks typedefs
        typedef boost::function<void (const NODE & ,const NODE &)>        MergeNodeCallBackType;
        typedef boost::function<void (const EDGE & ,const EDGE &)>        MergeEdgeCallBackType;
        typedef boost::function<void (const EDGE &)>                      EraseEdgeCallBackType;

        NewMergeGraphCallbacks(){}


        template<class OBJ,class F>
        void registerMergeNodeCallBack(OBJ & obj,F  f){
            MergeNodeCallBackType internalF ;
            internalF = boost::bind(boost::mem_fn(f), &obj , _1,_2);
            mergeNodeCallbacks_.push_back(internalF);
        }
        template<class OBJ,class F>
        void registerMergeEdgeCallBack(OBJ & obj,F  f){
            MergeEdgeCallBackType internalF ;
            internalF = boost::bind(boost::mem_fn(f), &obj , _1,_2);
            mergeEdgeCallbacks_.push_back(internalF);
        }
        template<class OBJ,class F>
        void registerEraseEdgeCallBack(OBJ & obj,F  f){
            EraseEdgeCallBackType internalF ;
            internalF = boost::bind(boost::mem_fn(f), &obj , _1);
            eraseEdgeCallbacks_.push_back(internalF);
        }


        void registerMergeNodeCallBack(MergeNodeCallBackType  f){
            mergeNodeCallbacks_.push_back(f);
        }
        void registerMergeEdgeCallBack(MergeEdgeCallBackType  f){
            mergeEdgeCallbacks_.push_back(f);
        }
        void registerEraseEdgeCallBack(EraseEdgeCallBackType  f){
            eraseEdgeCallbacks_.push_back(f);
        }

    protected:
        void callMergeNodeCallbacks(const NODE & a,const NODE & b){
            for(size_t i=0;i<mergeNodeCallbacks_.size();++i)
                mergeNodeCallbacks_[i](a,b);
        }
        void callMergeEdgeCallbacks(const EDGE & a,const EDGE & b){
            for(size_t i=0;i<mergeEdgeCallbacks_.size();++i)
                mergeEdgeCallbacks_[i](a,b);
        }
        void callEraseEdgeCallbacks(const EDGE & a){
            for(size_t i=0;i<eraseEdgeCallbacks_.size();++i)
                eraseEdgeCallbacks_[i](a);
        }
    private:

        // callback vectors
        std::vector<MergeNodeCallBackType> mergeNodeCallbacks_;
        std::vector<MergeEdgeCallBackType> mergeEdgeCallbacks_;
        std::vector<EraseEdgeCallBackType> eraseEdgeCallbacks_;
};


} // end namespace vigra



#endif //VIGRA_MERGE_GRAPH_CALLBACKS_HXX