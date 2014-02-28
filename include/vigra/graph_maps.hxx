/************************************************************************/
/*                                                                      */
/*   Copyright 2014 by Ullrich Koethe  and Thorsten Beier               */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/


#ifndef VIGRA_GRAPH_MAPS
#define VIGRA_GRAPH_MAPS

/*vigra*/
#include <vigra/multi_array.hxx>
#include <vigra/graph_generalization.hxx>
#include <vigra/graphs.hxx>


namespace vigra{

template<class T,class KEY,class REF,class CREF>
class DenseReferenceMap
: public MultiArray<1,T>
{
public:
    typedef KEY  Key;
    typedef T    Value;
    typedef REF  Reference;
    typedef CREF ConstReference;

    typedef Key             key_type;
    typedef Value           value_type;
    typedef Reference       reference;
    typedef ConstReference  const_reference;
    typedef boost::read_write_property_map_tag category;


    typedef typename MultiArray<1,T>::difference_type Shape1Type;

    DenseReferenceMap()
    : MultiArray<1,T>(){
    }
    DenseReferenceMap(const size_t maxKey)
    : MultiArray<1,T>(Shape1Type(maxKey+1)){
    }
    DenseReferenceMap(const size_t maxKey,ConstReference  value)
    : MultiArray<1,T>(Shape1Type(maxKey+1),value){
    }




    ConstReference operator[](const KEY & key)const{
        //return this->operator[](key.id());
        return MultiArray<1,T>::operator()(key.id());
    }
    Reference operator[](const KEY & key){
        return MultiArray<1,T>::operator()(key.id());
    }

    size_t size()const{
        return this->shape(0);
    }
protected:
    void assign(const size_t maxKey){
        this->reshape(Shape1Type(maxKey+1));
    }
private:
    // NONE
};


template<class GRAPH,class ITEM,class T,class REF,class CREF>
class DenseGraphItemReferenceMap
: public DenseReferenceMap<T,ITEM,REF,CREF>
{
    typedef GRAPH Graph;
    typedef ITEM  Item;
    typedef DenseReferenceMap<T,ITEM,REF,CREF> DenseReferenceMapType;
    typedef GraphItemHelper<Graph,ITEM> ItemHelper;
    typedef typename ItemHelper::ItemIt ItemIt;

public:
    DenseGraphItemReferenceMap()
    :   DenseReferenceMapType(){

    }
    DenseGraphItemReferenceMap(const Graph & g)
    :   DenseReferenceMapType(ItemHelper::itemNum(g)==0 ? 0: ItemHelper::maxItemId(g) ){

    }
    DenseGraphItemReferenceMap(const Graph & g,typename DenseReferenceMapType::ConstReference value)
    :   DenseReferenceMapType(ItemHelper::itemNum(g)==0 ? 0: ItemHelper::maxItemId(g)){

    }
    void assign(const Graph & g){
        DenseReferenceMapType::assign(ItemHelper::itemNum(g)==0 ? 0: ItemHelper::maxItemId(g));
    }
};


template<class GRAPH,class T,class REF= T & ,class CREF = const T & >
class DenseNodeReferenceMap
: public DenseGraphItemReferenceMap<GRAPH,typename GRAPH::Node,T,REF,CREF>
{
    typedef typename GRAPH::Node Node;
    typedef DenseGraphItemReferenceMap<GRAPH,Node,T,REF,CREF> DenseGraphItemReferenceMapType;
    public:
        DenseNodeReferenceMap()
        : DenseGraphItemReferenceMapType(){
        }
        DenseNodeReferenceMap(const GRAPH & g)
        : DenseGraphItemReferenceMapType(g){
        }
        DenseNodeReferenceMap(const GRAPH & g,typename DenseGraphItemReferenceMapType::ConstReference value)
        : DenseGraphItemReferenceMapType(g,value){
        }
};

template<class GRAPH,class T,class REF= T & ,class CREF = const T & >
class DenseEdgeReferenceMap
: public DenseGraphItemReferenceMap<GRAPH,typename GRAPH::Edge,T,REF,CREF>
{
    typedef typename GRAPH::Edge Edge;
    typedef DenseGraphItemReferenceMap<GRAPH,Edge,T,REF,CREF> DenseGraphItemReferenceMapType;
    public:
        DenseEdgeReferenceMap()
        : DenseGraphItemReferenceMapType(){
        }
        DenseEdgeReferenceMap(const GRAPH & g)
        : DenseGraphItemReferenceMapType(g){
        }
        DenseEdgeReferenceMap(const GRAPH & g,typename DenseGraphItemReferenceMapType::ConstReference value)
        : DenseGraphItemReferenceMapType(g,value){
        }
};

template<class GRAPH,class T,class REF= T & ,class CREF = const T & >
class DenseArcReferenceMap
: public DenseGraphItemReferenceMap<GRAPH,typename GRAPH::Arc,T,REF,CREF>
{
    typedef typename GRAPH::Arc Arc;
    typedef DenseGraphItemReferenceMap<GRAPH,Arc,T,REF,CREF> DenseGraphItemReferenceMapType;
    public:
        DenseArcReferenceMap()
        : DenseGraphItemReferenceMapType(){
        }
        DenseArcReferenceMap(const GRAPH & g)
        : DenseGraphItemReferenceMapType(g){
        }
        DenseArcReferenceMap(const GRAPH & g,typename DenseGraphItemReferenceMapType::ConstReference value)
        : DenseGraphItemReferenceMapType(g,value){
        }
};


template<class G,class NODE_MAP,class FUNCTOR,class RESULT>
class OnTheFlyEdgeMap{

public:
    typedef G  Graph;
    typedef typename Graph::Node Node;
    typedef NODE_MAP  NodeMap;
    typedef typename  Graph::Edge      Key;
    typedef RESULT   Value;
    typedef RESULT   ConstReference;

    typedef Key             key_type;
    typedef Value           value_type;
    typedef ConstReference  const_reference;

    typedef boost::readable_property_map_tag category;

    OnTheFlyEdgeMap(const Graph & graph,const NodeMap & nodeMap,FUNCTOR & f)
    :   graph_(graph),
        nodeMap_(nodeMap),
        f_(f){
    }

    ConstReference operator[](const Key & key){
        const Node u(graph_.u(key));
        const Node v(graph_.v(key));
        return f_(nodeMap_[u],nodeMap_[v]);
    }

    ConstReference operator[](const Key & key)const{
        const Node u(graph_.u(key));
        const Node v(graph_.v(key));
        return f_(nodeMap_[u],nodeMap_[v]);
    }
private:

    const Graph & graph_;
    const NodeMap & nodeMap_;
    FUNCTOR & f_;
};

template<class G,class EDGE_MAP_A,class EDGE_MAP_B,class FUNCTOR,class RESULT>
class BinaryOpEdgeMap{
public:
    typedef G  Graph;
    typedef typename Graph::Edge Key;
    typedef RESULT   Value;
    typedef RESULT   ConstReference;

    typedef Key             key_type;
    typedef Value           value_type;
    typedef ConstReference  const_reference;

    typedef boost::readable_property_map_tag category;

    BinaryOpEdgeMap(const Graph & graph,const EDGE_MAP_A & edgeMapA,const EDGE_MAP_B & edgeMapB,FUNCTOR & f)
    :   graph_(graph),
        edgeMapA_(edgeMapA),
        edgeMapB_(edgeMapB),
        f_(f){
    }
    ConstReference operator[](const Key & key){
        return f_(edgeMapA_[key],edgeMapB_[key]);
    }
    ConstReference operator[](const Key & key)const{
        return f_(edgeMapA_[key],edgeMapB_[key]);
    }
private:

    const Graph & graph_;
    const EDGE_MAP_A & edgeMapA_;
    const EDGE_MAP_B & edgeMapB_;
    FUNCTOR & f_;
};




} // end namespace vigra

#endif // VIGRA_GRAPH_MAPS