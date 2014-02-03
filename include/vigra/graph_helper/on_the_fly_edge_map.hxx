/************************************************************************/
/*                                                                      */
/*                 Copyright 2011 by Ullrich Koethe                     */
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

#ifndef VIGRA_ON_THE_FLY_EDGE_MAP_HXX
#define VIGRA_ON_THE_FLY_EDGE_MAP_HXX

/*vigra*/
//#include <vigra/array_vector.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/graph_generalization.hxx>
#include <vigra/graphs.hxx>

namespace vigra{


template<class G,class NODE_MAP>
class OnTheFlyEdgeMap{

public:
    typedef G  Graph;
    typedef typename Graph::Node Node;
    typedef NODE_MAP  NodeMap;
    typedef typename  Graph::Edge      Key;
    typedef typename  NodeMap::Value   Value;
    typedef typename  NodeMap::Value   ConstReference;

    typedef Key             key_type;
    typedef Value           value_type;
    typedef ConstReference  const_reference;

    typedef boost::readable_property_map_tag category;

    OnTheFlyEdgeMap(const Graph & graph,const NodeMap & nodeMap)
    :   graph_(graph),
        nodeMap_(nodeMap){
    }

    ConstReference operator[](const Key & key){
        const Node u(graph_.u(key));
        const Node v(graph_.v(key));
        return (nodeMap_[u]+nodeMap_[v])/static_cast<Value>(2.0);
    }

    ConstReference operator[](const Key & key)const{
        const Node u(graph_.u(key));
        const Node v(graph_.v(key));
        return (nodeMap_[u]+nodeMap_[v])/static_cast<Value>(2.0);
    }
private:

    const Graph & graph_;
    const NodeMap & nodeMap_;

};







} // end namespace vigra


namespace boost{


    template<class G,class NODE_MAP>
    typename NODE_MAP::value_type
    get(const vigra::OnTheFlyEdgeMap<G,NODE_MAP> & pmap, const typename G::Edge & key){
        return pmap[key];
    }



    // property_traits class template
    template <class G,class NODE_MAP >
    struct property_traits< vigra::OnTheFlyEdgeMap<G,NODE_MAP> >
    {
        typedef vigra::OnTheFlyEdgeMap<G,NODE_MAP> PropMap;
        typedef typename PropMap::key_type    key_type;
        typedef typename PropMap::value_type  value_type; 
        //typedef typename PropMap::reference   reference;
        typedef typename PropMap::category    category;
    };

}


#endif //VIGRA_ON_THE_FLY_EDGE_MAP_HXX 