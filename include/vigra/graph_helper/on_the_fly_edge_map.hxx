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



} // end namespace vigra




#endif //VIGRA_ON_THE_FLY_EDGE_MAP_HXX 