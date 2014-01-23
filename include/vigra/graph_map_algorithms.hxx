/************************************************************************/
/*                                                                      */
/*     Copyright 2011-2012 Stefan Schmidt and Ullrich Koethe            */
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

/**
 * This header provides definitions of graph-related types
 * and optionally provides a gateway to popular graph libraries
 * (for now, BGL is supported).
 */

#ifndef VIGRA_GRAPH_MAP_ALGORITHMS_HXX
#define VIGRA_GRAPH_MAP_ALGORITHMS_HXX


#include <vigra/graphs.hxx>
#include <vigra/graph_generalization.hxx>
#include <vigra/multi_gridgraph.hxx>

#include <boost/iterator/transform_iterator.hpp>

namespace vigra{

    template<class GRAPH,class MAP>
    class EdgeMapIteratorHelper{
    public:
        typedef typename GraphMapTypeTraits<MAP>::Reference      Reference;
        typedef typename GraphMapTypeTraits<MAP>::ConstReference ConstReference;
        typedef typename GraphMapTypeTraits<MAP>::Value          Value;
    private:
        struct Transform{
            
            Transform(MAP & map)
            : map_(&map){
            }
            template<class ITEM>
            Reference operator()(const ITEM & item)const{
                return  map_->operator[](item);
            }
            mutable MAP * map_;
        };  
        struct ConstTransform{
            
            ConstTransform(const MAP & map)
            : map_(&map){
            }
            template<class ITEM>
            ConstReference operator()(const ITEM & item)const{
                return  map_->operator[](item);
            }
            const MAP * map_;
        }; 
    public:
        typedef  boost::transform_iterator< Transform,      typename GRAPH::EdgeIt,Reference      ,Value> iterator;
        typedef  boost::transform_iterator< ConstTransform, typename GRAPH::EdgeIt,ConstReference ,Value> const_iterator;
        static iterator
        begin(const GRAPH & g, MAP & map){
            Transform f(map);
            typename  GRAPH::EdgeIt iter = GraphIteratorAccessor<GRAPH>::edgesBegin(g);
            return iterator(iter,f);
        }
        static iterator
        end(const GRAPH & g, MAP & map){
            Transform f(map);
            typename  GRAPH::EdgeIt iter = GraphIteratorAccessor<GRAPH>::edgesEnd(g);
            return iterator(iter,f);
        }
        static const_iterator
        begin(const GRAPH & g, const MAP & map){
            ConstTransform f(map);
            typename  GRAPH::EdgeIt iter = GraphIteratorAccessor<GRAPH>::edgesBegin(g);
            return const_iterator(iter,f);
        }
        static const_iterator
        end(const GRAPH & g,const MAP & map){
            ConstTransform f(map);
            typename  GRAPH::EdgeIt iter = GraphIteratorAccessor<GRAPH>::edgesEnd(g);
            return const_iterator(iter,f);
        }
    private:

    };


    template<class GRAPH,class MAP>
    class NodeMapIteratorHelper{
    public:
        typedef typename GraphMapTypeTraits<MAP>::Reference      Reference;
        typedef typename GraphMapTypeTraits<MAP>::ConstReference ConstReference;
        typedef typename GraphMapTypeTraits<MAP>::Value          Value;
    private:
        struct Transform{
            
            Transform(MAP & map)
            : map_(&map){
            }
            template<class ITEM>
            Reference operator()(const ITEM & item)const{
                return  map_->operator[](item);
            }
            mutable MAP * map_;
        };  
        struct ConstTransform{
            
            ConstTransform(const MAP & map)
            : map_(&map){
            }
            template<class ITEM>
            ConstReference operator()(const ITEM & item)const{
                return  map_->operator[](item);
            }
            const MAP * map_;
        }; 
    public:
        typedef  boost::transform_iterator< Transform,      typename GRAPH::NodeIt,Reference      ,Value> iterator;
        typedef  boost::transform_iterator< ConstTransform, typename GRAPH::NodeIt,ConstReference ,Value> const_iterator;
        static iterator
        begin(const GRAPH & g, MAP & map){
            Transform f(map);
            typename  GRAPH::NodeIt iter = GraphIteratorAccessor<GRAPH>::nodesBegin(g);
            return iterator(iter,f);
        }
        static iterator
        end(const GRAPH & g, MAP & map){
            Transform f(map);
            typename  GRAPH::NodeIt iter = GraphIteratorAccessor<GRAPH>::nodesEnd(g);
            return iterator(iter,f);
        }
        static const_iterator
        begin(const GRAPH & g, const MAP & map){
            ConstTransform f(map);
            typename  GRAPH::NodeIt iter = GraphIteratorAccessor<GRAPH>::nodesBegin(g);
            return const_iterator(iter,f);
        }
        static const_iterator
        end(const GRAPH & g,const MAP & map){
            ConstTransform f(map);
            typename  GRAPH::NodeIt iter = GraphIteratorAccessor<GRAPH>::nodesEnd(g);
            return const_iterator(iter,f);
        }
    private:

    };


} // namespace vigra

#endif // VIGRA_GRAPH_MAP_ALGORITHMS_HXX
