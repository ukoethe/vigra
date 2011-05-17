/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by  Ullrich Koethe and Rahul Nair         */
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
#ifndef RF_VISITORS_HXX
#define RF_VISITORS_HXX

#ifdef HasHDF5
# include "vigra/hdf5impex.hxx"
#endif // HasHDF5
#include <vigra/windows.h>
#include <iostream>
#include <iomanip>
#include <vigra/timing.hxx>

namespace vigra
{
namespace rf
{
/** \addtogroup MachineLearning Machine Learning
**/
//@{

/**
  	This namespace contains all classes and methods related to extracting information during 
	learning of the random forest. All Visitors share the same interface defined in 
	visitors::VisitorBase. The member methods are invoked at certain points of the main code in 
	the order they were supplied.
   	
	For the Random Forest the  Visitor concept is implemented as a statically linked list 
	(Using templates). Each Visitor object is encapsulated in a detail::VisitorNode object. The 
	VisitorNode object calls the Next Visitor after one of its visit() methods have terminated.
	
	To simplify usage create_visitor() factory methods are supplied.
	Use the create_visitor() method to supply visitor objects to the RandomForest::learn() method.
	It is possible to supply more than one visitor. They will then be invoked in serial order.

	The calculated information are stored as public data members of the class. - see documentation
	of the individual visitors
	
	While creating a new visitor the new class should therefore publicly inherit from this class 
	(i.e.: see visitors::OOB_Error).

	\code

 	  typedef xxx feature_t \\ replace xxx with whichever type
 	  typedef yyy label_t   \\ meme chose. 
 	  MultiArrayView<2, feature_t> f = get_some_features();
 	  MultiArrayView<2, label_t>   l = get_some_labels();
 	  RandomForest<> rf()
	
	  //calculate OOB Error
	  visitors::OOB_Error oob_v;
	  //calculate Variable Importance
	  visitors::VariableImportanceVisitor varimp_v;

 	  double oob_error = rf.learn(f, l, visitors::create_visitor(oob_v, varimp_v);
	  //the data can be found in the attributes of oob_v and varimp_v now
      
	\endcode
*/
namespace visitors
{
    
    
/** Base Class from which all Visitors derive. Can be used as a template to create new 
 * Visitors.
 */
class VisitorBase
{
    public:
    bool active_;   
    bool is_active()
    {
        return active_;
    }

    bool has_value()
    {
        return false;
    }

    VisitorBase()
        : active_(true)
    {}

    void deactivate()
    {
        active_ = false;
    }
    void activate()
    {
        active_ = true;
    }
    
    /** do something after the the Split has decided how to process the Region
     * (Stack entry)
     *
     * \param tree      reference to the tree that is currently being learned
     * \param split     reference to the split object
     * \param parent    current stack entry  which was used to decide the split
     * \param leftChild left stack entry that will be pushed
     * \param rightChild
     *                  right stack entry that will be pushed.
     * \param features  features matrix
     * \param labels    label matrix
     * \sa RF_Traits::StackEntry_t
     */
    template<class Tree, class Split, class Region, class Feature_t, class Label_t>
    void visit_after_split( Tree          & tree, 
                            Split         & split,
                            Region        & parent,
                            Region        & leftChild,
                            Region        & rightChild,
                            Feature_t     & features,
                            Label_t       & labels)
    {}
    
    /** do something after each tree has been learned
     *
     * \param rf        reference to the random forest object that called this
     *                  visitor
     * \param pr        reference to the preprocessor that processed the input
     * \param sm        reference to the sampler object
     * \param st        reference to the first stack entry
     * \param index     index of current tree
     */
    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {}
    
    /** do something after all trees have been learned
     *
     * \param rf        reference to the random forest object that called this
     *                  visitor
     * \param pr        reference to the preprocessor that processed the input
     */
    template<class RF, class PR>
    void visit_at_end(RF const & rf, PR const & pr)
    {}
	
    /** do something before learning starts 
     *
     * \param rf        reference to the random forest object that called this
     *                  visitor
     * \param pr        reference to the Processor class used.
     */
    template<class RF, class PR>
    void visit_at_beginning(RF const & rf, PR const & pr)
    {}
    /** do some thing while traversing tree after it has been learned 
     *  (external nodes)
     *
     * \param tr        reference to the tree object that called this visitor
     * \param index     index in the topology_ array we currently are at
     * \param node_t    type of node we have (will be e_.... - )
     * \param weight    Node weight of current node. 
     * \sa  NodeTags;
     *
     * you can create the node by using a switch on node_tag and using the 
     * corresponding Node objects. Or - if you do not care about the type 
     * use the Nodebase class.
     */
    template<class TR, class IntT, class TopT,class Feat>
    void visit_external_node(TR & tr, IntT index, TopT node_t,Feat & features)
    {}
    
    /** do something when visiting a internal node after it has been learned
     *
     * \sa visit_external_node
     */
    template<class TR, class IntT, class TopT,class Feat>
    void visit_internal_node(TR & tr, IntT index, TopT node_t,Feat & features)
    {}

    /** return a double value.  The value of the first 
     * visitor encountered that has a return value is returned with the
     * RandomForest::learn() method - or -1.0 if no return value visitor
     * existed. This functionality basically only exists so that the 
     * OOB - visitor can return the oob error rate like in the old version 
     * of the random forest.
     */
    double return_val()
    {
        return -1.0;
    }
};


/** Last Visitor that should be called to stop the recursion.
 */
class StopVisiting: public VisitorBase
{
    public:
    bool has_value()
    {
        return true;
    }
    double return_val()
    {
        return -1.0;
    }
};
namespace detail
{
/** Container elements of the statically linked Visitor list.
 *
 * use the create_visitor() factory functions to create visitors up to size 10;
 *
 */
template <class Visitor, class Next = StopVisiting>
class VisitorNode
{
    public:
    
    StopVisiting    stop_;
    Next            next_;
    Visitor &       visitor_;   
    VisitorNode(Visitor & visitor, Next & next) 
    : 
        next_(next), visitor_(visitor)
    {}

    VisitorNode(Visitor &  visitor) 
    : 
        next_(stop_), visitor_(visitor)
    {}

    template<class Tree, class Split, class Region, class Feature_t, class Label_t>
    void visit_after_split( Tree          & tree, 
                            Split         & split,
                            Region        & parent,
                            Region        & leftChild,
                            Region        & rightChild,
                            Feature_t     & features,
                            Label_t       & labels)
    {
        if(visitor_.is_active())
            visitor_.visit_after_split(tree, split, 
                                       parent, leftChild, rightChild,
                                       features, labels);
        next_.visit_after_split(tree, split, parent, leftChild, rightChild,
                                features, labels);
    }

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        if(visitor_.is_active())
            visitor_.visit_after_tree(rf, pr, sm, st, index);
        next_.visit_after_tree(rf, pr, sm, st, index);
    }

    template<class RF, class PR>
    void visit_at_beginning(RF & rf, PR & pr)
    {
        if(visitor_.is_active())
            visitor_.visit_at_beginning(rf, pr);
        next_.visit_at_beginning(rf, pr);
    }
    template<class RF, class PR>
    void visit_at_end(RF & rf, PR & pr)
    {
        if(visitor_.is_active())
            visitor_.visit_at_end(rf, pr);
        next_.visit_at_end(rf, pr);
    }
    
    template<class TR, class IntT, class TopT,class Feat>
    void visit_external_node(TR & tr, IntT & index, TopT & node_t,Feat & features)
    {
        if(visitor_.is_active())
            visitor_.visit_external_node(tr, index, node_t,features);
        next_.visit_external_node(tr, index, node_t,features);
    }
    template<class TR, class IntT, class TopT,class Feat>
    void visit_internal_node(TR & tr, IntT & index, TopT & node_t,Feat & features)
    {
        if(visitor_.is_active())
            visitor_.visit_internal_node(tr, index, node_t,features);
        next_.visit_internal_node(tr, index, node_t,features);
    }

    double return_val()
    {
        if(visitor_.is_active() && visitor_.has_value())
            return visitor_.return_val();
        return next_.return_val();
    }
};

} //namespace detail

//////////////////////////////////////////////////////////////////////////////
//  Visitor Factory function up to 10 visitors                              //
//////////////////////////////////////////////////////////////////////////////

/** factory method to to be used with RandomForest::learn()
 */
template<class A>
detail::VisitorNode<A>
create_visitor(A & a)
{
   typedef detail::VisitorNode<A> _0_t;
   _0_t _0(a);
   return _0;
}


/** factory method to to be used with RandomForest::learn()
 */
template<class A, class B>
detail::VisitorNode<A, detail::VisitorNode<B> >
create_visitor(A & a, B & b)
{
   typedef detail::VisitorNode<B> _1_t;
   _1_t _1(b);
   typedef detail::VisitorNode<A, _1_t> _0_t;
   _0_t _0(a, _1);
   return _0;
}


/** factory method to to be used with RandomForest::learn()
 */
template<class A, class B, class C>
detail::VisitorNode<A, detail::VisitorNode<B, detail::VisitorNode<C> > >
create_visitor(A & a, B & b, C & c)
{
   typedef detail::VisitorNode<C> _2_t;
   _2_t _2(c);
   typedef detail::VisitorNode<B, _2_t> _1_t;
   _1_t _1(b, _2);
   typedef detail::VisitorNode<A, _1_t> _0_t;
   _0_t _0(a, _1);
   return _0;
}


/** factory method to to be used with RandomForest::learn()
 */
template<class A, class B, class C, class D>
detail::VisitorNode<A, detail::VisitorNode<B, detail::VisitorNode<C, 
    detail::VisitorNode<D> > > >
create_visitor(A & a, B & b, C & c, D & d)
{
   typedef detail::VisitorNode<D> _3_t;
   _3_t _3(d);
   typedef detail::VisitorNode<C, _3_t> _2_t;
   _2_t _2(c, _3);
   typedef detail::VisitorNode<B, _2_t> _1_t;
   _1_t _1(b, _2);
   typedef detail::VisitorNode<A, _1_t> _0_t;
   _0_t _0(a, _1);
   return _0;
}


/** factory method to to be used with RandomForest::learn()
 */
template<class A, class B, class C, class D, class E>
detail::VisitorNode<A, detail::VisitorNode<B, detail::VisitorNode<C, 
    detail::VisitorNode<D, detail::VisitorNode<E> > > > >
create_visitor(A & a, B & b, C & c, 
               D & d, E & e)
{
   typedef detail::VisitorNode<E> _4_t;
   _4_t _4(e);
   typedef detail::VisitorNode<D, _4_t> _3_t;
   _3_t _3(d, _4);
   typedef detail::VisitorNode<C, _3_t> _2_t;
   _2_t _2(c, _3);
   typedef detail::VisitorNode<B, _2_t> _1_t;
   _1_t _1(b, _2);
   typedef detail::VisitorNode<A, _1_t> _0_t;
   _0_t _0(a, _1);
   return _0;
}


/** factory method to to be used with RandomForest::learn()
 */
template<class A, class B, class C, class D, class E,
         class F>
detail::VisitorNode<A, detail::VisitorNode<B, detail::VisitorNode<C, 
    detail::VisitorNode<D, detail::VisitorNode<E, detail::VisitorNode<F> > > > > >
create_visitor(A & a, B & b, C & c, 
               D & d, E & e, F & f)
{
   typedef detail::VisitorNode<F> _5_t;
   _5_t _5(f);
   typedef detail::VisitorNode<E, _5_t> _4_t;
   _4_t _4(e, _5);
   typedef detail::VisitorNode<D, _4_t> _3_t;
   _3_t _3(d, _4);
   typedef detail::VisitorNode<C, _3_t> _2_t;
   _2_t _2(c, _3);
   typedef detail::VisitorNode<B, _2_t> _1_t;
   _1_t _1(b, _2);
   typedef detail::VisitorNode<A, _1_t> _0_t;
   _0_t _0(a, _1);
   return _0;
}


/** factory method to to be used with RandomForest::learn()
 */
template<class A, class B, class C, class D, class E,
         class F, class G>
detail::VisitorNode<A, detail::VisitorNode<B, detail::VisitorNode<C, 
    detail::VisitorNode<D, detail::VisitorNode<E, detail::VisitorNode<F, 
    detail::VisitorNode<G> > > > > > >
create_visitor(A & a, B & b, C & c, 
               D & d, E & e, F & f, G & g)
{
   typedef detail::VisitorNode<G> _6_t;
   _6_t _6(g);
   typedef detail::VisitorNode<F, _6_t> _5_t;
   _5_t _5(f, _6);
   typedef detail::VisitorNode<E, _5_t> _4_t;
   _4_t _4(e, _5);
   typedef detail::VisitorNode<D, _4_t> _3_t;
   _3_t _3(d, _4);
   typedef detail::VisitorNode<C, _3_t> _2_t;
   _2_t _2(c, _3);
   typedef detail::VisitorNode<B, _2_t> _1_t;
   _1_t _1(b, _2);
   typedef detail::VisitorNode<A, _1_t> _0_t;
   _0_t _0(a, _1);
   return _0;
}


/** factory method to to be used with RandomForest::learn()
 */
template<class A, class B, class C, class D, class E,
         class F, class G, class H>
detail::VisitorNode<A, detail::VisitorNode<B, detail::VisitorNode<C, 
    detail::VisitorNode<D, detail::VisitorNode<E, detail::VisitorNode<F, 
    detail::VisitorNode<G, detail::VisitorNode<H> > > > > > > >
create_visitor(A & a, B & b, C & c, 
               D & d, E & e, F & f, 
               G & g, H & h)
{
   typedef detail::VisitorNode<H> _7_t;
   _7_t _7(h);
   typedef detail::VisitorNode<G, _7_t> _6_t;
   _6_t _6(g, _7);
   typedef detail::VisitorNode<F, _6_t> _5_t;
   _5_t _5(f, _6);
   typedef detail::VisitorNode<E, _5_t> _4_t;
   _4_t _4(e, _5);
   typedef detail::VisitorNode<D, _4_t> _3_t;
   _3_t _3(d, _4);
   typedef detail::VisitorNode<C, _3_t> _2_t;
   _2_t _2(c, _3);
   typedef detail::VisitorNode<B, _2_t> _1_t;
   _1_t _1(b, _2);
   typedef detail::VisitorNode<A, _1_t> _0_t;
   _0_t _0(a, _1);
   return _0;
}


/** factory method to to be used with RandomForest::learn()
 */
template<class A, class B, class C, class D, class E,
         class F, class G, class H, class I>
detail::VisitorNode<A, detail::VisitorNode<B, detail::VisitorNode<C, 
    detail::VisitorNode<D, detail::VisitorNode<E, detail::VisitorNode<F, 
    detail::VisitorNode<G, detail::VisitorNode<H, detail::VisitorNode<I> > > > > > > > >
create_visitor(A & a, B & b, C & c, 
               D & d, E & e, F & f, 
               G & g, H & h, I & i)
{
   typedef detail::VisitorNode<I> _8_t;
   _8_t _8(i);
   typedef detail::VisitorNode<H, _8_t> _7_t;
   _7_t _7(h, _8);
   typedef detail::VisitorNode<G, _7_t> _6_t;
   _6_t _6(g, _7);
   typedef detail::VisitorNode<F, _6_t> _5_t;
   _5_t _5(f, _6);
   typedef detail::VisitorNode<E, _5_t> _4_t;
   _4_t _4(e, _5);
   typedef detail::VisitorNode<D, _4_t> _3_t;
   _3_t _3(d, _4);
   typedef detail::VisitorNode<C, _3_t> _2_t;
   _2_t _2(c, _3);
   typedef detail::VisitorNode<B, _2_t> _1_t;
   _1_t _1(b, _2);
   typedef detail::VisitorNode<A, _1_t> _0_t;
   _0_t _0(a, _1);
   return _0;
}

/** factory method to to be used with RandomForest::learn()
 */
template<class A, class B, class C, class D, class E,
         class F, class G, class H, class I, class J>
detail::VisitorNode<A, detail::VisitorNode<B, detail::VisitorNode<C, 
    detail::VisitorNode<D, detail::VisitorNode<E, detail::VisitorNode<F, 
    detail::VisitorNode<G, detail::VisitorNode<H, detail::VisitorNode<I,
    detail::VisitorNode<J> > > > > > > > > >
create_visitor(A & a, B & b, C & c, 
               D & d, E & e, F & f, 
               G & g, H & h, I & i,
               J & j)
{
   typedef detail::VisitorNode<J> _9_t;
   _9_t _9(j);
   typedef detail::VisitorNode<I, _9_t> _8_t;
   _8_t _8(i, _9);
   typedef detail::VisitorNode<H, _8_t> _7_t;
   _7_t _7(h, _8);
   typedef detail::VisitorNode<G, _7_t> _6_t;
   _6_t _6(g, _7);
   typedef detail::VisitorNode<F, _6_t> _5_t;
   _5_t _5(f, _6);
   typedef detail::VisitorNode<E, _5_t> _4_t;
   _4_t _4(e, _5);
   typedef detail::VisitorNode<D, _4_t> _3_t;
   _3_t _3(d, _4);
   typedef detail::VisitorNode<C, _3_t> _2_t;
   _2_t _2(c, _3);
   typedef detail::VisitorNode<B, _2_t> _1_t;
   _1_t _1(b, _2);
   typedef detail::VisitorNode<A, _1_t> _0_t;
   _0_t _0(a, _1);
   return _0;
}

//////////////////////////////////////////////////////////////////////////////
// Visitors of communal interest.                                           //
//////////////////////////////////////////////////////////////////////////////


/** Visitor to gain information, later needed for online learning.
 */

class OnlineLearnVisitor: public VisitorBase
{
public:
    //Set if we adjust thresholds
    bool adjust_thresholds;
    //Current tree id
    int tree_id;
    //Last node id for finding parent
    int last_node_id;
    //Need to now the label for interior node visiting
    vigra::Int32 current_label;
    //marginal distribution for interior nodes
	//
	OnlineLearnVisitor():
		adjust_thresholds(false), tree_id(0), last_node_id(0), current_label(0)
	{}
    struct MarginalDistribution
    {
        ArrayVector<Int32> leftCounts;
        Int32 leftTotalCounts;
        ArrayVector<Int32> rightCounts;
        Int32 rightTotalCounts;
        double gap_left;
        double gap_right;
    };
    typedef ArrayVector<vigra::Int32> IndexList;

    //All information for one tree
    struct TreeOnlineInformation
    {
        std::vector<MarginalDistribution> mag_distributions;
        std::vector<IndexList> index_lists;
        //map for linear index of mag_distiributions
        std::map<int,int> interior_to_index;
        //map for linear index of index_lists
        std::map<int,int> exterior_to_index;
    };

    //All trees
    std::vector<TreeOnlineInformation> trees_online_information;

    /** Initilize, set the number of trees
     */
    template<class RF,class PR>
    void visit_at_beginning(RF & rf,const PR & pr)
    {
        tree_id=0;
        trees_online_information.resize(rf.options_.tree_count_);
    }

    /** Reset a tree
     */
    void reset_tree(int tree_id)
    {
        trees_online_information[tree_id].mag_distributions.clear();
        trees_online_information[tree_id].index_lists.clear();
        trees_online_information[tree_id].interior_to_index.clear();
        trees_online_information[tree_id].exterior_to_index.clear();
    }

    /** simply increase the tree count
    */
    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        tree_id++;
    }
	
    template<class Tree, class Split, class Region, class Feature_t, class Label_t>
    void visit_after_split( Tree  	      & tree, 
			    Split         & split,
                            Region       & parent,
                            Region        & leftChild,
                            Region        & rightChild,
                            Feature_t     & features,
                            Label_t       & labels)
    {
        int linear_index;
        int addr=tree.topology_.size();
        if(split.createNode().typeID() == i_ThresholdNode)
        {
            if(adjust_thresholds)
            {
                //Store marginal distribution
                linear_index=trees_online_information[tree_id].mag_distributions.size();
                trees_online_information[tree_id].interior_to_index[addr]=linear_index;
                trees_online_information[tree_id].mag_distributions.push_back(MarginalDistribution());

                trees_online_information[tree_id].mag_distributions.back().leftCounts=leftChild.classCounts_;
                trees_online_information[tree_id].mag_distributions.back().rightCounts=rightChild.classCounts_;

                trees_online_information[tree_id].mag_distributions.back().leftTotalCounts=leftChild.size_;
                trees_online_information[tree_id].mag_distributions.back().rightTotalCounts=rightChild.size_;
                //Store the gap
                double gap_left,gap_right;
                int i;
                gap_left=features(leftChild[0],split.bestSplitColumn());
                for(i=1;i<leftChild.size();++i)
                    if(features(leftChild[i],split.bestSplitColumn())>gap_left)
                        gap_left=features(leftChild[i],split.bestSplitColumn());
                gap_right=features(rightChild[0],split.bestSplitColumn());
                for(i=1;i<rightChild.size();++i)
                    if(features(rightChild[i],split.bestSplitColumn())<gap_right)
                        gap_right=features(rightChild[i],split.bestSplitColumn());
                trees_online_information[tree_id].mag_distributions.back().gap_left=gap_left;
                trees_online_information[tree_id].mag_distributions.back().gap_right=gap_right;
            }
        }
        else
        {
            //Store index list
            linear_index=trees_online_information[tree_id].index_lists.size();
            trees_online_information[tree_id].exterior_to_index[addr]=linear_index;

            trees_online_information[tree_id].index_lists.push_back(IndexList());

            trees_online_information[tree_id].index_lists.back().resize(parent.size_,0);
            std::copy(parent.begin_,parent.end_,trees_online_information[tree_id].index_lists.back().begin());
        }
    }
    void add_to_index_list(int tree,int node,int index)
    {
        if(!this->active_)
            return;
        TreeOnlineInformation &ti=trees_online_information[tree];
        ti.index_lists[ti.exterior_to_index[node]].push_back(index);
    }
    void move_exterior_node(int src_tree,int src_index,int dst_tree,int dst_index)
    {
        if(!this->active_)
            return;
        trees_online_information[dst_tree].exterior_to_index[dst_index]=trees_online_information[src_tree].exterior_to_index[src_index];
        trees_online_information[src_tree].exterior_to_index.erase(src_index);
    }
    /** do something when visiting a internal node during getToLeaf
     *
     * remember as last node id, for finding the parent of the last external node
     * also: adjust class counts and borders
     */
    template<class TR, class IntT, class TopT,class Feat>
        void visit_internal_node(TR & tr, IntT index, TopT node_t,Feat & features)
        {
            last_node_id=index;
            if(adjust_thresholds)
            {
                vigra_assert(node_t==i_ThresholdNode,"We can only visit threshold nodes");
                //Check if we are in the gap
                double value=features(0, Node<i_ThresholdNode>(tr.topology_,tr.parameters_,index).column());
                TreeOnlineInformation &ti=trees_online_information[tree_id];
                MarginalDistribution &m=ti.mag_distributions[ti.interior_to_index[index]];
                if(value>m.gap_left && value<m.gap_right)
                {
                    //Check which site we want to go
                    if(m.leftCounts[current_label]/double(m.leftTotalCounts)>m.rightCounts[current_label]/double(m.rightTotalCounts))
                    {
                        //We want to go left
                        m.gap_left=value;
                    }
                    else
                    {
                        //We want to go right
                        m.gap_right=value;
                    }
                    Node<i_ThresholdNode>(tr.topology_,tr.parameters_,index).threshold()=(m.gap_right+m.gap_left)/2.0;
                }
                //Adjust class counts
                if(value>Node<i_ThresholdNode>(tr.topology_,tr.parameters_,index).threshold())
                {
                    ++m.rightTotalCounts;
                    ++m.rightCounts[current_label];
                }
                else
                {
                    ++m.leftTotalCounts;
                    ++m.rightCounts[current_label];
                }
            }
        }
    /** do something when visiting a extern node during getToLeaf
     * 
     * Store the new index!
     */
};

//////////////////////////////////////////////////////////////////////////////
// Out of Bag Error estimates                                               //
//////////////////////////////////////////////////////////////////////////////


/** Visitor that calculates the oob error of each individual randomized
 * decision tree. 
 *
 * After training a tree, all those samples that are OOB for this particular tree
 * are put down the tree and the error estimated. 
 * the per tree oob error is the average of the individual error estimates. 
 * (oobError = average error of one randomized tree)
 * Note: This is Not the OOB - Error estimate suggested by Breiman (See OOB_Error 
 * visitor)
 */
class OOB_PerTreeError:public VisitorBase
{
public:
	/** Average error of one randomized decision tree
	 */
    double oobError;

    int totalOobCount;
    ArrayVector<int> oobCount,oobErrorCount;

    OOB_PerTreeError()
    : oobError(0.0),
      totalOobCount(0)
    {}


    bool has_value()
    {
        return true;
    }


    /** does the basic calculation per tree*/
    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(    RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        //do the first time called.
        if(int(oobCount.size()) != rf.ext_param_.row_count_)
        {
            oobCount.resize(rf.ext_param_.row_count_, 0);
            oobErrorCount.resize(rf.ext_param_.row_count_, 0);
        }
        // go through the samples
        for(int l = 0; l < rf.ext_param_.row_count_; ++l)
        {
            // if the lth sample is oob...
            if(!sm.is_used()[l])
            {
                ++oobCount[l];
                if(     rf.tree(index)
                            .predictLabel(rowVector(pr.features(), l)) 
                    !=  pr.response()(l,0))
                {
                    ++oobErrorCount[l];
                }
            }

        }
    }

    /** Does the normalisation
     */
    template<class RF, class PR>
    void visit_at_end(RF & rf, PR & pr)
    {
        // do some normalisation
        for(int l=0; l < (int)rf.ext_param_.row_count_; ++l)
        {
            if(oobCount[l])
            {
                oobError += double(oobErrorCount[l]) / oobCount[l];
                ++totalOobCount;
            }
        } 
        oobError/=totalOobCount;
    }
    
};

/** Visitor that calculates the oob error of the ensemble
 *  This rate should be used to estimate the crossvalidation 
 *  error rate.
 *  Here each sample is put down those trees, for which this sample
 *  is OOB i.e. if sample #1 is  OOB for trees 1, 3 and 5 we calculate
 *  the output using the ensemble consisting only of trees 1 3 and 5. 
 *
 *  Using normal bagged sampling each sample is OOB for approx. 33% of trees
 *  The error rate obtained as such therefore corresponds to crossvalidation
 *  rate obtained using a ensemble containing 33% of the trees.
 */
class OOB_Error : public VisitorBase
{
    typedef MultiArrayShape<2>::type Shp;
    int class_count;
    bool is_weighted;
    MultiArray<2,double> tmp_prob;
    public:

    MultiArray<2, double>       prob_oob; 
	/** Ensemble oob error rate
	 */
    double                      oob_breiman;

    MultiArray<2, double>       oobCount;
    ArrayVector< int>           indices; 
    OOB_Error() : VisitorBase(), oob_breiman(0.0) {}
#ifdef HasHDF5
    void save(std::string filen, std::string pathn)
    {
        if(*(pathn.end()-1) != '/')
            pathn += "/";
        const char* filename = filen.c_str();
        MultiArray<2, double> temp(Shp(1,1), 0.0); 
        temp[0] = oob_breiman;
        writeHDF5(filename, (pathn + "breiman_error").c_str(), temp);
    }
#endif
    // negative value if sample was ib, number indicates how often.
    //  value >=0  if sample was oob, 0 means fail 1, corrrect

    template<class RF, class PR>
    void visit_at_beginning(RF & rf, PR & pr)
    {
        class_count = rf.class_count();
        tmp_prob.reshape(Shp(1, class_count), 0); 
        prob_oob.reshape(Shp(rf.ext_param().row_count_,class_count), 0);
        is_weighted = rf.options().predict_weighted_;
        indices.resize(rf.ext_param().row_count_);
        if(int(oobCount.size()) != rf.ext_param_.row_count_)
        {
            oobCount.reshape(Shp(rf.ext_param_.row_count_, 1), 0);
        }
        for(int ii = 0; ii < rf.ext_param().row_count_; ++ii)
        {
            indices[ii] = ii;
        }
    }

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        // go through the samples
        int total_oob =0;
        int wrong_oob =0;
        if(rf.ext_param_.actual_msample_ < pr.features().shape(0)- 10000)
        {
            ArrayVector<int> oob_indices;
            ArrayVector<int> cts(class_count, 0);
            std::random_shuffle(indices.begin(), indices.end());
            for(int ii = 0; ii < rf.ext_param_.row_count_; ++ii)
            {
                if(!sm.is_used()[indices[ii]] && cts[pr.response()(indices[ii], 0)] < 40000)
                {
                    oob_indices.push_back(indices[ii]);
                    ++cts[pr.response()(indices[ii], 0)];
                }
            }
            for(int ll = 0; ll < oob_indices.size(); ++ll)
            {
                // update number of trees in which current sample is oob
                ++oobCount[oob_indices[ll]];

                // update number of oob samples in this tree.
                ++total_oob; 
                // get the predicted votes ---> tmp_prob;
                int pos =  rf.tree(index).getToLeaf(rowVector(pr.features(),oob_indices[ll]));
                Node<e_ConstProbNode> node ( rf.tree(index).topology_, 
                                                    rf.tree(index).parameters_,
                                                    pos);
                tmp_prob.init(0); 
                for(int ii = 0; ii < class_count; ++ii)
                {
                    tmp_prob[ii] = node.prob_begin()[ii];
                }
                if(is_weighted)
                {
                    for(int ii = 0; ii < class_count; ++ii)
                        tmp_prob[ii] = tmp_prob[ii] * (*(node.prob_begin()-1));
                }
                rowVector(prob_oob, oob_indices[ll]) += tmp_prob;
                int label = argMax(tmp_prob); 
                
            }
        }else
        {
            for(int ll = 0; ll < rf.ext_param_.row_count_; ++ll)
            {
                // if the lth sample is oob...
                if(!sm.is_used()[ll])
                {
                    // update number of trees in which current sample is oob
                    ++oobCount[ll];

                    // update number of oob samples in this tree.
                    ++total_oob; 
                    // get the predicted votes ---> tmp_prob;
                    int pos =  rf.tree(index).getToLeaf(rowVector(pr.features(),ll));
                    Node<e_ConstProbNode> node ( rf.tree(index).topology_, 
                                                        rf.tree(index).parameters_,
                                                        pos);
                    tmp_prob.init(0); 
                    for(int ii = 0; ii < class_count; ++ii)
                    {
                        tmp_prob[ii] = node.prob_begin()[ii];
                    }
                    if(is_weighted)
                    {
                        for(int ii = 0; ii < class_count; ++ii)
                            tmp_prob[ii] = tmp_prob[ii] * (*(node.prob_begin()-1));
                    }
                    rowVector(prob_oob, ll) += tmp_prob;
                    int label = argMax(tmp_prob); 
                    
                }
            }
        }
        // go through the ib samples; 
    }

    /** Normalise variable importance after the number of trees is known.
     */
    template<class RF, class PR>
    void visit_at_end(RF & rf, PR & pr)
    {
        // ullis original metric and breiman style stuff
        int totalOobCount =0;
        int breimanstyle = 0;
        for(int ll=0; ll < (int)rf.ext_param_.row_count_; ++ll)
        {
            if(oobCount[ll])
            {
                if(argMax(rowVector(prob_oob, ll)) != pr.response()(ll, 0))
                   ++breimanstyle;
                ++totalOobCount;
            }
        }
        oob_breiman = double(breimanstyle)/totalOobCount; 
    }
};


/** Visitor that calculates different OOB error statistics
 */
class CompleteOOBInfo : public VisitorBase
{
    typedef MultiArrayShape<2>::type Shp;
    int class_count;
    bool is_weighted;
    MultiArray<2,double> tmp_prob;
    public:

    /** OOB Error rate of each individual tree
	 */
    MultiArray<2, double>       oob_per_tree;
	/** Mean of oob_per_tree
	 */
    double                      oob_mean;
	/**Standard deviation of oob_per_tree
	 */
    double                      oob_std;
    
    MultiArray<2, double>       prob_oob; 
	/** Ensemble OOB error
	 *
	 * \sa OOB_Error
	 */
    double                      oob_breiman;

    MultiArray<2, double>       oobCount;
    MultiArray<2, double>       oobErrorCount;
	/** Per Tree OOB error calculated as in OOB_PerTreeError
	 * (Ulli's version)
	 */
    double                      oob_per_tree2;

	/**Column containing the development of the Ensemble
	 * error rate with increasing number of trees
	 */
    MultiArray<2, double>       breiman_per_tree;
	/** 4 dimensional array containing the development of confusion matrices 
	 * with number of trees - can be used to estimate ROC curves etc.
	 *
	 * oobroc_per_tree(ii,jj,kk,ll) 
	 * corresponds true label = ii 
	 * predicted label = jj
	 * confusion matrix after ll trees
	 *
	 * explaination of third index:
	 *
	 * Two class case:
	 * kk = 0 - (treeCount-1)
	 *         Threshold is on Probability for class 0  is kk/(treeCount-1);
	 * More classes:
	 * kk = 0. Threshold on probability set by argMax of the probability array.
	 */
    MultiArray<4, double>       oobroc_per_tree;
    
    CompleteOOBInfo() : VisitorBase(), oob_mean(0), oob_std(0), oob_per_tree2(0)  {}

	/** save to HDF5 file
	 */
#ifdef HasHDF5
    void save(std::string filen, std::string pathn)
    {
        if(*(pathn.end()-1) != '/')
            pathn += "/";
        const char* filename = filen.c_str();
        MultiArray<2, double> temp(Shp(1,1), 0.0); 
        writeHDF5(filename, (pathn + "oob_per_tree").c_str(), oob_per_tree);
        writeHDF5(filename, (pathn + "oobroc_per_tree").c_str(), oobroc_per_tree);
        writeHDF5(filename, (pathn + "breiman_per_tree").c_str(), breiman_per_tree);
        temp[0] = oob_mean;
        writeHDF5(filename, (pathn + "per_tree_error").c_str(), temp);
        temp[0] = oob_std;
        writeHDF5(filename, (pathn + "per_tree_error_std").c_str(), temp);
        temp[0] = oob_breiman;
        writeHDF5(filename, (pathn + "breiman_error").c_str(), temp);
        temp[0] = oob_per_tree2;
        writeHDF5(filename, (pathn + "ulli_error").c_str(), temp);
    }
#endif
    // negative value if sample was ib, number indicates how often.
    //  value >=0  if sample was oob, 0 means fail 1, corrrect

    template<class RF, class PR>
    void visit_at_beginning(RF & rf, PR & pr)
    {
        class_count = rf.class_count();
        if(class_count == 2)
            oobroc_per_tree.reshape(MultiArrayShape<4>::type(2,2,rf.tree_count(), rf.tree_count()));
        else
            oobroc_per_tree.reshape(MultiArrayShape<4>::type(rf.class_count(),rf.class_count(),1, rf.tree_count()));
        tmp_prob.reshape(Shp(1, class_count), 0); 
        prob_oob.reshape(Shp(rf.ext_param().row_count_,class_count), 0);
        is_weighted = rf.options().predict_weighted_;
        oob_per_tree.reshape(Shp(1, rf.tree_count()), 0);
        breiman_per_tree.reshape(Shp(1, rf.tree_count()), 0);
        //do the first time called.
        if(int(oobCount.size()) != rf.ext_param_.row_count_)
        {
            oobCount.reshape(Shp(rf.ext_param_.row_count_, 1), 0);
            oobErrorCount.reshape(Shp(rf.ext_param_.row_count_,1), 0);
        }
    }

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        // go through the samples
        int total_oob =0;
        int wrong_oob =0;
        for(int ll = 0; ll < rf.ext_param_.row_count_; ++ll)
        {
            // if the lth sample is oob...
            if(!sm.is_used()[ll])
            {
                // update number of trees in which current sample is oob
                ++oobCount[ll];

                // update number of oob samples in this tree.
                ++total_oob; 
                // get the predicted votes ---> tmp_prob;
                int pos =  rf.tree(index).getToLeaf(rowVector(pr.features(),ll));
                Node<e_ConstProbNode> node ( rf.tree(index).topology_, 
                                                    rf.tree(index).parameters_,
                                                    pos);
                tmp_prob.init(0); 
                for(int ii = 0; ii < class_count; ++ii)
                {
                    tmp_prob[ii] = node.prob_begin()[ii];
                }
                if(is_weighted)
                {
                    for(int ii = 0; ii < class_count; ++ii)
                        tmp_prob[ii] = tmp_prob[ii] * (*(node.prob_begin()-1));
                }
                rowVector(prob_oob, ll) += tmp_prob;
                int label = argMax(tmp_prob); 
                
                if(label != pr.response()(ll, 0))
                {
                    // update number of wrong oob samples in this tree.
                    ++wrong_oob;
                    // update number of trees in which current sample is wrong oob
                    ++oobErrorCount[ll];
                }
            }
        }
        int breimanstyle = 0;
        int totalOobCount = 0;
        for(int ll=0; ll < (int)rf.ext_param_.row_count_; ++ll)
        {
            if(oobCount[ll])
            {
                if(argMax(rowVector(prob_oob, ll)) != pr.response()(ll, 0))
                   ++breimanstyle;
                ++totalOobCount;
                if(oobroc_per_tree.shape(2) == 1)
                {
                    oobroc_per_tree(pr.response()(ll,0), argMax(rowVector(prob_oob, ll)),0 ,index)++;
                }
            }
        }
        if(oobroc_per_tree.shape(2) == 1)
            oobroc_per_tree.bindOuter(index)/=totalOobCount;
        if(oobroc_per_tree.shape(2) > 1)
        {
            MultiArrayView<3, double> current_roc 
                    = oobroc_per_tree.bindOuter(index);
            for(int gg = 0; gg < current_roc.shape(2); ++gg)
            {
                for(int ll=0; ll < (int)rf.ext_param_.row_count_; ++ll)
                {
                    if(oobCount[ll])
                    {
                        int pred = prob_oob(ll, 1) > (double(gg)/double(current_roc.shape(2)))?
                                        1 : 0; 
                        current_roc(pr.response()(ll, 0), pred, gg)+= 1; 
                    }
                }
                current_roc.bindOuter(gg)/= totalOobCount;
            }
        }
        breiman_per_tree[index] = double(breimanstyle)/double(totalOobCount);
        oob_per_tree[index] = double(wrong_oob)/double(total_oob);
        // go through the ib samples; 
    }

    /** Normalise variable importance after the number of trees is known.
     */
    template<class RF, class PR>
    void visit_at_end(RF & rf, PR & pr)
    {
        // ullis original metric and breiman style stuff
        oob_per_tree2 = 0; 
        int totalOobCount =0;
        int breimanstyle = 0;
        for(int ll=0; ll < (int)rf.ext_param_.row_count_; ++ll)
        {
            if(oobCount[ll])
            {
                if(argMax(rowVector(prob_oob, ll)) != pr.response()(ll, 0))
                   ++breimanstyle;
                oob_per_tree2 += double(oobErrorCount[ll]) / oobCount[ll];
                ++totalOobCount;
            }
        }
        oob_per_tree2 /= totalOobCount; 
        oob_breiman = double(breimanstyle)/totalOobCount; 
        // mean error of each tree
        MultiArrayView<2, double> mean(Shp(1,1), &oob_mean);
        MultiArrayView<2, double> stdDev(Shp(1,1), &oob_std);
        rowStatistics(oob_per_tree, mean, stdDev);
    }
};

/** calculate variable importance while learning.
 */
class VariableImportanceVisitor : public VisitorBase
{
    public:

    /** This Array has the same entries as the R - random forest variable
     *  importance.
	 *  Matrix is   featureCount by (classCount +2)
	 *  variable_importance_(ii,jj) is the variable importance measure of 
	 *  the ii-th variable according to:
	 *  jj = 0 - (classCount-1)
	 *  	classwise permutation importance 
	 *  jj = rowCount(variable_importance_) -2
	 *  	permutation importance
	 *  jj = rowCount(variable_importance_) -1
	 *  	gini decrease importance.
	 *
	 *  permutation importance:
	 *  The difference between the fraction of OOB samples classified correctly
	 *  before and after permuting (randomizing) the ii-th column is calculated.
	 *  The ii-th column is permuted rep_cnt times.
	 *
	 *  class wise permutation importance:
	 *  same as permutation importance. We only look at those OOB samples whose 
	 *  response corresponds to class jj.
	 *
	 *  gini decrease importance:
	 *  row ii corresponds to the sum of all gini decreases induced by variable ii 
	 *  in each node of the random forest.
     */
    MultiArray<2, double>       variable_importance_;
    int                         repetition_count_;
    bool                        in_place_;

#ifdef HasHDF5
    void save(std::string filename, std::string prefix)
    {
        prefix = "variable_importance_" + prefix;
        writeHDF5(filename.c_str(), 
                        prefix.c_str(), 
                        variable_importance_);
    }
#endif
	/** Constructor
	 * \param rep_cnt (defautl: 10) how often should 
	 * the permutation take place. Set to 1 to make calculation faster (but
	 * possibly more instable)
	 */
    VariableImportanceVisitor(int rep_cnt = 10) 
    :   repetition_count_(rep_cnt)

    {}

    /** calculates impurity decrease based variable importance after every
     * split.  
     */
    template<class Tree, class Split, class Region, class Feature_t, class Label_t>
    void visit_after_split( Tree          & tree, 
                            Split         & split,
                            Region        & parent,
                            Region        & leftChild,
                            Region        & rightChild,
                            Feature_t     & features,
                            Label_t       & labels)
    {
        //resize to right size when called the first time
        
        Int32 const  class_count = tree.ext_param_.class_count_;
        Int32 const  column_count = tree.ext_param_.column_count_;
        if(variable_importance_.size() == 0)
        {
            
            variable_importance_
                .reshape(MultiArrayShape<2>::type(column_count, 
                                                 class_count+2));
        }

        if(split.createNode().typeID() == i_ThresholdNode)
        {
            Node<i_ThresholdNode> node(split.createNode());
            variable_importance_(node.column(),class_count+1) 
                += split.region_gini_ - split.minGini();
        }
    }

    /**compute permutation based var imp. 
     * (Only an Array of size oob_sample_count x 1 is created.
     *  - apposed to oob_sample_count x feature_count in the other method.
     * 
     * \sa FieldProxy
     */
    template<class RF, class PR, class SM, class ST>
    void after_tree_ip_impl(RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        typedef MultiArrayShape<2>::type Shp_t;
        Int32                   column_count = rf.ext_param_.column_count_;
        Int32                   class_count  = rf.ext_param_.class_count_;  
        
		/* This solution saves memory uptake but not multithreading
		 * compatible
		 */
        // remove the const cast on the features (yep , I know what I am 
        // doing here.) data is not destroyed.
        //typename PR::Feature_t & features 
        //    = const_cast<typename PR::Feature_t &>(pr.features());

		typename PR::FeatureWithMemory_t features = pr.features();

        //find the oob indices of current tree. 
        ArrayVector<Int32>      oob_indices;
        ArrayVector<Int32>::iterator
                                iter;
        for(int ii = 0; ii < rf.ext_param_.row_count_; ++ii)
            if(!sm.is_used()[ii])
                oob_indices.push_back(ii);

        //create space to back up a column      
        std::vector<double>     backup_column;

        // Random foo
#ifdef CLASSIFIER_TEST
        RandomMT19937           random(1);
#else 
        RandomMT19937           random(RandomSeed);
#endif
        UniformIntRandomFunctor<RandomMT19937>  
                                randint(random);


        //make some space for the results
        MultiArray<2, double>
                    oob_right(Shp_t(1, class_count + 1)); 
        MultiArray<2, double>
                    perm_oob_right (Shp_t(1, class_count + 1)); 
            
        
        // get the oob success rate with the original samples
        for(iter = oob_indices.begin(); 
            iter != oob_indices.end(); 
            ++iter)
        {
            if(rf.tree(index)
                    .predictLabel(rowVector(features, *iter)) 
                ==  pr.response()(*iter, 0))
            {
                //per class
                ++oob_right[pr.response()(*iter,0)];
                //total
                ++oob_right[class_count];
            }
        }
        //get the oob rate after permuting the ii'th dimension.
        for(int ii = 0; ii < column_count; ++ii)
        {
            perm_oob_right.init(0.0); 
            //make backup of orinal column
            backup_column.clear();
            for(iter = oob_indices.begin(); 
                iter != oob_indices.end(); 
                ++iter)
            {
                backup_column.push_back(features(*iter,ii));
            }
            
            //get the oob rate after permuting the ii'th dimension.
            for(int rr = 0; rr < repetition_count_; ++rr)
            {               
                //permute dimension. 
                int n = oob_indices.size();
                for(int jj = 1; jj < n; ++jj)
                    std::swap(features(oob_indices[jj], ii), 
                              features(oob_indices[randint(jj+1)], ii));

                //get the oob sucess rate after permuting
                for(iter = oob_indices.begin(); 
                    iter != oob_indices.end(); 
                    ++iter)
                {
                    if(rf.tree(index)
                            .predictLabel(rowVector(features, *iter)) 
                        ==  pr.response()(*iter, 0))
                    {
                        //per class
                        ++perm_oob_right[pr.response()(*iter, 0)];
                        //total
                        ++perm_oob_right[class_count];
                    }
                }
            }
            
            
            //normalise and add to the variable_importance array.
            perm_oob_right  /=  repetition_count_;
            perm_oob_right -=oob_right;
            perm_oob_right *= -1;
            perm_oob_right      /=  oob_indices.size();
            variable_importance_
                .subarray(Shp_t(ii,0), 
                          Shp_t(ii+1,class_count+1)) += perm_oob_right;
            //copy back permuted dimension
            for(int jj = 0; jj < int(oob_indices.size()); ++jj)
                features(oob_indices[jj], ii) = backup_column[jj];
        }
    }

    /** calculate permutation based impurity after every tree has been 
     * learned  default behaviour is that this happens out of place.
     * If you have very big data sets and want to avoid copying of data 
     * set the in_place_ flag to true. 
     */
    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
            after_tree_ip_impl(rf, pr, sm, st, index);
    }

    /** Normalise variable importance after the number of trees is known.
     */
    template<class RF, class PR>
    void visit_at_end(RF & rf, PR & pr)
    {
        variable_importance_ /= rf.trees_.size();
    }
};

/** Verbose output
 */
class RandomForestProgressVisitor : public VisitorBase {
    public:
    RandomForestProgressVisitor() : VisitorBase() {}

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(RF& rf, PR & pr,  SM & sm, ST & st, int index){
        if(index != rf.options().tree_count_-1) {
            std::cout << "\r[" << std::setw(10) << (index+1)/static_cast<double>(rf.options().tree_count_)*100 << "%]"
                      << " (" << index+1 << " of " << rf.options().tree_count_ << ") done" << std::flush;
        }
        else {
            std::cout << "\r[" << std::setw(10) << 100.0 << "%]" << std::endl;
        }
    }
    
    template<class RF, class PR>
    void visit_at_end(RF const & rf, PR const & pr) {
        std::string a = TOCS;
        std::cout << "all " << rf.options().tree_count_ << " trees have been learned in " << a  << std::endl;
    }
    
    template<class RF, class PR>
    void visit_at_beginning(RF const & rf, PR const & pr) {
        TIC;
        std::cout << "growing random forest, which will have " << rf.options().tree_count_ << " trees" << std::endl;
    }
    
    private:
    USETICTOC;
};


/** Computes Correlation/Similarity Matrix of features while learning
 * random forest.
 */
class CorrelationVisitor : public VisitorBase
{
    public:
	/** gini_missc(ii, jj) describes how well variable jj can describe a partition
	 * created on variable ii(when variable ii was chosen)
	 */ 
    MultiArray<2, double>   gini_missc;
    MultiArray<2, int>      tmp_labels;
	/** additional noise features. 
	 */
    MultiArray<2, double>   noise;
    MultiArray<2, double>   noise_l;
	/** how well can a noise column describe a partition created on variable ii.
	 */
    MultiArray<2, double>   corr_noise;
    MultiArray<2, double>   corr_l;

	/** Similarity Matrix
	 * 
	 * (numberOfFeatures + 1) by (number Of Features + 1) Matrix
	 * gini_missc 
	 *  - row normalized by the number of times the column was chosen
	 *  - mean of corr_noise subtracted
	 *  - and symmetrised. 
	 *          
	 */
    MultiArray<2, double>   similarity;
	/** Distance Matrix 1-similarity
	 */
    MultiArray<2, double>   distance;
    ArrayVector<int>        tmp_cc;
	
	/** How often was variable ii chosen
	 */
    ArrayVector<int>        numChoices;
    typedef BestGiniOfColumn<GiniCriterion> ColumnDecisionFunctor;
    BestGiniOfColumn<GiniCriterion>         bgfunc;
    void save(std::string file, std::string prefix)
    {
		/*
        std::string tmp;
#define VAR_WRITE(NAME) \
        tmp = #NAME;\
        tmp += "_";\
        tmp += prefix;\
        vigra::writeToHDF5File(file.c_str(), tmp.c_str(), NAME);
        VAR_WRITE(gini_missc);
        VAR_WRITE(corr_noise);
        VAR_WRITE(distance);
        VAR_WRITE(similarity);
        vigra::writeToHDF5File(file.c_str(), "nChoices", MultiArrayView<2, int>(MultiArrayShape<2>::type(numChoices.size(),1), numChoices.data()));
#undef VAR_WRITE
*/
    }
    template<class RF, class PR>
    void visit_at_beginning(RF const & rf, PR  & pr)
    {
        typedef MultiArrayShape<2>::type Shp;
        int n = rf.ext_param_.column_count_;
        gini_missc.reshape(Shp(n +1,n+ 1));
        corr_noise.reshape(Shp(n + 1, 10));
        corr_l.reshape(Shp(n +1, 10));

        noise.reshape(Shp(pr.features().shape(0), 10));
        noise_l.reshape(Shp(pr.features().shape(0), 10));
        RandomMT19937 random(RandomSeed);
        for(int ii = 0; ii < noise.size(); ++ii)
        {
            noise[ii]   = random.uniform53();
            noise_l[ii] = random.uniform53()  > 0.5;
        }
        bgfunc = ColumnDecisionFunctor( rf.ext_param_);
        tmp_labels.reshape(pr.response().shape()); 
        tmp_cc.resize(2);
        numChoices.resize(n+1);
        // look at allaxes
    }
    template<class RF, class PR>
    void visit_at_end(RF const & rf, PR const & pr)
    {
        typedef MultiArrayShape<2>::type Shp;
        similarity.reshape(gini_missc.shape());
        similarity = gini_missc;;
        MultiArray<2, double> mean_noise(Shp(corr_noise.shape(0), 1));
        rowStatistics(corr_noise, mean_noise);
        mean_noise/= MultiArrayView<2, int>(mean_noise.shape(), numChoices.data());        
        int rC = similarity.shape(0);
        for(int jj = 0; jj < rC-1; ++jj)
        {
            rowVector(similarity, jj) /= numChoices[jj];
            rowVector(similarity, jj) -= mean_noise(jj, 0);
        }
        for(int jj = 0; jj < rC; ++jj)
        {
            similarity(rC -1, jj) /= numChoices[jj];
        }
        rowVector(similarity, rC -  1) -= mean_noise(rC-1, 0);
        similarity = abs(similarity);
        FindMinMax<double> minmax;
        inspectMultiArray(srcMultiArrayRange(similarity), minmax);
        
        for(int jj = 0; jj < rC; ++jj)
            similarity(jj, jj) = minmax.max;
        
        similarity.subarray(Shp(0,0), Shp(rC-1, rC-1)) 
            += similarity.subarray(Shp(0,0), Shp(rC-1, rC-1)).transpose();
        similarity.subarray(Shp(0,0), Shp(rC-1, rC-1))/= 2;  
        columnVector(similarity, rC-1) = rowVector(similarity, rC-1).transpose();
        for(int jj = 0; jj < rC; ++jj)
            similarity(jj, jj) = 0;
        
        FindMinMax<double> minmax2;
        inspectMultiArray(srcMultiArrayRange(similarity), minmax2);
        for(int jj = 0; jj < rC; ++jj)
            similarity(jj, jj) = minmax2.max;
        distance.reshape(gini_missc.shape(), minmax2.max);
        distance -= similarity; 
    }

    template<class Tree, class Split, class Region, class Feature_t, class Label_t>
    void visit_after_split( Tree          & tree, 
                            Split         & split,
                            Region        & parent,
                            Region        & leftChild,
                            Region        & rightChild,
                            Feature_t     & features,
                            Label_t       & labels)
    {
        if(split.createNode().typeID() == i_ThresholdNode)
        {
            double wgini;
            tmp_cc.init(0); 
            for(int ii = 0; ii < parent.size(); ++ii)
            {
                tmp_labels[parent[ii]] 
                    = (features(parent[ii], split.bestSplitColumn()) < split.bestSplitThreshold());
                ++tmp_cc[tmp_labels[parent[ii]]];
            }
            double region_gini = bgfunc.loss_of_region(tmp_labels, 
                                                       parent.begin(),
                                                       parent.end(),
                                                       tmp_cc);

            int n = split.bestSplitColumn(); 
            ++numChoices[n];
            ++(*(numChoices.end()-1));
            //this functor does all the work
            for(int k = 0; k < features.shape(1); ++k)
            {
                bgfunc(columnVector(features, k),
                       tmp_labels, 
                       parent.begin(), parent.end(), 
                       tmp_cc);
                wgini = (region_gini - bgfunc.min_gini_);
                gini_missc(n, k) 
                    += wgini;
            }
            for(int k = 0; k < 10; ++k)
            {
                bgfunc(columnVector(noise, k),
                       tmp_labels, 
                       parent.begin(), parent.end(), 
                       tmp_cc);
                wgini = (region_gini - bgfunc.min_gini_);
                corr_noise(n, k) 
                    += wgini;
            }
            
            for(int k = 0; k < 10; ++k)
            {
                bgfunc(columnVector(noise_l, k),
                       tmp_labels, 
                       parent.begin(), parent.end(), 
                       tmp_cc);
                wgini = (region_gini - bgfunc.min_gini_);
                corr_l(n, k) 
                    += wgini;
            }
            bgfunc(labels, tmp_labels, parent.begin(), parent.end(),tmp_cc);
            wgini = (region_gini - bgfunc.min_gini_);
            gini_missc(n, columnCount(gini_missc)-1) 
                += wgini;
            
            region_gini = split.region_gini_;
#if 1 
            Node<i_ThresholdNode> node(split.createNode());
            gini_missc(rowCount(gini_missc)-1, 
                                  node.column()) 
                 +=split.region_gini_ - split.minGini();
#endif
            for(int k = 0; k < 10; ++k)
            {
                split.bgfunc(columnVector(noise, k),
                             labels, 
                             parent.begin(), parent.end(), 
                             parent.classCounts());
                corr_noise(rowCount(gini_missc)-1, 
                           k) 
                     += wgini;
            }
#if 0
            for(int k = 0; k < tree.ext_param_.actual_mtry_; ++k)
            {
                wgini = region_gini - split.min_gini_[k];
                
                gini_missc(rowCount(gini_missc)-1, 
                                      split.splitColumns[k]) 
                     += wgini;
            }
            
            for(int k=tree.ext_param_.actual_mtry_; k<features.shape(1); ++k)
            {
                split.bgfunc(columnVector(features, split.splitColumns[k]),
                             labels, 
                             parent.begin(), parent.end(), 
                             parent.classCounts());
                wgini = region_gini - split.bgfunc.min_gini_;
                gini_missc(rowCount(gini_missc)-1, 
                                      split.splitColumns[k]) += wgini;
            }
#endif
            // remember to partition the data according to the best.
                gini_missc(rowCount(gini_missc)-1, 
                           columnCount(gini_missc)-1) 
                     += region_gini;
                SortSamplesByDimensions<Feature_t> 
                sorter(features, split.bestSplitColumn(), split.bestSplitThreshold());
            std::partition(parent.begin(), parent.end(), sorter);
        }
    }
};


} // namespace visitors
} // namespace rf
} // namespace vigra

//@}
#endif // RF_VISITORS_HXX
