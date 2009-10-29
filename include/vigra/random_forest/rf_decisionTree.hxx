/************************************************************************/
/*                                                                      */
/*        Copyright 2008-2009 by  Ullrich Koethe and Rahul Nair         */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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

#ifndef VIGRA_RANDOM_FOREST_DT_HXX
#define VIGRA_RANDOM_FOREST_DT_HXX

#include <algorithm>
#include <map>
#include <numeric>
#include <iostream>
#include "vigra/multi_array.hxx"
#include "vigra/mathutil.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/sized_int.hxx"
#include "vigra/matrix.hxx"
#include "vigra/random.hxx"
#include "vigra/functorexpression.hxx"
#include <vector>

#include "rf_common.hxx"
#include "rf_nodeproxy.hxx"
namespace vigra
{

namespace detail
{
/** decisiontree classifier. 
 * todo FINALLY DECIDE TO USE CAMEL CASE OR UNDERSCORES !!!!!!
 *
 * This class is actually meant to be used in conjunction with the 
 * Random Forest Classifier 
 * - My suggestion would be to use the RandomForest classifier with 
 *   following parameters instead of directly using this 
 *   class (Preprocessing default values etc is handled in there):
 *
 * \code
 * 		RandomForest decisionTree(RF_Traits::Options_t()
 * 									.features_per_node(RF_ALL)
 * 									.tree_count(1)			  );
 * \endcode
 * 
 * \todo remove the classCount and featurecount from the topology
 * 		 array. Pass ext_param_ to the nodes!
 * \todo Use relative addressing of nodes?
 */
class DecisionTree
{
	/**\todo make private?*/
  public:
	
	/** value type of container array. use whenever referencing it
	 */
    typedef Int32 TreeInt;

    ArrayVector<TreeInt>  topology_;
    ArrayVector<double>   parameters_;

    RF_Traits::ProblemSpec_t ext_param_;
    unsigned int classCount_;


  public:

	/** \Brief Create tree with parameters */
    DecisionTree(RF_Traits::ProblemSpec_t ext_param)
    :
        ext_param_(ext_param),
        classCount_(ext_param.class_count_)
    {}

    /**clears all memory used.
	 */
    void reset(unsigned int classCount = 0)
    {
        if(classCount)
            classCount_ = classCount;
        topology_.clear();
        parameters_.clear();
    }


    /** learn a Tree
	 *
	 * \tparam 	StackEntry_t The Stackentry containing Node/StackEntry_t 
	 * 			Information used during learing. Each Split functor has a 
	 * 			Stack entry associated with it (Split_t::StackEntry_t)
	 * \sa RandomForest::learn()
	 */
    template <  class U, class C,
                class U2, class C2,
                class StackEntry_t,
                class Stop_t,
                class Split_t,
                class Visitor_t,
                class Random_t >
    void learn(     MultiArrayView<2, U, C> const      & features,
                    MultiArrayView<2, U2, C2> const    & labels,
                    StackEntry_t const &                 stack_entry,
                    Split_t                        		 split,
                    Stop_t              				 stop,
                    Visitor_t &                          visitor,
                    Random_t &                           randint);

    /** is a node a Leaf Node? */
    inline bool isLeafNode(TreeInt in)
    {
        return (in & LeafNodeTag) == LeafNodeTag;
    }

    // TODO if DoNothingClass is not wegoptimiert then two version of this function have to be defined.
    // goToLeaf is used for predict, leafID and searchDepth functions
    template<class U, class C>
    TreeInt getToLeaf(MultiArrayView<2, U, C> const & features)
	{
		TreeInt index = 2;
        while(!isLeafNode(topology_[index]))
		{

            switch(topology_[index])
            {
                case i_ThresholdNode:
                {
                    Node<i_ThresholdNode> node(topology_, parameters_, index);
                    index = node.next(features);
                    break;
                }
                case i_HyperplaneNode:
                {
                    Node<i_HyperplaneNode> node(topology_, parameters_, index);
                    index = node.next(features);
                    break;
                }
                case i_HypersphereNode:
                {
                    Node<i_HypersphereNode> node(topology_, parameters_, index);
                    index = node.next(features);
                    break;
                }
                default:
                    vigra_fail("DecisionTree::getToLeaf() : encountered unknown Node Type");
            }
        }
        return index;
    }

    //template< class DoAtNodeFunctor, class DoAtLeafFunctor>
    //void TraverseTree(DoAtNodeFunctor doAtNode, DoAtLeafFunctor doAtLeaf);
    //TODO tree traversing iterator <- do we really need this?



    

    template <class U, class C>
    ArrayVector<double>::iterator
    predict(MultiArrayView<2, U, C> const & features)
    {
        TreeInt nodeindex = getToLeaf(features);
        switch(topology_[nodeindex])
        {
            case e_ConstProbNode:
            {
//                Node<e_ConstProbNode> da(topology_, parameters_,nodeindex);
////                std::cout << ">>>>>>>>>>>>>>>>>>>>PROB: " << Node<e_ConstProbNode>(topology_, parameters_,nodeindex).prob_begin()[0] << " " << Node<e_ConstProbNode>(topology_, parameters_,nodeindex).prob_begin()[1] << std::endl;
//                for(int ii = 0; ii < da.parameters_size(); ++ii)
//                {
//                    std::cout << da.parameters_begin()[ii];
//                }
                return Node<e_ConstProbNode>(topology_, parameters_,nodeindex).prob_begin();
            }
            	break;

            case e_LogRegProbNode:
                // TODO: uncomment this once class below is complete
                //return ExteriorNode<e_LogRegProbNode>(topology_, parameters_,nodeindex).getProbabilities();
			default:
                vigra_fail("DecisionTree::predict() : encountered unknown Node Type");
        }
        return ArrayVector<double>::iterator();
    }



    template <class U, class C>
    Int32 predictLabel(MultiArrayView<2, U, C> const & features)
    {
        ArrayVector<double>::const_iterator weights = predict(features);
        return argMax(weights, weights+classCount_) - weights;
    }


};


template <  class U, class C,
            class U2, class C2,
            class StackEntry_t,
            class Stop_t,
            class Split_t,
            class Visitor_t,
            class Random_t>
void DecisionTree::learn(   MultiArrayView<2, U, C> const     	& features,
                            MultiArrayView<2, U2, C2> const   	& labels,
                            StackEntry_t const &                  stack_entry,
                            Split_t                        		  split,
                            Stop_t               				  stop,
                            Visitor_t &                           visitor,
                            Random_t &                            randint)
{
    this->reset();
    topology_.reserve(256);
    parameters_.reserve(256);
    topology_.push_back(features.shape(1));
    topology_.push_back(classCount_);


    std::vector<StackEntry_t> stack;
    stack.reserve(128);
    ArrayVector<StackEntry_t> child_stack_entry(2, stack_entry);
    stack.push_back(stack_entry);


    while(!stack.empty())
    {

        // Take an element of the stack.
        StackEntry_t top = stack.back();
        stack.pop_back();

        // Make sure no data from the last round has remained in Pipeline;
        child_stack_entry[0].reset();
        child_stack_entry[1].reset();
        split.reset();



        //Either the StoppingCriterion decides that the split should 
		//produce a Terminal Node or the Split itself decides what 
		//kind of node to make
        TreeInt NodeID;
        if(stop(top))
            NodeID = split.makeTerminalNode(features, 
											labels, 
											top, 
											randint);
        else
            NodeID = split.findBestSplit(features, 
										 labels, 
										 top, 
										 child_stack_entry, 
										 randint);

		// do some visiting yawn - just added this comment as eye candy
		// (looks odd otherwise with my syntax highlighting....
        visitor.visit_after_split(*this, split, top, 
								  child_stack_entry[0], 
								  child_stack_entry[1]);

        // Update the Child entries of the parent
        // Using InteriorNodeBase because exact parameter form not needed.
		// look at the Node base before getting scared.
        if(top.leftParent != StackEntry_t::DecisionTreeNoParent)
            NodeBase(topology_, 
					 parameters_, 
					 top.leftParent).child(0) = topology_.size();
        else if(top.rightParent != StackEntry_t::DecisionTreeNoParent)
            NodeBase(topology_, 
					 parameters_, 
					 top.rightParent).child(1) = topology_.size();

        // Supply the split functor with the Node type it requires.
		// set the address to which the children of this node should point 
		// to and push back children onto stack
       	child_stack_entry[0].leftParent = topology_.size();
     	child_stack_entry[1].rightParent = topology_.size();    
		
        NodeBase(split.createNode(), topology_, parameters_ );

		if(!isLeafNode(NodeID))
        {

			stack.push_back(child_stack_entry[0]);
            stack.push_back(child_stack_entry[1]);
        }

    }
}

} //namespace detail

} //namespace vigra

#endif //VIGRA_RANDOM_FOREST_DT_HXX
