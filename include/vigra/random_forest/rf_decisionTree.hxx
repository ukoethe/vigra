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

#ifndef VIGRA_RANDOM_FOREST_DT_HXX
#define VIGRA_RANDOM_FOREST_DT_HXX

#include <algorithm>
#include <map>
#include <numeric>
#include "vigra/multi_array.hxx"
#include "vigra/mathutil.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/sized_int.hxx"
#include "vigra/matrix.hxx"
#include "vigra/random.hxx"
#include "vigra/functorexpression.hxx"
#include <vector>

#include "rf_common.hxx"
#include "rf_visitors.hxx"
#include "rf_nodeproxy.hxx"
namespace vigra
{

namespace detail
{
 // todo FINALLY DECIDE TO USE CAMEL CASE OR UNDERSCORES !!!!!!
/** decisiontree classifier. 
 *
 * This class is actually meant to be used in conjunction with the 
 * Random Forest Classifier 
 * - My suggestion would be to use the RandomForest classifier with 
 *   following parameters instead of directly using this 
 *   class (Preprocessing default values etc is handled in there):
 *
 * \code
 *      RandomForest decisionTree(RF_Traits::Options_t()
 *                                  .features_per_node(RF_ALL)
 *                                  .tree_count(1)            );
 * \endcode
 * 
 * \todo remove the classCount and featurecount from the topology
 *       array. Pass ext_param_ to the nodes!
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

    ProblemSpec<> ext_param_;
    unsigned int classCount_;


  public:
    /** \brief Create tree with parameters */
    template<class T>
    DecisionTree(ProblemSpec<T> ext_param)
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
     * \tparam  StackEntry_t The Stackentry containing Node/StackEntry_t 
     *          Information used during learing. Each Split functor has a 
     *          Stack entry associated with it (Split_t::StackEntry_t)
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
                    Split_t                              split,
                    Stop_t                               stop,
                    Visitor_t &                          visitor,
                    Random_t &                           randint);
    template <  class U, class C,
             class U2, class C2,
             class StackEntry_t,
             class Stop_t,
             class Split_t,
             class Visitor_t,
             class Random_t>
    void continueLearn(   MultiArrayView<2, U, C> const       & features,
                          MultiArrayView<2, U2, C2> const     & labels,
                          StackEntry_t const &                  stack_entry,
                          Split_t                               split,
                          Stop_t                                stop,
                          Visitor_t &                           visitor,
                          Random_t &                            randint,
                          //an index to which the last created exterior node will be moved (because it is not used anymore)
                          int                                   garbaged_child=-1);

    /** is a node a Leaf Node? */
    inline bool isLeafNode(TreeInt in) const
    {
        return (in & LeafNodeTag) == LeafNodeTag;
    }

    /** data driven traversal from root to leaf
     *
     * traverse through tree with data given in features. Use Visitors to 
     * collect statistics along the way. 
     */
    template<class U, class C, class Visitor_t>
    TreeInt getToLeaf(MultiArrayView<2, U, C> const & features, 
                      Visitor_t  & visitor) const
    {
        TreeInt index = 2;
        while(!isLeafNode(topology_[index]))
        {
            visitor.visit_internal_node(*this, index, topology_[index],features);
            switch(topology_[index])
            {
                case i_ThresholdNode:
                {
                    Node<i_ThresholdNode> 
                                node(topology_, parameters_, index);
                    index = node.next(features);
                    break;
                }
                case i_HyperplaneNode:
                {
                    Node<i_HyperplaneNode> 
                                node(topology_, parameters_, index);
                    index = node.next(features);
                    break;
                }
                case i_HypersphereNode:
                {
                    Node<i_HypersphereNode> 
                                node(topology_, parameters_, index);
                    index = node.next(features);
                    break;
                }
#if 0 
                // for quick prototyping! has to be implemented.
                case i_VirtualNode:
                {
                    Node<i_VirtualNode> 
                                node(topology_, parameters, index);
                    index = node.next(features);
                }
#endif
                default:
                    vigra_fail("DecisionTree::getToLeaf():"
                               "encountered unknown internal Node Type");
            }
        }
        visitor.visit_external_node(*this, index, topology_[index],features);
        return index;
    }
    /** traverse tree to get statistics
     *
     * Tree is traversed in order the Nodes are in memory (i.e. if no 
     * relearning//pruning scheme is utilized this will be pre order)
     */
    template<class Visitor_t>
    void traverse_mem_order(Visitor_t visitor) const
    {
        TreeInt index = 2;
        Int32 ii = 0;
        while(index < topology_.size())
        {
            if(isLeafNode(topology_[index]))
            {
                visitor
                    .visit_external_node(*this, index, topology_[index]);
            }
            else
            {
                visitor
                    ._internal_node(*this, index, topology_[index]);
            }
        }
    }

    template<class Visitor_t>
    void traverse_post_order(Visitor_t visitor,  TreeInt start = 2) const
    {
        typedef TinyVector<double, 2> Entry; 
        std::vector<Entry > stack;
        std::vector<double> result_stack;
        stack.push_back(Entry(2, 0));
        int addr; 
        while(!stack.empty())
        {
            addr = stack.back()[0];
            NodeBase node(topology_, parameters_, stack.back()[0]);
            if(stack.back()[1] == 1)
            {
                stack.pop_back();
                double leftRes = result_stack.back();
                double rightRes = result_stack.back();
                result_stack.pop_back();
                result_stack.pop_back();
                result_stack.push_back(rightRes+ leftRes);
                visitor.visit_internal_node(*this, 
                                            addr, 
                                            node.typeID(), 
                                            rightRes+leftRes);
            }
            else
            {
                if(isLeafNode(node.typeID()))
                {
                    visitor.visit_external_node(*this, 
                                                addr, 
                                                node.typeID(), 
                                                node.weights());
                    stack.pop_back();
                    result_stack.push_back(node.weights());
                }
                else
                {
                    stack.back()[1] = 1; 
                    stack.push_back(Entry(node.child(0), 0));
                    stack.push_back(Entry(node.child(1), 0));
                }
                    
            }
        }
    }

    /** same thing as above, without any visitors */
    template<class U, class C>
    TreeInt getToLeaf(MultiArrayView<2, U, C> const & features) const
    {
        rf::StopVisiting stop;
        return getToLeaf(features, stop);
    }


    template <class U, class C>
    ArrayVector<double>::iterator
    predict(MultiArrayView<2, U, C> const & features) const
    {
        TreeInt nodeindex = getToLeaf(features);
        switch(topology_[nodeindex])
        {
            case e_ConstProbNode:
                return Node<e_ConstProbNode>(topology_, 
                                             parameters_,
                                             nodeindex).prob_begin();
                break;
#if 0 
            //first make the Logistic regression stuff...
            case e_LogRegProbNode:
                return Node<e_LogRegProbNode>(topology_, 
                                              parameters_,
                                              nodeindex).prob_begin();
#endif            
            default:
                vigra_fail("DecisionTree::predict() :"
                           " encountered unknown external Node Type");
        }
        return ArrayVector<double>::iterator();
    }



    template <class U, class C>
    Int32 predictLabel(MultiArrayView<2, U, C> const & features) const
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
void DecisionTree::learn(   MultiArrayView<2, U, C> const       & features,
                            MultiArrayView<2, U2, C2> const     & labels,
                            StackEntry_t const &                  stack_entry,
                            Split_t                               split,
                            Stop_t                                stop,
                            Visitor_t &                           visitor,
                            Random_t &                            randint)
{
    this->reset();
    topology_.reserve(256);
    parameters_.reserve(256);
    topology_.push_back(features.shape(1));
    topology_.push_back(classCount_);
    continueLearn(features,labels,stack_entry,split,stop,visitor,randint);
}

template <  class U, class C,
            class U2, class C2,
            class StackEntry_t,
            class Stop_t,
            class Split_t,
            class Visitor_t,
            class Random_t>
void DecisionTree::continueLearn(   MultiArrayView<2, U, C> const       & features,
                            MultiArrayView<2, U2, C2> const     & labels,
                            StackEntry_t const &                  stack_entry,
                            Split_t                               split,
                            Stop_t                                stop,
                            Visitor_t &                           visitor,
                            Random_t &                            randint,
                            //an index to which the last created exterior node will be moved (because it is not used anymore)
                            int                                   garbaged_child)
{
    std::vector<StackEntry_t> stack;
    stack.reserve(128);
    ArrayVector<StackEntry_t> child_stack_entry(2, stack_entry);
    stack.push_back(stack_entry);
    size_t last_node_pos = 0;
    StackEntry_t top=stack.back();

    while(!stack.empty())
    {

        // Take an element of the stack. Obvious ain't it?
        top = stack.back();
        stack.pop_back();

        // Make sure no data from the last round has remained in Pipeline;
        child_stack_entry[0].reset();
        child_stack_entry[1].reset();
        split.reset();


        //Either the Stopping criterion decides that the split should 
        //produce a Terminal Node or the Split itself decides what 
        //kind of node to make
        TreeInt NodeID;
        
        if(stop(top))
            NodeID = split.makeTerminalNode(features, 
                                            labels, 
                                            top, 
                                            randint);
        else
        {
            //TIC;
            NodeID = split.findBestSplit(features, 
                                         labels, 
                                         top, 
                                         child_stack_entry, 
                                         randint);
            //std::cerr << TOC <<" " << NodeID << ";" <<std::endl;
        }

        // do some visiting yawn - just added this comment as eye candy
        // (looks odd otherwise with my syntax highlighting....
        visitor.visit_after_split(*this, split, top, 
                                  child_stack_entry[0], 
                                  child_stack_entry[1],
                                  features, 
                                  labels);


        // Update the Child entries of the parent
        // Using InteriorNodeBase because exact parameter form not needed.
        // look at the Node base before getting scared.
        last_node_pos = topology_.size();
        if(top.leftParent != StackEntry_t::DecisionTreeNoParent)
        {
            NodeBase(topology_, 
                     parameters_, 
                     top.leftParent).child(0) = last_node_pos;
        }
        else if(top.rightParent != StackEntry_t::DecisionTreeNoParent)
        {
            NodeBase(topology_, 
                     parameters_, 
                     top.rightParent).child(1) = last_node_pos;
        }


        // Supply the split functor with the Node type it requires.
        // set the address to which the children of this node should point 
        // to and push back children onto stack
        if(!isLeafNode(NodeID))
        {
            child_stack_entry[0].leftParent = topology_.size();
            child_stack_entry[1].rightParent = topology_.size();    
            child_stack_entry[0].rightParent = -1;
            child_stack_entry[1].leftParent = -1;
            stack.push_back(child_stack_entry[0]);
            stack.push_back(child_stack_entry[1]);
        }

        //copy the newly created node form the split functor to the
        //decision tree.
        NodeBase(split.createNode(), topology_, parameters_ );
    }
    if(garbaged_child!=-1)
    {
        Node<e_ConstProbNode>(topology_,parameters_,garbaged_child).copy(Node<e_ConstProbNode>(topology_,parameters_,last_node_pos));

        int last_parameter_size = Node<e_ConstProbNode>(topology_,parameters_,garbaged_child).parameters_size();
        topology_.resize(last_node_pos);
        parameters_.resize(parameters_.size() - last_parameter_size);
    
        if(top.leftParent != StackEntry_t::DecisionTreeNoParent)
            NodeBase(topology_, 
                     parameters_, 
                     top.leftParent).child(0) = garbaged_child;
        else if(top.rightParent != StackEntry_t::DecisionTreeNoParent)
            NodeBase(topology_, 
                     parameters_, 
                     top.rightParent).child(1) = garbaged_child;
    }
}

} //namespace detail

} //namespace vigra

#endif //VIGRA_RANDOM_FOREST_DT_HXX
