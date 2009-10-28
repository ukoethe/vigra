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

#include "vigra/rf_common.hxx"
#include "vigra/rf_nodeproxy.hxx"
namespace vigra
{


class DecisionTree
{
  public:
    typedef Int32 TreeInt;

    ArrayVector<TreeInt>  topology_;
    ArrayVector<double>   parameters_;

    RF_Traits::ProblemSpec_t ext_param_;
    unsigned int classCount_;


  public:

    DecisionTree(RF_Traits::ProblemSpec_t ext_param)
    :
        ext_param_(ext_param),
        classCount_(ext_param.class_count_)
    {}

    //clears all memory used.
    void reset(unsigned int classCount = 0)
    {
        if(classCount)
            classCount_ = classCount;
        topology_.clear();
        parameters_.clear();
    }


    // learn functions
    template <  class U, class C,
                class U2, class C2,
                class Region,
                class EarlyStopingPredicate,
                class SplitFunctor,
                class Visitor,
                class Random >
    void learn(     MultiArrayView<2, U, C> const      & features,
                    MultiArrayView<2, U2, C2> const    & labels,
                    Region const &                      region,
                    SplitFunctor                        split,
                    EarlyStopingPredicate               earlyStoppingCriterion,
                    Visitor &                           visitor,
                    Random &                            randint);

    // TODO isLeafNode() function does not really belong to DecisionTreeClass - rather
    // to nodeTags enum ... but defining a struct makes longer typing time
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



    //TODO Add the traversing visitors

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

    //todo get depth functor
    /*template <class U, class C>
    Int32 search_depth(MultiArrayView<2, U, C> const & features) const
    {
        GetDepthFunctor func;
        getToLeaf(features, func);
        return func.depth;
    }*/



    // TODO reimplement print and printStatistics with treeiterator.

};


template <  class U, class C,
            class U2, class C2,
            class Region,
            class EarlyStopingPredicate,
            class SplitFunctor,
            class Visitor,
            class Random>
void DecisionTree::learn(   MultiArrayView<2, U, C> const     & features,
                            MultiArrayView<2, U2, C2> const             & labels,
                            Region const &                      region,
                            SplitFunctor                        split,
                            EarlyStopingPredicate               earlyStoppingCriterion,
                            Visitor &                           visitor,
                            Random &                            randint)
{
    this->reset();
    topology_.reserve(256);
    parameters_.reserve(256);
    topology_.push_back(features.shape(1));
    topology_.push_back(classCount_);


    std::vector<Region> stack;
    stack.reserve(128);
    ArrayVector<Region> childStatistics(2, region);
    stack.push_back(region);


    while(!stack.empty())
    {

        // Take an element of the stack.
        Region top = stack.back();
        stack.pop_back();

        // Make sure no data from the last round has remained in Pipeline;
        childStatistics[0].reset();
        childStatistics[1].reset();
        split.reset();
        //Either the StoppingCriterion decides that the split should produce a Terminal Node
        //Or the Split itself decides what kind of node to make


        Int32 NodeID;
        if(earlyStoppingCriterion(top))
            NodeID = split.makeTerminalNode(features, labels, top, randint);
        else
            NodeID = split.findBestSplit(features, labels, top, childStatistics, randint);


        // TODO visitor
        visitor.visitAfterSplit(split, top, childStatistics[0], childStatistics[1]);

        // Update the Child entries of the parent
        // Using InteriorNodeBase because exact parameter form not needed.


        if(top.leftParent != Region::DecisionTreeNoParent)
            NodeBase(topology_, parameters_, top.leftParent).child(0) = topology_.size();
        else if(top.rightParent != Region::DecisionTreeNoParent)
            NodeBase(topology_, parameters_, top.rightParent).child(1) = topology_.size();

        // Supply the split functor with the Node type it requires.
        childStatistics[0].leftParent = topology_.size();
        childStatistics[1].rightParent = topology_.size();

        NodeBase(split.createNode(), topology_, parameters_ );
        if(!isLeafNode(NodeID))
        {
            stack.push_back(childStatistics[0]);
            stack.push_back(childStatistics[1]);
        }
    }
}


} //namespace vigra

#endif //VIGRA_RANDOM_FOREST_DT_HXX
