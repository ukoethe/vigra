/************************************************************************/
/*                                                                      */
/*        Copyright 2014-2015 by Ullrich Koethe and Philip Schill       */
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
#ifndef VIGRA_RANDOM_FOREST_VISITORS_HXX
#define VIGRA_RANDOM_FOREST_VISITORS_HXX

#include <vector>
#include <memory>


namespace vigra
{

/**
 * @brief Base class from which all random forest visitors derive.
 * 
 * @details
 * Due to the parallel training, we cannot simply use a single visitor for all trees.
 * Instead, each tree gets a copy of the original visitor.
 * 
 * The random forest training with visitors looks as follows:
 * - Do the random forest preprocessing (translate labels to 0, 1, 2, ...).
 * - Call visit_at_beginning() on the original visitor.
 * - For each tree:
 * - - Copy the original visitor and give the copy to the tree.
 * - - Do the preprocessing (create the bootstrap sample, assign weights to the data points, ...).
 * - - Call visit_before_tree() on the visitor copy.
 * - - Do the node splitting until the tree is fully trained.
 * - - Call visit_after_tree() on the visitor copy.
 * - Call visit_at_end (which gets a vector with pointers to the visitor copies) on the original visitor.
 */
class RFVisitorBase
{
public:

    RFVisitorBase()
        :
        active_(true)
    {}

    /**
     * @brief Do something before training starts.
     */
    void visit_before_training()
    {}

    /**
     * @brief Do something after all trees have been learned.
     * 
     * @param v vector with pointers to the visitor copies
     * @param rf the trained random forest
     */
    template <typename VISITORS, typename RF, typename FEATURES, typename LABELS>
    void visit_after_training(VISITORS & v, RF & rf, FEATURES & features, LABELS & labels)
    {}

    /**
     * @brief Do something before a tree has been learned.
     * 
     * @param weights the actual instance weights (after bootstrap sampling and class weights)
     */
    template <typename FEATURES, typename LABELS, typename WEIGHTS>
    void visit_before_tree(FEATURES & features, LABELS & labels, WEIGHTS & weights)
    {}

    /**
     * @brief Do something after a tree has been learned.
     */
    template <typename RF, typename FEATURES, typename LABELS, typename WEIGHTS>
    void visit_after_tree(RF & rf,
                          FEATURES & features,
                          LABELS & labels,
                          WEIGHTS & weights)
    {}

    /**
     * @brief Do something after the split was made.
     */
    template <typename FEATURES,
              typename LABELS,
              typename WEIGHTS,
              typename SCORER,
              typename ITER>
    void visit_after_split(FEATURES & features,
                           LABELS & labels,
                           WEIGHTS & weights,
                           SCORER & scorer,
                           ITER begin,
                           ITER split,
                           ITER end)
    {}

    /**
     * @brief Return whether the visitor is active or not.
     */
    bool is_active() const
    {
        return active_;
    }

    /**
     * @brief Activate the visitor.
     */
    void activate()
    {
        active_ = true;
    }

    /**
     * @brief Deactivate the visitor.
     */
    void deactivate()
    {
        active_ = false;
    }

private:

    bool active_;

};

/////////////////////////////////////////////
//          The concrete visitors          //
/////////////////////////////////////////////

/**
 * @brief Compute the out of bag error.
 * 
 * After training, each data point is put down those trees for which it is OOB.
 * Using bootstrap sampling, each data point is OOB for approximately 37% of
 * the trees.
 */
class OOBError : public RFVisitorBase
{
public:

    /**
     * Save whether a data point is in-bag (weight > 0) or out-of-bag (weight == 0).
     */
    template <typename FEATURES, typename LABELS, typename WEIGHTS>
    void visit_before_tree(
            FEATURES & features,
            LABELS & labels,
            WEIGHTS & weights
    ){
        double const EPS = 1e-20;
        bool found = false;

        // Save the in-bags.
        is_in_bag_.resize(weights.size(), true);
        for (size_t i = 0; i < weights.size(); ++i)
        {
            if (std::abs(weights[i]) < EPS)
            {
                is_in_bag_[i] = false;
                found = true;
            }
        }

        if (!found)
            throw std::runtime_error("OOBError::visit_after_tree(): The tree has no out-of-bags.");
    }

    /**
     * Compute the out-of-bag error.
     */
    template <typename VISITORS, typename RF, typename FEATURES, typename LABELS>
    void visit_after_training(
            VISITORS & visitors,
            RF & rf,
            FEATURES & features,
            LABELS & labels
    ){
        // Check the input sizes.
        vigra_precondition(rf.num_trees() > 0, "OOBError::visit_after_training(): Number of trees must be greater than zero after training.");
        vigra_precondition(visitors.size() == rf.num_trees(), "OOBError::visit_after_training(): Number of visitors must be equal to number of trees.");
        size_t const num_instances = features.shape()[0];//visitors[0]->is_in_bag_.size();
        size_t const num_features = features.shape()[1];
        for (auto vptr : visitors)
            vigra_precondition(vptr->is_in_bag_.size() == num_instances, "OOBError::visit_after_training(): Some visitors have the wrong number of data points.");

        // Get a prediction for each data point using only the trees where it is out of bag.
        typedef typename std::remove_const<LABELS>::type Labels;
        Labels pred(Shape1(1));
        oob_err_ = 0.0;
        for (size_t i = 0; i < num_instances; ++i)
        {
            // Get the indices of the trees where the data points is out of bag.
            std::vector<size_t> tree_indices;
            for (size_t k = 0; k < visitors.size(); ++k)
                if (!visitors[k]->is_in_bag_[i])
                    tree_indices.push_back(k);

            // Get the prediction using the above trees.
            auto const sub_features = features.subarray(Shape2(i, 0), Shape2(i+1, num_features));
            rf.predict(sub_features, pred, 1, tree_indices);
            if (pred(0) != labels(i))
                oob_err_ += 1.0;
        }
        oob_err_ /= num_instances;
    }

    /**
     * the out-of-bag error
     */
    double oob_err_;

private:
    std::vector<bool> is_in_bag_; // whether a data point is in-bag or out-of-bag
};



/**
 * @brief Compute the variable importance.
 */
class VariableImportance : public RFVisitorBase
{
public:

    template <typename FEATURES,
              typename LABELS,
              typename WEIGHTS,
              typename SCORER,
              typename ITER>
    void visit_after_split(FEATURES & features,
                           LABELS & labels,
                           WEIGHTS & weights,
                           SCORER & scorer,
                           ITER begin,
                           ITER split,
                           ITER end)
    {



    }

    /**
     * the variable importance
     */
    MultiArray<2, double> variable_importance_;

    /**
     * how often the permutation takes place
     */
    int repetition_count_;

};



/////////////////////////////////////////////
//            The visitor chain            //
/////////////////////////////////////////////

/**
 * @brief The default visitor node (= "do nothing").
 */
class RFStopVisiting : public RFVisitorBase
{};

namespace detail
{

/**
 * @brief Container elements of the statically linked visitor list. Use the create_visitor() functions to create visitors up to size 10.
 */
template <typename VISITOR, typename NEXT = RFStopVisiting, bool CPY = false>
class RFVisitorNode
{
public:

    typedef VISITOR Visitor;
    typedef NEXT Next;

    typename std::conditional<CPY, Visitor, Visitor &>::type visitor_;
    Next next_;
    
    RFVisitorNode(Visitor & visitor, Next next)
        :
        visitor_(visitor),
        next_(next)
    {}

    explicit RFVisitorNode(Visitor & visitor)
        :
        visitor_(visitor),
        next_(RFStopVisiting())
    {}

    explicit RFVisitorNode(RFVisitorNode<Visitor, Next, !CPY> & other)
        :
        visitor_(other.visitor_),
        next_(other.next_)
    {}

    explicit RFVisitorNode(RFVisitorNode<Visitor, Next, !CPY> const & other)
        :
        visitor_(other.visitor_),
        next_(other.next_)
    {}

    void visit_before_training()
    {
        if (visitor_.is_active())
            visitor_.visit_before_training();
        next_.visit_before_training();
    }

    template <typename VISITORS, typename RF, typename FEATURES, typename LABELS>
    void visit_after_training(VISITORS & v, RF & rf, FEATURES & features, LABELS & labels)
    {
        typedef typename VISITORS::value_type VisitorNodeType;
        typedef typename VisitorNodeType::Visitor VisitorType;
        typedef typename VisitorNodeType::Next NextType;
        
        // We want to call the visit_after_training function of the concrete visitor (e. g. OOBError).
        // Since v is a vector of visitor nodes (and not a vector of concrete visitors), we have to
        // extract the concrete visitors.
        // A vector cannot hold references, so we use pointers instead.
        if (visitor_.is_active())
        {
            std::vector<VisitorType*> visitors;
            for (auto & x : v)
                visitors.push_back(&x.visitor_);
            visitor_.visit_after_training(visitors, rf, features, labels);
        }

        // Remove the concrete visitors that we just visited.
        std::vector<NextType> nexts;
        for (auto & x : v)
            nexts.push_back(x.next_);

        // Call the next visitor node in the chain.
        next_.visit_after_training(nexts, rf, features, labels);
    }

    template <typename FEATURES, typename LABELS, typename WEIGHTS>
    void visit_before_tree(FEATURES & features, LABELS & labels, WEIGHTS & weights)
    {
        if (visitor_.is_active())
            visitor_.visit_before_tree(features, labels, weights);
        next_.visit_before_tree(features, labels, weights);
    }

    template <typename RF, typename FEATURES, typename LABELS, typename WEIGHTS>
    void visit_after_tree(RF & rf,
                          FEATURES & features,
                          LABELS & labels,
                          WEIGHTS & weights)
    {
        if (visitor_.is_active())
            visitor_.visit_after_tree(rf, features, labels, weights);
        next_.visit_after_tree(rf, features, labels, weights);
    }

    template <typename FEATURES,
              typename LABELS,
              typename WEIGHTS,
              typename SCORER,
              typename ITER>
    void visit_after_split(FEATURES & features,
                           LABELS & labels,
                           WEIGHTS & weights,
                           SCORER & scorer,
                           ITER begin,
                           ITER split,
                           ITER end)
    {
        if (visitor_.is_active())
            visitor_.visit_after_split(features, labels, weights, scorer, begin, split, end);
        next_.visit_after_split(features, labels, weights, scorer, begin, split, end);
    }

};

} // namespace detail

/**
 * The VisitorCopy can be used to set the copy argument of the given visitor chain to true.
 */
template <typename VISITOR>
struct VisitorCopy
{
    typedef detail::RFVisitorNode<typename VISITOR::Visitor, typename VisitorCopy<typename VISITOR::Next>::type, true> type;
};

template <>
struct VisitorCopy<RFStopVisiting>
{
    typedef RFStopVisiting type;
};



//////////////////////////////////////////////////////////
//   Visitor factory functions for up to 10 visitors.   //
//      FIXME: This should be a variadic template.      //
//////////////////////////////////////////////////////////

template<typename A>
detail::RFVisitorNode<A>
create_visitor(A & a)
{
    typedef detail::RFVisitorNode<A> _0_t;
    _0_t _0(a);
    return _0;
}

template<typename A, typename B>
detail::RFVisitorNode<A, detail::RFVisitorNode<B> >
create_visitor(A & a, B & b)
{
    typedef detail::RFVisitorNode<B> _1_t;
    _1_t _1(b);
    typedef detail::RFVisitorNode<A, _1_t> _0_t;
    _0_t _0(a, _1);
    return _0;
}

template<typename A, typename B, typename C>
detail::RFVisitorNode<A, detail::RFVisitorNode<B, detail::RFVisitorNode<C> > >
create_visitor(A & a, B & b, C & c)
{
    typedef detail::RFVisitorNode<C> _2_t;
    _2_t _2(c);
    typedef detail::RFVisitorNode<B, _2_t> _1_t;
    _1_t _1(b, _2);
    typedef detail::RFVisitorNode<A, _1_t> _0_t;
    _0_t _0(a, _1);
    return _0;
}

template<typename A, typename B, typename C, typename D>
detail::RFVisitorNode<A, detail::RFVisitorNode<B, detail::RFVisitorNode<C, 
    detail::RFVisitorNode<D> > > >
create_visitor(A & a, B & b, C & c, D & d)
{
    typedef detail::RFVisitorNode<D> _3_t;
    _3_t _3(d);
    typedef detail::RFVisitorNode<C, _3_t> _2_t;
    _2_t _2(c, _3);
    typedef detail::RFVisitorNode<B, _2_t> _1_t;
    _1_t _1(b, _2);
    typedef detail::RFVisitorNode<A, _1_t> _0_t;
    _0_t _0(a, _1);
    return _0;
}

template<typename A, typename B, typename C, typename D, typename E>
detail::RFVisitorNode<A, detail::RFVisitorNode<B, detail::RFVisitorNode<C, 
    detail::RFVisitorNode<D, detail::RFVisitorNode<E> > > > >
create_visitor(A & a, B & b, C & c, D & d, E & e)
{
    typedef detail::RFVisitorNode<E> _4_t;
    _4_t _4(e);
    typedef detail::RFVisitorNode<D, _4_t> _3_t;
    _3_t _3(d, _4);
    typedef detail::RFVisitorNode<C, _3_t> _2_t;
    _2_t _2(c, _3);
    typedef detail::RFVisitorNode<B, _2_t> _1_t;
    _1_t _1(b, _2);
    typedef detail::RFVisitorNode<A, _1_t> _0_t;
    _0_t _0(a, _1);
    return _0;
}

template<typename A, typename B, typename C, typename D, typename E,
         typename F>
detail::RFVisitorNode<A, detail::RFVisitorNode<B, detail::RFVisitorNode<C, 
    detail::RFVisitorNode<D, detail::RFVisitorNode<E, detail::RFVisitorNode<F> > > > > >
create_visitor(A & a, B & b, C & c, D & d, E & e, F & f)
{
    typedef detail::RFVisitorNode<F> _5_t;
    _5_t _5(f);
    typedef detail::RFVisitorNode<E, _5_t> _4_t;
    _4_t _4(e, _5);
    typedef detail::RFVisitorNode<D, _4_t> _3_t;
    _3_t _3(d, _4);
    typedef detail::RFVisitorNode<C, _3_t> _2_t;
    _2_t _2(c, _3);
    typedef detail::RFVisitorNode<B, _2_t> _1_t;
    _1_t _1(b, _2);
    typedef detail::RFVisitorNode<A, _1_t> _0_t;
    _0_t _0(a, _1);
    return _0;
}

template<typename A, typename B, typename C, typename D, typename E,
         typename F, typename G>
detail::RFVisitorNode<A, detail::RFVisitorNode<B, detail::RFVisitorNode<C, 
    detail::RFVisitorNode<D, detail::RFVisitorNode<E, detail::RFVisitorNode<F,
    detail::RFVisitorNode<G> > > > > > >
create_visitor(A & a, B & b, C & c, D & d, E & e, F & f, G & g)
{
    typedef detail::RFVisitorNode<G> _6_t;
    _6_t _6(g);
    typedef detail::RFVisitorNode<F, _6_t> _5_t;
    _5_t _5(f, _6);
    typedef detail::RFVisitorNode<E, _5_t> _4_t;
    _4_t _4(e, _5);
    typedef detail::RFVisitorNode<D, _4_t> _3_t;
    _3_t _3(d, _4);
    typedef detail::RFVisitorNode<C, _3_t> _2_t;
    _2_t _2(c, _3);
    typedef detail::RFVisitorNode<B, _2_t> _1_t;
    _1_t _1(b, _2);
    typedef detail::RFVisitorNode<A, _1_t> _0_t;
    _0_t _0(a, _1);
    return _0;
}

template<typename A, typename B, typename C, typename D, typename E,
         typename F, typename G, typename H>
detail::RFVisitorNode<A, detail::RFVisitorNode<B, detail::RFVisitorNode<C, 
    detail::RFVisitorNode<D, detail::RFVisitorNode<E, detail::RFVisitorNode<F,
    detail::RFVisitorNode<G, detail::RFVisitorNode<H> > > > > > > >
create_visitor(A & a, B & b, C & c, D & d, E & e, F & f, G & g, H & h)
{
    typedef detail::RFVisitorNode<H> _7_t;
    _7_t _7(h);
    typedef detail::RFVisitorNode<G, _7_t> _6_t;
    _6_t _6(g, _7);
    typedef detail::RFVisitorNode<F, _6_t> _5_t;
    _5_t _5(f, _6);
    typedef detail::RFVisitorNode<E, _5_t> _4_t;
    _4_t _4(e, _5);
    typedef detail::RFVisitorNode<D, _4_t> _3_t;
    _3_t _3(d, _4);
    typedef detail::RFVisitorNode<C, _3_t> _2_t;
    _2_t _2(c, _3);
    typedef detail::RFVisitorNode<B, _2_t> _1_t;
    _1_t _1(b, _2);
    typedef detail::RFVisitorNode<A, _1_t> _0_t;
    _0_t _0(a, _1);
    return _0;
}

template<typename A, typename B, typename C, typename D, typename E,
         typename F, typename G, typename H, typename I>
detail::RFVisitorNode<A, detail::RFVisitorNode<B, detail::RFVisitorNode<C, 
    detail::RFVisitorNode<D, detail::RFVisitorNode<E, detail::RFVisitorNode<F,
    detail::RFVisitorNode<G, detail::RFVisitorNode<H, detail::RFVisitorNode<I> > > > > > > > >
create_visitor(A & a, B & b, C & c, D & d, E & e, F & f, G & g, H & h, I & i)
{
    typedef detail::RFVisitorNode<I> _8_t;
    _8_t _8(i);
    typedef detail::RFVisitorNode<H, _8_t> _7_t;
    _7_t _7(h, _8);
    typedef detail::RFVisitorNode<G, _7_t> _6_t;
    _6_t _6(g, _7);
    typedef detail::RFVisitorNode<F, _6_t> _5_t;
    _5_t _5(f, _6);
    typedef detail::RFVisitorNode<E, _5_t> _4_t;
    _4_t _4(e, _5);
    typedef detail::RFVisitorNode<D, _4_t> _3_t;
    _3_t _3(d, _4);
    typedef detail::RFVisitorNode<C, _3_t> _2_t;
    _2_t _2(c, _3);
    typedef detail::RFVisitorNode<B, _2_t> _1_t;
    _1_t _1(b, _2);
    typedef detail::RFVisitorNode<A, _1_t> _0_t;
    _0_t _0(a, _1);
    return _0;
}

template<typename A, typename B, typename C, typename D, typename E,
         typename F, typename G, typename H, typename I, typename J>
detail::RFVisitorNode<A, detail::RFVisitorNode<B, detail::RFVisitorNode<C, 
    detail::RFVisitorNode<D, detail::RFVisitorNode<E, detail::RFVisitorNode<F,
    detail::RFVisitorNode<G, detail::RFVisitorNode<H, detail::RFVisitorNode<I,
    detail::RFVisitorNode<J> > > > > > > > > >
create_visitor(A & a, B & b, C & c, D & d, E & e, F & f, G & g, H & h, I & i,
               J & j)
{
    typedef detail::RFVisitorNode<J> _9_t;
    _9_t _9(j);
    typedef detail::RFVisitorNode<I, _9_t> _8_t;
    _8_t _8(i, _9);
    typedef detail::RFVisitorNode<H, _8_t> _7_t;
    _7_t _7(h, _8);
    typedef detail::RFVisitorNode<G, _7_t> _6_t;
    _6_t _6(g, _7);
    typedef detail::RFVisitorNode<F, _6_t> _5_t;
    _5_t _5(f, _6);
    typedef detail::RFVisitorNode<E, _5_t> _4_t;
    _4_t _4(e, _5);
    typedef detail::RFVisitorNode<D, _4_t> _3_t;
    _3_t _3(d, _4);
    typedef detail::RFVisitorNode<C, _3_t> _2_t;
    _2_t _2(c, _3);
    typedef detail::RFVisitorNode<B, _2_t> _1_t;
    _1_t _1(b, _2);
    typedef detail::RFVisitorNode<A, _1_t> _0_t;
    _0_t _0(a, _1);
    return _0;
}







} // namespace vigra

#endif
