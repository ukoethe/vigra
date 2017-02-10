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
#ifndef VIGRA_RF3_VISITORS_HXX
#define VIGRA_RF3_VISITORS_HXX

#include <vector>
#include <memory>
#include "../multi_array.hxx"
#include "../multi_shape.hxx"
#include <typeinfo>


namespace vigra
{
namespace rf3
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
    void visit_after_training(VISITORS &, RF &, const FEATURES &, const LABELS &)
    {}

    /**
     * @brief Do something before a tree has been learned.
     * 
     * @param weights the actual instance weights (after bootstrap sampling and class weights)
     */
    template <typename TREE, typename FEATURES, typename LABELS, typename WEIGHTS>
    void visit_before_tree(TREE &, FEATURES &, LABELS &, WEIGHTS &)
    {}

    /**
     * @brief Do something after a tree has been learned.
     */
    template <typename RF, typename FEATURES, typename LABELS, typename WEIGHTS>
    void visit_after_tree(RF &,
                          FEATURES &,
                          LABELS &,
                          WEIGHTS &)
    {}

    /**
     * @brief Do something after the split was made.
     */
    template <typename TREE,
              typename FEATURES,
              typename LABELS,
              typename WEIGHTS,
              typename SCORER,
              typename ITER>
    void visit_after_split(TREE &,
                           FEATURES &,
                           LABELS &,
                           WEIGHTS &,
                           SCORER &,
                           ITER,
                           ITER,
                           ITER)
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
    template <typename TREE, typename FEATURES, typename LABELS, typename WEIGHTS>
    void visit_before_tree(
            TREE & /*tree*/,
            FEATURES & /*features*/,
            LABELS & /*labels*/,
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
            throw std::runtime_error("OOBError::visit_before_tree(): The tree has no out-of-bags.");
    }

    /**
     * Compute the out-of-bag error.
     */
    template <typename VISITORS, typename RF, typename FEATURES, typename LABELS>
    void visit_after_training(
            VISITORS & visitors,
            RF & rf,
            const FEATURES & features,
            const LABELS & labels
    ){
        // Check the input sizes.
        vigra_precondition(rf.num_trees() > 0, "OOBError::visit_after_training(): Number of trees must be greater than zero after training.");
        vigra_precondition(visitors.size() == rf.num_trees(), "OOBError::visit_after_training(): Number of visitors must be equal to number of trees.");
        size_t const num_instances = features.shape()[0];
        auto const num_features = features.shape()[1];
        for (auto vptr : visitors)
            vigra_precondition(vptr->is_in_bag_.size() == num_instances, "OOBError::visit_after_training(): Some visitors have the wrong number of data points.");

        // Get a prediction for each data point using only the trees where it is out of bag.
        typedef typename std::remove_const<LABELS>::type Labels;
        Labels pred(Shape1(1));
        oob_err_ = 0.0;
        for (size_t i = 0; i < (size_t)num_instances; ++i)
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

    VariableImportance(size_t repetition_count = 10)
        :
        repetition_count_(repetition_count)
    {}

    /**
     * Resize the variable importance array and store in-bag / out-of-bag information.
     */
    template <typename TREE, typename FEATURES, typename LABELS, typename WEIGHTS>
    void visit_before_tree(
            TREE & tree,
            FEATURES & features,
            LABELS & /*labels*/,
            WEIGHTS & weights
    ){
        // Resize the variable importance array.
        // The shape differs from the shape of the actual output, since the single trees
        // only store the impurity decrease without the permutation importances.
        auto const num_features = features.shape()[1];
        variable_importance_.reshape(Shape2(num_features, tree.num_classes()+2), 0.0);

        // Save the in-bags.
        double const EPS = 1e-20;
        bool found = false;
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
            throw std::runtime_error("VariableImportance::visit_before_tree(): The tree has no out-of-bags.");
    }

    /**
     * Calculate the impurity decrease based variable importance after every split.
     */
    template <typename TREE,
              typename FEATURES,
              typename LABELS,
              typename WEIGHTS,
              typename SCORER,
              typename ITER>
    void visit_after_split(TREE & tree,
                           FEATURES & /*features*/,
                           LABELS & labels,
                           WEIGHTS & weights,
                           SCORER & scorer,
                           ITER begin,
                           ITER /*split*/,
                           ITER end)
    {
        // Update the impurity decrease.
        typename SCORER::Functor functor;
        auto const region_impurity = functor.region_score(labels, weights, begin, end);
        auto const split_impurity = scorer.best_score_;
        variable_importance_(scorer.best_dim_, tree.num_classes()+1) += region_impurity - split_impurity;
    }

    /**
     * Compute the permutation importance.
     */
    template <typename RF, typename FEATURES, typename LABELS, typename WEIGHTS>
    void visit_after_tree(RF & rf,
                          const FEATURES & features,
                          const LABELS & labels,
                          WEIGHTS & /*weights*/)
    {
        // Non-const types of features and labels.
        typedef typename std::remove_const<FEATURES>::type Features;
        typedef typename std::remove_const<LABELS>::type Labels;
        
        typedef typename Features::value_type FeatureType;

        auto const num_features = features.shape()[1];

        // For the permutation importance, the features must be permuted (obviously).
        // This cannot be done on the original feature matrix, since it would interfere
        // with other threads in concurrent training. Therefore, we have to make a copy.
        Features feats;
        Labels labs;
        copy_out_of_bags(features, labels, feats, labs);
        auto const num_oobs = feats.shape()[0];

        // Compute (standard and class-wise) out-of-bag success rate with the original sample.
        MultiArray<1, double> oob_right(Shape1(rf.num_classes()+1), 0.0);
        vigra::MultiArray<1,int> pred( (Shape1(num_oobs)) );
        rf.predict(feats, pred, 1);
        for (size_t i = 0; i < (size_t)labs.size(); ++i)
        {
            if (labs(i) == pred(i))
            {
                oob_right(labs(i)) += 1.0; // per class
                oob_right(rf.num_classes()) += 1.0; // total
            }
        }

        // Get out-of-bag success rate after permuting the j'th dimension.
        UniformIntRandomFunctor<MersenneTwister> randint;
        for (size_t j = 0; j < (size_t)num_features; ++j)
        {
            MultiArray<1, FeatureType> backup(( Shape1(num_oobs) ));
            backup = feats.template bind<1>(j);
            MultiArray<2, double> perm_oob_right(Shape2(1, rf.num_classes()+1), 0.0);

            for (size_t k = 0; k < repetition_count_; ++k)
            {
                // Permute the current dimension.
                for (int ii = num_oobs-1; ii >= 1; --ii)
                    std::swap(feats(ii, j), feats(randint(ii+1), j));

                // Get the out-of-bag success rate after permuting.
                rf.predict(feats, pred, 1);
                for (size_t i = 0; i < (size_t)labs.size(); ++i)
                {
                    if (labs(i) == pred(i))
                    {
                        perm_oob_right(0, labs(i)) += 1.0; // per class
                        perm_oob_right(0, rf.num_classes()) += 1.0; // total
                    }
                }
            }

            // Normalize and add to the importance matrix.
            perm_oob_right /= repetition_count_;
            perm_oob_right.bind<0>(0) -= oob_right;
            perm_oob_right *= -1;
            perm_oob_right /= num_oobs;
            variable_importance_.subarray(Shape2(j, 0), Shape2(j+1, rf.num_classes()+1)) += perm_oob_right;

            // Copy back the permuted dimension.
            feats.template bind<1>(j) = backup;
        }
    }

    /**
     * Accumulate the variable importances from the single trees.
     */
    template <typename VISITORS, typename RF, typename FEATURES, typename LABELS>
    void visit_after_training(
            VISITORS & visitors,
            RF & rf,
            const FEATURES & features,
            const LABELS & /*labels*/
    ){
        vigra_precondition(rf.num_trees() > 0, "VariableImportance::visit_after_training(): Number of trees must be greater than zero after training.");
        vigra_precondition(visitors.size() == rf.num_trees(), "VariableImportance::visit_after_training(): Number of visitors must be equal to number of trees.");

        // Sum the variable importances from the single trees.
        auto const num_features = features.shape()[1];
        variable_importance_.reshape(Shape2(num_features, rf.num_classes()+2), 0.0);
        for (auto vptr : visitors)
        {
            vigra_precondition(vptr->variable_importance_.shape() == variable_importance_.shape(),
                               "VariableImportance::visit_after_training(): Shape mismatch.");
            variable_importance_ += vptr->variable_importance_;
        }

        // Normalize the variable importance.
        variable_importance_ /= rf.num_trees();
    }

    /**
     * This Array has the same entries as the R - random forest variable
     * importance.
     * Matrix is   featureCount by (classCount +2)
     * variable_importance_(ii,jj) is the variable importance measure of 
     * the ii-th variable according to:
     * jj = 0 - (classCount-1)
     *     classwise permutation importance 
     * jj = rowCount(variable_importance_) -2
     *     permutation importance
     * jj = rowCount(variable_importance_) -1
     *     gini decrease importance.
     *     
     * permutation importance:
     * The difference between the fraction of OOB samples classified correctly
     * before and after permuting (randomizing) the ii-th column is calculated.
     * The ii-th column is permuted rep_cnt times.
     *
     * class wise permutation importance:
     * same as permutation importance. We only look at those OOB samples whose 
     * response corresponds to class jj.
     *
     * gini decrease importance:
     * row ii corresponds to the sum of all gini decreases induced by variable ii 
     * in each node of the random forest.
     */
    MultiArray<2, double> variable_importance_;

    /**
     * how often the permutation takes place
     */
    size_t repetition_count_;

private:

    /**
     * Copy the out-of-bag features and labels.
     */
    template <typename F0, typename L0, typename F1, typename L1>
    void copy_out_of_bags(
            F0 const & features_in,
            L0 const & labels_in,
            F1 & features_out,
            L1 & labels_out
    ) const {
        auto const num_instances = features_in.shape()[0];
        auto const num_features = features_in.shape()[1];

        // Count the out-of-bags.
        size_t num_oobs = 0;
        for (auto x : is_in_bag_)
            if (!x)
                ++num_oobs;

        // Copy the out-of-bags.
        features_out.reshape(Shape2(num_oobs, num_features));
        labels_out.reshape(Shape1(num_oobs));
        size_t current = 0;
        for (size_t i = 0; i < (size_t)num_instances; ++i)
        {
            if (!is_in_bag_[i])
            {
                auto const src = features_in.template bind<0>(i);
                auto out = features_out.template bind<0>(current);
                out = src;
                labels_out(current) = labels_in(i);
                ++current;
            }
        }
    }

    std::vector<bool> is_in_bag_; // whether a data point is in-bag or out-of-bag
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
    void visit_after_training(VISITORS & v, RF & rf, const FEATURES & features, const LABELS & labels)
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

    template <typename TREE, typename FEATURES, typename LABELS, typename WEIGHTS>
    void visit_before_tree(TREE & tree, FEATURES & features, LABELS & labels, WEIGHTS & weights)
    {
        if (visitor_.is_active())
            visitor_.visit_before_tree(tree, features, labels, weights);
        next_.visit_before_tree(tree, features, labels, weights);
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

    template <typename TREE,
              typename FEATURES,
              typename LABELS,
              typename WEIGHTS,
              typename SCORER,
              typename ITER>
    void visit_after_split(TREE & tree,
                           FEATURES & features,
                           LABELS & labels,
                           WEIGHTS & weights,
                           SCORER & scorer,
                           ITER begin,
                           ITER split,
                           ITER end)
    {
        if (visitor_.is_active())
            visitor_.visit_after_split(tree, features, labels, weights, scorer, begin, split, end);
        next_.visit_after_split(tree, features, labels, weights, scorer, begin, split, end);
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



} // namespace rf3
} // namespace vigra

#endif
