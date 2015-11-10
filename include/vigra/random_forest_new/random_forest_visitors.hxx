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
 * - Call visit_at_end (which gets a vector with the visitor copies) on the original visitor.
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
     * @param v vector with the visitor copies
     * @param rf the trained random forest
     */
    template <typename VISITORS, typename RF>
    void visit_after_training(VISITORS & v, RF & rf)
    {}

    /**
     * @brief Do something before a tree has been learned.
     */
    void visit_before_tree()
    {}

    /**
     * @brief Do something after a tree has been learned.
     */
    void visit_after_tree()
    {}

    /**
     * @brief Do something after the split was made.
     */
    template <typename FEATURES,
              typename LABELS,
              typename WEIGHTS,
              typename SPLIT,
              typename RANGE>
    void visit_after_split(FEATURES & features,
                           LABELS & labels,
                           WEIGHTS & weights,
                           SPLIT & split,
                           RANGE & parent,
                           RANGE & left_child,
                           RANGE & right_child)
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
 * @brief The default visitor (= "do nothing").
 */
class RFStopVisiting : public RFVisitorBase
{

};



/**
 * @brief Compute the out of bag error.
 */
class OOBError : public RFVisitorBase
{
public:



};



/////////////////////////////////////////////
//            The visitor chain            //
/////////////////////////////////////////////

namespace detail
{

/**
 * @brief Container elements of the statically linked visitor list. Use the create_visitor() functions to create visitors up to size 10.
 */
template <typename Visitor, typename Next = RFStopVisiting>
class RFVisitorNode
{
public:

    Visitor & visitor_;
    Next next_; // This is no reference, because Next is either RFStopVisiting or RFVisitorNode.
    
    RFVisitorNode(Visitor & visitor, Next next)
        :
        visitor_(visitor),
        next_(next)
    {}

    RFVisitorNode(Visitor & visitor)
        :
        visitor_(visitor),
        next_(RFStopVisiting())
    {}

    void visit_before_training()
    {
        if (visitor_.is_active())
            visitor_.visit_before_training();
        next_.visit_before_training();
    }

    template <typename VISITORS, typename RF>
    void visit_after_training(VISITORS & v, RF & rf)
    {
        if (visitor_.is_active())
            visitor_.visit_after_training(v, rf);
        next_.visit_after_training(v, rf);
    }

    void visit_before_tree()
    {
        if (visitor_.is_active())
            visitor_.visit_before_tree();
        next_.visit_before_tree();
    }

    void visit_after_tree()
    {
        if (visitor_.is_active())
            visitor_.visit_after_tree();
        next_.visit_after_tree();
    }

};

} // namespace detail



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
