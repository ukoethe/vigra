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
#ifndef RF_VISITORS_HXX
#define RF_VISITORS_HXX

namespace vigra
{
template<class Tag>
class RandomForest;

/** Last Visitor that should be called to stop the recursion.
 */
class StopVisiting
{
    public:
    template<class Tree, class Split, class Region>
    void visit_after_split( Tree 	      & tree, 
						   	Split    	  & split,
                            Region        & parent,
                            Region        & leftChild,
                            Region        & rightChild)
    {}
    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(    RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {}
    template<class RF, class PR>
    void visit_at_end(RF & rf, PR & pr)
    {}
	template<class TR, class IntT, class TopT>
	void visit_external_node(TR & tr, IntT index, TopT node_t)
	{}
	template<class TR, class IntT, class TopT>
	void visit_internal_node(TR & tr, IntT index, TopT node_t)
	{}
    double return_val()
    {
		return -1.0;
	}
};


/** Base Class from which all Visitors derive
 */
template <class Next>
class VisitorBase
{
    public:
    StopVisiting	a;
    Next & 			next;

    VisitorBase(Next & next_) : next(next_)
    {}

    VisitorBase() : next(a)
    {}

    template<class Tree, class Split, class Region>
    void visit_after_split( Tree 	      & tree, 
						   	Split         & split,
                            Region        & parent,
                            Region        & leftChild,
                            Region        & rightChild)
    {
        next.visit_after_split(tree, split, parent, leftChild, rightChild);
    }

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        next.visit_after_tree(rf, pr, sm, st, index);
    }

    template<class RF, class PR>
    void visit_at_end(RF & rf, PR & pr)
    {
        next.visit_at_end(rf, pr);
    }
	
	template<class TR, class IntT, class TopT>
	void visit_external_node(TR & tr, IntT index, TopT node_t)
	{
		next.visit_external_node();
	}
	template<class TR, class IntT, class TopT>
	void visit_internal_node(TR & tr, IntT index, TopT node_t)
	{
		next.visit_internal_node();
	}

    double return_val()
    {
        return next.return_val();
    }
};




/** Template to make new visitors 
 *
 * if you don't need a function just delete it. - The base class will take 
 * care of the rest. - leave the constructors
 *
 */
#define Visitor_Template Your_Name_Here
template <class Next = StopVisiting>
class Visitor_Template : public VisitorBase<Next>
{
    public:
	typedef VisitorBase<Next> Base_t;

	/** construct with next visitor object add options or whatever. Just see
	 * to it that the base constructor is called - otherwise -kabumm
	 */
    Visitor_Template(Next & next_) : Base_t(next_)
    {}

	/** construct as last visitor in the chain - use this constructor only
	 * if Next is of type StopVisiting
	 */
    Visitor_Template() : Base_t()
    {}

	/** do something after the the Split has decided how to process the Region
	 * (Stack entry)
	 *
	 * \param tree 		reference to the tree that is currently being learned
	 * \param split 	reference to the split object
	 * \param parent 	current stack entry  which was used to decide the split
	 * \param leftChild	left stack entry that will be pushed
	 * \param rightChild
	 * 					right stack entry that will be pushed.
	 * \sa RF_Traits::StackEntry_t
	 */
    template<class Tree, class Split, class Region>
    void visit_after_split( Tree 	      & tree, 
						   	Split         & split,
                            Region        & parent,
                            Region        & leftChild,
                            Region        & rightChild)
    {
		Base_t::visit_after_split(tree, split, 
								  parent, leftChild, rightChild);
    }

	/** do something after each tree has been learned
	 *
	 * \param rf 		reference to the random forest object that called this
	 * 					visitor
	 * \param pr 		reference to the preprocessor that processed the input
	 * \param sm		reference to the sampler object
	 * \param st		reference to the first stack entry
	 * \param index 	index of current tree
	 */
    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
		Base_t::visit_after_tree(rf, pr, sm, st, index);
    }
	
	/** do something after all trees have been learned
	 *
	 * \param rf		reference to the random forest object that called this
	 * 					visitor
	 * \param pr		reference to the preprocessor that processed the input
	 */
    template<class RF, class PR>
    void visit_at_end(RF & rf, PR & pr)
    {
		Base_t::visit_at_end(rf, pr);
    }
	
	/** do some thing while traversing tree after it has been learned 
	 * 	(external nodes)
	 *
	 * \param tr 		reference to the tree object that called this visitor
	 * \param index		index in the topology_ array we currently are at
	 * \param node_tag	type of node we have (will be e_.... - )
	 * \sa 	NodeTags;
	 *
	 * you can create the node by using a switch on node_tag and using the 
	 * corresponding Node objects. Or - if you do not care about the type 
	 * use the Nodebase class.
	 */
	template<class TR, class IntT, class TopT>
	void visit_external_node(TR & tr, IntT index, TopT node_t)
	{
		Base_t::visit_external_node();
	}

	/** do something when visiting a internal node after it has been learned
	 *
	 * \sa visit_external_node
	 */
	template<class TR, class IntT, class TopT>
	void visit_internal_node(TR & tr, IntT index, TopT node_t)
	{
		Base_t::visit_internal_node();
	}

	/** return a double value.  The value of the first 
	 * visitor encountered that has a return value is returned with the
	 * RandomForest::learn() method - or -1.0 if no return value visitor
	 * existed. This functionality basically only exists so that the 
	 * OOB - visitor can return the oob error rate like in the old version 
	 * of the random forest.
	 */
    double return_val()
    {
		// if nothing to return just ignore this func
		// return do_some_foo();
		return Base_t::return_val();
    }
};
#undef Visitor_Template



//////////////////////////////////////////////////////////////////////////////
// Visitors of communal interest. Do not spam this file with stuff          //
// nobody wants.												 	        //
//////////////////////////////////////////////////////////////////////////////





/** Visitor that calculates the oob error of the random forest. 
 * this is the default visitor used. 
 *
 * To bored to comment each line of this class - trust me it works.
 */
template <class Next = StopVisiting>
class OOB_Visitor:public VisitorBase<Next>
{
public:
    typedef VisitorBase<Next> BT;
    double oobError;
    int totalOobCount;
    ArrayVector<int> oobCount,oobErrorCount;

    OOB_Visitor()
    : oobError(0.0),
      totalOobCount(0)
    {}
    OOB_Visitor(Next & next_)
    : BT(next_),
        oobError(0.0),
      totalOobCount(0)
    {}

	/** does the basic calculation per tree*/
    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(    RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
		//do the first time called.
        if(oobCount.size() != rf.ext_param_.row_count_)
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
                if(		rf.tree(index)
							.predictLabel(rowVector(pr.features(), l)) 
					!= 	pr.response()(l,0))
				{
                    ++oobErrorCount[l];
				}
            }

        }
        BT::next.visit_after_tree(rf, pr, sm, st, index);
    }

	/** Does the normalisation
	 */
    template<class RF, class PR>
    void visit_at_end(RF & rf, PR & pr)
    {
		// do some normalisation
        for(int l=0; l < (int)rf.ext_param_.row_count_; ++l)
        if(oobCount[l])
        {
            oobError += double(oobErrorCount[l]) / oobCount[l];
            ++totalOobCount;
        }
        BT::visit_at_end(rf, pr);
    }
	
	//returns value of the learn function. 
    double return_val()
    {
        return oobError/totalOobCount;
    }
};


} // namespace vigra
#endif // RF_VISITORS_HXX
