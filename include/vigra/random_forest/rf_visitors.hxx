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
#if 0
/** non intrusive implementation of the static visitor lists
 */
template <class Visitor, class Next = StopVisiting>
class VisitorNode
{
    public:
	
    StopVisiting	stop_;
    Next & 			next_;
	Visitor 		visitor_;	
    VisitorBase(Visitor & visitor, Next & next) 
	: 
		nexti_(next), visitor_(visitor)
    {}

    VisitorBase(Visitor &  visitor) 
	: 
		next(stop_), visitor_(visitor)
    {}

    template<class Tree, class Split, class Region>
    void visit_after_split( Tree 	      & tree, 
						   	Split         & split,
                            Region        & parent,
                            Region        & leftChild,
                            Region        & rightChild)
    {
		if(visitor_.is_active())
			visitor_.visit_after_split(tree, split, 
									   parent, leftChild, rightChild);
        next_.visit_after_split(tree, split, parent, leftChild, rightChild);
    }

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
		if(visitor_.is_active())
        	visitor_.visit_after_tree(rf, pr, sm, st, index);
        next_.visit_after_tree(rf, pr, sm, st, index);
    }

    template<class RF, class PR>
    void visit_at_end(RF & rf, PR & pr)
    {
		if(visitor_.is_active())
        	visitor_.visit_at_end(rf, pr);
        next_.visit_at_end(rf, pr);
    }
	
	template<class TR, class IntT, class TopT>
	void visit_external_node(TR & tr, IntT & index, TopT & node_t)
	{
		if(visitor_.is_active())
			visitor_.visit_external_node(tr, index, node_t);
		next_.visit_external_node(tr, index, node_t);
	}
	template<class TR, class IntT, class TopT>
	void visit_internal_node(TR & tr, IntT & index, TopT & node_t)
	{
		if(visitor_.is_active())
			visitor_.visit_internal_node(tr, index, node_t);
		next_.visit_internal_node(tr, index, node_t);
	}

    double return_val()
    {
		if(visitor_.is_active() && visitor_.has_value())
			return visitor_.return_val();
		return next_.return_val();
    }
};
#endif 
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



/** FieldProxy class 
 *
 * This helper class is there to be able to make non contiguous data 
 * appear as contiguous data for usage with stl algorithm.
 *
 * The idea is to make an array of pointers to the desired memory 
 * (e.g.: each pointer points to some arbitrary point in a multiArray)
 * 
 * Instead of having to use the *x syntax the assignment and the cast to 
 * the underlying data type and swap of FieldProxy were overloaded.
 *
 * so when random_shuffle permutes the elements in an Array of FieldProxys
 * not the FieldProxy objects but the values the are pointing to are 
 * permuted
 */
template<class T>
struct FieldProxy
{

    T* something;
    FieldProxy(T& in): something(&in){};

    FieldProxy(): something(0){};
    FieldProxy& operator= (const T rhs)
    {
        *something = rhs;
        return *this;
    }
    FieldProxy& operator= (const FieldProxy rhs)
    {
        *something = *(rhs.something);
        return *this;
    }

    operator T() const
    {
        return *something;
    }


};
} //namespace vigra

namespace std
{
/** really don't know why swap must be defined explicitly but fine
 *	has been tested and works. 
 */
template<class T>
void swap(vigra::FieldProxy<T> & rhs, vigra::FieldProxy<T> & lhs)
{
    T temp= *(lhs.something);
    *(lhs.something) = *(rhs.something);
    *(rhs.something) = temp;
}
/** calculate variable importance while learning.
 */

}// namespace std
namespace vigra
{
template <class Next = StopVisiting>
class VariableImportanceVisitor : public VisitorBase<Next>
{
    public:
	typedef VisitorBase<Next> 	Base_t;

	/** This Array has the same entries as the R - random forest variable
	 *  importance
	 */
	MultiArray<2, double> 		variable_importance_;
	int 						repetition_count_;
	bool						in_place_;
	/** construct with next visitor object add options or whatever. Just see
	 * to it that the base constructor is called - otherwise -kabumm
	 */
    VariableImportanceVisitor(Next & next_, 
							  int rep_cnt = 10, 
							  bool in_pl = false) 
	: 	Base_t(next_),
		repetition_count_(rep_cnt),
		in_place_(in_pl)
    {}

	/** construct as last visitor in the chain - use this constructor only
	 * if Next is of type StopVisiting
	 */
    VariableImportanceVisitor(int repetition_count = 10, 
							  bool  in_place = false) 
	: 	Base_t(),
		repetition_count_(repetition_count),
		in_place_(in_place)
    {}

	/** calculates impurity decrease based variable importance after every
	 * split.  
	 */
    template<class Tree, class Split, class Region>
    void visit_after_split( Tree 	      & tree, 
						   	Split         & split,
                            Region        & parent,
                            Region        & leftChild,
                            Region        & rightChild)
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
				+= split.giniDecrease;
		}
		Base_t::visit_after_split(tree, split, 
								  parent, leftChild, rightChild);


    }

	/** copies oob samples and calculates the permutation based var imp
	 * 	with the new arrays.
	 */
    template<class RF, class PR, class SM, class ST>
    void after_tree_oop_impl(RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
		typedef MultiArrayShape<2>::type Shp_t;
		Int32 					column_count = rf.ext_param_.column_count_;
	   	Int32 					class_count  = rf.ext_param_.class_count_;	


		//find the oob indices of current tree. 
		ArrayVector<Int32> 		oob_indices;
		for(int ii = 0; ii < rf.ext_param_.row_count_; ++ii)
			if(!sm.is_used()[ii])
				oob_indices.push_back(ii);
		//Copy the oob samples. This is done twice to be able to restore
		//a column after permutation.
		MultiArray<2, double> oob_samples(Shp_t(oob_indices.size(),
												column_count));
		MultiArray<2, double> oob_labels(Shp_t(oob_indices.size(),
											   column_count));

		for(int ii = 0; ii < oob_indices.size(); ++ii)
		{
			rowVector(oob_samples, ii) = rowVector(pr.features(),
												   oob_indices[ii]);
			oob_labels[ii] 			   = pr.response()(oob_indices[ii], 0); 
		}
		MultiArray<2, double> perm_oob_samples = oob_samples;	


		// More Initialising Foo
#ifdef CLASSIFIER_TEST
		RandomMT19937			random(1);
#endif
#ifndef CLASSIFIER_TEST
		RandomMT19937  			random(RandomSeed);
#endif
		UniformIntRandomFunctor<RandomMT19937> 	
								randint(random);


		//make some space for the results
		MultiArray<2, double>
					oob_right(Shp_t(1, class_count + 1)); 
		MultiArray<2, double>
					perm_oob_right (Shp_t(1, class_count + 1)); 
			
		
		// get the oob success rate with the original samples
		for(int jj = 0; jj < oob_indices.size(); ++jj)
		{
			if(rf.tree(index).predictLabel(rowVector(oob_samples, jj)) 
				== 	oob_labels[jj])
			{
				//per class
				++oob_right[oob_labels[jj]];
				//total
				++oob_right[class_count];
			}
		}
		
		//get the oob rate after permuting the ii'th dimension.
		for(int ii = 0; ii < column_count; ++ii)
		{
			perm_oob_right.init(0.0);
			// the whole thing is repeated repetition_count_ times to 
			// stabilize the result
			for(int rr = 0; rr < repetition_count_; ++rr)
			{
				//permute dimension. 
				std::random_shuffle(columnVector(perm_oob_samples, ii).data(), 
									columnVector(perm_oob_samples, ii).data()
													   	 + oob_indices.size(),
									randint);

				//get the oob sucess rate after permuting
				for(int jj = 0; jj < oob_indices.size(); ++jj)
				{
					if(rf.tree(index)
					   	 	.predictLabel(rowVector(perm_oob_samples, jj)) 
						== 	oob_labels[jj])
					{
						//per class
						++perm_oob_right[oob_labels[jj]];
						//total
						++perm_oob_right[class_count];
					}
				}
			}


			//normalise and add to the variable_importance array.
			perm_oob_right 	/= 	repetition_count_;
			perm_oob_right -=oob_right;
			perm_oob_right *= -1;
			perm_oob_right 		/= 	oob_indices.size();
			variable_importance_
				.subarray(Shp_t(ii,0), 
						  Shp_t(ii+1,class_count+1)) += perm_oob_right;
			//copy back permuted dimension
			columnVector(perm_oob_samples, ii) = columnVector(oob_samples, ii);
		}
    }


	/** same thing as above, using less additional system memory. 
	 * (Only an Array of size oob_sample_count x 1 is created.
	 *  - apposed to oob_sample_count x feature_count in the other method.
	 * 
	 * \sa FieldProxy
	 */
    template<class RF, class PR, class SM, class ST>
    void after_tree_ip_impl(RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
		typedef MultiArrayShape<2>::type Shp_t;
		Int32 					column_count = rf.ext_param_.column_count_;
	   	Int32 					class_count  = rf.ext_param_.class_count_;	
		
		// remove the const cast on the features (yep , I know what I am 
		// doing here.) data is not destroyed.
		typename PR::Feature_t & features 
			= const_cast<typename PR::Feature_t &>(pr.features());

		//find the oob indices of current tree. 
		ArrayVector<Int32> 		oob_indices;
		ArrayVector<Int32>::iterator
								iter;
		for(int ii = 0; ii < rf.ext_param_.row_count_; ++ii)
			if(!sm.is_used()[ii])
				oob_indices.push_back(ii);

		//create space to back up a column		
    	std::vector<double> 	backup_column;
    	std::vector<FieldProxy<typename PR::Feature_t::value_type> >    
								original_column;

		// Random foo
#ifdef CLASSIFIER_TEST
		RandomMT19937			random(1);
#endif
#ifndef CLASSIFIER_TEST
		RandomMT19937  			random(RandomSeed);
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
				== 	pr.response()(*iter, 0))
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
            original_column.clear();
			for(iter = oob_indices.begin(); 
				iter != oob_indices.end(); 
				++iter)
            {
                backup_column.push_back(features(*iter,ii));
                original_column.push_back(features(*iter,ii));
            }
			
			//get the oob rate after permuting the ii'th dimension.
			for(int rr = 0; rr < repetition_count_; ++rr)
			{				
				//permute dimension. 
				std::random_shuffle(original_column.begin(),
									original_column.end(),
									randint);

				//get the oob sucess rate after permuting
				for(iter = oob_indices.begin(); 
					iter != oob_indices.end(); 
					++iter)
				{
					if(rf.tree(index)
					   	 	.predictLabel(rowVector(features, *iter)) 
						== 	pr.response()(*iter, 0))
					{
						//per class
						++perm_oob_right[pr.response()(*iter, 0)];
						//total
						++perm_oob_right[class_count];
					}
				}
			}
			
			
			//normalise and add to the variable_importance array.
			perm_oob_right 	/= 	repetition_count_;
			perm_oob_right -=oob_right;
			perm_oob_right *= -1;
			perm_oob_right 		/= 	oob_indices.size();
			variable_importance_
				.subarray(Shp_t(ii,0), 
						  Shp_t(ii+1,class_count+1)) += perm_oob_right;
			//copy back permuted dimension
            std::copy(backup_column.begin(), 
					  backup_column.end(), 
					  original_column.begin());
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
		if(in_place_)
			after_tree_ip_impl(rf, pr, sm, st, index);
		else
			after_tree_oop_impl(rf, pr, sm, st, index);

		Base_t::visit_after_tree(rf, pr, sm, st, index);
    }

	/** Normalise variable importance after the number of trees is known.
	 */
    template<class RF, class PR>
    void visit_at_end(RF & rf, PR & pr)
    {
		variable_importance_ /= rf.trees_.size();
		Base_t::visit_at_end(rf, pr);
    }
};

} // namespace vigra
#endif // RF_VISITORS_HXX
