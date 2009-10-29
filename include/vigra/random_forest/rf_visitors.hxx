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

#include <fstream>
#include <sstream>
namespace vigra
{
template<class Tag>
class RandomForest;


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
	{

	}
	template<class TR, class IntT, class TopT>
	void visit_internal_node(TR & tr, IntT index, TopT node_t)
	{

	}
    double return_val()
    {
        return 0.0;
    }
};

template <class Next>
class VisitorBase
{
    public:
    StopVisiting a;
    Next & next;
    VisitorBase(Next & next_)
    : next(next_)
    {}

    VisitorBase()
    : next(a)
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
    void visit_after_tree(    RF& rf, PR & pr,  SM & sm, ST & st, int index)
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



template <class Next = StopVisiting>
class TestVisitor: public VisitorBase<Next>
{
    typedef VisitorBase<Next> BT;



    public:
    std::ofstream fout;
    TestVisitor(Next & next_, std::string output = std::string("RandomForestNodeTest.log"))
    :
        BT(next_),
        fout(output.c_str())
    {    }

    //TODO split must be const
    template<class Tree, class Split, class Region>
    void visit_after_split( Tree 	      & tree, 
						   	Split         & split,
                            Region        & parent,
                            Region        & leftChild,
                            Region        & rightChild)
    {
        if(split.createNode().typeID() == i_ThresholdNode)
        {
            std::ostringstream s1;
            s1.precision(10);
            s1.setf(std::ios::fixed,std::ios::floatfield);
            s1 << split.minGini;
            fout   << "minGini: " << s1.str().substr(0, s1.str().size()-4) << " - " << split.bestSplitColumn << std::endl;
//            fout   << "Threshold: " << *(split.node_.parameters_begin() +1) << std::endl;
            fout   << "Region.size: "      << parent.size() << std::endl;
            fout   << "LeftChild.size: "  << leftChild.size() << std::endl;
            fout   << "LeftChild.ClassCounts: ";
            for(int ii = 0; ii < split.classCount(); ++ii)
                fout << leftChild.classCounts()[ii] << " " ;
            fout   << std::endl;
            fout   << "RightChild.size: "   << rightChild.size() << std::endl;
            fout   << "RightChild.ClassCounts: ";
            for(int ii = 0; ii < split.classCount(); ++ii)
                fout << rightChild.classCounts()[ii] << " " ;
            fout   << std::endl;
            fout   << std::endl;
        }
        BT::visit_after_split(tree, split, parent, leftChild, rightChild);
    }

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(    RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        fout << std::endl << std::endl << "Tree Number: " << index << " finished." << std::endl << std::endl;
		BT::visit_after_tree(rf, pr, sm, st, index);
    }

    ~TestVisitor()
    {
        fout.close();
    }
};

template <class Next = StopVisiting>
class SetTestVisitor: public VisitorBase<Next>
{
    typedef VisitorBase<Next> BT;



    public:
    std::ostringstream sout;
    std::set<std::string> treesset;

    SetTestVisitor(Next & next_)
    :
        BT(next_)
    {    }

    template<class Tree, class Split, class Region>
    void visit_after_split( Tree 	      & tree, 
						   	Split         & split,
                            Region        & parent,
                            Region        & leftChild,
                            Region        & rightChild)
    {
        if(split.createNode().typeID() == i_ThresholdNode)
        {
            sout   << "minGini: " << split.minGini << " - " << split.bestSplitColumn << std::endl;
            sout   << "Region.size: "      << parent.size() << std::endl;
            sout   << "LeftChild.size: "  << leftChild.size() << std::endl;
            sout   << "LeftChild.ClassCounts: ";
            for(int ii = 0; ii < split.classCount(); ++ii)
                sout << leftChild.classCounts()[ii] << " " ;
            sout   << std::endl;
            sout   << "RightChild.size: "   << rightChild.size() << std::endl;
            sout   << "RightChild.ClassCounts: ";
            for(int ii = 0; ii < split.classCount(); ++ii)
                sout << rightChild.classCounts()[ii] << " " ;
            sout   << std::endl;
            sout   << std::endl;
        }
        BT::visit_after_split(tree, split, parent, leftChild, rightChild);
    }

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(    RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        treesset.insert(sout.str());
        sout.str(std::string());
		BT::visit_after_tree(rf, pr, sm, st, index);
    }

    template<class RF, class PR>
    void visit_at_end(RF & rf, PR & pr)
    {
        std::ofstream fout("setTest.log");
        std::set<std::string>::iterator iter;
        int k = 0;
        for(iter = treesset.begin(); iter != treesset.end(); ++iter)
        {
            fout << *iter;
            fout << std::endl << std::endl << "Tree Number: " << k << " finished." << std::endl << std::endl;
            ++k;
        }
        fout.close();
        BT::visit_at_end(rf, pr);
	}
};

template <class T1, class C1, class T2, class C2, class Next = StopVisiting>
class AllOutputVisitor: public VisitorBase<Next>
{
    typedef VisitorBase<Next> BT;
    MultiArrayView<2, T1, C1> * features;
    MultiArrayView<2, T2, C2> * labels;


    public:
    std::ofstream fout;
    AllOutputVisitor(Next & next_, std::string output = std::string("RandomForestNodeTest.log"))
    :
        BT(next_),
        fout(output.c_str())
    {    }

    //TODO split must be const
    template<class Tree, class Split, class Region>
    void visit_after_split( Tree 	      & tree, 
						   	Split         & split,
                            Region        & parent,
                            Region        & leftChild,
                            Region        & rightChild)
    {
        if(split.createNode().typeID() == i_ThresholdNode)
        {
            fout   << "minGini: " << split.minGini << " - " << split.bestSplitColumn << std::endl;
            fout   << "Threshold: " << *(split.node_.parameters_begin() +1) << std::endl;
            fout   << "Region.size: "      << parent.size() << std::endl;
            for(int ii = 0 ; ii < parent.size(); ++ii)
            {
                fout << "(" << parent[ii] << ", " << (*labels)[parent[ii]] << ", " << (*features)(parent[ii],split.bestSplitColumn) << ") ";
                if(ii%4 == 0)fout << std::endl;
            }
            fout << std::endl;
            fout   << "LeftChild.size: "  << leftChild.size() << std::endl;
            fout   << "LeftChild.ClassCounts: ";
            for(int ii = 0; ii < split.classCount(); ++ii)
                fout << leftChild.classCounts()[ii] << " " ;
            fout   << std::endl;
            for(int ii = 0 ; ii < leftChild.size(); ++ii)
            {
                fout << "(" << leftChild[ii] << ", " << (*labels)[leftChild[ii]] << ", " << (*features)(leftChild[ii],split.bestSplitColumn) << ") ";
                if(ii%4 == 0)fout << std::endl;
            }
            fout << std::endl;
            fout   << "RightChild.size: "   << rightChild.size() << std::endl;
            fout   << "RightChild.ClassCounts: ";
            for(int ii = 0; ii < split.classCount(); ++ii)
                fout << rightChild.classCounts()[ii] << " " ;
            fout   << std::endl;
            for(int ii = 0 ; ii < rightChild.size(); ++ii)
            {
            fout << "(" << rightChild[ii] << ", " << (*labels)[rightChild[ii]] << ", " << (*features)(rightChild[ii],split.bestSplitColumn) << ") ";
                if(ii%4 == 0)fout << std::endl;
            }
            fout   << std::endl;
        }
        BT::visit_after_split(tree, split, parent, leftChild, rightChild);
    }

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(    RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        fout << std::endl << std::endl << "Tree Number: " << index << " finished." << std::endl << std::endl;
		BT::visit_after_tree(rf, pr, sm, st, index);
	}

    void setDataSource(MultiArrayView<2, T1, C1> * feat, MultiArrayView<2, T2, C2> * label)
    {
        features = feat;
        labels = label;
    }


    ~AllOutputVisitor()
    {
        fout.close();
    }
};

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


    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(    RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        if(oobCount.size() != rf.ext_param_.row_count_)
        {
            oobCount.resize(rf.ext_param_.row_count_, 0);
            oobErrorCount.resize(rf.ext_param_.row_count_, 0);
        }

        for(int l = 0; l < rf.ext_param_.row_count_; ++l)
        {

            if(!sm.is_used()[l])
            {
//                std::cout << "Predicting Label: " << l << std::endl;
                ++oobCount[l];
                if(rf.tree(index).predictLabel(rowVector(pr.features(), l)) != pr.response()(l,0))
                    ++oobErrorCount[l];
            }

        }
        BT::next.visit_after_tree(rf, pr, sm, st, index);
    }

    template<class RF, class PR>
    void visit_at_end(RF & rf, PR & pr)
    {
        for(int l=0; l < (int)rf.ext_param_.row_count_; ++l)
        if(oobCount[l])
        {
            oobError += double(oobErrorCount[l]) / oobCount[l];
            ++totalOobCount;
        }
        BT::visit_at_end(rf, pr);
    }

    double return_val()
    {
        return oobError/totalOobCount;
    }
};


#if 0
template <class Next = StopVisiting>
class oobVisitor: public VisitorBase<Next>
{
    typedef VisitorBase<Next> BT;



    public:
    oobVisitor(Next & next_)
    :
        BT(next_)
    {    }

    //TODO split must be const
    template<class Split, class Region>
    void visit_after_split(   Split    & split,
                            Region   & parent,
                            Region & leftChild,
                            Region & rightChild)
    {

        BT::next.visit_after_split(split, parent, leftChild, rightChild);
    }

    void visitAfterTree(    ArrayVector<DecisionTree> trees,
                            int k                        )
    {

    }

    void visitAtEnd()
    {

    }

};

template <class Next = StopVisiting>
class rankingVisitor: public VisitorBase<Next>
{
    typedef VisitorBase<Next> BT;



    public:
    oobVisitor(Next & next_)
    :
        BT(next_)
    {    }

    //TODO split must be const
    template<class Split, class Region>
    void visit_after_split(   Split    & split,
                            Region   & parent,
                            Region & leftChild,
                            Region & rightChild)
    {

        BT::next.visit_after_split(split, parent, leftChild, rightChild);
    }

    void visitAfterTree(    ArrayVector<DecisionTree> trees,
                            int k                        )
    {

    }

    void visitAtEnd()
    {

    }

};

template <class Next = StopVisiting>
class isInformativeVisitor: public VisitorBase<Next>
{
    typedef VisitorBase<Next> BT;



    public:
    oobVisitor(Next & next_)
    :
        BT(next_)
    {    }

    //TODO split must be const
    template<class Split, class Region>
    void visit_after_split(   Split    & split,
                            Region   & parent,
                            Region & leftChild,
                            Region & rightChild)
    {

        BT::next.visit_after_split(split, parent, leftChild, rightChild);
    }

    void visitAfterTree(    ArrayVector<DecisionTree> trees,
                            int k                        )
    {

    }

    void visitAtEnd()
    {

    }

};

template <class Next = StopVisiting>
class nodeToStackMappingVisitor: public VisitorBase<Next>
{
    typedef VisitorBase<Next> BT;



    public:
    oobVisitor(Next & next_)
    :
        BT(next_)
    {    }

    //TODO split must be const
    template<class Split, class Region>
    void visit_after_split(   Split    & split,
                            Region   & parent,
                            Region & leftChild,
                            Region & rightChild)
    {

        BT::next.visit_after_split(split, parent, leftChild, rightChild);
    }

    void visitAfterTree(    ArrayVector<DecisionTree> trees,
                            int k                        )
    {

    }

    void visitAtEnd()
    {

    }

};

template <class Next = StopVisiting>
class variableVisitor: public VisitorBase<Next>
{
    typedef VisitorBase<Next> BT;



    public:
    oobVisitor(Next & next_)
    :
        BT(next_)
    {    }

    //TODO split must be const
    template<class Split, class Region>
    void visit_after_split(   Split    & split,
                            Region   & parent,
                            Region & leftChild,
                            Region & rightChild)
    {

        BT::next.visit_after_split(split, parent, leftChild, rightChild);
    }

    void visitAfterTree(    ArrayVector<DecisionTree> trees,
                            int k                        )
    {

    }

    void visitAtEnd()
    {

    }

};
#endif // unmaked visitors
} // namespace vigra
#endif // RF_VISITORS_HXX
