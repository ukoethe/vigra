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
#ifndef VIGRA_RF_TEST_VISITOR_HXX
#define VIGRA_RF_TEST_VISITOR_HXX

#include <vigra/random_forest.hxx>
#include <fstream>
#include <sstream>
namespace vigra
{

class TestVisitor: public VisitorBase
{

    public:
    std::ofstream fout;
    TestVisitor(std::string output = std::string("RandomForestNodeTest.log"))
    :
        fout(output.c_str())
    {    }

    //TODO split must be const
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
            std::ostringstream s1;
            s1.precision(10);
            s1.setf(std::ios::fixed,std::ios::floatfield);
            s1 << split.minGini();
            fout   << "minGini: " << s1.str().substr(0, s1.str().size()-4) << " - " << split.bestSplitColumn() << std::endl;
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
    }

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(    RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        fout << std::endl << std::endl << "Tree Number: " << index << " finished." << std::endl << std::endl;
    }

    ~TestVisitor()
    {
        fout.close();
    }
};

class SetTestVisitor: public VisitorBase
{
    public:
    std::ostringstream sout;
    std::set<std::string> treesset;


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
            sout   << "minGini: " << split.minGini() << " - " << split.bestSplitColumn() << std::endl;
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
    }

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(    RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        treesset.insert(sout.str());
        sout.str(std::string());
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
	}
};

template <class T1, class C1, class T2, class C2>
class AllOutputVisitor: public VisitorBase
{
    MultiArrayView<2, T1, C1> * features;
    MultiArrayView<2, T2, C2> * labels;


    public:
    std::ofstream fout;
    AllOutputVisitor(std::string output = std::string("RandomForestNodeTest.log"))
    :
        fout(output.c_str())
    {    }

    //TODO split must be const
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
    }

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(    RF& rf, PR & pr,  SM & sm, ST & st, int index)
    {
        fout << std::endl << std::endl << "Tree Number: " << index << " finished." << std::endl << std::endl;
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


} // namespace vigra

#endif //VIGRA_RF_TEST_VISITOR_HXX
