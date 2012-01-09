#ifndef RANDOMFORESTPROGRESSVISITOR_HXX
#define RANDOMFORESTPROGRESSVISITOR_HXX


#include <iostream>
#include <iomanip>

#include <vigra/random_forest.hxx>
//#include <vigra/random_forest_hdf5_impex.hxx>

#include <vigra/timing.hxx>

#ifdef _OPENMP
  #include <omp.h>
#endif

struct MatlabRandomForestTreeInfo
{
    unsigned int maxDepthLeft;
    unsigned int maxDepthRight; //depth to right and left of first node
};

namespace vigra
{
namespace rf
{
namespace visitors
{
class MatlabRandomForestProgressVisitor : public VisitorBase {
private:
    volatile int numProcessedTrees;
    bool         computeMaxTreeDepth;

    template<class RF>
    unsigned int maxTreeDepth( RF& rf, int treeIndex, bool leftChild )
    {
        typedef vigra::detail::DecisionTree            DecisionTree_t;
        const DecisionTree_t &tree = rf.trees_[treeIndex];

        typedef typename DecisionTree_t::TreeInt    TreeInt;

        unsigned int maxDepth = 0;

        std::vector<TreeInt> toExplore;
        std::vector<unsigned int> theDepth;

        int firstChild = 1;
        if (leftChild)
            firstChild = 0;

        bool isFirstRound = true;

        toExplore.push_back(2);
        theDepth.push_back(1);

        while( !toExplore.empty() )
        {
            unsigned int thisDepth = theDepth.back();
            TreeInt index = toExplore.back();

            theDepth.pop_back();
            toExplore.pop_back();

            // leaf node?
            if ( (tree.topology_[index] & LeafNodeTag) == LeafNodeTag ) {
                if (maxDepth < thisDepth)
                    maxDepth = thisDepth;
                continue;
            }

            {
                switch(tree.topology_[index])
                {
                    case i_ThresholdNode:
                    {
                        Node<i_ThresholdNode>
                                    node(tree.topology_, tree.parameters_, index);

                        if (isFirstRound) {
                            toExplore.push_back(node.child(firstChild));
                            theDepth.push_back(thisDepth+1);
                        } else {
                            toExplore.push_back( node.child(0) );
                            toExplore.push_back( node.child(1) );

                            theDepth.push_back( thisDepth + 1 );
                            theDepth.push_back( thisDepth + 1 );
                        }
                        break;
                    }
                    case i_HyperplaneNode:
                    {
                        Node<i_HyperplaneNode>
                                    node(tree.topology_, tree.parameters_, index);
                        if (isFirstRound) {
                            toExplore.push_back(node.child(firstChild));
                            theDepth.push_back(thisDepth+1);
                        } else {
                            toExplore.push_back( node.child(0) );
                            toExplore.push_back( node.child(1) );

                            theDepth.push_back( thisDepth + 1 );
                            theDepth.push_back( thisDepth + 1 );
                        }
                        break;
                    }
                    case i_HypersphereNode:
                    {
                        Node<i_HypersphereNode>
                                    node(tree.topology_, tree.parameters_, index);

                        if (isFirstRound) {
                            toExplore.push_back(node.child(firstChild));
                            theDepth.push_back(thisDepth+1);
                        } else {
                            toExplore.push_back( node.child(0) );
                            toExplore.push_back( node.child(1) );

                            theDepth.push_back( thisDepth + 1 );
                            theDepth.push_back( thisDepth + 1 );
                        }

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

            isFirstRound = false;
        }

        return maxDepth;
    }

public:
    void setComputeMaxTreeDepth(bool yes) { computeMaxTreeDepth = yes; }
    MatlabRandomForestProgressVisitor() : VisitorBase() { numProcessedTrees = 0; computeMaxTreeDepth = false; }
    std::vector<MatlabRandomForestTreeInfo>  treeInfo;

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(RF& rf, PR & pr,  SM & sm, ST & st, int treeIndex){
        #pragma omp critical
        {
            if (treeInfo.size() != rf.options().tree_count_)
                treeInfo.resize( rf.options().tree_count_ );

            if (computeMaxTreeDepth)
            {
                treeInfo.at(numProcessedTrees).maxDepthLeft = maxTreeDepth( rf, treeIndex, true );
                treeInfo.at(numProcessedTrees).maxDepthRight = maxTreeDepth( rf, treeIndex, false );
            }
            numProcessedTrees++;
        }
	
	int totalTrees = rf.options().tree_count_;
	
#ifdef _OPENMP
        if (omp_get_thread_num() == 0) {	// main thread?
#endif
        mexPrintf("%d of %d learned (%.1f%%)\n", numProcessedTrees, totalTrees, numProcessedTrees * 100.0 / totalTrees); 
#ifdef _OPENMP
        mexEvalString("drawnow");
      }
#endif
    }
    
    template<class RF, class PR>
    void visit_at_end(RF const & rf, PR const & pr) {
        mexPrintf("all trees learned\n");
    }
    
    template<class RF, class PR>
    void visit_at_beginning(RF const & rf, PR const & pr) {
        mexPrintf("growing random forest, which will have %d trees\n", rf.options().tree_count_);
    }
};
}
}
}
#endif //RANDOMFORESTPROGRESSVISITOR_HXX
