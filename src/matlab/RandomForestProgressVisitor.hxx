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

namespace vigra
{
namespace rf
{
namespace visitors
{
class MatlabRandomForestProgressVisitor : public VisitorBase {
    private:
    volatile int numProcessedTrees;
      
    public:
    MatlabRandomForestProgressVisitor() : VisitorBase() { numProcessedTrees = 0; }

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(RF& rf, PR & pr,  SM & sm, ST & st, int index){
        #pragma omp critical
        numProcessedTrees++;
	
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
