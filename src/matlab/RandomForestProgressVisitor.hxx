#ifndef RANDOMFORESTPROGRESSVISITOR_HXX
#define RANDOMFORESTPROGRESSVISITOR_HXX


#include <iostream>
#include <iomanip>

#include <vigra/random_forest.hxx>
//#include <vigra/random_forest_hdf5_impex.hxx>

#include <vigra/timing.hxx>
namespace vigra
{
namespace rf
{
namespace visitors
{
class MatlabRandomForestProgressVisitor : public VisitorBase {
    public:
    MatlabRandomForestProgressVisitor() : VisitorBase() {}

    template<class RF, class PR, class SM, class ST>
    void visit_after_tree(RF& rf, PR & pr,  SM & sm, ST & st, int index){
        mexPrintf("%d of %d learned\n", index, rf.options().tree_count_); 
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
