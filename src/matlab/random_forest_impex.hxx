/************************************************************************/
/*                                                                      */
/*                  Copyright 2008 by Ullrich Koethe                    */
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


#ifndef VIGRA_RANDOM_FOREST_IMPEX_HXX
#define VIGRA_RANDOM_FOREST_IMPEX_HXX

#include <set>
#include <memory>
#include <vigra/matlab.hxx>
#include <vigra/random_forest.hxx>

namespace vigra {
namespace matlab {

template <class T>
std::auto_ptr<RandomForest<T> >
importRandomForest(ConstCellArray cells)
{
    // read RF parameters
    MultiArrayView<1, UInt32> parameters = getArray<UInt32>(cells[0]);
    if(parameters.size() != 3)
        mexErrMsgTxt("Parameter array must have size 3.");
    UInt32 labelCount = parameters(0);
    UInt32 numberOfFeatures = parameters(1);
    UInt32 Ntree = parameters(2);

    // read array of possible class labels
    MultiArrayView<1, double> classLabels = getArray<double>(cells[1]);
    if(labelCount != classLabels.size())
        mexErrMsgTxt("Class label array has wrong size.");

    ArrayVector<ArrayVector<Int32> >  trees;
    ArrayVector<ArrayVector<double> > weights;

    // for all decision trees
    for(UInt32 k=0; k<Ntree; ++k)
    {
        // read int tree array
        MultiArrayView<1, Int32> tree = getArray<Int32>(cells[2*k+2]);
        trees.push_back(ArrayVector<Int32>(tree.traverser_begin(), tree.traverser_end()));

        // read double weight/threshold array
        MultiArrayView<1, double> weight = getArray<double>(cells[2*k+3]);
        weights.push_back(ArrayVector<double>(weight.traverser_begin(), weight.traverser_end()));
    }
    return std::auto_ptr<RandomForest<T> >(
                           new RandomForest<T>(classLabels.traverser_begin(), classLabels.traverser_end(),
                               Ntree, numberOfFeatures, trees.begin(), weights.begin()));
}

template <class T>
void
exportRandomForest(RandomForest<T> const & rf, CellArray cells)
{
    // write RF parameters
    int parameterCount = 3;
    MultiArrayView<1, UInt32> parameters = createArray<UInt32>(parameterCount, cells[0]);
    parameters(0) = rf.labelCount();
    parameters(1) = rf.featureCount();
    parameters(2) = rf.treeCount();

    // write array of possible class labels
    MultiArrayView<1, double> classLabels = createArray<double>(rf.labelCount(), cells[1]);
    for(unsigned int k =0; k<rf.labelCount(); ++k)
        classLabels(k) = rf.classes_[k];

    // for all decision trees
    for(unsigned int k=0; k<rf.treeCount(); ++k)
    {
        // write int tree array
        MultiArrayView<1, Int32> tree =
            createArray<Int32>(rf.trees_[k].tree_.size(), cells[2*k+2]);
        for(unsigned int i =0; i<rf.trees_[k].tree_.size(); ++i)
            tree(i) = rf.trees_[k].tree_[i];

        // write double weight/threshold array
        MultiArrayView<1, double> weights =
            createArray<double>(rf.trees_[k].terminalWeights_.size(), cells[2*k+3]);
        for(unsigned int i =0; i<rf.trees_[k].terminalWeights_.size(); ++i)
            weights(i) = rf.trees_[k].terminalWeights_[i];
    }
}

}} // namespace vigra::matlab

#endif // VIGRA_RANDOM_FOREST_IMPEX_HXX
