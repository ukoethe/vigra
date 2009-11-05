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

template<class T>
void importRandomForest(vigra::RandomForest<T> &rf,ConstCellArray cells)
{
    // read RF parameters
    MultiArrayView<1, UInt32> e_param = getArray<UInt32>(cells[0]);
	rf.ext_param_.unserialize(e_param.data(), e_param.data()+ e_param.size());

    MultiArrayView<1, UInt32> opt = getArray<UInt32>(cells[1]);
	rf.options_.unserialize(opt.data(), opt.data()+ opt.size());


	rf.trees_.resize(rf.options_.tree_count_, rf.ext_param_);
    // for all decision trees
    for(UInt32 k=0; k<rf.options_.tree_count_; ++k)
    {
		
        // read int tree array
        MultiArrayView<1, Int32> tree = getArray<Int32>(cells[2*k+2]);
        rf.tree(k).topology_.resize(tree.size());
		std::copy(tree.traverser_begin(), tree.traverser_end(),
				  rf.tree(k).topology_.begin());

		
        MultiArrayView<1, double> weight = getArray<double>(cells[2*k+3]);
        rf.tree(k).parameters_.resize(weight.size());
		std::copy(weight.traverser_begin(), weight.traverser_end(),
				  rf.tree(k).parameters_.begin());
    }
}

template <class T>
void
exportRandomForest(RandomForest<T> const & rf, CellArray cells)
{
    // write RF parameters
    int parameterCount = rf.ext_param_.serialized_size();
    MultiArrayView<1, UInt32> parameters = createArray<UInt32>(parameterCount, cells[0]);
	rf.ext_param_.serialize(parameters.data(), parameters.data()+parameterCount);


    int optCount = rf.options_.serialized_size();
    MultiArrayView<1, UInt32> opt = createArray<UInt32>(optCount, cells[0]);
	rf.options_.serialize(opt.data(), opt.data() + optCount);

    // for all decision trees
    for(int k=0; k<rf.options_.tree_count_; ++k)
    {
        // write int topology array
        MultiArrayView<1, Int32> tree =
            createArray<Int32>(rf.tree(k).topology_.size(), cells[2*k+2]);
		std::copy(rf.tree(k).topology_.begin(),
				  rf.tree(k).topology_.end(),
				  tree.data());

        // write double parameters array
        MultiArrayView<1, double> weights =
            createArray<double>(rf.tree(k).parameters_.size(), cells[2*k+3]);
		std::copy(rf.tree(k).parameters_.begin(),
				  rf.tree(k).parameters_.end(),
				  weights.data());
    }
}

}} // namespace vigra::matlab

#endif // VIGRA_RANDOM_FOREST_IMPEX_HXX
