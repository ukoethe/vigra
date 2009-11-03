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
#ifndef VIGRA_RANDOM_FOREST_REGION_HXX
#define VIGRA_RANDOM_FOREST_REGION_HXX
#include <algorithm>
#include <map>
#include <numeric>
#include <iostream>
#include "vigra/mathutil.hxx"
#include "vigra/array_vector.hxx"
#include "vigra/sized_int.hxx"
#include "vigra/matrix.hxx"
#include "vigra/random.hxx"
#include "vigra/functorexpression.hxx"



namespace vigra
{


/** Standard Stackentry used to Build a Tree. Contains Information 
 * About the current region being split
 */
template <class Iter>
class DT_StackEntry
{
  public:
    typedef Iter IndexIterator;
    // tree specific stuff
    enum  ParentTag
    {
        DecisionTreeNoParent = -1
    };

	/** Address of left and Right parent in the topology container
	 */
    Int32 									leftParent;
	Int32									rightParent;
	/** rule associated with current node
	 */
    ArrayVector<std::pair<Int32, double> > 	rule;


    // RegionSpecificStuff
    ArrayVector<Int32>  					classCounts_;
    ArrayVector<double> 					weightedClassCounts_;
	bool 									classCountsIsValid;
	bool									weightedClassCountsIsValid;
	IndexIterator            				begin_,  end_;
    int                   					size_; 
	IndexIterator            				oob_begin_, oob_end_;
	int										oob_size_;

    Int32 depth()
    {
        return rule.size();
    }

    void setRange(IndexIterator s, IndexIterator e)
    {
        begin_      = s;
        end_        = e;
        size_       = e-s;
    }
    void set_oob_range(IndexIterator s, IndexIterator e)
    {
        oob_begin_   = s;
        oob_end_     = e;
        oob_size_       = e-s;
    }

    void reset()
    {
        classCountsIsValid = false;
    }

    bool  isPure()
    {
        int num = 0;

        for(int ii = 0; ii < (int)classCounts().size(); ++ii)
        {
            num += classCounts()[ii] > 0;
        }
        return num <= 1;
    }

    int&  operator[](int i)
    {
        return *(begin_+i);
    }

    IndexIterator & begin()
    {
        return begin_;
    }

    IndexIterator & end()
    {
        return end_;
    }
    IndexIterator & oob_begin()
    {
        return oob_begin_;
    }

    IndexIterator & oob_end()
    {
        return oob_end_;
    }
    ArrayVector<Int32> & classCounts()
    {
        return classCounts_;
    }
    ArrayVector<Int32> & weightedClassCounts()
    {
        return classCounts_;
    }
    bool  classCountsValid(bool u)
    {
        classCountsIsValid = u;
        return classCountsIsValid;

    }

    void classCounts(ArrayVector<Int32> in);

    DT_StackEntry( IndexIterator i, IndexIterator e,
                        int classCount,
                        Int32 lp = DecisionTreeNoParent,
                        Int32 rp = DecisionTreeNoParent)
    :
        leftParent(lp),
        rightParent(rp),
        classCounts_(classCount, 0u),
        classCountsIsValid(false),
        begin_(i),
        end_(e),
        size_(e-i)
    {}

	
    Int32 size()const
    {
        return size_;
    }


    Int32 oob_size()const
    {
        return oob_size_;
    }

};


}
//namespace vigra

#endif // VIGRA_RANDOM_FOREST_REGION_HXX
