/************************************************************************/
/*                                                                      */
/*     Copyright 2013-2014 by Martin Bidlingmaier and Ullrich Koethe    */
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

#ifndef VIGRA_BLOCKWISE_ALGORITHMS_TEST_UTILS_HXX
#define VIGRA_BLOCKWISE_ALGORITHMS_TEST_UTILS_HXX

#include <cstdlib>
#include <vector>

template <class Iterator1,class Iterator2>
bool equivalentLabels(Iterator1 begin1, Iterator1 end1,
                      Iterator2 begin2, Iterator2 end2)
{
    using namespace std;

    if(end1 - begin1 != end2 - begin2)
        return false;
    
    typedef vector<int> LabelMap;
    LabelMap left_to_right; // from range 1 to range 2
    LabelMap right_to_left;
    for( ; begin1 != end1; ++begin1, ++begin2)
    {
        if(left_to_right.size() <= *begin1)
            left_to_right.resize(*begin1 + 1, -1); // "-1" means unmapped
        if(right_to_left.size() <= *begin2)
            right_to_left.resize(*begin2 + 1, -1);

        if(left_to_right[*begin1] == -1) // unmapped -> map it
            left_to_right[*begin1] = *begin2;
        else if(left_to_right[*begin1] != *begin2) // already mapped to different value -> not equivalent labels
            return false;

        if(right_to_left[*begin2] == -1)
            right_to_left[*begin2] = *begin1;
        else if(right_to_left[*begin2] != *begin1)
            return false;
    }
    
    return true;
}

template <class Iterator>
void fillRandom(Iterator begin, Iterator end, int maximum)
{
    using namespace std;

    for( ; begin != end; ++begin)
        *begin = rand() % maximum;
}

#endif // VIGRA_BLOCKWISE_ALGORITHMS_TEST_UTILS_HXX
