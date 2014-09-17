#ifndef BLOCKWISE_ALGORITHMS_TEST_UTILS_
#define BLOCKWISE_ALGORITHMS_TEST_UTILS_

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

#endif

