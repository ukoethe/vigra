#ifndef MULTI_ARRAY_RANGE_HXX_
#define MULTI_ARRAY_RANGE_HXX_

#include <vigra/range/strided_ptr_range.hxx>

namespace vigra
{

template <class ArrayType>
StridedPtrRange<typename ArrayType::value_type, typename ArrayType::difference_type> multi_range(ArrayType& array)
{
    return { array.data(), array.stride(), array.shape() };
}

}

#endif

