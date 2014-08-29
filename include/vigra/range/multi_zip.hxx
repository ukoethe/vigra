#ifndef MULTI_ZIP_HXX_
#define MULTI_ZIP_HXX_

#include <vigra/range/zip.hxx>

namespace vigra
{

//namespace zipped_multi_range_detail
//{
//
//template <class First, class Second, unsigned int N>
//class ZippedMultiRangeImpl
//  : public ZippedRange<First, Second>
//{
//private:
//    typedef ZippedRange<First, Second> base_type;
//public:
//    enum{ dimension = First::dimension };
//    typedef ZippedMultiRangeImpl<typename First::value_type, typename Second::value_type, N - 1> value_type;
//    
//
////    using base_type::base_type;
//
//    value_type front() const
//    {
//        return value_type(base_type::front());
//    }
//    value_type back() const
//    {
//        return value_type(base_type::back());
//    }
//};
//
//template <class First, class Second>
//class ZippedMultiRangeImpl<First, Second, 1>
//  : public ZippedRange<First, Second>
//{
//private:
//    typedef ZippedRange<First, Second> base_type;
//public:
//    enum{ dimension = 1 };
//
////    using base_type::base_type;
//};
//
//
//}

//template <class First, class Second>
//using ZippedMultiRange = zipped_multi_range_detail::ZippedMultiRangeImpl<First, Second, First::dimension>;

template <class First, class Second, unsigned int N=First::dimension>
class ZippedMultiRange
  : public ZippedRange<First, Second>
{
private:
    typedef ZippedRange<First, Second> base_type;
public:
    enum{ dimension = First::dimension };
    typedef ZippedMultiRange<typename First::value_type, typename Second::value_type, N - 1> value_type;
    

//    using base_type::base_type;

    ZippedMultiRange(First const & f, Second const & s)
    : base_type(f, s)
    {}
    
    ZippedMultiRange(std::pair<First, Second> const & p)
    : base_type(p)
    {}

    value_type front() const
    {
        return value_type(base_type::front());
    }
    value_type back() const
    {
        return value_type(base_type::back());
    }
};

template <class First, class Second>
class ZippedMultiRange<First, Second, 1>
  : public ZippedRange<First, Second>
{
private:
    typedef ZippedRange<First, Second> base_type;
public:
    enum{ dimension = 1 };

//    using base_type::base_type;
    ZippedMultiRange(First const & f, Second const & s)
    : base_type(f, s)
    {}
    
    ZippedMultiRange(std::pair<First, Second> const & p)
    : base_type(p)
    {}

};

template <class First, class Second>
ZippedMultiRange<First, Second> multi_zip(First first, Second second)
{
    return ZippedMultiRange<First, Second>(std::move(first), std::move(second));
}

}

#endif

