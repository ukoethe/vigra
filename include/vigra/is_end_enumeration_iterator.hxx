#ifndef VIGRA_IS_END_ENUMERATION_ITERATOR_HXX
#define VIGRA_IS_END_ENUMERATION_ITERATOR_HXX

#include <vigra/graphs.hxx>
#include <boost/iterator/iterator_facade.hpp>

namespace vigra{



    template<class FROM_COUNTER>
    class EnumerationIterator
    : public boost::iterator_facade<
        EnumerationIterator<FROM_COUNTER>,
        const FROM_COUNTER,
        //boost::forward_traversal_tag 
        boost::random_access_traversal_tag
    >
    {
    public:

        EnumerationIterator(const lemon::Invalid  invalid)
        :   current_(0),
            size_(0){
        }


        EnumerationIterator()
        :   current_(0),
            size_(0){
        }

        EnumerationIterator(const EnumerationIterator & other)
        :   current_(other.current_),
            size_(other.size_){
        }




        explicit EnumerationIterator(const ptrdiff_t current , const ptrdiff_t size)
        :   current_(current),
            size_(size){
        }

        bool isEnd()const{
            return current_>=size_;
        }

        bool isBegin()const{
            return current_ == 0 && size_!=0;
        }


    private:
        friend class boost::iterator_core_access;

        void increment() {
            ++current_;
        }
        
        void derement() {
            --current_;
        }

        void advance(const ptrdiff_t n){
            current_+=n;
        }
        ptrdiff_t distance_to(const EnumerationIterator & other)const{
            return other.current_-current_;
        }
        
        bool equal(const EnumerationIterator & other) const{
            return current_ == other.current_ && size_==other.size_;
        }

        const FROM_COUNTER & dereference() const { 
            val_ =  FROM_COUNTER(current_);
            return val_;
        }




        ptrdiff_t current_;
        ptrdiff_t size_;

        mutable FROM_COUNTER val_;
    };


    template<class FROM_COUNTER>
    inline bool operator == (const EnumerationIterator<FROM_COUNTER> & iter,const lemon::Invalid & iv){
        return iter.isEnd(); 
    }
    template<class FROM_COUNTER>
    inline bool operator == (const lemon::Invalid & iv , const EnumerationIterator<FROM_COUNTER> & iter){
       return iter.isEnd(); 
    }

    template<class FROM_COUNTER>
    inline bool operator != (const EnumerationIterator<FROM_COUNTER> & iter,const lemon::Invalid & iv){
       return !iter.isEnd(); 
    }
    template<class FROM_COUNTER>
    inline bool operator != (const lemon::Invalid & iv , const EnumerationIterator<FROM_COUNTER> & iter){
        return !iter.isEnd(); 
    }



} // end namespace vigra






#endif //VIGRA_IS_END_ENUMERATION_ITERATOR_HXX