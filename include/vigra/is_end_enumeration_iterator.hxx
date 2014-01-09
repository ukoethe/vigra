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
            size_(0),
            offset_(0){
        }


        EnumerationIterator()
        :   current_(0),
            size_(0),
            offset_(0){
        }

        EnumerationIterator(const EnumerationIterator & other)
        :   current_(other.current_),
            size_(other.size_),
            offset_(other.offset_){
        }




        explicit EnumerationIterator(const ptrdiff_t current , const ptrdiff_t size, const ptrdiff_t offset = 0)
        :   current_(current),
            size_(size),
            offset_(offset){
        }

        bool isEnd()const{
            return current_>=size_;
        }

        bool isBegin()const{
            return current_ == 0;// && size_!=0;
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
            if(current_>size_)
                current_=size_;
        }
        ptrdiff_t distance_to(const EnumerationIterator & other)const{

            if(isEnd() && other.isEnd()){
                return 0;
            }
            else if(!isEnd() && other.isEnd()){
                return size_-current_;
            }
            else if(isEnd() && !other.isEnd()){
               return -1*(other.size_-other.current_);
            }
            else{
                return other.current_-current_;
            }

            
        }
        
        bool equal(const EnumerationIterator & other) const{
            return   (isEnd() && other.isEnd() ) || (current_ == other.current_ && size_==other.size_);
        }

        const FROM_COUNTER & dereference() const { 
            val_ =  FROM_COUNTER(current_+offset_);
            return val_;
        }




        ptrdiff_t current_;
        ptrdiff_t size_;
        ptrdiff_t offset_;

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