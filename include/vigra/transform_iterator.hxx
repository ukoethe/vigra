


namespace vigra{


    template <typename T>
    class proxy_holder {
        T t;
    public:
        proxy_holder(const T& t) : t(t) {}
        T* operator ->() const { return &t; }
    };



    struct UseDefault{
    };


    template<class T>
    struct IdentityResultType{
        typedef T result_type;
    };



    template <
        class UnaryFunction, 
        class Iterator
    >
    class TransformIterator{

    public:
        
        typedef typename UnaryFunction::result_type function_result_type;
        typedef typename UnqualifiedType<function_result_type>::type  value_type;
        typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;


        








        TransformIterator(const Iterator & iter, const UnaryFunction & f)
        :   iter_(iter),
            f_(f){
        }


        #define TRANSFORMITERATOR_CP_OP_GEN(OP)\
        bool operator OP (const Iterator & rhs)const{\
            return iter_ OP rhs.iter_;\
        }

        TRANSFORMITERATOR_CP_OP_GEN(==);
        TRANSFORMITERATOR_CP_OP_GEN(!=);
        TRANSFORMITERATOR_CP_OP_GEN(<);
        TRANSFORMITERATOR_CP_OP_GEN(<=);
        TRANSFORMITERATOR_CP_OP_GEN(>);
        TRANSFORMITERATOR_CP_OP_GEN(>=);

        #undef TRANSFORMITERATOR_CP_OP_GEN

        TransformIterator & operator ++ (){
            ++iter_;
            return *this;
        }
        TransformIterator & operator -- (){
            --iter_;
            return *this;
        }
        TransformIterator operator ++ (int) const{
            TransformIterator res(*this);
            ++res.iter_;
            return res;
        }
        TransformIterator operator -- (int)const{
            TransformIterator res(*this);
            --res.iter_;
            return res;
        }
        TransformIterator & operator+=( const std::ptrdiff_t i ){
            iter_ += i;
            return *this;
        } 
        TransformIterator & operator-=( const std::ptrdiff_t i ){
            iter_ -= i;
            return *this;
        }
        TransformIterator operator+( const std::ptrdiff_t i )const{
            TransformIterator res(*this);
            res += i;
            return res;
        } 
        TransformIterator operator-( const std::ptrdiff_t i )const{
            TransformIterator res(*this);
            res -= i;
            return res;
        } 

        std::ptrdiff_t operator - (const TransformIterator rhs) const{
            return iter_ - rhs.iter_;
        }



    private:
        Iterator iter_;
        UnaryFunction f_;
    };
}
