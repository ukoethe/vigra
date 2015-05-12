


namespace vigra{




    template<class T>
    class TransformIterValProxy{
    public:
        typedef const T & reference;
        typedef const T * pointer;
        typedef T value_type;

        reference getRef(const T & functionReturn){
            t_ = functionReturn;
            return t_;
        }
        pointer getPr(const T & functionReturn){
            t_ = functionReturn;
            return &t_;
        }
    private:
        T t_;
    };


    template<class T>
    class TransformIterValProxy<const T &>{
        typedef const T & reference;
        typedef const T * pointer;
        typedef T value_type;

        reference getRef(const T & functionReturn){
            t_ = functionReturn;
            return t_;
        }
        pointer getPr(const T & functionReturn){
            t_ = functionReturn;
            return &t_;
        }
    private:
        T t_;
    };

    template<class T>
    class TransformIterValProxy<T &>{
        typedef T & reference;
        typedef T * pointer;
        typedef T value_type;

        reference getRef(const T & functionReturn){
            t_ = functionReturn;
            return t_;
        }
        pointer getPr(const T & functionReturn){
            t_ = functionReturn;
            return &t_;
        }
    private:
        T t_;
    };


    template <
        class UnaryFunction, 
        class Iterator
    >
    class TransformIterator{

    public:
        
        typedef typename UnaryFunction::result_type function_result_type;
        typedef TransformIterValProxy<function_result_type> RetHelper;

        typedef typename RetHelper::value_type  value_type;
        typedef typename RetHelper::reference   reference;
        typedef typename RetHelper::pointer     pointer;

        typedef typename std::iterator_traits<Iterator>::difference_type    difference_type;
        typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;

        TransformIterator(const Iterator & iter = Iterator(), const UnaryFunction & f = UnaryFunction())
        :   iter_(iter),
            f_(f){
        }

        reference  operator * () const{
            return retHelper_.getRef(f_(*iter_)); 
        } 

        reference  operator[](const difference_type i) const{
            return retHelper_.getRef(f_(iter_[i])); 
        }

        pointer  operator -> () const{
            return retHelper_.getRef(f_(*iter_)); 
        }


        #define TRANSFORMITERATOR_CP_OP_GEN(OP)\
        bool operator OP (const TransformIterator & rhs)const{\
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
        TransformIterator & operator+=( const difference_type i ){
            iter_ += i;
            return *this;
        } 
        TransformIterator & operator-=( const difference_type i ){
            iter_ -= i;
            return *this;
        }
        TransformIterator operator+( const difference_type i )const{
            TransformIterator res(*this);
            res += i;
            return res;
        } 
        TransformIterator operator-( const difference_type i )const{
            TransformIterator res(*this);
            res -= i;
            return res;
        } 

        difference_type operator - (const TransformIterator rhs) const{
            return iter_ - rhs.iter_;
        }

    protected:
        const Iterator & baseIterator()const{
            return iter_;
        }
        const UnaryFunction & unaryFunction()const{
            return f_;
        }
    private:
        Iterator iter_;
        UnaryFunction f_;
        mutable RetHelper retHelper_;
    };


    template <
        class UnaryFunction, 
        class Iterator
    >
    class EndAwareTransformIterator
    : public TransformIterator<UnaryFunction, Iterator>
    {
    public:
        EndAwareTransformIterator(const Iterator & iter = Iterator(), const UnaryFunction & f = UnaryFunction())
        :   TransformIterator<UnaryFunction, Iterator>(iter,f){
        }

        EndAwareTransformIterator getEndIterator()const{
            return EndAwareTransformIterator(this->baseIterator().getEndIterator(),
                                             this->unaryFunction());
        }
    private:    

    };


}
