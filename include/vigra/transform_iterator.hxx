


namespace vigra{


    template <typename T>
    class proxy_holder {
      T t;
    public:
      proxy_holder(const T& t) : t(t) {}
      T* operator ->() const { return &t; }
    };



    template<class REF>
    class ReturnValueHelp{
        vigra::If<>
    }



    template <class UnaryFunction, class Iterator, class Reference, class Value>
    class TransformIterator{
    public:
        TransformIterator(const Iterator & iter, const UnaryFunction & f)
        :   iter_(iter),
            f_(f){
        }


        #define TRANSFORMITERATOR_CP_OP_GEN(OP)\
        bool operator OP (const Iterator & rhs)const{\
            return iter_ OPrhs.iter_;\
        }

        TRANSFORMITERATOR_CP_OP_GEN(==);
        TRANSFORMITERATOR_CP_OP_GEN(!=);
        TRANSFORMITERATOR_CP_OP_GEN(<);
        TRANSFORMITERATOR_CP_OP_GEN(<=);
        TRANSFORMITERATOR_CP_OP_GEN(>);
        TRANSFORMITERATOR_CP_OP_GEN(>=);

        #undef TRANSFORMITERATOR_CP_OP_GEN

        Iterator & operator ++ (){
            return ++iter;
        }
        Iterator operator ++ (int){
            return iter++;
        }
        Iterator & operator -- (){
            return --iter;
        }
        Iterator operator -- (int){
            return iter--;
        }

    private:
        Iterator iter_;
        UnaryFunction f_;
        Value val_;
    }
}
