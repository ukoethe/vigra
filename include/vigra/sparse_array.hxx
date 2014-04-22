/************************************************************************/
/*                                                                      */
/*               Copyright 2014 by Thorsten Beier                       */
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

#ifndef VIGRA_SPARSE_ARRAY
#define VIGRA_SPARSE_ARRAY


#include <map>


namespace vigra{


    template<class SA>
    class SparseReturnProxy{
    public:
        //typedef typename SA::MapType MapType;
        //typedef typename SA::ConstMapIter ConstMapIter;
        //typedef typename SA::MapIter MapIter;
        typedef typename SA::value_type value_type;
        typedef value_type & reference;
        typedef typename SA::coordinate_value_pair coordinate_value_pair;


        SparseReturnProxy(SA * sa = NULL,const size_t index = 0):
        sa_(sa),
        index_(index){

        }


        // write access
        // - convert from value_type to proxy
        SparseReturnProxy & operator=( const value_type & value ) {
            sa_->writableRef(index_)=value;
            return *this;
        }

        // read access
        // - convert from proxy to value_type
        operator value_type()const{
            return sa_->read(index_);
        }
        // macro to make ++ , ++(int) , --, --(int)
        // +=  and  -=  are also generated
        #define VIGRA_SPARSE_MAKE_OP_MACRO(OP,OP_PRE,OP_POST,SIG) \
        SparseReturnProxy<SA> & operator OP (SIG){ \
            reference val = sa_->writableRef(index_); \
            OP_PRE val OP_POST; \
            return *this; \
        }  
        VIGRA_SPARSE_MAKE_OP_MACRO(++,++, , void)
        VIGRA_SPARSE_MAKE_OP_MACRO(++,  ,++,int)
        VIGRA_SPARSE_MAKE_OP_MACRO(--,--, , void)
        VIGRA_SPARSE_MAKE_OP_MACRO(--,  ,--,int)
        VIGRA_SPARSE_MAKE_OP_MACRO(+=, ,+=value, const value_type value)
        VIGRA_SPARSE_MAKE_OP_MACRO(-=, ,-=value, const value_type value)

        #undef VIGRA_SPARSE_MAKE_OP_MACRO




        #if 0
        SparseReturnProxy<SA> & operator ++ (){
            MapIter iter = sa_->find(index_);
            if(iter==sa_->storageEnd()){
                value_type val(sa_->zeroValue());
                ++val;
                if(!sa_->eqComp_(val,sa_->zeroValue()))
                    sa_->write(index_, val);
            }
            else
                ++iter->second;
            return *this;
        }
        SparseReturnProxy<SA> & operator ++ (int){
            MapIter iter = sa_->storage_.find(index_);
            if(iter==sa_->storage_.end()){
                value_type val(sa_->zeroValue());
                val++;
                if(!sa_->eqComp_(val,sa_->zeroValue()))
                    sa_->write(index_, val);
            }
            else
                iter->second++;
            return *this;
        }


        SparseReturnProxy<SA> & operator += (const value_type & value){
            MapIter iter = sa_->storage_.find(index_);
            if(iter==sa_->storage_.end()){
                value_type val(sa_->zeroValue());
                val+=value;
                if(!sa_->eqComp_(val,sa_->zeroValue()))
                    sa_->write(index_, val);
            }
            else
                iter->second+=value;
            return *this;
        }
        SparseReturnProxy<SA> & operator -= (const value_type & value){
            MapIter iter = sa_->storage_.find(index_);
            if(iter==sa_->storage_.end()){
                value_type val(sa_->zeroValue());
                val-=value;
                if(!sa_->eqComp_(val,sa_->zeroValue()))
                    sa_->write(index_, val);
            }
            else
                iter->second-=value;
            return *this;
        }
        SparseReturnProxy<SA> & operator *= (const value_type & value){
            MapIter iter = sa_->storage_.find(index_);
            if(iter==sa_->storage_.end()){
                value_type val(sa_->zeroValue());
                val*=value;
                if(!sa_->eqComp_(val,sa_->zeroValue()))
                    sa_->write(index_, val);
            }
            else
                iter->second*=value;
            return *this;
        }

        SparseReturnProxy<SA> & operator /= (const value_type & value){
            MapIter iter = sa_->storage_.find(index_);
            if(iter==sa_->storage_.end()){
                value_type val(sa_->zeroValue());
                val/=value;
                if(!sa_->eqComp_(val,sa_->zeroValue()))
                    sa_->write(index_, val);
            }
            else
                iter->second/=value;
            return *this;
        }

        #endif

    private:
        SA * sa_;
        size_t index_;
    };



    template<class SPARSE_VECTOR, class T,class C, class EQUAL_COMP>
    class SparseVectorCrtpBase{
    private:
        friend class SparseReturnProxy< SPARSE_VECTOR >;
        typedef SparseVectorCrtpBase<SPARSE_VECTOR,T,C,EQUAL_COMP> SelfType;
        typedef SPARSE_VECTOR Child;
    public:
        typedef SparseReturnProxy<Child> Proxy;
        typedef T value_type;
        typedef const T & const_reference;
        typedef std::pair<T,C> coordinate_value_pair;


        Proxy operator()(const size_t index){
            return Proxy(getPtr(),index);
        }

        Proxy operator[](const size_t index){
            return Proxy(getPtr(),index);
        }


        const_reference operator()(const size_t index)const{
            return getRef().read(index);
        }

        const_reference operator[](const size_t index)const{
            return getRef().read(index);
        }

        const SelfType & asConst()const{
            return *this;
        }   




    private:
        const SPARSE_VECTOR * getPtr()const{
            return static_cast<SPARSE_VECTOR const *>(this);
        }
        SPARSE_VECTOR * getPtr(){
            return static_cast<SPARSE_VECTOR *>(this);
        }

        const SPARSE_VECTOR & getRef()const{
            return *static_cast<SPARSE_VECTOR const *>(this);
        }
        SPARSE_VECTOR & getRef(){
            return *static_cast<SPARSE_VECTOR *>(this);
        }
    };


    template<class T ,class EQUAL_COMP = std::equal_to<T> >
    class SparseMapVector

    : public SparseVectorCrtpBase< SparseMapVector< T,EQUAL_COMP > ,
                                   T,size_t,EQUAL_COMP >

    {


    private:
        typedef SparseMapVector<T,EQUAL_COMP> SelfType;
        typedef SparseReturnProxy<SelfType> Proxy;
        typedef EQUAL_COMP EqualCompare;
        friend class SparseReturnProxy<SparseMapVector<T> >;

        typedef std::map<size_t,T> MapType;
        typedef typename MapType::const_iterator ConstMapIter;
        typedef typename MapType::iterator MapIter;
        
    public:
        typedef std::pair<size_t,T> coordinate_value_pair;
        typedef size_t coordinate_type;
        typedef T value_type;
        typedef const T & const_reference;

        SparseMapVector(const size_t size=0, const T & zeroValue  = T(0),
                        const EqualCompare & eqComp = EqualCompare())
        :   size_(size),
            storage_(),
            zeroValue_(zeroValue),
            eqComp_(eqComp)
        {

        }

        template<class INDEX_VALUE_PAIR_ITER>
        SparseMapVector(const size_t size,INDEX_VALUE_PAIR_ITER indexValuePairBegin,
                        INDEX_VALUE_PAIR_ITER indexValuePairEnd,
                        const T & zeroValue  = T(0),
                        const EqualCompare & eqComp = EqualCompare() )
        :   size_(size),
            storage_(indexValuePairBegin,indexValuePairEnd),
            zeroValue_(zeroValue),
            eqComp_(eqComp)
        {

        }


        template<class ITER_INDEX , class ITER_VALUE>
        SparseMapVector(const size_t size, ITER_INDEX indexBegin,
                        ITER_INDEX indexEnd,ITER_VALUE valueBegin,
                        const T & zeroValue  = T(0),
                        const EqualCompare & eqComp = EqualCompare())
        :   size_(size),
            storage_(),
            zeroValue_(zeroValue),
            eqComp_(eqComp)
        {
            const size_t nVals = std::distance(indexBegin,indexEnd);
            std::vector<coordinate_value_pair> indexValuePairVec(nVals);
            for(size_t i=0;i<nVals;++i){
                coordinate_value_pair & mp = indexValuePairVec[i];
                mp.first=*indexBegin;
                mp.second=*valueBegin;
                ++indexBegin;
                ++valueBegin;
            }
            MapType tmp(indexValuePairVec.begin(),indexValuePairVec.end());
            storage_.swap(tmp);
        }



        const T & read(const size_t index)const{
            const ConstMapIter iter = storage_.find(index);
            return iter==storage_.end() ? zeroValue_ : iter->second;
        }

        T & writableRef(const size_t index){
           return storage_.insert(coordinate_value_pair(index,zeroValue())).first->second;
        }


        //template<class FUNCTOR>
        //void manipulate(const size_t index ,const,FUNCTOR){
        //    // http://stackoverflow.com/questions/97050/stdmap-insert-or-stdmap-find
        //    MapIter lb = storage_.lower_bound(index);            
        //    if(lb != storage_.end() && !(storage_.key_comp()(k, lb->first))){
        //        // key already exists
        //        // update lb->second 
        //        fFound(lb->second);
        //    }
        //    else
        //    {
        //        // the key does not exist in the map
        //        // add it to the map
        //        mymap.insert(lb, MapType::value_type(k, v));    // Use lb as a hint to insert,
        //                                                        // so it can avoid another lookup
        //    }
        //}
//
        //void write(const size_t index,const T & value){
        //    std::pair<MapIter,bool> ret=storage_.insert(coordinate_value_pair(index,value));
        //    if(ret.second==false){
        //        ret.first->second=value;
        //    }
        //}
//
        size_t size()const{
            return size_;
        }

        const_reference zeroValue()const{
            return zeroValue_;
        }


        void  swap(SelfType & other){
            storage_.swap(other.storage_);
            std::swap(size_,other.size_);
            std::swap(zeroValue_,other.zeroValue_);
            std::swap(eqComp_,other.eqComp_);
        }

    private:

        MapIter find(const size_t index){
            return storage_.find(index);
        }
        ConstMapIter find(const size_t index)const{
            return storage_.find(index);
        }


        MapIter storageEnd(){
            return storage_.end();
        }
        ConstMapIter storageEnd()const{
            return storage_.end();
        }


        size_t size_; 
        MapType storage_;
        T zeroValue_;
        EqualCompare eqComp_;
    };

} // namespace vigra

#endif /* VIGRA_SPARSE_ARRAY */
