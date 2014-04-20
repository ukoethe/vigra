/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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
#ifndef VIGRA_ITERATORFACADE_HXX
#define VIGRA_ITERATORFACADE_HXX


#include "metaprogramming.hxx"

namespace vigra {

// facade needs to make this class 
// a friend class
class IteratorFacadeCoreAccess{
public:
    template<class F>
    static bool equal(const F & fa,const F & fb){
        return fa.equal(fb);
    }

    template<class F,class REFERENCE>
    static REFERENCE dereference(const F & f){
        return f.dereference();
    }

    template<class F>
    static void increment(F & f){
        f.increment();
    }
};


// see boost iterator facade 
template<class FACADE,class VALUE_TYPE,bool IS_CONST = true>
class ForwardIteratorFacade{
private:
    
public:

    typedef std::forward_iterator_tag iterator_category;
    typedef typename UnqualifiedType<VALUE_TYPE>::type value_type;
    typedef typename IfBool<IS_CONST, value_type const * , value_type *>::type  pointer;
    typedef typename IfBool<IS_CONST, const value_type  & , value_type &>::type  reference;
    typedef std::ptrdiff_t difference_type;
    

    FACADE & operator++()
    {
        IteratorFacadeCoreAccess::increment(getF());
        return getF();
    }

    FACADE operator++(int)
    {
        FACADE res(getF());
        IteratorFacadeCoreAccess::increment(res);
        return res;
    }

    bool operator ==(const FACADE & f)const{
        return IteratorFacadeCoreAccess::equal(getF(),f);
    }
    bool operator !=(const FACADE & f)const{
        return !IteratorFacadeCoreAccess::equal(getF(),f);
    }

    reference operator*()const{
        return IteratorFacadeCoreAccess:: template dereference<FACADE,reference>(getF());
    }

    pointer operator->()const{
        return *IteratorFacadeCoreAccess:: template dereference<FACADE,reference>(getF());
    }

private:


    const FACADE & getF()const{
        return *static_cast<FACADE const *>(this);
    }
    FACADE & getF(){
        return *static_cast<FACADE *>(this);
    }
};

} // namespace vigra


#endif /* VIGRA_ITERATORFACADE_HXX */
