/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2000 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/
 
#ifndef VIGRA_FUNCTOREXPRESSION_HXX 
#define VIGRA_FUNCTOREXPRESSION_HXX 

#include <vigra/numerictraits.hxx>

#if !defined(NO_PARTIAL_TEMPLATE_SPECIALIZATION)

/************************************************************/
/*                                                          */
/*                 unary functor base template              */
/*                                                          */
/************************************************************/

template <class Operation, class T1>
struct ResultTraits1
{
    typedef T1 Res;
};

template <class Operation, class T1, class T2>
struct ResultTraits2
{
    typedef typename PromoteTraits<T1, T2>::Promote Res;
};

template <class Operation, class T1, class T2, class T3>
struct ResultTraits3
{
    typedef typename PromoteTraits<T1, T2>::Promote P1;
    typedef typename PromoteTraits<P1, T3>::Promote Res;
};

template <class EXPR>
struct UnaryFunctor
{
    UnaryFunctor(EXPR const & e)
    : expr_(e)
    {}
    
    template <class T1>
    typename ResultTraits1<EXPR, T1>::Res 
    operator()(T1 const & v) const
    {
        return expr_(v);
    }
    
    template <class T1, class T2>
    typename ResultTraits2<EXPR, T1, T2>::Res 
    operator()(T1 const & v1, T2 const & v2) const
    {
        return expr_(v1, v2);
    }
    
    template <class T1, class T2, class T3>
    typename ResultTraits3<EXPR, T1, T2, T3>::Res 
    operator()(T1 const & v1, T2 const & v2, T3 const & v3) const
    {
        return expr_(v1, v2, v3);
    }
  
  protected:  
    EXPR expr_;
};

template <class Expr, class T1>
struct ResultTraits1<UnaryFunctor<Expr>, T1>
{
    typedef typename ResultTraits1<Expr, T1>::Res Res;
};

template <class Expr, class T1, class T2>
struct ResultTraits2<UnaryFunctor<Expr>, T1, T2>
{
    typedef typename ResultTraits2<Expr, T1, T2>::Res Res;
};

template <class Expr, class T1, class T2, class T3>
struct ResultTraits3<UnaryFunctor<Expr>, T1, T2, T3>
{
    typedef typename ResultTraits3<Expr, T1, T2, T3>::Res Res;
};

/************************************************************/
/*                                                          */
/*                 unary functors for arguments             */
/*                                                          */
/************************************************************/

struct ArgumentFunctor1; 
struct ArgumentFunctor2;
struct ArgumentFunctor3;

template <>
struct UnaryFunctor<ArgumentFunctor1>
{
    UnaryFunctor()
    {}
    
    template <class T1>
    T1 operator()(T1 const & v1) const
    {
        return v1;
    }
    
    template <class T1, class T2>
    T1 operator()(T1 const & v1, T2 const &) const
    {
        return v1;
    }
    
    template <class T1, class T2, class T3>
    T1 operator()(T1 const & v1, T2 const &, T3 const &) const
    {
        return v1;
    }
};

template <class T1>
struct ResultTraits1<UnaryFunctor<ArgumentFunctor1>, T1>
{
    typedef T1 Res;
};

template <class T1, class T2>
struct ResultTraits2<UnaryFunctor<ArgumentFunctor1>, T1, T2>
{
    typedef T1 Res;
};

template <class T1, class T2, class T3>
struct ResultTraits3<UnaryFunctor<ArgumentFunctor1>, T1, T2, T3>
{
    typedef T1 Res;
};

/************************************************************/

UnaryFunctor<ArgumentFunctor1> 
Arg1()
{
    return UnaryFunctor<ArgumentFunctor1>();
}

/************************************************************/

template <>
struct UnaryFunctor<ArgumentFunctor2>
{
    UnaryFunctor()
    {}
    
    template <class T1, class T2>
    T2 operator()(T1 const &, T2 const & v2) const
    {
        return v2;
    }
    
    template <class T1, class T2, class T3>
    T2 operator()(T1 const &, T2 const & v2, T3 const &) const
    {
        return v2;
    }
};

template <class T1, class T2>
struct ResultTraits2<UnaryFunctor<ArgumentFunctor2>, T1, T2>
{
    typedef T2 Res;
};

template <class T1, class T2, class T3>
struct ResultTraits3<UnaryFunctor<ArgumentFunctor2>, T1, T2, T3>
{
    typedef T2 Res;
};

/************************************************************/

UnaryFunctor<ArgumentFunctor2> 
Arg2()
{
    return UnaryFunctor<ArgumentFunctor2>();
}

/************************************************************/

template <>
struct UnaryFunctor<ArgumentFunctor3>
{
    UnaryFunctor()
    {}
    
    template <class T1, class T2, class T3>
    T3 operator()(T1 const &, T2 const &, T3 const & v3) const
    {
        return v3;
    }
};

template <class T1, class T2, class T3>
struct ResultTraits3<UnaryFunctor<ArgumentFunctor3>, T1, T2, T3>
{
    typedef T3 Res;
};

/************************************************************/

UnaryFunctor<ArgumentFunctor3> 
Arg3()
{
    return UnaryFunctor<ArgumentFunctor3>();
}

/************************************************************/
/*                                                          */
/*                    constant parameters                   */
/*                                                          */
/************************************************************/

template <class T>
struct ParameterFunctor
{
    ParameterFunctor(T v)
    : value_(v)
    {}
    
    template <class U1>
    T operator()(U1 const &) const
    {
        return value_;
    }
    
    template <class U1, class U2>
    T operator()(U1 const &, U2 const &) const
    {
        return value_;
    }
    
    template <class U1, class U2, class U3>
    T operator()(U1 const &, U2 const &, U3 const &) const
    {
        return value_;
    }
    
  protected:
    T value_;
};

template <class T, class T1>
struct ResultTraits1<ParameterFunctor<T>, T1>
{
    typedef T Res;
};

template <class T, class T1, class T2>
struct ResultTraits2<ParameterFunctor<T>, T1, T2>
{
    typedef T Res;
};

template <class T, class T1, class T2, class T3>
struct ResultTraits3<ParameterFunctor<T>, T1, T2, T3>
{
    typedef T Res;
};

template <class T>
UnaryFunctor<ParameterFunctor<T> >
Param(T const & v)
{
    ParameterFunctor<T> fv(v);
    return UnaryFunctor<ParameterFunctor<T> >(fv);
}

/************************************************************/
/*                                                          */
/*                unary analyser base template              */
/*                                                          */
/************************************************************/


template <class EXPR>
struct UnaryAnalyser
{
    UnaryAnalyser(EXPR const & e)
    : expr_(e)
    {}
    
    template <class T1>
    void operator()(T1 const & v) const
    {
        expr_(v);
    }
    
    template <class T1, class T2>
    void operator()(T1 const & v1, T2 const & v2) const
    {
        expr_(v1, v2);
    }
    
    template <class T1, class T2, class T3>
    void operator()(T1 const & v1, T2 const & v2, T3 const & v3) const
    {
        expr_(v1, v2, v3);
    }
  protected:
  
    EXPR expr_;
};

/************************************************************/
/*                                                          */
/*                     variable assignment                  */
/*                                                          */
/************************************************************/

template <class T>
struct VarFunctor;

template <class T>
struct UnaryFunctor<VarFunctor<T> >;

/************************************************************/

#define makeAssignmentFunctor(name, op) \
    template <class V, class EXPR> \
    struct AssignmentFunctor_##name \
    { \
        AssignmentFunctor_##name(UnaryFunctor<VarFunctor<V> > v,  \
                                 UnaryFunctor<EXPR> const & e) \
        : value_(v.value_), expr_(e) \
        {} \
         \
        template <class T1>  \
        void operator()(T1 const & v1) const \
        { \
            const_cast<V &>(value_) op expr_(v1); \
        } \
         \
        template <class T1, class T2>  \
        void operator()(T1 const & v1, T2 const & v2) const \
        { \
            const_cast<V &>(value_) op expr_(v1, v2); \
        } \
         \
        template <class T1, class T2, class T3>  \
        void operator()(T1 const & v1, T2 const & v2, T3 const & v3) const \
        { \
            const_cast<V &>(value_) op expr_(v1, v2, v3); \
        } \
         \
      private: \
        V & value_; \
        UnaryFunctor<EXPR> expr_; \
    }; 

/************************************************************/

makeAssignmentFunctor(assign, =);
makeAssignmentFunctor(add, +=);
makeAssignmentFunctor(subtract, -=);
makeAssignmentFunctor(multiply, *=);
makeAssignmentFunctor(divide, /=);

/************************************************************/
/*                                                          */
/*                          variables                       */
/*                                                          */
/************************************************************/

template <class T>
struct UnaryFunctor<VarFunctor<T> >
{
    UnaryFunctor(T & v)
    : value_(v)
    {}
        
    template <class EXPR>
    UnaryAnalyser< AssignmentFunctor_assign<T, UnaryFunctor<EXPR> > >
    operator=(UnaryFunctor<EXPR> const & e)
    {
        AssignmentFunctor_assign<T, UnaryFunctor<EXPR> > va(*this, e);
        return UnaryAnalyser< AssignmentFunctor_assign<T, UnaryFunctor<EXPR> > >(va);
    }
    
    template <class EXPR>
    UnaryAnalyser< AssignmentFunctor_add<T, UnaryFunctor<EXPR> > >
    operator+=(UnaryFunctor<EXPR> const & e)
    {
        AssignmentFunctor_add<T, UnaryFunctor<EXPR> > va(*this, e);
        return UnaryAnalyser< AssignmentFunctor_add<T, UnaryFunctor<EXPR> > >(va);
    }
    
    template <class EXPR>
    UnaryAnalyser< AssignmentFunctor_subtract<T, UnaryFunctor<EXPR> > >
    operator-=(UnaryFunctor<EXPR> const & e)
    {
        AssignmentFunctor_subtract<T, UnaryFunctor<EXPR> > va(*this, e);
        return UnaryAnalyser< AssignmentFunctor_subtract<T, UnaryFunctor<EXPR> > >(va);
    }
    
    template <class EXPR>
    UnaryAnalyser< AssignmentFunctor_multiply<T, UnaryFunctor<EXPR> > >
    operator*=(UnaryFunctor<EXPR> const & e)
    {
        AssignmentFunctor_multiply<T, UnaryFunctor<EXPR> > va(*this, e);
        return UnaryAnalyser< AssignmentFunctor_multiply<T, UnaryFunctor<EXPR> > >(va);
    }
    
    template <class EXPR>
    UnaryAnalyser< AssignmentFunctor_divide<T, UnaryFunctor<EXPR> > >
    operator/=(UnaryFunctor<EXPR> const & e)
    {
        AssignmentFunctor_divide<T, UnaryFunctor<EXPR> > va(*this, e);
        return UnaryAnalyser< AssignmentFunctor_divide<T, UnaryFunctor<EXPR> > >(va);
    }
    
    template <class U1>
    T operator()(U1 const &) const
    {
        return value_;
    }
    
    template <class U1, class U2>
    T operator()(U1 const &, U2 const &) const
    {
        return value_;
    }
    
    template <class U1, class U2, class U3>
    T operator()(U1 const &, U2 const &, U3 const &) const
    {
        return value_;
    }
    
    T & value_;
};

template <class T, class T1>
struct ResultTraits1<UnaryFunctor<VarFunctor<T> >, T1>
{
    typedef T Res;
};

template <class T, class T1, class T2>
struct ResultTraits2<UnaryFunctor<VarFunctor<T> >, T1, T2>
{
    typedef T Res;
};

template <class T, class T1, class T2, class T3>
struct ResultTraits3<UnaryFunctor<VarFunctor<T> >, T1, T2, T3>
{
    typedef T Res;
};

template <class T>
UnaryFunctor<VarFunctor<T> >
Var(T & v)
{
    return UnaryFunctor<VarFunctor<T> >(v);
}

/************************************************************/
/*                                                          */
/*                          if then                         */
/*                                                          */
/************************************************************/

template <class EXPR1, class EXPR2>
struct IfThenFunctor
{
    IfThenFunctor(EXPR1 const & e1, EXPR2 const & e2)
    : expr1_(e1), expr2_(e2)
    {}
    
    template <class T> 
    void operator()(T const & v1) const 
    {
        if( expr1_(v1) ) expr2_(v1);
    }

    template <class T1, class T2> 
    void operator()(T1 const & v1, T2 const & v2) const 
    {
        if( expr1_(v1, v2) ) expr2_(v1, v2);
    }

    template <class T1, class T2, class T3> 
    void operator()(T1 const & v1, T2 const & v2, T3 const & v3) const 
    {
        if( expr1_(v1, v2, v3) ) expr2_(v1, v2, v3);
    }
    
  private:
  
    EXPR1 expr1_;
    EXPR2 expr2_;
};

template <class EXPR1, class EXPR2>
UnaryAnalyser<IfThenFunctor<UnaryFunctor<EXPR1>, 
                            UnaryAnalyser<EXPR2> > >
ifThen(UnaryFunctor<EXPR1> const & e1, 
       UnaryAnalyser<EXPR2> const & e2)
{
    IfThenFunctor<UnaryFunctor<EXPR1>, 
                  UnaryAnalyser<EXPR2> > p(e1, e2);
    return UnaryAnalyser<IfThenFunctor<UnaryFunctor<EXPR1>, 
                                       UnaryAnalyser<EXPR2> > >(p);
}

/************************************************************/
/*                                                          */
/*                         if then else                     */
/*                                                          */
/************************************************************/

template <class EXPR1, class EXPR2, class EXPR3>
struct IfThenElseFunctor;

template <class EXPR1, class EXPR2, class EXPR3, class T1>
struct ResultTraits1<IfThenElseFunctor<EXPR1, EXPR2, EXPR3>, T1>
{
    typedef typename ResultTraits1<EXPR2, T1>::Res R2;
    typedef typename ResultTraits1<EXPR3, T1>::Res R3;
    typedef typename PromoteTraits<R2, R3>::Promote Res;
};

template <class EXPR1, class EXPR2, class EXPR3, class T1, class T2>
struct ResultTraits2<IfThenElseFunctor<EXPR1, EXPR2, EXPR3>, T1, T2>
{
    typedef typename ResultTraits2<EXPR2, T1, T2>::Res R2;
    typedef typename ResultTraits2<EXPR3, T1, T2>::Res R3;
    typedef typename PromoteTraits<R2, R3>::Promote Res;
};

template <class EXPR1, class EXPR2, class EXPR3, class T1, class T2, class T3>
struct ResultTraits3<IfThenElseFunctor<EXPR1, EXPR2, EXPR3>, T1, T2, T3>
{
    typedef typename ResultTraits3<EXPR2, T1, T2, T3>::Res R2;
    typedef typename ResultTraits3<EXPR3, T1, T2, T3>::Res R3;
    typedef typename PromoteTraits<R2, R3>::Promote Res;
};

template <class EXPR1, class EXPR2, class EXPR3>
struct IfThenElseFunctor
{
    IfThenElseFunctor(EXPR1 const & e1, EXPR2 const & e2, EXPR3 const & e3)
    : expr1_(e1), expr2_(e2), expr3_(e3)
    {}
    
    template <class T> 
    typename ResultTraits1<IfThenElseFunctor, T>::Res 
    operator()(T const & v1) const 
    {
        typename 
            ResultTraits1<IfThenElseFunctor, T>::Res 
            r2(expr2_(v1, v2, v3));
        typename 
            ResultTraits1<IfThenElseFunctor, T>::Res 
            r3(expr3_(v1, v2, v3));
        return expr1_(v1, v2, v3) ? r2 : r3;
    }

    template <class T1, class T2> 
    typename ResultTraits2<IfThenElseFunctor, T1, T2>::Res 
    operator()(T1 const & v1, T2 const & v2) const 
    {
        typename 
            ResultTraits2<IfThenElseFunctor, T1, T2>::Res 
            r2(expr2_(v1, v2, v3));
        typename 
            ResultTraits2<IfThenElseFunctor, T1, T2>::Res 
            r3(expr3_(v1, v2, v3));
        return expr1_(v1, v2, v3) ? r2 : r3;
    }

    template <class T1, class T2, class T3> 
    typename ResultTraits3<IfThenElseFunctor, T1, T2, T3>::Res 
    operator()(T1 const & v1, T2 const & v2, T3 const & v3) const 
    {
        typename 
            ResultTraits3<IfThenElseFunctor, T1, T2, T3>::Res 
            r2(expr2_(v1, v2, v3));
        typename 
            ResultTraits3<IfThenElseFunctor, T1, T2, T3>::Res 
            r3(expr3_(v1, v2, v3));
        return expr1_(v1, v2, v3) ? r2 : r3;
    }
    
  private:
  
    EXPR1 expr1_;
    EXPR2 expr2_;
    EXPR3 expr3_;
};

template <class EXPR1, class EXPR2, class EXPR3>
UnaryFunctor<IfThenElseFunctor<UnaryFunctor<EXPR1>, 
                               UnaryFunctor<EXPR2>, 
                               UnaryFunctor<EXPR3> > >
ifThenElse(UnaryFunctor<EXPR1> const & e1, 
           UnaryFunctor<EXPR2> const & e2, 
           UnaryFunctor<EXPR3> const & e3)
{
    IfThenElseFunctor<UnaryFunctor<EXPR1>, 
                      UnaryFunctor<EXPR2>, 
                      UnaryFunctor<EXPR3> > p(e1, e2, e3);
    return UnaryFunctor<IfThenElseFunctor<UnaryFunctor<EXPR1>, 
                                          UnaryFunctor<EXPR2>, 
                                          UnaryFunctor<EXPR3> > >(p);
}

/************************************************************/
/*                                                          */
/*                functors for unary functions              */
/*                                                          */
/************************************************************/

#define makeFunctorUnaryFunction(function) \
    template <class EXPR> \
    struct Functor_##function; \
    \
    template <class EXPR, class T1> \
    struct ResultTraits1<Functor_##function<EXPR>, T1> \
    { \
        typedef typename ResultTraits1<EXPR, T1>::Res R1; \
        typedef typename NumericTraits<R1>::RealPromote Res; \
    }; \
    \
    template <class EXPR, class T1, class T2> \
    struct ResultTraits2<Functor_##function<EXPR>, T1, T2> \
    { \
        typedef typename ResultTraits2<EXPR, T1, T2>::Res R1; \
        typedef typename NumericTraits<R1>::RealPromote Res; \
    }; \
    \
    template <class EXPR, class T1, class T2, class T3> \
    struct ResultTraits3<Functor_##function<EXPR>, T1, T2, T3> \
    { \
        typedef typename ResultTraits3<EXPR, T1, T2, T3>::Res R1; \
        typedef typename NumericTraits<R1>::RealPromote Res; \
    }; \
    \
    template <class EXPR> \
    struct Functor_##function \
    { \
        Functor_##function(EXPR const & e) \
        : expr_(e) \
        {} \
         \
        template <class T> \
        typename ResultTraits1<Functor_##function, T>::Res \
        operator()(T const & v1) const \
        { \
            return function(expr_(v1)); \
        } \
         \
        template <class T1, class T2> \
        typename ResultTraits2<Functor_##function, T1, T2>::Res \
        operator()(T1 const & v1, T2 const & v2) const \
        { \
            return function(expr_(v1, v2)); \
        } \
         \
        template <class T1, class T2, class T3> \
        typename ResultTraits3<Functor_##function, T1, T2, T3>::Res \
        operator()(T1 const & v1, T2 const & v2, T3 const & v3) const \
        { \
            return function(expr_(v1, v2, v3)); \
        } \
         \
      protected: \
       \
        EXPR expr_; \
    }; \
     \
    template <class EXPR> \
    UnaryFunctor<Functor_##function<UnaryFunctor<EXPR> > > \
    function(UnaryFunctor<EXPR> const & e) \
    { \
        Functor_##function<UnaryFunctor<EXPR> > p(e); \
        return UnaryFunctor<Functor_##function<UnaryFunctor<EXPR> > >(p); \
    }

/************************************************************/

makeFunctorUnaryFunction(sqrt);
makeFunctorUnaryFunction(exp);
makeFunctorUnaryFunction(log);
makeFunctorUnaryFunction(log10);
makeFunctorUnaryFunction(sin);
makeFunctorUnaryFunction(asin);
makeFunctorUnaryFunction(cos);
makeFunctorUnaryFunction(acos);
makeFunctorUnaryFunction(tan);
makeFunctorUnaryFunction(atan);
makeFunctorUnaryFunction(abs);
makeFunctorUnaryFunction(floor);
makeFunctorUnaryFunction(ceil);
makeFunctorUnaryFunction(rint);

/************************************************************/
/*                                                          */
/*                functors for unary operators              */
/*                                                          */
/************************************************************/

#define makeFunctorUnaryOperator(name, op) \
    template <class EXPR> \
    struct Functor_##name; \
    \
    template <class EXPR, class T1> \
    struct ResultTraits1<Functor_##name<EXPR>, T1> \
    { \
        typedef typename ResultTraits1<EXPR, T1>::Res Res; \
    }; \
    \
    template <class EXPR, class T1, class T2> \
    struct ResultTraits2<Functor_##name<EXPR>, T1, T2> \
    { \
        typedef typename ResultTraits2<EXPR, T1, T2>::Res Res; \
    }; \
    \
    template <class EXPR, class T1, class T2, class T3> \
    struct ResultTraits3<Functor_##name<EXPR>, T1, T2, T3> \
    { \
        typedef typename ResultTraits3<EXPR, T1, T2, T3>::Res Res; \
    }; \
    \
    template <class EXPR> \
    struct Functor_##name \
    { \
        Functor_##name(EXPR const & e) \
        : expr_(e) \
        {} \
         \
        template <class T> \
        typename ResultTraits1<Functor_##name, T>::Res \
        operator()(T const & v1) const \
        { \
            return op expr_(v1); \
        } \
         \
        template <class T1, class T2> \
        typename ResultTraits2<Functor_##name, T1, T2>::Res \
        operator()(T1 const & v1, T2 const & v2) const \
        { \
            return op expr_(v1, v2); \
        } \
         \
        template <class T1, class T2, class T3> \
        typename ResultTraits3<Functor_##name, T1, T2, T3>::Res \
        operator()(T1 const & v1, T2 const & v2, T3 const & v3) const \
        { \
            return op expr_(v1, v2, v3); \
        } \
      protected: \
       \
        EXPR expr_; \
    }; \
     \
    template <class EXPR> \
    UnaryFunctor<Functor_##name<UnaryFunctor<EXPR> > > \
    operator op(UnaryFunctor<EXPR> const & e) \
    { \
        Functor_##name<UnaryFunctor<EXPR> > p(e); \
        return UnaryFunctor<Functor_##name<UnaryFunctor<EXPR> > >(p); \
    }


/************************************************************/

makeFunctorUnaryOperator(minus, -);
makeFunctorUnaryOperator(negate, !);
makeFunctorUnaryOperator(bitNegate, ~);

/************************************************************/
/*                                                          */
/*               functors for binary functions              */
/*                                                          */
/************************************************************/

#define makeFunctorBinaryFunction(function) \
    template <class EXPR1, class EXPR2> \
    struct Functor_##function; \
    \
    template <class EXPR1, class EXPR2, class T1> \
    struct ResultTraits1<Functor_##function<EXPR1, EXPR2>, T1> \
    { \
        typedef typename ResultTraits1<EXPR1, T1>::Res R1; \
        typedef typename ResultTraits1<EXPR2, T1>::Res R2; \
        typedef typename PromoteTraits<R1, R2>::Promote R3; \
        typedef typename NumericTraits<R3>::RealPromote Res; \
    }; \
    \
    template <class EXPR1, class EXPR2, class T1, class T2> \
    struct ResultTraits2<Functor_##function<EXPR1, EXPR2>, T1, T2> \
    { \
        typedef typename ResultTraits2<EXPR1, T1, T2>::Res R1; \
        typedef typename ResultTraits2<EXPR2, T1, T2>::Res R2; \
        typedef typename PromoteTraits<R1, R2>::Promote R3; \
        typedef typename NumericTraits<R3>::RealPromote Res; \
    }; \
    \
    template <class EXPR1, class EXPR2, class T1, class T2, class T3> \
    struct ResultTraits3<Functor_##function<EXPR1, EXPR2>, T1, T2, T3> \
    { \
        typedef typename ResultTraits3<EXPR1, T1, T2, T3>::Res R1; \
        typedef typename ResultTraits3<EXPR2, T1, T2, T3>::Res R2; \
        typedef typename PromoteTraits<R1, R2>::Promote R3; \
        typedef typename NumericTraits<R3>::RealPromote Res; \
    }; \
    \
    template <class EXPR1, class EXPR2> \
    struct Functor_##function \
    { \
        Functor_##function(EXPR1 const & e1, EXPR2 const & e2) \
        : expr1_(e1), expr2_(e2) \
        {} \
         \
        template <class T> \
        typename ResultTraits1<Functor_##function, T>::Res \
        operator()(T const & v1) const \
        { \
            return function(expr1_(v1), expr2_(v1)); \
        } \
         \
        template <class T1, class T2> \
        typename ResultTraits2<Functor_##function, T1, T2>::Res \
        operator()(T1 const & v1, T2 const & v2) const \
        { \
            return function(expr1_(v1, v2), expr2_(v1, v2)); \
        } \
         \
        template <class T1, class T2, class T3> \
        typename ResultTraits3<Functor_##function, T1, T2, T3>::Res \
        operator()(T1 const & v1, T2 const & v2, T3 const & v3) const \
        { \
            return function(expr1_(v1, v2, v3), expr2_(v1, v2, v3)); \
        } \
         \
      private: \
         \
        EXPR1 expr1_; \
        EXPR2 expr2_; \
    }; \
     \
    template <class EXPR1, class EXPR2> \
    UnaryFunctor<Functor_##function<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > > \
    function(UnaryFunctor<EXPR1> const & e1, UnaryFunctor<EXPR2> const & e2) \
    { \
        Functor_##function<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > p(e1, e2); \
        return UnaryFunctor<Functor_##function<UnaryFunctor<EXPR1>,  \
                                        UnaryFunctor<EXPR2> > >(p); \
    }

/************************************************************/

makeFunctorBinaryFunction(pow);
makeFunctorBinaryFunction(atan2);
makeFunctorBinaryFunction(fmod);

/************************************************************/

#define makeFunctorMinMax(name, op) \
    template <class EXPR1, class EXPR2> \
    struct Functor_##name; \
    \
    template <class EXPR1, class EXPR2, class T1> \
    struct ResultTraits1<Functor_##name<EXPR1, EXPR2>, T1> \
    { \
        typedef typename ResultTraits1<EXPR1, T1>::Res R1; \
        typedef typename ResultTraits1<EXPR2, T1>::Res R2; \
        typedef typename PromoteTraits<R1, R2>::Promote Res; \
    }; \
    \
    template <class EXPR1, class EXPR2, class T1, class T2> \
    struct ResultTraits2<Functor_##name<EXPR1, EXPR2>, T1, T2> \
    { \
        typedef typename ResultTraits2<EXPR1, T1, T2>::Res R1; \
        typedef typename ResultTraits2<EXPR2, T1, T2>::Res R2; \
        typedef typename PromoteTraits<R1, R2>::Promote Res; \
    }; \
    \
    template <class EXPR1, class EXPR2, class T1, class T2, class T3> \
    struct ResultTraits3<Functor_##name<EXPR1, EXPR2>, T1, T2, T3> \
    { \
        typedef typename ResultTraits3<EXPR1, T1, T2, T3>::Res R1; \
        typedef typename ResultTraits3<EXPR2, T1, T2, T3>::Res R2; \
        typedef typename PromoteTraits<R1, R2>::Promote Res; \
    }; \
    \
    template <class EXPR1, class EXPR2> \
    struct Functor_##name \
    { \
        Functor_##name(EXPR1 const & e1, EXPR2 const & e2) \
        : expr1_(e1), expr2_(e2) \
        {} \
         \
        template <class T> \
        typename ResultTraits1<Functor_##name, T>::Res \
        operator()(T const & v1) const \
        { \
            typename \
            ResultTraits1<Functor_##name<EXPR1, EXPR2>, T>::R1 r1(expr1_(v1)); \
            typename \
            ResultTraits1<Functor_##name<EXPR1, EXPR2>, T>::R2 r2(expr2_(v1)); \
            return (r1 op r2) ? r1 : r2; \
        } \
         \
        template <class T1, class T2> \
        typename ResultTraits2<Functor_##name, T1, T2>::Res \
        operator()(T1 const & v1, T2 const & v2) const \
        { \
            typename \
            ResultTraits2<Functor_##name<EXPR1, EXPR2>, T1, T2>::R1 r1(expr1_(v1, v2)); \
            typename \
            ResultTraits2<Functor_##name<EXPR1, EXPR2>, T1, T2>::R2 r2(expr2_(v1, v2)); \
            return (r1 op r2) ? r1 : r2; \
        } \
         \
        template <class T1, class T2, class T3> \
        typename ResultTraits3<Functor_##name, T1, T2, T3>::Res \
        operator()(T1 const & v1, T2 const & v2, T3 const & v3) const \
        { \
            typename \
            ResultTraits3<Functor_##name<EXPR1, EXPR2>, T1, T2, T3>::R1 r1(expr1_(v1, v2, v3)); \
            typename \
            ResultTraits3<Functor_##name<EXPR1, EXPR2>, T1, T2, T3>::R2 r2(expr2_(v1, v2, v3)); \
            return (r1 op r2) ? r1 : r2; \
        } \
         \
      private: \
         \
        EXPR1 expr1_; \
        EXPR2 expr2_; \
    }; \
     \
    template <class EXPR1, class EXPR2> \
    UnaryFunctor<Functor_##name<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > > \
    name(UnaryFunctor<EXPR1> const & e1, UnaryFunctor<EXPR2> const & e2) \
    { \
        Functor_##name<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > p(e1, e2); \
        return UnaryFunctor<Functor_##name<UnaryFunctor<EXPR1>,  \
                                        UnaryFunctor<EXPR2> > >(p); \
    }

makeFunctorMinMax(min, <);
makeFunctorMinMax(max, >);

/************************************************************/
/*                                                          */
/*               functors for binary operators              */
/*                                                          */
/************************************************************/

#define makeFunctorBinaryOperator(name, op) \
    template <class EXPR1, class EXPR2> \
    struct Functor_##name; \
    \
    template <class EXPR1, class EXPR2, class T1> \
    struct ResultTraits1<Functor_##name<EXPR1, EXPR2>, T1> \
    { \
        typedef typename ResultTraits1<EXPR1, T1>::Res R1; \
        typedef typename ResultTraits1<EXPR2, T1>::Res R2; \
        typedef typename PromoteTraits<R1, R2>::Promote Res; \
    }; \
    \
    template <class EXPR1, class EXPR2, class T1, class T2> \
    struct ResultTraits2<Functor_##name<EXPR1, EXPR2>, T1, T2> \
    { \
        typedef typename ResultTraits2<EXPR1, T1, T2>::Res R1; \
        typedef typename ResultTraits2<EXPR2, T1, T2>::Res R2; \
        typedef typename PromoteTraits<R1, R2>::Promote Res; \
    }; \
    \
    template <class EXPR1, class EXPR2, class T1, class T2, class T3> \
    struct ResultTraits3<Functor_##name<EXPR1, EXPR2>, T1, T2, T3> \
    { \
        typedef typename ResultTraits3<EXPR1, T1, T2, T3>::Res R1; \
        typedef typename ResultTraits3<EXPR2, T1, T2, T3>::Res R2; \
        typedef typename PromoteTraits<R1, R2>::Promote Res; \
    }; \
    \
    template <class EXPR1, class EXPR2> \
    struct Functor_##name \
    { \
        Functor_##name(EXPR1 const & e1, EXPR2 const & e2) \
        : expr1_(e1), expr2_(e2) \
        {} \
         \
        template <class T> \
        typename ResultTraits1<Functor_##name, T>::Res \
        operator()(T const & v1) const \
        { \
            return expr1_(v1) op expr2_(v1); \
        } \
         \
        template <class T1, class T2> \
        typename ResultTraits2<Functor_##name, T1, T2>::Res \
        operator()(T1 const & v1, T2 const & v2) const \
        { \
            return expr1_(v1, v2) op expr2_(v1, v2); \
        } \
         \
        template <class T1, class T2, class T3> \
        typename ResultTraits3<Functor_##name, T1, T2, T3>::Res \
        operator()(T1 const & v1, T2 const & v2, T3 const & v3) const \
        { \
            return expr1_(v1, v2, v3) op expr2_(v1, v2, v3); \
        } \
         \
      private: \
         \
        EXPR1 expr1_; \
        EXPR2 expr2_; \
    }; \
     \
    template <class EXPR1, class EXPR2> \
    UnaryFunctor<Functor_##name<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > > \
    operator op(UnaryFunctor<EXPR1> const & e1, UnaryFunctor<EXPR2> const & e2) \
    { \
        Functor_##name<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > p(e1, e2); \
        return UnaryFunctor<Functor_##name<UnaryFunctor<EXPR1>,  \
                                        UnaryFunctor<EXPR2> > >(p); \
    }

/************************************************************/

makeFunctorBinaryOperator(add, +);
makeFunctorBinaryOperator(subtract, -);
makeFunctorBinaryOperator(multiply, *);
makeFunctorBinaryOperator(divide, /);
makeFunctorBinaryOperator(bitAnd, &);
makeFunctorBinaryOperator(bitOr, |);
makeFunctorBinaryOperator(bitXor, ^);

/************************************************************/

#define makeFunctorBinaryOperatorBool(name, op) \
    template <class EXPR1, class EXPR2> \
    struct Functor_##name; \
    \
    template <class EXPR1, class EXPR2, class T1> \
    struct ResultTraits1<Functor_##name<EXPR1, EXPR2>, T1> \
    { \
        typedef bool Res; \
    }; \
    \
    template <class EXPR1, class EXPR2, class T1, class T2> \
    struct ResultTraits2<Functor_##name<EXPR1, EXPR2>, T1, T2> \
    { \
        typedef bool Res; \
    }; \
    \
    template <class EXPR1, class EXPR2, class T1, class T2, class T3> \
    struct ResultTraits3<Functor_##name<EXPR1, EXPR2>, T1, T2, T3> \
    { \
        typedef bool Res; \
    }; \
    \
    template <class EXPR1, class EXPR2> \
    struct Functor_##name \
    { \
        Functor_##name(EXPR1 const & e1, EXPR2 const & e2) \
        : expr1_(e1), expr2_(e2) \
        {} \
         \
        template <class T> \
        bool operator()(T const & v1) const \
        { \
            return expr1_(v1) op expr2_(v1); \
        } \
         \
        template <class T1, class T2> \
        bool operator()(T1 const & v1, T2 const & v2) const \
        { \
            return expr1_(v1, v2) op expr2_(v1, v2); \
        } \
         \
        template <class T1, class T2, class T3> \
        bool operator()(T1 const & v1, T2 const & v2, T3 const & v3) const \
        { \
            return expr1_(v1, v2, v3) op expr2_(v1, v2, v3); \
        } \
         \
      private: \
         \
        EXPR1 expr1_; \
        EXPR2 expr2_; \
    }; \
     \
    template <class EXPR1, class EXPR2> \
    UnaryFunctor<Functor_##name<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > > \
    operator op(UnaryFunctor<EXPR1> const & e1, UnaryFunctor<EXPR2> const & e2) \
    { \
        Functor_##name<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > p(e1, e2); \
        return UnaryFunctor<Functor_##name<UnaryFunctor<EXPR1>,  \
                                        UnaryFunctor<EXPR2> > >(p); \
    }

/************************************************************/

makeFunctorBinaryOperatorBool(equals, ==);
makeFunctorBinaryOperatorBool(differs, !=);
makeFunctorBinaryOperatorBool(less, <);
makeFunctorBinaryOperatorBool(lessEqual, <=);
makeFunctorBinaryOperatorBool(greater, >);
makeFunctorBinaryOperatorBool(greaterEqual, >=);
makeFunctorBinaryOperatorBool(and, &&);
makeFunctorBinaryOperatorBool(or, ||);

/************************************************************/
/*                                                          */
/*                         unary apply                      */
/*                                                          */
/************************************************************/

template <class EXPR, class RES, class ARG>
struct UnaryFctPtrFunctor
{
    UnaryFctPtrFunctor(EXPR const & e, RES (*fct)(ARG))
    : expr_(e), f_(fct)
    {}
    
    template <class T> 
    RES operator()(T const & v1) const 
    {
        return f_(expr_(v1));
    }
    
    template <class T1, class T2> 
    RES operator()(T1 const & v1, T2 const & v2) const 
    {
        return f_(expr_(v1, v2));
    }
    
    template <class T1, class T2, class T3> 
    RES operator()(T1 const & v1, T2 const & v2, T3 const & v3) const 
    {
        return f_(expr_(v1, v2, v3));
    }
  protected:
  
    EXPR expr_;
    RES (*f_)(ARG);
};

template <class EXPR, class RES, class ARG, class T1>
struct ResultTraits1<UnaryFctPtrFunctor<EXPR, RES, ARG>, T1>
{
    typedef RES Res;
};

template <class EXPR, class RES, class ARG, class T1, class T2>
struct ResultTraits2<UnaryFctPtrFunctor<EXPR, RES, ARG>, T1, T2>
{
    typedef RES Res;
};

template <class EXPR, class RES, class ARG, class T1, class T2, class T3>
struct ResultTraits3<UnaryFctPtrFunctor<EXPR, RES, ARG>, T1, T2, T3>
{
    typedef RES Res;
};

template <class EXPR, class RES, class ARG>
UnaryFunctor<UnaryFctPtrFunctor<UnaryFunctor<EXPR>, RES, ARG> >
applyFct(RES (*f)(ARG), UnaryFunctor<EXPR> const & e)
{
    UnaryFctPtrFunctor<UnaryFunctor<EXPR>, RES, ARG> p(e, f);
    return UnaryFunctor<UnaryFctPtrFunctor<UnaryFunctor<EXPR>, RES, ARG> >(p);
}

/************************************************************/
/*                                                          */
/*                        binary apply                      */
/*                                                          */
/************************************************************/

template <class EXPR1, class EXPR2, class RES, class ARG1, class ARG2>
struct BinaryFctPtrFunctor
{
    BinaryFctPtrFunctor(EXPR1 const & e1, EXPR2 const & e2, 
                        RES (*f)(ARG1, ARG2))
    : expr1_(e1), expr2_(e2), f_(f)
    {}
    
    template <class T> 
    RES operator()(T const & v1) const 
    {
        return f_(expr1_(v1), expr2_(v1));
    }
    
    template <class T1, class T2> 
    RES operator()(T1 const & v1, T2 const & v2) const 
    {
        return f_(expr1_(v1, v2), expr2_(v1, v2));
    }
    
    template <class T1, class T2, class T3> 
    RES operator()(T1 const & v1, T2 const & v2, T3 const & v3) const 
    {
        return f_(expr1_(v1, v2, v3), expr2_(v1, v2, v3));
    }
  protected:
  
    EXPR1 expr1_;
    EXPR2 expr2_;
    RES (*f_)(ARG1, ARG2);
};

template <class EXPR1, class EXPR2, class RES, class ARG1, class ARG2, 
          class T1>
struct ResultTraits1<BinaryFctPtrFunctor<EXPR1, EXPR2, RES, ARG1, ARG2>, T1>
{
    typedef RES Res;
};

template <class EXPR1, class EXPR2, class RES, class ARG1, class ARG2, 
          class T1, class T2>
struct ResultTraits2<BinaryFctPtrFunctor<EXPR1, EXPR2, RES, ARG1, ARG2>, T1, T2>
{
    typedef RES Res;
};

template <class EXPR1, class EXPR2, class RES, class ARG1, class ARG2, 
          class T1, class T2, class T3>
struct ResultTraits3<BinaryFctPtrFunctor<EXPR1, EXPR2, RES, ARG1, ARG2>, T1, T2, T3>
{
    typedef RES Res;
};

template <class EXPR1, class EXPR2, class RES, class ARG1, class ARG2>
UnaryFunctor<BinaryFctPtrFunctor<UnaryFunctor<EXPR1>, 
                                 UnaryFunctor<EXPR2>, 
                                 RES, ARG1, ARG2> >
applyFct(RES (*f)(ARG1, ARG2), UnaryFunctor<EXPR1> const & e1, 
         UnaryFunctor<EXPR2> const & e2)
{
    BinaryFctPtrFunctor<UnaryFunctor<EXPR1>, 
                        UnaryFunctor<EXPR2>, 
                        RES, ARG1, ARG2>  p(e1, e2, f);
    return UnaryFunctor<BinaryFctPtrFunctor<UnaryFunctor<EXPR1>, 
                                            UnaryFunctor<EXPR2>, 
                                            RES, ARG1, ARG2> >(p);
}

/************************************************************/
/*                                                          */
/*                      comma operator                      */
/*                                                          */
/************************************************************/

template <class EXPR1, class EXPR2>
struct CommaFunctor
{
    CommaFunctor(EXPR1 const & e1, EXPR2 const & e2)
    : expr1_(e1), expr2_(e2)
    {}
    
    template <class T> 
    typename ResultTraits1<EXPR2, T>::Res 
    operator()(T const & v1) const 
    {
        expr1_(v1);
        return expr2_(v1);
    }
    
    template <class T1, class T2> 
    typename ResultTraits2<EXPR2, T1, T2>::Res 
    operator()(T1 const & v1, T2 const & v2) const 
    {
        expr1_(v1, v2);
        return expr2_(v1, v2);
    }
    
    template <class T1, class T2, class T3> 
    typename ResultTraits3<EXPR2, T1, T2, T3>::Res 
    operator()(T1 const & v1, T2 const & v2, T3 const & v3) const 
    {
        expr1_(v1, v2, v3);
        return expr2_(v1, v2, v3);
    }
    
  protected:
  
    EXPR1 expr1_;
    EXPR2 expr2_;
};

template <class Expr1, class Expr2, class T1>
struct ResultTraits1<CommaFunctor<Expr1, Expr2>, T1>
{
    typedef typename ResultTraits1<Expr2, T1>::Res Res;
};

template <class Expr1, class Expr2, class T1, class T2>
struct ResultTraits2<CommaFunctor<Expr1, Expr2>, T1, T2>
{
    typedef typename ResultTraits2<Expr2, T1, T2>::Res Res;
};

template <class Expr1, class Expr2, class T1, class T2, class T3>
struct ResultTraits3<CommaFunctor<Expr1, Expr2>, T1, T2, T3>
{
    typedef typename ResultTraits3<Expr2, T1, T2, T3>::Res Res;
};

template <class EXPR1, class EXPR2>
UnaryFunctor<CommaFunctor<UnaryAnalyser<EXPR1>, 
                            UnaryFunctor<EXPR2> > >
operator,(UnaryAnalyser<EXPR1> const & e1, 
          UnaryFunctor<EXPR2> const & e2)
{
    CommaFunctor<UnaryAnalyser<EXPR1>, 
                            UnaryFunctor<EXPR2> >  p(e1, e2);
    return UnaryFunctor<CommaFunctor<UnaryAnalyser<EXPR1>, 
                            UnaryFunctor<EXPR2> > >(p);
}

/************************************************************/

template <class EXPR1, class EXPR2>
struct CommaAnalyser
{
    CommaAnalyser(EXPR1 const & e1, EXPR2 const & e2)
    : expr1_(e1), expr2_(e2)
    {}
    
    template <class T> 
    void operator()(T const & v1) const 
    {
        expr1_(v1);
        expr2_(v1);
    }
    
    template <class T1, class T2> 
    void operator()(T1 const & v1, T2 const & v2) const 
    {
        expr1_(v1, v2);
        expr2_(v1, v2);
    }
    
    template <class T1, class T2, class T3> 
    void operator()(T1 const & v1, T2 const & v2, T3 const & v3) const 
    {
        expr1_(v1, v2, v3);
        expr2_(v1, v2, v3);
    }
    
  protected:
  
    EXPR1 expr1_;
    EXPR2 expr2_;
};

template <class EXPR1, class EXPR2>
UnaryAnalyser<CommaAnalyser<UnaryAnalyser<EXPR1>, 
                            UnaryAnalyser<EXPR2> > >
operator,(UnaryAnalyser<EXPR1> const & e1, 
          UnaryAnalyser<EXPR2> const & e2)
{
    CommaAnalyser<UnaryAnalyser<EXPR1>, 
                            UnaryAnalyser<EXPR2> >  p(e1, e2);
    return UnaryAnalyser<CommaAnalyser<UnaryAnalyser<EXPR1>, 
                            UnaryAnalyser<EXPR2> > >(p);
}

#endif /* NO_PARTIAL_TEMPLATE_SPECIALIZATION */

#endif /* VIGRA_FUNCTOREXPRESSION_HXX  */
