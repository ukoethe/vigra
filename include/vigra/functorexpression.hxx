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


/** @heading Functor Creation      

    Include-File:
    \URL[vigra/functorexpression.hxx]{../include/vigra/functorexpression.hxx}\\
    Namespace: vigra::functor
    
    {\bf Note:} This functionality is not available under Microsoft Visual C++, 
    because support for partial template specialization is required.

    {\bf Motivation}
    
    Many generic algorithms are made more flexible by means of functors
    which define part of the algorithms' behavior according to the
    needs of a specific situation. For example, we can apply an exponential
    to each pixel by passing a pointer to the #exp# function 
    to #transformImage()#:
    
    \begin{verbatim}
    vigra::FImage src(w,h), dest(w,h);
    ... // fill src
    
    vigra::transformImage(srcImageRange(src), destImage(dest), &exp);    
    \end{verbatim}
    
    However, this only works for simple operations. If we wanted to 
    apply the exponential to a scaled pixel value (i.e. we want to execute
    #exp(-beta*v)#), we first need to implement a new functor:
    
    \begin{verbatim}
    struct Exponential
    {
        Exponential(double b)
        : beta(b)
        {}
        
        template <class PixelType>
        PixelType operator()(PixelType const& v) const
        {
            return exp(-beta*v);
        }
        
        double beta;
    };
    \end{verbatim}
    
    This functor would be used like this:
    
    \begin{verbatim}
    double beta =  ...;
    vigra::transformImage(srcImageRange(src), destImage(dest), 
                   Exponential(beta));    
    \end{verbatim}
    
    However, this approach has some disadvantages:
    
    \begin{itemize}
    
    \item Writing a functor is more work then simply programm the loop
          directly, i.e. non-generically. Programmers will tend to
          avoid generic constructs, if they require so much writing. 
    \item Often, functors are only needed for a single expression. 
          It is not desirable to get into the trouble of introducing 
          and documenting a new class if that class is used only once.
    \item Functors cannot be implemented directly at the point of use.
          Thus, to find out exactly what a functor is doing, one needs
          to look somewhere else. This complicates use and maintainance
          ot generic code.
    
    \end{itemize}
    
    Therefore, it is necessary to provide a means to generate functors on 
    the fly where they are needed. The C++ standard library contains so called
    "functor combinators" that allow to construct complicated functors from 
    simpler ones. The above problem "apply #exp(-beta*v)# to every pixel"
    would be solved like this:
    
    \begin{verbatim}
    float beta = ...;
    
    vigra::transformImage(srcImageRange(src), destImage(dest), 
                   std::compose1(std::ptr_fun(exp),
                                 std::bind1st(std::multiplies<float>(), -beta)));
    \end{verbatim}
 
    I won't go into details on how this works. Suffice it to say that
    this technique requires a functional programming style that is unfamiliar
    to many programmers, and thus leads to code that is difficult to 
    understand. Moreover, this technique has some limitations that prevent 
    certain expressions from being implementable this way. Therefore, VIGRA
    provides a better and simpler means to create functors on the fly.
    
    {\bf Automatic Functor Creation}
    
    Automatic functor creation in VIGRA is based on a technique called
    \URL[Expression Templates]{http://extreme.indiana.edu/~tveldhui/papers/Expression-Templates/exprtmpl.html}.
    This means that C++ operators are
    overloaded so that they don't execute the specified operation directly, 
    but instead produce a functor which will later calculate the result.
    This technique has the big advantage that the familiar operator notation
    can be used, while all the flexibility of generic programming is preserved.
    Unfortunately, it requires partial template specialization, so these capabilities
    are not available on compilers that dont support this C++ feature
    (in particular, on Microsoft Visual C++).
    
    The above problem "apply #exp(-beta*v)# to every pixel" will be solved
    like this:
    
    \begin{verbatim}
    using namespace vigra::functor;
    
    float beta = ...;
    
    transformImage(srcImageRange(src), destImage(dest), 
                   exp(Param(-beta)*Arg1()));
    \end{verbatim}
    
    Here, four expression templates have been used to create the desired
    functor:
    
    \begin{description}
    
    \item[#Param(-beta):#] creates a functor that represents a 
         constant (#-beta# in this case)
         
    \item[#Arg1():#] represents the first argument of the expression (i.e.
         the pixels of image #src# in the example). Likewise, #Arg2()# and
         #Arg3()# are defined to represent more arguments. These are needed
         for algorithms that have multiple input images, such as
         \Ref{combineTwoImages}() and \Ref{combineThreeImages}().
         
    \item[* (multiplication):] creates a functor that returns the product of
         its arguments. Likewise, the other C++ operators (i.e. 
         #+, -, *, /, %, ==, !=, <, <=, >, >=, &&, ||, &, |, ^, !, ~#) 
         are overloaded.
    
    \item[#exp():#] creates a functor that takes the exponential of its 
        argument. Likewise, the other algebraic functions
        (i.e. #sqrt, exp, log, log10, sin, asin, cos, acos, tan, 
        atan, abs, floor, ceil, rint, pow, atan2, fmod, min, max#) 
        are overloaded.
.
    
    \end{description}
    
    We will explain additional capabilities of the functor creation mechanism 
    by means of examples.
    
    The same argument can be used several times in the expression. 
    For example, to calculate the gradient magnitude from the components
    of the gradient vector, you may write:
    
    \begin{verbatim}
    using namespace vigra::functor;
    
    vigra::FImage gradient_x(w,h), gradient_y(w,h), magnitude(w,h);
    ... // calculate gradient_x and gradient_y
    
    combineTwoImages(srcImageRange(gradient_x), srcImage(gradient_y),
                     destImage(magnitude),
                     sqrt(Arg1()*Arg1() + Arg2()*Arg2()));
    \end{verbatim}
    
    It is also possible to build other functions into functor expressions. Suppose 
    you want to apply #my_complicated_function()# to the sum of two images:
    
    \begin{verbatim}
    using namespace vigra::functor;
    
    vigra::FImage src1(w,h), src2(w,h), dest(w,h);
    
    double my_complicated_function(double);
    
    combineTwoImages(srcImageRange(src1), srcImage(src2), destImage(dest),
                     applyFct(&my_complicated_function, Arg1()+Arg2()));    
    \end{verbatim}
    
    [Note that the arguments of the wrapped function are passed as additional
    arguments to #applyFct()#]
    
    You can implement conditional expression by means of the #ifThenElse()# 
    functor. It corresponds to the "? :" operator that cannot be overloaded.
    #ifThenElse()# can be used, for example, to threshold an image:
    
    \begin{verbatim}
    using namespace vigra::functor;
    
    vigra::FImage src(w,h), thresholded(w,h);
    ...// fill src
    
    float threshold = ...;
    
    transformImage(srcImageRange(src), destImage(thresholded),
                   ifThenElse(Arg1() < Param(threshold),
                              Param(0.0),    // yes branch
                              Param(1.0))    // no  branch
                  );
    \end{verbatim}

    You can use the #Var()# functor to assign values to a variable 
    (#=, +=, -=, *=, /=#&nbsp; are suported). For example, the average gray
    value of the image is calculated like this:
    
    \begin{verbatim}
    using namespace vigra::functor;
    
    vigra::FImage src(w,h);
    ...// fill src
    
    double sum = 0.0;
    
    inspectImage(srcImageRange(src), Var(sum) += Arg1());
    
    std::cout << "Average: " << (sum / (w*h)) << std::endl;
    \end{verbatim}
    
    For use in \Ref{inspectImage}() and its relatives, there is a second
    conditional functor #ifThen()# that emulates the #if()# statement
    and does not return a value. Using #ifThen()#, we can calculate the size
    of an image region:
    
    \begin{verbatim}
    using namespace vigra::functor;
    
    vigra::IImage label_image(w,h);
    ...// mark regions by labels in label_image
    
    int region_label = ...; // the region we want to inspect
    int size = 0;
    
    inspectImage(srcImageRange(label_image),
                 ifThen(Arg1() == Param(region_label),
                        Var(size) += Param(1)));
                        
    std::cout << "Size of region " << region_label << ": " << size << std::endl;
    \end{verbatim}
    
    Often, we want to execute several commands in one functor. This can be done
    by means of the overloaded #operator,()# ("operator comma"). Expressions
    seperated by a comma will be executed in succession. We can thus 
    simultaneously find the size and the average gray value of a region:
    
    \begin{verbatim}
    using namespace vigra::functor;
    
    vigra::FImage src(w,h);
    vigra::IImage label_image(w,h);
    ...// segment src and mark regions in label_image
    
    int region_label = ...; // the region we want to inspect
    int size = 0;
    double sum = 0.0;
    
    inspectTwoImages(srcImageRange(src), srcImage(label_image),
                     ifThen(Arg2() == Param(region_label),
                     (
                        Var(size) += Param(1), // the comma operator is invoked
                        Var(sum) += Arg1()
                     )));

    std::cout << "Region " << region_label << ": size = " << size << 
                                              ", average = " << sum / size << std::endl;
    \end{verbatim}
    
    [Note that the list of comma-separated expressions must be enclosed in parentheses.]
    
    A comma separated list of expressions can also be applied in the context of
    \Ref{transformImage}() and its cousins. Here, a general rule of C++ applies: The 
    return value of a comma expression is the value of its last subexpression.
    For example, we can initialize an image so that each pixel contains its 
    address in scan order:
    
    \begin{verbatim}
    using namespace vigra::functor;
    
    vigra::IImage img(w,h);
    
    int count = -1;
    
    initImage(destImageRange(img),
              (
                  Var(count) += 1,  
                  Var(count)     // this is the result of the comma expression
              ));
    \end{verbatim}
    
    Further information about how this mechanism works can be found in
    \URL[this paper]{documents/FunctorFactory.ps} (sorry, slightly out of date).
    
    @memo Expression templates to automate functor creation. 
*/

#include <cmath>
#include <vigra/numerictraits.hxx>

#if !defined(NO_PARTIAL_TEMPLATE_SPECIALIZATION)


namespace vigra {

namespace functor {

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
    T1 const & operator()(T1 const & v1) const
    {
        return v1;
    }
    
    template <class T1, class T2>
    T1 const & operator()(T1 const & v1, T2 const &) const
    {
        return v1;
    }
    
    template <class T1, class T2, class T3>
    T1 const & operator()(T1 const & v1, T2 const &, T3 const &) const
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

inline
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
    T2 const & operator()(T1 const &, T2 const & v2) const
    {
        return v2;
    }
    
    template <class T1, class T2, class T3>
    T2 const & operator()(T1 const &, T2 const & v2, T3 const &) const
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

inline
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
    T3 const & operator()(T1 const &, T2 const &, T3 const & v3) const
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

inline
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
    T const & operator()(U1 const &) const
    {
        return value_;
    }
    
    template <class U1, class U2>
    T const & operator()(U1 const &, U2 const &) const
    {
        return value_;
    }
    
    template <class U1, class U2, class U3>
    T const & operator()(U1 const &, U2 const &, U3 const &) const
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
        V & operator()(T1 const & v1) const \
        { \
            return const_cast<V &>(value_) op expr_(v1); \
        } \
         \
        template <class T1, class T2>  \
        V & operator()(T1 const & v1, T2 const & v2) const \
        { \
            return const_cast<V &>(value_) op expr_(v1, v2); \
        } \
         \
        template <class T1, class T2, class T3>  \
        V & operator()(T1 const & v1, T2 const & v2, T3 const & v3) const \
        { \
            return const_cast<V &>(value_) op expr_(v1, v2, v3); \
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
    T const & operator()(U1 const &) const
    {
        return value_;
    }
    
    template <class U1, class U2>
    T const & operator()(U1 const &, U2 const &) const
    {
        return value_;
    }
    
    template <class U1, class U2, class U3>
    T const & operator()(U1 const &, U2 const &, U3 const &) const
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
            r2(expr2_(v1));
        typename 
            ResultTraits1<IfThenElseFunctor, T>::Res 
            r3(expr3_(v1));
        return expr1_(v1) ? r2 : r3;
    }

    template <class T1, class T2> 
    typename ResultTraits2<IfThenElseFunctor, T1, T2>::Res 
    operator()(T1 const & v1, T2 const & v2) const 
    {
        typename 
            ResultTraits2<IfThenElseFunctor, T1, T2>::Res 
            r2(expr2_(v1, v2));
        typename 
            ResultTraits2<IfThenElseFunctor, T1, T2>::Res 
            r3(expr3_(v1, v2));
        return expr1_(v1, v2) ? r2 : r3;
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
    using std::function; \
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
    using std::function; \
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
makeFunctorBinaryOperator(modulo, %);
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

} // namespace functor

} // namespace vigra

#endif /* VIGRA_FUNCTOREXPRESSION_HXX  */

