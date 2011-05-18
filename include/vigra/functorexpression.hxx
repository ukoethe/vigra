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
 
#ifndef VIGRA_FUNCTOREXPRESSION_HXX 
#define VIGRA_FUNCTOREXPRESSION_HXX 


/** \page FunctorExpressions Functor Expressions  

    Simple automatic functor creation by means of expression templates
    (also known as a "lambda library").    

    <b>\#include</b> \<vigra/functorexpression.hxx\><br>
    Namespace: vigra::functor
    
    <b> Motivation</b>
    
    Many generic algorithms are made more flexible by means of functors
    which define part of the algorithms' behavior according to the
    needs of a specific situation. For example, we can apply an exponential
    to each pixel by passing a pointer to the <TT>exp</TT> function 
    to <TT>transformImage()</TT>:
    
    \code
    vigra::FImage src(w,h), dest(w,h);
    ... // fill src
    
    vigra::transformImage(srcImageRange(src), destImage(dest), &exp);    
    \endcode
    
    However, this only works for simple operations. If we wanted to 
    apply the exponential to a scaled pixel value (i.e. we want to execute
    <TT>exp(-beta*v)</TT>), we first need to implement a new functor:
    
    \code
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
    \endcode
    
    This functor would be used like this:
    
    \code
    double beta =  ...;
    vigra::transformImage(srcImageRange(src), destImage(dest), 
                   Exponential(beta));    
    \endcode
    
    However, this approach has some disadvantages:
    
    <UL>
    
    <li> Writing a functor is more work then simply programm the loop
          directly, i.e. non-generically. Programmers will tend to
          avoid generic constructs, if they require so much writing. 
    <li> Often, functors are only needed for a single expression. 
          It is not desirable to get into the trouble of introducing 
          and documenting a new class if that class is used only once.
    <li> Functors cannot be implemented directly at the point of use.
          Thus, to find out exactly what a functor is doing, one needs
          to look somewhere else. This complicates use and maintainance
          ot generic code.
    
    </UL>
    
    Therefore, it is necessary to provide a means to generate functors on 
    the fly where they are needed. The C++ standard library contains so called
    "functor combinators" that allow to construct complicated functors from 
    simpler ones. The above problem "apply <TT>exp(-beta*v)</TT> to every pixel"
    would be solved like this:
    
    \code
    float beta = ...;
    
    vigra::transformImage(srcImageRange(src), destImage(dest), 
                   std::compose1(std::ptr_fun(exp),
                                 std::bind1st(std::multiplies<float>(), -beta)));
    \endcode
 
    I won't go into details on how this works. Suffice it to say that
    this technique requires a functional programming style that is unfamiliar
    to many programmers, and thus leads to code that is difficult to 
    understand. Moreover, this technique has some limitations that prevent 
    certain expressions from being implementable this way. Therefore, VIGRA
    provides a better and simpler means to create functors on the fly.
    
    <b> Automatic Functor Creation</b>
    
    Automatic functor creation in VIGRA is based on a technique called
    <a href="http://extreme.indiana.edu/~tveldhui/papers/Expression-Templates/exprtmpl.html">Expression Templates</a>.
    This means that C++ operators are
    overloaded so that they don't execute the specified operation directly, 
    but instead produce a functor which will later calculate the result.
    This technique has the big advantage that the familiar operator notation
    can be used, while all the flexibility of generic programming is preserved.
    
    The above problem "apply <TT>exp(-beta*v)</TT> to every pixel" will be solved
    like this:
    
    \code
    using namespace vigra::functor;
    
    float beta = ...;
    
    transformImage(srcImageRange(src), destImage(dest), 
                   exp(Param(-beta)*Arg1()));
    \endcode
    
    Here, four expression templates have been used to create the desired
    functor:
    
    <DL>
    
    <DT><b><TT>Param(-beta):</TT></b><DD> creates a functor that represents a 
         constant (<TT>-beta</TT> in this case)
         
    <DT><b><TT>Arg1():</TT></b><DD> represents the first argument of the expression (i.e.
         the pixels of image <TT>src</TT> in the example). Likewise, <TT>Arg2()</TT> and
         <TT>Arg3()</TT> are defined to represent more arguments. These are needed
         for algorithms that have multiple input images, such as
         \ref combineTwoImages() and \ref combineThreeImages().
         
    <DT><b>* (multiplication):</b><DD> creates a functor that returns the product of
         its arguments. Likewise, the other C++ operators (i.e. 
         <TT>+, -, *, /, %, ==, !=, <, <=, >, >=, &&, ||, &, |, ^, !, ~</TT>) 
         are overloaded.
    
    <DT><b><TT>exp():</TT></b><DD> creates a functor that takes the exponential of its 
        argument. Likewise, the other algebraic functions
        (i.e. <TT>sq, sqrt, exp, log, log10, sin, asin, cos, acos, tan, 
        atan, abs, floor, ceil, pow, atan2, fmod, min, max</TT>) 
        are overloaded.
    
    </DL>
    
    We will explain additional capabilities of the functor creation mechanism 
    by means of examples.
    
    The same argument can be used several times in the expression. 
    For example, to calculate the gradient magnitude from the components
    of the gradient vector, you may write:
    
    \code
    using namespace vigra::functor;
    
    vigra::FImage gradient_x(w,h), gradient_y(w,h), magnitude(w,h);
    ... // calculate gradient_x and gradient_y
    
    combineTwoImages(srcImageRange(gradient_x), srcImage(gradient_y),
                     destImage(magnitude),
                     sqrt(Arg1()*Arg1() + Arg2()*Arg2()));
    \endcode
    
    It is also possible to build other functions into functor expressions. Suppose 
    you want to apply <TT>my_complicated_function()</TT> to the sum of two images:
    
    \code
    using namespace vigra::functor;
    
    vigra::FImage src1(w,h), src2(w,h), dest(w,h);
    
    double my_complicated_function(double);
    
    combineTwoImages(srcImageRange(src1), srcImage(src2), destImage(dest),
                     applyFct(&my_complicated_function, Arg1()+Arg2()));    
    \endcode
    
    [Note that the arguments of the wrapped function are passed as additional
    arguments to <TT>applyFct()</TT>]
    
    You can implement conditional expression by means of the <TT>ifThenElse()</TT> 
    functor. It corresponds to the "? :" operator that cannot be overloaded.
    <TT>ifThenElse()</TT> can be used, for example, to threshold an image:
    
    \code
    using namespace vigra::functor;
    
    vigra::FImage src(w,h), thresholded(w,h);
    ...// fill src
    
    float threshold = ...;
    
    transformImage(srcImageRange(src), destImage(thresholded),
                   ifThenElse(Arg1() < Param(threshold),
                              Param(0.0),    // yes branch
                              Param(1.0))    // no  branch
                  );
    \endcode

    You can use the <TT>Var()</TT> functor to assign values to a variable 
    (<TT>=, +=, -=, *=, /=</TT>&nbsp; are suported). For example, the average gray
    value of the image is calculated like this:
    
    \code
    using namespace vigra::functor;
    
    vigra::FImage src(w,h);
    ...// fill src
    
    double sum = 0.0;
    
    inspectImage(srcImageRange(src), Var(sum) += Arg1());
    
    std::cout << "Average: " << (sum / (w*h)) << std::endl;
    \endcode
    
    For use in \ref inspectImage() and its relatives, there is a second
    conditional functor <TT>ifThen()</TT> that emulates the <TT>if()</TT> statement
    and does not return a value. Using <TT>ifThen()</TT>, we can calculate the size
    of an image region:
    
    \code
    using namespace vigra::functor;
    
    vigra::IImage label_image(w,h);
    ...// mark regions by labels in label_image
    
    int region_label = ...; // the region we want to inspect
    int size = 0;
    
    inspectImage(srcImageRange(label_image),
                 ifThen(Arg1() == Param(region_label),
                        Var(size) += Param(1)));
                        
    std::cout << "Size of region " << region_label << ": " << size << std::endl;
    \endcode
    
    Often, we want to execute several commands in one functor. This can be done
    by means of the overloaded <TT>operator,()</TT> ("operator comma"). Expressions
    seperated by a comma will be executed in succession. We can thus 
    simultaneously find the size and the average gray value of a region:
    
    \code
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
    \endcode
    
    [Note that the list of comma-separated expressions must be enclosed in parentheses.]
    
    A comma separated list of expressions can also be applied in the context of
    \ref transformImage() and its cousins. Here, a general rule of C++ applies: The 
    return value of a comma expression is the value of its last subexpression.
    For example, we can initialize an image so that each pixel contains its 
    address in scan order:
    
    \code
    using namespace vigra::functor;
    
    vigra::IImage img(w,h);
    
    int count = -1;
    
    initImageWithFunctor(destImageRange(img),
                         (
                              Var(count) += Param(1),  
                              Var(count)     // this is the result of the comma expression
                         ));
    \endcode
    
    Further information about how this mechanism works can be found in
    <a href="documents/FunctorFactory.ps">this paper</a> (sorry, slightly out of date).
*/

#ifndef DOXYGEN

#if !defined(NO_PARTIAL_TEMPLATE_SPECIALIZATION)

#include <cmath>
#include "numerictraits.hxx"
#include "mathutil.hxx"
#include "functortraits.hxx"


namespace vigra {

namespace functor {

/************************************************************/
/*                                                          */
/*                 unary functor base template              */
/*                                                          */
/************************************************************/


struct ErrorType;

template <class Operation>
struct ResultTraits0;

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
    
//    typename ResultTraits0<EXPR>::Res 
    typename ResultTraits0<EXPR>::Res 
    operator()() const
    {
        return expr_();
    }
    
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
  
  private:
    UnaryFunctor & operator=(UnaryFunctor const &); // not implemented
};

template <class Expr>
struct ResultTraits0<UnaryFunctor<Expr> >
{
    typedef typename ResultTraits0<Expr>::Res Res;
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
  
  private:
    UnaryFunctor & operator=(UnaryFunctor const &); // not implemented
};

typedef UnaryFunctor<ArgumentFunctor1> Identity;

template <>
struct ResultTraits0<UnaryFunctor<ArgumentFunctor1> >
{
    typedef ErrorType Res;
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
  
  private:
    UnaryFunctor & operator=(UnaryFunctor const &); // not implemented
};

template <>
struct ResultTraits0<UnaryFunctor<ArgumentFunctor2> >
{
    typedef ErrorType Res;
};

template <class T1>
struct ResultTraits1<UnaryFunctor<ArgumentFunctor2>, T1>
{
    typedef ErrorType Res;
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
  
  private:
    UnaryFunctor & operator=(UnaryFunctor const &); // not implemented
};

template <>
struct ResultTraits0<UnaryFunctor<ArgumentFunctor3> >
{
    typedef ErrorType Res;
};

template <class T1>
struct ResultTraits1<UnaryFunctor<ArgumentFunctor3>, T1>
{
    typedef ErrorType Res;
};

template <class T1, class T2>
struct ResultTraits2<UnaryFunctor<ArgumentFunctor3>, T1, T2>
{
    typedef ErrorType Res;
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
    
    T const & operator()() const
    {
        return value_;
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
    
  protected:
    T value_;
  
  private:
    ParameterFunctor & operator=(ParameterFunctor const &); // not implemented
};

template <class T>
struct ResultTraits0<ParameterFunctor<T> >
{
    typedef T Res;
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
inline UnaryFunctor<ParameterFunctor<T> >
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
class UnaryAnalyser
{
  public:
    UnaryAnalyser(EXPR const & e)
    : expr_(e)
    {}
    
    void operator()() const
    {
        expr_();
    }
    
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
  
  private:
    UnaryAnalyser & operator=(UnaryAnalyser const &); // not implemented
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

#define MAKE_ASSIGNMENT_FUNCTOR(name, op) \
    template <class V, class EXPR> \
    struct AssignmentFunctor_##name \
    { \
        AssignmentFunctor_##name(UnaryFunctor<VarFunctor<V> > v,  \
                                 UnaryFunctor<EXPR> const & e) \
        : value_(v.value_), expr_(e) \
        {} \
         \
        V & operator()() const \
        { \
            const_cast<V &>(value_) op expr_(); \
            return const_cast<V &>(value_); \
        } \
         \
        template <class T1>  \
        V & operator()(T1 const & v1) const \
        { \
            const_cast<V &>(value_) op expr_(v1); \
            return const_cast<V &>(value_); \
        } \
         \
        template <class T1, class T2>  \
        V & operator()(T1 const & v1, T2 const & v2) const \
        { \
            const_cast<V &>(value_) op expr_(v1, v2); \
            return const_cast<V &>(value_); \
        } \
         \
        template <class T1, class T2, class T3>  \
        V & operator()(T1 const & v1, T2 const & v2, T3 const & v3) const \
        { \
            const_cast<V &>(value_) op expr_(v1, v2, v3); \
            return const_cast<V &>(value_); \
        } \
         \
      private: \
        V & value_; \
        UnaryFunctor<EXPR> expr_; \
        \
        AssignmentFunctor_##name & operator=(AssignmentFunctor_##name const &);\
    }; 

/************************************************************/

MAKE_ASSIGNMENT_FUNCTOR(assign, =)
MAKE_ASSIGNMENT_FUNCTOR(add, +=)
MAKE_ASSIGNMENT_FUNCTOR(subtract, -=)
MAKE_ASSIGNMENT_FUNCTOR(multiply, *=)
MAKE_ASSIGNMENT_FUNCTOR(divide, /=)

#undef MAKE_ASSIGNMENT_FUNCTOR

/************************************************************/
/*                                                          */
/*                          variables                       */
/*                                                          */
/************************************************************/

template <class T>
struct UnaryFunctor<VarFunctor<T> >
{
    typedef T Res;
    
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
    
    T const & operator()() const
    {
        return value_;
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
  
  private:
    UnaryFunctor & operator=(UnaryFunctor const &); // not implemented
};

template <class T>
struct ResultTraits0<UnaryFunctor<VarFunctor<T> > >
{
    typedef T Res;
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
inline UnaryFunctor<VarFunctor<T> >
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
    typedef void Res;
    
    IfThenFunctor(EXPR1 const & e1, EXPR2 const & e2)
    : expr1_(e1), expr2_(e2)
    {}
    
    void operator()() const 
    {
        if( expr1_() ) expr2_();
    }

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
  
  private:
    IfThenFunctor & operator=(IfThenFunctor const &); // not implemented
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

template <class EXPR1, class EXPR2, class EXPR3>
struct ResultTraits0<IfThenElseFunctor<EXPR1, EXPR2, EXPR3> >
{
    typedef typename ResultTraits0<EXPR2>::Res R2;
    typedef typename ResultTraits0<EXPR3>::Res R3;
    typedef typename PromoteTraits<R2, R3>::Promote Res;
};

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
    
    typename ResultTraits0<IfThenElseFunctor>::Res 
    operator()() const 
    {
        if(expr1_())
        {
            return typename ResultTraits0<IfThenElseFunctor>::Res(expr2_());
        }
        else
        {
            return typename ResultTraits0<IfThenElseFunctor>::Res(expr3_());
        }
    }

    template <class T> 
    typename ResultTraits1<IfThenElseFunctor, T>::Res 
    operator()(T const & v1) const 
    {
        if(expr1_(v1))
        {
            return typename ResultTraits1<IfThenElseFunctor, T>::Res(expr2_(v1));
        }
        else
        {
            return typename ResultTraits1<IfThenElseFunctor, T>::Res(expr3_(v1));
        }
    }

    template <class T1, class T2> 
    typename ResultTraits2<IfThenElseFunctor, T1, T2>::Res 
    operator()(T1 const & v1, T2 const & v2) const 
    {
        if(expr1_(v1, v2))
        {
            return typename ResultTraits2<IfThenElseFunctor, T1, T2>::Res(expr2_(v1, v2));
        }
        else
        {
            return typename ResultTraits2<IfThenElseFunctor, T1, T2>::Res(expr3_(v1, v2));
        }
    }

    template <class T1, class T2, class T3> 
    typename ResultTraits3<IfThenElseFunctor, T1, T2, T3>::Res 
    operator()(T1 const & v1, T2 const & v2, T3 const & v3) const 
    {
        if(expr1_(v1, v2, v3))
        {
            return typename ResultTraits3<IfThenElseFunctor, T1, T2, T3>::Res(expr2_(v1, v2, v3));
        }
        else
        {
            return typename ResultTraits3<IfThenElseFunctor, T1, T2, T3>::Res(expr3_(v1, v2, v3));
        }
    }
    
  private:
  
    EXPR1 expr1_;
    EXPR2 expr2_;
    EXPR3 expr3_;
  
    IfThenElseFunctor & operator=(IfThenElseFunctor const &); // not implemented
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

#define MAKE_FUNCTOR_UNARY_FUNCTION(function, namespc, traitsClass, traitsValue) \
    using ::namespc::function; \
    template <class EXPR> \
    struct Functor_##function; \
    \
    template <class EXPR> \
    struct ResultTraits0<Functor_##function<EXPR> > \
    { \
        typedef typename ResultTraits0<EXPR>::Res R1; \
        typedef typename traitsClass<R1>::traitsValue Res; \
    }; \
    \
    template <class EXPR, class T1> \
    struct ResultTraits1<Functor_##function<EXPR>, T1> \
    { \
        typedef typename ResultTraits1<EXPR, T1>::Res R1; \
        typedef typename traitsClass<R1>::traitsValue Res; \
    }; \
    \
    template <class EXPR, class T1, class T2> \
    struct ResultTraits2<Functor_##function<EXPR>, T1, T2> \
    { \
        typedef typename ResultTraits2<EXPR, T1, T2>::Res R1; \
        typedef typename traitsClass<R1>::traitsValue Res; \
    }; \
    \
    template <class EXPR, class T1, class T2, class T3> \
    struct ResultTraits3<Functor_##function<EXPR>, T1, T2, T3> \
    { \
        typedef typename ResultTraits3<EXPR, T1, T2, T3>::Res R1; \
        typedef typename traitsClass<R1>::traitsValue Res; \
    }; \
    \
    template <class EXPR> \
    struct Functor_##function \
    { \
        Functor_##function(EXPR const & e) \
        : expr_(e) \
        {} \
         \
        typename ResultTraits0<Functor_##function>::Res \
        operator()() const \
        { \
            return function(expr_()); \
        } \
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
       \
      private: \
        Functor_##function & operator=(Functor_##function const &); \
    }; \
     \
    template <class EXPR> \
    inline UnaryFunctor<Functor_##function<UnaryFunctor<EXPR> > > \
    function(UnaryFunctor<EXPR> const & e) \
    { \
        Functor_##function<UnaryFunctor<EXPR> > p(e); \
        return UnaryFunctor<Functor_##function<UnaryFunctor<EXPR> > >(p); \
    }

/************************************************************/

MAKE_FUNCTOR_UNARY_FUNCTION(sq, vigra, NumericTraits, RealPromote)
MAKE_FUNCTOR_UNARY_FUNCTION(sqrt, std, NumericTraits, RealPromote)
MAKE_FUNCTOR_UNARY_FUNCTION(exp, std, NumericTraits, RealPromote)
MAKE_FUNCTOR_UNARY_FUNCTION(log, std, NumericTraits, RealPromote)
MAKE_FUNCTOR_UNARY_FUNCTION(log10, std, NumericTraits, RealPromote)
MAKE_FUNCTOR_UNARY_FUNCTION(sin, std, NumericTraits, RealPromote)
MAKE_FUNCTOR_UNARY_FUNCTION(asin, std, NumericTraits, RealPromote)
MAKE_FUNCTOR_UNARY_FUNCTION(cos, std, NumericTraits, RealPromote)
MAKE_FUNCTOR_UNARY_FUNCTION(acos, std, NumericTraits, RealPromote)
MAKE_FUNCTOR_UNARY_FUNCTION(tan, std, NumericTraits, RealPromote)
MAKE_FUNCTOR_UNARY_FUNCTION(atan, std, NumericTraits, RealPromote)
MAKE_FUNCTOR_UNARY_FUNCTION(floor, std, NumericTraits, RealPromote)
MAKE_FUNCTOR_UNARY_FUNCTION(ceil, std, NumericTraits, RealPromote)
MAKE_FUNCTOR_UNARY_FUNCTION(abs, vigra, NumericTraits, RealPromote)
MAKE_FUNCTOR_UNARY_FUNCTION(norm, vigra, NormTraits, NormType)
MAKE_FUNCTOR_UNARY_FUNCTION(squaredNorm, vigra, NormTraits, SquaredNormType)

#undef MAKE_FUNCTOR_UNARY_FUNCTION

/************************************************************/
/*                                                          */
/*                functors for unary operators              */
/*                                                          */
/************************************************************/

#define MAKE_FUNCTOR_UNARY_OPERATOR(name, op) \
    template <class EXPR> \
    struct Functor_##name; \
    \
    template <class EXPR> \
    struct ResultTraits0<Functor_##name<EXPR> > \
    { \
        typedef typename ResultTraits0<EXPR>::Res Res; \
    }; \
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
        typename ResultTraits0<Functor_##name>::Res \
        operator()() const \
        { \
            return op expr_(); \
        } \
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
       \
      private: \
        Functor_##name & operator=(Functor_##name const &);\
    }; \
     \
    template <class EXPR> \
    inline UnaryFunctor<Functor_##name<UnaryFunctor<EXPR> > > \
    operator op(UnaryFunctor<EXPR> const & e) \
    { \
        Functor_##name<UnaryFunctor<EXPR> > p(e); \
        return UnaryFunctor<Functor_##name<UnaryFunctor<EXPR> > >(p); \
    }


/************************************************************/

MAKE_FUNCTOR_UNARY_OPERATOR(minus, -)
MAKE_FUNCTOR_UNARY_OPERATOR(negate, !)
MAKE_FUNCTOR_UNARY_OPERATOR(bitNegate, ~)

#undef MAKE_FUNCTOR_UNARY_OPERATOR

/************************************************************/
/*                                                          */
/*               functors for binary functions              */
/*                                                          */
/************************************************************/

#define MAKE_FUNCTOR_BINARY_FUNCTION(function) \
    using std::function; \
    template <class EXPR1, class EXPR2> \
    struct Functor_##function; \
    \
    template <class EXPR1, class EXPR2> \
    struct ResultTraits0<Functor_##function<EXPR1, EXPR2> > \
    { \
        typedef typename ResultTraits0<EXPR1>::Res R1; \
        typedef typename ResultTraits0<EXPR2>::Res R2; \
        typedef typename PromoteTraits<R1, R2>::Promote R3; \
        typedef typename NumericTraits<R3>::RealPromote Res; \
    }; \
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
        typename ResultTraits0<Functor_##function>::Res \
        operator()() const \
        { \
            return function(expr1_(), expr2_()); \
        } \
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
        \
        Functor_##function & operator=(Functor_##function const &); \
    }; \
     \
    template <class EXPR1, class EXPR2> \
    inline UnaryFunctor<Functor_##function<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > > \
    function(UnaryFunctor<EXPR1> const & e1, UnaryFunctor<EXPR2> const & e2) \
    { \
        Functor_##function<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > p(e1, e2); \
        return UnaryFunctor<Functor_##function<UnaryFunctor<EXPR1>,  \
                                        UnaryFunctor<EXPR2> > >(p); \
    }

/************************************************************/

MAKE_FUNCTOR_BINARY_FUNCTION(pow)
MAKE_FUNCTOR_BINARY_FUNCTION(atan2)
MAKE_FUNCTOR_BINARY_FUNCTION(fmod)

#undef MAKE_FUNCTOR_BINARY_FUNCTION

/************************************************************/

#define MAKE_FUNCTOR_MINMAX(name, op) \
    template <class EXPR1, class EXPR2> \
    struct Functor_##name; \
    \
    template <class EXPR1, class EXPR2> \
    struct ResultTraits0<Functor_##name<EXPR1, EXPR2> > \
    { \
        typedef typename ResultTraits0<EXPR1>::Res R1; \
        typedef typename ResultTraits0<EXPR2>::Res R2; \
        typedef typename PromoteTraits<R1, R2>::Promote Res; \
    }; \
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
        typename ResultTraits0<Functor_##name>::Res \
        operator()() const \
        { \
            typename \
            ResultTraits0<Functor_##name<EXPR1, EXPR2> >::R1 r1(expr1_()); \
            typename \
            ResultTraits0<Functor_##name<EXPR1, EXPR2> >::R2 r2(expr2_()); \
            return (r1 op r2) ? r1 : r2; \
        } \
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
        \
        Functor_##name & operator=(Functor_##name const &); \
    }; \
     \
    template <class EXPR1, class EXPR2> \
    inline UnaryFunctor<Functor_##name<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > > \
    name(UnaryFunctor<EXPR1> const & e1, UnaryFunctor<EXPR2> const & e2) \
    { \
        Functor_##name<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > p(e1, e2); \
        return UnaryFunctor<Functor_##name<UnaryFunctor<EXPR1>,  \
                                        UnaryFunctor<EXPR2> > >(p); \
    }

MAKE_FUNCTOR_MINMAX(min, <)
MAKE_FUNCTOR_MINMAX(max, >)

#undef MAKE_FUNCTOR_MINMAX

/************************************************************/
/*                                                          */
/*               functors for binary operators              */
/*                                                          */
/************************************************************/

#define MAKE_FUNCTOR_BINARY_OPERATOR(name, op) \
    template <class EXPR1, class EXPR2> \
    struct Functor_##name; \
    \
    template <class EXPR1, class EXPR2> \
    struct ResultTraits0<Functor_##name<EXPR1, EXPR2> > \
    { \
        typedef typename ResultTraits0<EXPR1>::Res R1; \
        typedef typename ResultTraits0<EXPR2>::Res R2; \
        typedef typename PromoteTraits<R1, R2>::Promote Res; \
    }; \
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
        typename ResultTraits0<Functor_##name>::Res \
        operator()() const \
        { \
            return expr1_() op expr2_(); \
        } \
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
        \
        Functor_##name & operator=(Functor_##name const &); \
    }; \
     \
    template <class EXPR1, class EXPR2> \
    inline UnaryFunctor<Functor_##name<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > > \
    operator op(UnaryFunctor<EXPR1> const & e1, UnaryFunctor<EXPR2> const & e2) \
    { \
        Functor_##name<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > p(e1, e2); \
        return UnaryFunctor<Functor_##name<UnaryFunctor<EXPR1>,  \
                                        UnaryFunctor<EXPR2> > >(p); \
    }

/************************************************************/

MAKE_FUNCTOR_BINARY_OPERATOR(add, +)
MAKE_FUNCTOR_BINARY_OPERATOR(subtract, -)
MAKE_FUNCTOR_BINARY_OPERATOR(multiply, *)
MAKE_FUNCTOR_BINARY_OPERATOR(divide, /)
MAKE_FUNCTOR_BINARY_OPERATOR(modulo, %)
MAKE_FUNCTOR_BINARY_OPERATOR(bitAnd, &)
MAKE_FUNCTOR_BINARY_OPERATOR(bitOr, |)
MAKE_FUNCTOR_BINARY_OPERATOR(bitXor, ^)

#undef MAKE_FUNCTOR_BINARY_OPERATOR

/************************************************************/

#define MAKE_FUNCTOR_BINARY_OPERATOR_BOOL(name, op) \
    template <class EXPR1, class EXPR2> \
    struct Functor_##name; \
    \
    template <class EXPR1, class EXPR2> \
    struct ResultTraits0<Functor_##name<EXPR1, EXPR2> > \
    { \
        typedef bool Res; \
    }; \
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
        bool operator()() const \
        { \
            return expr1_() op expr2_(); \
        } \
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
        \
        Functor_##name & operator=(Functor_##name const &); \
    }; \
     \
    template <class EXPR1, class EXPR2> \
    inline UnaryFunctor<Functor_##name<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > > \
    operator op(UnaryFunctor<EXPR1> const & e1, UnaryFunctor<EXPR2> const & e2) \
    { \
        Functor_##name<UnaryFunctor<EXPR1>, UnaryFunctor<EXPR2> > p(e1, e2); \
        return UnaryFunctor<Functor_##name<UnaryFunctor<EXPR1>,  \
                                        UnaryFunctor<EXPR2> > >(p); \
    }

/************************************************************/

MAKE_FUNCTOR_BINARY_OPERATOR_BOOL(equals, ==)
MAKE_FUNCTOR_BINARY_OPERATOR_BOOL(differs, !=)
MAKE_FUNCTOR_BINARY_OPERATOR_BOOL(less, <)
MAKE_FUNCTOR_BINARY_OPERATOR_BOOL(lessEqual, <=)
MAKE_FUNCTOR_BINARY_OPERATOR_BOOL(greater, >)
MAKE_FUNCTOR_BINARY_OPERATOR_BOOL(greaterEqual, >=)
MAKE_FUNCTOR_BINARY_OPERATOR_BOOL(and, &&)
MAKE_FUNCTOR_BINARY_OPERATOR_BOOL(or, ||)

#undef MAKE_FUNCTOR_BINARY_OPERATOR_BOOL

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
    
    RES operator()() const 
    {
        return f_(expr_());
    }
    
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
  
  private:
    UnaryFctPtrFunctor & operator=(UnaryFctPtrFunctor const &); // not implemented
};

template <class EXPR, class RES, class ARG>
struct ResultTraits0<UnaryFctPtrFunctor<EXPR, RES, ARG> >
{
    typedef RES Res;
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
inline UnaryFunctor<UnaryFctPtrFunctor<UnaryFunctor<EXPR>, RES, ARG> >
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
    
    RES operator()() const 
    {
        return f_(expr1_(), expr2_());
    }
    
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
  
  private:
    BinaryFctPtrFunctor & operator=(BinaryFctPtrFunctor const &); // not implemented
};

template <class EXPR1, class EXPR2, class RES, class ARG1, class ARG2>
struct ResultTraits0<BinaryFctPtrFunctor<EXPR1, EXPR2, RES, ARG1, ARG2> >
{
    typedef RES Res;
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
inline UnaryFunctor<BinaryFctPtrFunctor<UnaryFunctor<EXPR1>, 
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
    
    typename ResultTraits0<EXPR2>::Res 
    operator()() const 
    {
        expr1_();
        return expr2_();
    }
    
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
  
  private:
    CommaFunctor & operator=(CommaFunctor const &); // not implemented
};

template <class Expr1, class Expr2>
struct ResultTraits0<CommaFunctor<Expr1, Expr2> >
{
    typedef typename ResultTraits0<Expr2>::Res Res;
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
inline UnaryFunctor<CommaFunctor<UnaryAnalyser<EXPR1>, 
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
    
    void operator()() const 
    {
        expr1_();
        expr2_();
    }
    
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
  
  private:
    CommaAnalyser & operator=(CommaAnalyser const &); // not implemented
};

template <class EXPR1, class EXPR2>
inline UnaryAnalyser<CommaAnalyser<UnaryAnalyser<EXPR1>, 
                            UnaryAnalyser<EXPR2> > >
operator,(UnaryAnalyser<EXPR1> const & e1, 
          UnaryAnalyser<EXPR2> const & e2)
{
    CommaAnalyser<UnaryAnalyser<EXPR1>, 
                            UnaryAnalyser<EXPR2> >  p(e1, e2);
    return UnaryAnalyser<CommaAnalyser<UnaryAnalyser<EXPR1>, 
                            UnaryAnalyser<EXPR2> > >(p);
}

} // namespace functor

#if defined(__GNUC__) &&  __GNUC__ < 3
using functor::Arg1;
using functor::Arg2;
using functor::Arg3;
using functor::Param;
#endif

template <class T>
class FunctorTraits<functor::UnaryFunctor<T> >
: public FunctorTraitsBase<functor::UnaryFunctor<T> >
{
  public:
    typedef VigraTrueType isInitializer;
    typedef VigraTrueType isUnaryFunctor;
    typedef VigraTrueType isBinaryFunctor;
    typedef VigraTrueType isTernaryFunctor;
};

template <class T>
class FunctorTraits<functor::UnaryAnalyser<T> >
: public FunctorTraitsBase<functor::UnaryAnalyser<T> >
{
  public:
    typedef VigraTrueType isUnaryAnalyser;
    typedef VigraTrueType isBinaryAnalyser;
    typedef VigraTrueType isTernaryAnalyser;
};



} // namespace vigra

#endif /* NO_PARTIAL_TEMPLATE_SPECIALIZATION */

#endif // DOXYGEN

#endif /* VIGRA_FUNCTOREXPRESSION_HXX  */


