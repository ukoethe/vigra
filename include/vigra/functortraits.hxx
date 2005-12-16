/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2005 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de          or                  */
/*        vigra@kogs1.informatik.uni-hamburg.de                         */
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


#ifndef VIGRA_FUNCTORTRAITS_HXX
#define VIGRA_FUNCTORTRAITS_HXX

#include <functional>
#include <vigra/metaprogramming.hxx>

namespace vigra {

template <class T>
class FunctorTraitsBase
{
  public:
    typedef T type;
    
    typedef VigraFalseType isInitializer;
    
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};



/** \addtogroup Functors
*/
//@{
/** \brief Export associated information for a functor.

    The FunctorTraits class contains the following fields:

    \code
    template <class T>
    struct FunctorTraits
    {
        typedef T type;
        
        typedef ... isInitializer;

        typedef ... isUnaryFunctor;
        typedef ... isBinaryFunctor;
        typedef ... isTernaryFunctor;
        
        typedef ... isUnaryAnalyser;
        typedef ... isBinaryAnalyser;
        typedef ... isTernaryAnalyser;
    };
    \endcode

    Where the dots are either <tt>VigraTrueType</tt> or <tt>VigraFalseType</tt>
    depending on whether the functor supports the respective functionality or not.
    If a functor <tt>f<tt> is a model of these categories, it supports the following
    calls (<tt>v</tt> is a variable such that the result type of the functor
    calls can be converted into <tt>v</tt>'s type, and <tt>a1, a2, a3</tt> are
    variables convertible into the functor's argument types):
    
    <DL>
    <DT><b>Initializer</b>
        <DD> <tt>v = f()</tt> (used with initImageWithFunctor())
    <DT><b>UnaryFunctor</b>
        <DD> <tt>v = f(a1)</tt> (used with transformImage())
    <DT><b>BinaryFunctor</b>
        <DD> <tt>v = f(a1, a2)</tt> (used with combineTwoImages())
    <DT><b>TernaryFunctor</b>
        <DD> <tt>v = f(a1, a2, a3)</tt> (used with combineThreeImages())
    <DT><b>UnaryAnalyser</b>
        <DD> <tt>f(a1)</tt> (return type <tt>void>/tt>, used with inspectImage())
    <DT><b>BinaryAnalyser</b>
        <DD> <tt>f(a1, a2)</tt> (return type <tt>void>/tt>, used with inspectTwoImages())
    <DT><b>TernaryAnalyser</b>
        <DD> <tt>f(a1, a2, a3)</tt> (return type <tt>void>/tt>)
    </DL>
    
    It should be noted that the functor's argument and result types are not contained
    in the traits class: Since the function calls are often member template functions in 
    VIGRA, many functors do not have fixed argument types. Neither are the result
    types fixed in this case because they are computed (via a template meta-program)
    from the argument types.

    <b>\#include</b> "<a href="functortraits_8hxx-source.html">vigra/functortraits.hxx</a>"
    Namespace: vigra
*/
template <class T>
class FunctorTraits
: public FunctorTraitsBase<T>
{};

#define VIGRA_DEFINE_STL_FUNCTOR(name, unary, binary) \
template <class T> \
class FunctorTraits<name<T> > \
{ \
  public: \
    typedef T type; \
     \
    typedef VigraFalseType isInitializer; \
     \
    typedef unary          isUnaryFunctor; \
    typedef binary         isBinaryFunctor; \
    typedef VigraFalseType isTernaryFunctor; \
     \
    typedef VigraFalseType isUnaryAnalyser; \
    typedef VigraFalseType isBinaryAnalyser; \
    typedef VigraFalseType isTernaryAnalyser; \
};

// ???TODO: these should also be specialized for the ptr_fun and mem_fun_ptr wrappers
VIGRA_DEFINE_STL_FUNCTOR(std::plus, VigraFalseType, VigraTrueType)
VIGRA_DEFINE_STL_FUNCTOR(std::minus, VigraFalseType, VigraTrueType)
VIGRA_DEFINE_STL_FUNCTOR(std::multiplies, VigraFalseType, VigraTrueType)
VIGRA_DEFINE_STL_FUNCTOR(std::divides, VigraFalseType, VigraTrueType)
VIGRA_DEFINE_STL_FUNCTOR(std::modulus, VigraFalseType, VigraTrueType)
VIGRA_DEFINE_STL_FUNCTOR(std::equal_to, VigraFalseType, VigraTrueType)
VIGRA_DEFINE_STL_FUNCTOR(std::not_equal_to, VigraFalseType, VigraTrueType)
VIGRA_DEFINE_STL_FUNCTOR(std::greater, VigraFalseType, VigraTrueType)
VIGRA_DEFINE_STL_FUNCTOR(std::less, VigraFalseType, VigraTrueType)
VIGRA_DEFINE_STL_FUNCTOR(std::greater_equal, VigraFalseType, VigraTrueType)
VIGRA_DEFINE_STL_FUNCTOR(std::less_equal, VigraFalseType, VigraTrueType)
VIGRA_DEFINE_STL_FUNCTOR(std::logical_and, VigraFalseType, VigraTrueType)
VIGRA_DEFINE_STL_FUNCTOR(std::logical_or, VigraFalseType, VigraTrueType)
VIGRA_DEFINE_STL_FUNCTOR(std::binary_negate, VigraFalseType, VigraTrueType)

VIGRA_DEFINE_STL_FUNCTOR(std::negate, VigraTrueType, VigraFalseType)
VIGRA_DEFINE_STL_FUNCTOR(std::logical_not, VigraTrueType, VigraFalseType)
VIGRA_DEFINE_STL_FUNCTOR(std::unary_negate, VigraTrueType, VigraFalseType)
VIGRA_DEFINE_STL_FUNCTOR(std::binder1st, VigraTrueType, VigraFalseType)
VIGRA_DEFINE_STL_FUNCTOR(std::binder2nd, VigraTrueType, VigraFalseType)
#undef VIGRA_DEFINE_STL_FUNCTOR

template <class R>
class FunctorTraits<R (*)()>
{
  public:
    typedef R (*type)();
    
    typedef VigraTrueType  isInitializer;
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};

template <class R, class T>
class FunctorTraits<R (*)(T)>
{
  public:
    typedef R (*type)(T);
    
    typedef VigraFalseType isInitializer;
    typedef VigraTrueType  isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};

template <class R, class T1, class T2>
class FunctorTraits<R (*)(T1, T2)>
{
  public:
    typedef R (*type)(T1, T2);
    
    typedef VigraFalseType isInitializer;
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraTrueType  isBinaryFunctor;
    typedef VigraFalseType isTernaryFunctor;
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};

template <class R, class T1, class T2, class T3>
class FunctorTraits<R (*)(T1, T2, T3)>
{
  public:
    typedef R (*type)(T1, T2, T3);
    
    typedef VigraFalseType isInitializer;
    typedef VigraFalseType isUnaryFunctor;
    typedef VigraFalseType isBinaryFunctor;
    typedef VigraTrueType  isTernaryFunctor;
    typedef VigraFalseType isUnaryAnalyser;
    typedef VigraFalseType isBinaryAnalyser;
    typedef VigraFalseType isTernaryAnalyser;
};

//@}

} // namespace vigra

#endif // VIGRA_FUNCTORTRAITS_HXX
