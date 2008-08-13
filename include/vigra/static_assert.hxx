/************************************************************************/
/*                                                                      */
/*               Copyright 2004-2005 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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

#ifndef VIGRA_STATIC_ASSERT_HXX
#define VIGRA_STATIC_ASSERT_HXX

// based on the static assertion design in boost::mpl (see www.boost.org)

#define VIGRA_PREPROCESSOR_CONCATENATE(a, b) VIGRA_PREPROCESSOR_CONCATENATE_IMPL(a, b)
#define VIGRA_PREPROCESSOR_CONCATENATE_IMPL(a, b) a ## b

namespace vigra {

namespace staticAssert {

template <bool Predicate>
struct AssertBool;

template <>
struct AssertBool<true>
{
    typedef int type;
    typedef void * not_type;
};

template <>
struct AssertBool<false>
{
    typedef void * type;
    typedef int not_type;
};

template <class T>
struct Assert;

template <>
struct Assert<VigraTrueType>
{
    typedef int type;
    typedef void * not_type;
};

template <>
struct Assert<VigraFalseType>
{
    typedef void * type;
    typedef int not_type;
};

struct failure{};
struct success {};
inline int check( success ) { return 0; }

template< typename Predicate >
failure ************ (Predicate::************ 
      assertImpl( void (*)(Predicate), typename Predicate::not_type )
    );

template< typename Predicate >
success
assertImpl( void (*)(Predicate), typename Predicate::type );

/* Usage:

1. Define an assertion class, derived from vigra::staticAssert::Assert,
   whose name serves as an error message:
   
    template <int N>
    struct FixedPoint_overflow_error__More_than_31_bits_requested
    : vigra::staticAssert::AssertBool<(N < 32)>
    {};

2. Call VIGRA_STATIC_ASSERT() with the assertion class:

    template <int N>
    void test()
    {
        // signal error if N > 31
        VIGRA_STATIC_ASSERT((FixedPoint_overflow_error__More_than_31_bits_requested<N>));
    }
    
TODO: provide more assertion base classes for other (non boolean) types of tests
*/
#if !defined(__GNUC__) || __GNUC__ > 2
#define VIGRA_STATIC_ASSERT(Predicate) \
enum { \
    VIGRA_PREPROCESSOR_CONCATENATE(vigra_assertion_in_line_, __LINE__) = sizeof( \
         staticAssert::check( \
              staticAssert::assertImpl( (void (*) Predicate)0, 1 ) \
            ) \
        ) \
}
#else
#define VIGRA_STATIC_ASSERT(Predicate)
#endif

} // namespace staticAssert

} // namespace vigra

#endif // VIGRA_STATIC_ASSERT_HXX
