/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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

#ifndef VIGRA_TUPLE_HXX
#define VIGRA_TUPLE_HXX

#include <utility>    // for pair

namespace vigra {

/*! \page TupleTypes Tuple Types

    pair, triple, tuple4, tuple5

    <b>\#include</b> "<a href="utilities_8hxx-source.html">vigra/utilities.hxx</a>"<br>
    Namespace: vigra

    VIGRA defines tuple types \p vigra::triple, \p vigra::tuple4, \p vigra::tuple5.
    In addition, \p std::pair is imported into namespace vigra from the C++ standard
    library. All these types are defined similarly:

    <ul>

    <li> They are parameterized by the respective number of types. For each tuple,
    a constructor is defined that takes that many arguments, e.g.:
    \code
    template <class First, class Second, class Third>
    class Triple { ... };
    \endcode
    </li>
    <li> A number of \p typedef's tells the types stored in the tuple:

    \code
    typedef ... first_type;
    typedef ... second_type;
    typedef ... third_type;  // triple, tuple4, tuple5 only
    typedef ... forth_type;  // tuple4, tuple5 only
    typedef ... fifth_type;  // tuple5 only
    \endcode
    </li>
    <li> Items are stored in the following public attributes:

    \code

    first;
    second;
    third;  // triple, tuple4, tuple5 only
    forth;  // tuple4, tuple5 only
    fifth;  // tuple5 only

    \endcode
    </li>
    </ul>


*/

/********************************************************/
/*                                                      */
/*                          pair                        */
/*                                                      */
/********************************************************/

using std::pair;

/********************************************************/
/*                                                      */
/*                          triple                      */
/*                                                      */
/********************************************************/

template <class T1, class T2, class T3>
struct triple {
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;

    T1 first;
    T2 second;
    T3 third;
    triple() {}
    triple(const T1& a, const T2& b, const T3& c)
    : first(a), second(b), third(c) {}
};

template <class T1, class T2, class T3>
triple<T1,T2,T3> make_triple( T1 t1, T2 t2, T3 t3 )
{ return triple<T1,T2,T3>( t1, t2, t3 ); }

/********************************************************/
/*                                                      */
/*                          tuple4                      */
/*                                                      */
/********************************************************/

template <class T1, class T2, class T3, class T4>
struct tuple4 {
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;
    typedef T4 fourth_type;

    T1 first;
    T2 second;
    T3 third;
    T4 fourth;
    tuple4() {}
    tuple4(const T1& a, const T2& b, const T3& c, const T4& d)
    : first(a), second(b), third(c), fourth(d) {}
};

template <class T1, class T2, class T3, class T4>
tuple4<T1,T2,T3,T4> make_tuple4( T1 t1, T2 t2, T3 t3, T4 t4 )
{ return tuple4<T1,T2,T3,T4>( t1, t2, t3, t4 ); }

/********************************************************/
/*                                                      */
/*                          tuple5                      */
/*                                                      */
/********************************************************/

template <class T1, class T2, class T3, class T4, class T5>
struct tuple5 {
    typedef T1 first_type;
    typedef T2 second_type;
    typedef T3 third_type;
    typedef T4 fourth_type;
    typedef T5 fifth_type;

    T1 first;
    T2 second;
    T3 third;
    T4 fourth;
    T5 fifth;
    tuple5() {}
    tuple5(const T1& a, const T2& b, const T3& c, const T4& d, const T5& e)
    : first(a), second(b), third(c), fourth(d), fifth(e) {}
};

template <class T1, class T2, class T3, class T4, class T5>
tuple5<T1,T2,T3,T4,T5> make_tuple5( T1 t1, T2 t2, T3 t3, T4 t4, T5 t5 )
{ return tuple5<T1,T2,T3,T4,T5>( t1, t2, t3, t4, t5 ); }


} // namespace vigra



#endif /* VIGRA_TUPLE_HXX */
