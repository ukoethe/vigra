#ifndef VIGRA_TYPE_LISTS_HXX
#define VIGRA_TYPE_LISTS_HXX

#include <iostream>
#include <typeinfo>
#include <utility>
#include <algorithm>
#include "metaprogramming.hxx"
#include "bit_array.hxx"
#include "error.hxx"

namespace vigra {

// mask cl.exe shortcomings [begin]
#if defined(_MSC_VER)
#pragma warning( push )
#pragma warning( disable : 4503 )
#endif

namespace type_lists {

    struct nil; // end-of-list marker.
    template <class T>
    struct nil_t;
    
    // type lists of size >= 1.

    template <class A, class B = nil> struct cons
    {
        typedef A first;
        typedef B rest;
    };

    template <class X, class A, class B> struct if_nil
    {
        typedef B type;
    };
    template <class A, class B> struct if_nil <nil, A, B>
    {
        typedef A type;
    };

    // truncate type list L (using class NIL as ending marker)
    // at the first occurence of type X
    template <class X, class L, class NIL = nil> struct truncate
    {
        typedef cons<typename L::first,
                     typename truncate<X, typename L::rest, NIL>::type> type;
    };
    template <class L, class NIL> struct truncate<typename L::first, L, NIL>
    {
        typedef nil type; // do the actual truncation
    };
    template <class X, class NIL> struct truncate<X, NIL, NIL>
    {
        typedef nil type;
    };

    template <class NIL, class A = NIL, class B = NIL, class C = NIL,
                         class D = NIL, class E = NIL, class F = NIL,
                         class G = NIL, class H = NIL, class I = NIL,
                         class J = NIL, class K = NIL, class L = NIL,
                         class M = NIL, class N = NIL, class O = NIL,
                         class P = NIL, class Q = NIL, class R = NIL,
                         class S = NIL, class T = NIL, class U = NIL,
                         class V = NIL, class W = NIL, class X = NIL,
                         class Y = NIL, class Z = NIL>
    struct make_list_nil {
        typedef typename truncate<NIL, cons<A, cons<B, cons<C, cons<D, cons<E,
                                       cons<F, cons<G, cons<H, cons<I, cons<J,
                                       cons<K, cons<L, cons<M, cons<N, cons<O,
                                       cons<P, cons<Q, cons<R, cons<S, cons<T,
                                       cons<U, cons<V, cons<W, cons<X, cons<Y,
                                       cons<Z, NIL> > > > > > > > > > > > > > >
                                       > > > > > > > > > > >, NIL>::type type;
    };

    template <class A = nil, class B = nil, class C = nil, class D = nil,
              class E = nil, class F = nil, class G = nil, class H = nil,
              class I = nil, class J = nil, class K = nil, class L = nil,
              class M = nil, class N = nil, class O = nil, class P = nil,
              class Q = nil, class R = nil, class S = nil, class T = nil,
              class U = nil, class V = nil, class W = nil, class X = nil,
              class Y = nil, class Z = nil>

    struct make_list {
        typedef typename make_list_nil<nil, A, B, C, D, E, F, G, H, I,
                                            J, K, L, M, N, O, P, Q, R,
                                            S, T, U, V, W, X, Y, Z
                                      >::type type;
    };

    template <class T_, template<class> class A = nil_t,
                        template<class> class B = nil_t,
                        template<class> class C = nil_t,
                        template<class> class D = nil_t,
                        template<class> class E = nil_t,
                        template<class> class F = nil_t,
                        template<class> class G = nil_t,
                        template<class> class H = nil_t,
                        template<class> class I = nil_t,
                        template<class> class J = nil_t,
                        template<class> class K = nil_t,
                        template<class> class L = nil_t,
                        template<class> class M = nil_t,
                        template<class> class N = nil_t,
                        template<class> class O = nil_t,
                        template<class> class P = nil_t,
                        template<class> class Q = nil_t,
                        template<class> class R = nil_t,
                        template<class> class S = nil_t,
                        template<class> class T = nil_t,
                        template<class> class U = nil_t,
                        template<class> class V = nil_t,
                        template<class> class W = nil_t,
                        template<class> class X = nil_t,
                        template<class> class Y = nil_t,
                        template<class> class Z = nil_t>
    struct make_list_template {
        typedef typename make_list_nil<nil_t<T_>,
                                       A<T_>, B<T_>, C<T_>, D<T_>, E<T_>,
                                       F<T_>, G<T_>, H<T_>, I<T_>, J<T_>,
                                       K<T_>, L<T_>, M<T_>, N<T_>, O<T_>,
                                       P<T_>, Q<T_>, R<T_>, S<T_>, T<T_>,
                                       U<T_>, V<T_>, W<T_>, X<T_>, Y<T_>,
                                       Z<T_> >::type type;
    };

    // a means to partially compensate for the lack of templated typedefs.
    template <template<class, class> class BASE, class T_,
                                                template<class> class A = nil_t,
                                                template<class> class B = nil_t,
                                                template<class> class C = nil_t,
                                                template<class> class D = nil_t,
                                                template<class> class E = nil_t,
                                                template<class> class F = nil_t,
                                                template<class> class G = nil_t,
                                                template<class> class H = nil_t,
                                                template<class> class I = nil_t,
                                                template<class> class J = nil_t,
                                                template<class> class K = nil_t,
                                                template<class> class L = nil_t,
                                                template<class> class M = nil_t,
                                                template<class> class N = nil_t,
                                                template<class> class O = nil_t,
                                                template<class> class P = nil_t,
                                                template<class> class Q = nil_t,
                                                template<class> class R = nil_t,
                                                template<class> class S = nil_t,
                                                template<class> class T = nil_t,
                                                template<class> class U = nil_t,
                                                template<class> class V = nil_t,
                                                template<class> class W = nil_t,
                                                template<class> class X = nil_t,
                                                template<class> class Y = nil_t,
                                                template<class> class Z = nil_t>
    struct use_template_list
        : public BASE<T_, typename make_list_template<T_, A, B, C, D, E, F, G,
                                                          H, I, J, K, L, M, N,
                                                          O, P, Q, R, S, T, U,
                                                          V, W, X, Y, Z>::type>
    {};

    // use first and rest only when possible:
    template <class T>
    struct has_first_rest : public sfinae_test<T, has_first_rest>
    {
	    template <class U>
        has_first_rest(U*, typename U::first* = 0, typename U::rest* = 0);
    };
    template <bool P, class A>
    struct cond_cons_rest;
    template <class A>
    struct cond_cons_rest<false, A>
    {
        typedef void* type;
    };
    template <class A>
    struct cond_cons_rest<true, A>
    {
        typedef typename A::rest type;
    };
    // test if a type is a list in the above sense.
    template <class A> struct is_list
    {
        static const bool value = is_list<typename
                      cond_cons_rest<has_first_rest<A>::value, A>::type>::value;
    };
    template <> struct is_list<nil>
    {
        static const bool value = true;
    };
    template <> struct is_list<void*>
    {
        static const bool value = false;
    };

    template <class A> struct list_guard
    {
        typedef typename IfBool<is_list<A>::value, A, nil>::type type;
    };

    template <class A> struct size
    {
        static const unsigned of = size<typename A::rest>::of + 1;
    };
    template <> struct size<nil>
    {
        static const unsigned of = 0;
    };

    template <class X, class L> struct append
    {
        typedef cons<typename L::first,
                     typename append<X, typename L::rest>::type> type;
    };
    template <class X> struct append<X, nil>
    {
        typedef cons<X, nil> type;
    };
    template <> struct append<nil, nil>
    {
        typedef nil type;
    };

    template <class L, class R = nil> struct reverse
    {
        typedef typename reverse<typename L::rest,
                                 cons<typename L::first, R> >::type type;
    };
    template <class R> struct reverse<nil, R>
    {
        typedef R type;
    };

    template <template<class> class P, class Q, class L>
    struct max_value
    {
        static const bool is_nil = false;
        static const Q first_value = P<typename L::first>::value;
        typedef max_value<P, Q, typename L::rest> rest_type;
        static const Q rest_value = rest_type::value;
        static const bool gt = first_value > rest_value || rest_type::is_nil;
        static const Q value = gt * first_value + !gt * rest_value;
    };
    template <template<class> class P, class Q>
    struct max_value<P, Q, nil>
    {
        static const Q value = 0;
        static const bool is_nil = true;
    };

    // remove the all occurences of type X in type list L
    template <class X, class L> struct remove // recursion
    {
        typedef cons<typename L::first,
                     typename remove<X, typename L::rest>::type> type;
    };
    template <class L> struct remove<typename L::first, L> // actual removal
    {
        typedef typename remove<typename L::first, typename L::rest>::type type;
    };
    template <class X> struct remove<X, nil> // list end
    {
        typedef nil type;
    };

    // remove the all occurences of type list L where predicate P equals value
    template <template<class> class P, class L, bool value = true>
    struct remove_if
    {
        typedef typename
            IfBool<
                value == P<typename L::first>::value, typename
                remove_if<P, typename L::rest, value>::type,
                cons<typename
                    L::first, typename
                    remove_if<P, typename L::rest, value>::type
                >
            >::type type;
    };
    template <template<class> class P, bool value>
    struct remove_if<P, nil, value>
    {
        typedef nil type;
    };

    template <template<class> class P, class L>
    struct remove_if_not
    {
        typedef typename remove_if<P, L, false>::type type;
    };

    template <class X, class L> struct contains
    {
        static const bool value = contains<X, typename L::rest>::value;
    };
    template <class L> struct contains<typename L::first, L>
    {
        static const bool value = true;
    };
    template <class X> struct contains<X, nil>
    {
        static const bool value = false;
    };

    // simple, unstable merge
    template <class X, class L> struct merge
    {
        typedef typename L::first first;
        typedef typename
            merge<
                    typename IfBool<contains<first, X>::value,
                               X,
                               cons<first, X>
                          >::type,
                    typename L::rest
                 >::type type;
    };
    template <class X> struct merge<X, nil>
    {
        typedef X type;
    };

    // simple, unstable unique
    template <class L> struct unique
    {
        typedef typename merge<nil, L>::type type;
    };

    template <class T_, template<class> class A = nil_t,
                        template<class> class B = nil_t,
                        template<class> class C = nil_t,
                        template<class> class D = nil_t,
                        template<class> class E = nil_t,
                        template<class> class F = nil_t,
                        template<class> class G = nil_t,
                        template<class> class H = nil_t,
                        template<class> class I = nil_t,
                        template<class> class J = nil_t,
                        template<class> class K = nil_t,
                        template<class> class L = nil_t,
                        template<class> class M = nil_t,
                        template<class> class N = nil_t,
                        template<class> class O = nil_t,
                        template<class> class P = nil_t,
                        template<class> class Q = nil_t,
                        template<class> class R = nil_t,
                        template<class> class S = nil_t,
                        template<class> class T = nil_t,
                        template<class> class U = nil_t,
                        template<class> class V = nil_t,
                        template<class> class W = nil_t,
                        template<class> class X = nil_t,
                        template<class> class Y = nil_t,
                        template<class> class Z = nil_t>
    struct implies_template
    {
        typedef typename make_list_template<T_, A, B, C, D, E, F, G, H, I, J, K,
                                                L, M, N, O, P, Q, R, S, T, U, V,
                                                W, X, Y, Z>::type implies_types;
    };

    template <class T_, template<class> class A = nil_t,
                        template<class> class B = nil_t,
                        template<class> class C = nil_t,
                        template<class> class D = nil_t,
                        template<class> class E = nil_t,
                        template<class> class F = nil_t,
                        template<class> class G = nil_t,
                        template<class> class H = nil_t,
                        template<class> class I = nil_t,
                        template<class> class J = nil_t,
                        template<class> class K = nil_t,
                        template<class> class L = nil_t,
                        template<class> class M = nil_t,
                        template<class> class N = nil_t,
                        template<class> class O = nil_t,
                        template<class> class P = nil_t,
                        template<class> class Q = nil_t,
                        template<class> class R = nil_t,
                        template<class> class S = nil_t,
                        template<class> class T = nil_t,
                        template<class> class U = nil_t,
                        template<class> class V = nil_t,
                        template<class> class W = nil_t,
                        template<class> class X = nil_t,
                        template<class> class Y = nil_t,
                        template<class> class Z = nil_t>
    struct follows_template
    {
        typedef typename make_list_template<T_, A, B, C, D, E, F, G, H, I, J, K,
                                                L, M, N, O, P, Q, R, S, T, U, V,
                                                W, X, Y, Z>::type follows_types;
    };

    template <class T_, template<class> class A = nil_t,
                        template<class> class B = nil_t,
                        template<class> class C = nil_t,
                        template<class> class D = nil_t,
                        template<class> class E = nil_t,
                        template<class> class F = nil_t,
                        template<class> class G = nil_t,
                        template<class> class H = nil_t,
                        template<class> class I = nil_t,
                        template<class> class J = nil_t,
                        template<class> class K = nil_t,
                        template<class> class L = nil_t,
                        template<class> class M = nil_t,
                        template<class> class N = nil_t,
                        template<class> class O = nil_t,
                        template<class> class P = nil_t,
                        template<class> class Q = nil_t,
                        template<class> class R = nil_t,
                        template<class> class S = nil_t,
                        template<class> class T = nil_t,
                        template<class> class U = nil_t,
                        template<class> class V = nil_t,
                        template<class> class W = nil_t,
                        template<class> class X = nil_t,
                        template<class> class Y = nil_t,
                        template<class> class Z = nil_t>
    struct depends_on_template
    {
        typedef typename make_list_template<T_, A, B, C, D, E, F, G, H, I, J, K,
                                                L, M, N, O, P, Q, R, S, T, U, V,
                                                W, X, Y, Z>::type depends_on;
    };

    template <class T_u, template<class> class A = nil_t,
                         template<class> class B = nil_t,
                         template<class> class C = nil_t,
                         template<class> class D = nil_t,
                         template<class> class E = nil_t,
                         template<class> class F = nil_t,
                         template<class> class G = nil_t,
                         template<class> class H = nil_t,
                         template<class> class I = nil_t,
                         template<class> class J = nil_t,
                         template<class> class K = nil_t,
                         template<class> class L = nil_t,
                         template<class> class M = nil_t,
                         template<class> class N = nil_t,
                         template<class> class O = nil_t,
                         template<class> class P = nil_t,
                         template<class> class Q = nil_t,
                         template<class> class R = nil_t,
                         template<class> class S = nil_t,
                         template<class> class T = nil_t,
                         template<class> class U = nil_t,
                         template<class> class V = nil_t,
                         template<class> class W = nil_t,
                         template<class> class X = nil_t,
                         template<class> class Y = nil_t,
                         template<class> class Z = nil_t>
    struct uses_template
        : public implies_template<T_u, A, B, C, D, E, F, G, H, I, J, K, L, M,
                                       N, O, P, Q, R, S, T, U, V, W, X, Y, Z>,
          public depends_on_template<T_u, A, B, C, D, E, F, G, H, I, J, K, L, M,
                                       N, O, P, Q, R, S, T, U, V, W, X, Y, Z>
    {
        template <template<class> class A_ = nil_t,
                  template<class> class B_ = nil_t,
                  template<class> class C_ = nil_t,
                  template<class> class D_ = nil_t,
                  template<class> class E_ = nil_t,
                  template<class> class F_ = nil_t,
                  template<class> class G_ = nil_t,
                  template<class> class H_ = nil_t,
                  template<class> class I_ = nil_t,
                  template<class> class J_ = nil_t,
                  template<class> class K_ = nil_t,
                  template<class> class L_ = nil_t,
                  template<class> class M_ = nil_t,
                  template<class> class N_ = nil_t,
                  template<class> class O_ = nil_t,
                  template<class> class P_ = nil_t,
                  template<class> class Q_ = nil_t,
                  template<class> class R_ = nil_t,
                  template<class> class S_ = nil_t,
                  template<class> class T_ = nil_t,
                  template<class> class U_ = nil_t,
                  template<class> class V_ = nil_t,
                  template<class> class W_ = nil_t,
                  template<class> class X_ = nil_t,
                  template<class> class Y_ = nil_t,
                  template<class> class Z_ = nil_t>
        struct follows
            : public uses_template
            , public follows_template<T_u, A_, B_, C_, D_, E_, F_, G_, H_, I_,
                                           J_, K_, L_, M_, N_, O_, P_, Q_, R_,
                                           S_, T_, U_, V_, W_, X_, Y_, Z_> {};

        template <template<class> class A_ = nil_t,
                  template<class> class B_ = nil_t,
                  template<class> class C_ = nil_t,
                  template<class> class D_ = nil_t,
                  template<class> class E_ = nil_t,
                  template<class> class F_ = nil_t,
                  template<class> class G_ = nil_t,
                  template<class> class H_ = nil_t,
                  template<class> class I_ = nil_t,
                  template<class> class J_ = nil_t,
                  template<class> class K_ = nil_t,
                  template<class> class L_ = nil_t,
                  template<class> class M_ = nil_t,
                  template<class> class N_ = nil_t,
                  template<class> class O_ = nil_t,
                  template<class> class P_ = nil_t,
                  template<class> class Q_ = nil_t,
                  template<class> class R_ = nil_t,
                  template<class> class S_ = nil_t,
                  template<class> class T_ = nil_t,
                  template<class> class U_ = nil_t,
                  template<class> class V_ = nil_t,
                  template<class> class W_ = nil_t,
                  template<class> class X_ = nil_t,
                  template<class> class Y_ = nil_t,
                  template<class> class Z_ = nil_t>
        struct implies
            : public uses_template
        {
            typedef typename
                merge<typename
                    uses_template::implies_types, typename
                    implies_template<T_u, A_, B_, C_, D_, E_, F_, G_, H_, I_,
                                          J_, K_, L_, M_, N_, O_, P_, Q_, R_,
                                          S_, T_, U_, V_, W_, X_, Y_, Z_>
                        ::implies_types
                >::type
            implies_types;
        };
    };

    // for_all() helper class.
    template <template<class> class EXEC, class L> struct for_exec
    {
        template <class TX>
        static void all(TX & tx)
        {
            EXEC<typename L::first>::exec(tx);
            for_exec<EXEC, typename L::rest>::all(tx);
        }
    };
    template <template<class> class EXEC> struct for_exec<EXEC, nil>
    {
        template <class TX> static void all(TX &) {}
    };
    // for_all on type lists.
    // for all types T in the list L,
    // calls the static member function EXEC<T>::exec(TX & tx).
    template <class L, template<class> class EXEC, class TX>
    inline void for_all(TX & tx)
    {
        for_exec<EXEC, L>::all(tx);
    }

    template <class T>
    struct has_depends_on : public sfinae_test<T, has_depends_on>
    {
	    template <class U> has_depends_on(U*, typename U::depends_on* = 0);
    };
    template <class T>
    struct has_implies : public sfinae_test<T, has_implies>
    {
	    template <class U> has_implies(U*, typename U::implies_types* = 0);
    };
    template <class T>
    struct has_follows : public sfinae_test<T, has_follows>
    {
	    template <class U> has_follows(U*, typename U::follows_types* = 0);
    };

    // use empty list in case of lacking / faulty depends_on or implies_types:
    template <bool P, class T>
    struct depends_on_guard;
    template <class T>
    struct depends_on_guard<false, T>
    {
        typedef nil type;
    };
    template <class T>
    struct depends_on_guard<true, T>
    {
        typedef typename list_guard<typename T::depends_on>::type type;
    };
    template <class T>
    struct get_pure_depends_on
    {
        typedef typename depends_on_guard<has_depends_on<T>::value, T>::type
            type;
    };

    template <bool P, class T>
    struct follows_guard;
    template <class T>
    struct follows_guard<false, T>
    {
        typedef nil type;
    };
    template <class T>
    struct follows_guard<true, T>
    {
        typedef typename list_guard<typename T::follows_types>::type type;
    };
    template <class T>
    struct get_follows
    {
        typedef typename follows_guard<has_follows<T>::value, T>::type
            type;
    };

    template <class T>
    struct get_depends_on
    {
        typedef typename merge<typename get_pure_depends_on<T>::type,
                               typename get_follows<T>::type       >::type type;
    };

    template <bool P, class T>
    struct implies_guard;
    template <class T>
    struct implies_guard<false, T>
    {
        typedef nil type;
    };
    template <class T>
    struct implies_guard<true, T>
    {
        typedef typename list_guard<typename T::implies_types>::type type;
    };
    template <class T>
    struct get_implies
    {
        typedef typename implies_guard<has_implies<T>::value, T>::type
            type;
    };

    template <class L> struct implies_expand
    {
        typedef typename L::first first;
        typedef typename L::rest  rest;
        typedef
            cons<first, typename
                        merge<typename
                              implies_expand<rest>::type, typename
                              implies_expand<typename
                                  get_implies<first>
        ::type>::type>::type> type;
    };
    template <> struct implies_expand<nil>
    {
        typedef nil type;
    };

    // for_all with type list == T::depends_on (if any.)
    template <class T, template<class> class EXEC, class TX>
    inline void for_all_used(TX & tx)
    {
        for_all<typename get_pure_depends_on<T>::type, EXEC>(tx);
    }

    template <class X, class T>
    struct contains_dependent
    {
        static const bool value
            = contains<X, typename get_depends_on<T>::type>::value;
    };

    template <class X, class XL> struct is_independent_on
    {
        static const bool value
            = ChooseBool<contains_dependent<X, typename XL::first>,
                            VigraFalseType,
                            is_independent_on<X, typename XL::rest>
                        >::value;
    };
    template <class X> struct is_independent_on<X, nil>
    {
        static const bool value = true;
    };

    template <class XL, class YL = XL> struct get_independent
    {
        typedef typename YL::first YL_first;
        typedef typename
        IfBool<is_independent_on<YL_first, XL>::value,
                          YL_first,
                          typename get_independent<XL, typename YL::rest>::type
              >::type type;
    };
    template <class XL> struct get_independent<XL, nil>
    {
        typedef nil type;
    };

    // the output is a list of types in reverse order, starting with the
    // most depedent types.
    template <class L> struct topo_sort
    {
        typedef typename get_independent<L>::type indep;
        typedef typename
        if_nil<indep,
                  nil,
                  cons<indep,
                      typename topo_sort<typename remove<indep, L>::type>::type
                      >
              >::type type;
    };
    template <> struct topo_sort<nil>
    {
        typedef nil type;
    };

    // Topological sort -- the input is a list of types (see below),
    // each of which may, optionally, have an embedded typedef 'depends_on'
    // set to a singly-linked-list of types declared
    // using vigra::type_lists::cons, such as
    // cons<type_a, cons<type_b, cons<type_c> > >
    // (a one-parameter cons will add the trailing nil automatically),
    // -- the output is a list of types with increasing dependence,
    // starting with the indepedent types.
    // Types that should be said lists -- but are in fact not -- are silently
    // replaced by empty lists.

    template <class L> struct topological_sort
    {
        typedef typename
            reverse<typename
                topo_sort<typename
                    unique<typename
                        list_guard<L>
        ::type>::type>::type>::type type;
    };

    template <class L> struct topological_sort_expanded
    {
        typedef typename
            topological_sort<typename
                implies_expand<L>
        ::type>::type type;
    };

    template <class V, unsigned pos = 0>
    class cond_val : public V
    {
        typedef V load_type;
    protected:
        bool is_set_;
    public:
        cond_val() : is_set_(false) {}
        const V & val() const { return *this; }
              V & val()       { return *this; }

        template <class TUPLE>
        bool is_set(const TUPLE &) const
        {
            return is_set_;
        }

        template <class TUPLE>
        void set(TUPLE &)
        {
            is_set_ = true;
        }
        template <class TUPLE>
        void set(const V & v, TUPLE & tuple)
        {
            set(tuple);
            val() = v;
        }
        friend std::ostream & operator<<(std::ostream & os, const cond_val & x)
        {
            if (x.is_set_)
                os << x.val(x);
            else
                os << "<nil>";
            return os;
        }
    };

    template <class V, unsigned pos>
    class bit_cond_val : public V
    {
        typedef V load_type;
    public:
        const V & val() const { return *this; }
              V & val()       { return *this; }

        template <class TUPLE>
        bool is_set(const TUPLE & tuple) const
        {
            return tuple.template is_bit_set<pos>();
        }

        template <class TUPLE>
        void set(TUPLE & tuple)
        {
            tuple.template set_bit<pos>();
        }
        template <class TUPLE>
        void set(const V & v, TUPLE & tuple)
        {
            set(tuple);
            val() = v;
        }
    };

    // no templated virtual functions in C++ ...
    //
    // simple_member_dispatch: sample base and polymorphic adapter classes
    // for cond_virtual_tuple / etc.
    // -- the names 'member_base_type' and 'load_type' of the templates,
    // and also their parameter types, are fixed.
    // -- note that the member_base_type template of any "member_dispatch"
    // class can not directly serve as a base class for any member type that
    // is used within the tuple, since the signatures of most of the
    // virtual functions that are usually needed necessarily require the
    // type of the tuple itself.
    //
    struct simple_member_dispatch
    {
        template <class ACX, class T>
        struct member_base_type
        {
            virtual void operator()() {}
            virtual void operator()(ACX &, const T &) {}
            virtual void operator()(const ACX &, const ACX &, const ACX &) {}
            virtual void call(ACX &, const T &, unsigned) {}
            virtual ~member_base_type() {}
        };
        template <class ACX, class T, class V>
        struct load_type : public member_base_type<ACX, T>, public V
        {
            load_type() {}
            load_type(const V & v) : V(v) {}
            void operator()()
            {
                V::operator()();
            }
            void operator()(ACX & x, const T & v)
            {
                V::operator()(x, v);
            }
            void operator()(const ACX & z, const ACX & x, const ACX & y)
            {
                V::operator()(z, x, y);
            }
            void call(ACX & x, const T & v, unsigned n)
            {
                V::call(x, v, n);
            }
        };
    };

    // polymorphic (conditional) tuple entry, modelled after cond_val
    template <class ACX, class T,  class Z, class V>
    class tuple_entry
    {
        typedef typename Z::template load_type<ACX, T, V> load_type;
        typedef load_type*                                ptr_type;

    protected:
        ptr_type p;
        
    public:
        tuple_entry() : p(0) {}
        template <class TUPLE>
        bool is_set(const TUPLE &) const { return p != 0; }

    protected:
        void make_load()
        {
            if (!p)
                p = new load_type;
        }
        void assign(const V & v)
        {
            ptr_type tmp = new load_type(v);
            delete p;
            p = tmp;
        }
        void check_pointer() const
        {
            if (!p)
                vigra_fail("tuple_entry::operator V &: unused tuple entry "
                           "type V = [" + std::string(typeid(V).name()) + "], "
                           "use set() to create an entry.");
        }
    public:
        operator const V & () const
        {
            check_pointer();
            return *p;
        }
        operator V & ()
        {
            check_pointer();
            return *p;
        }

        template <class TUPLE> // not neccearily identical to ACX
        void set(TUPLE & tuple)
        {
            make_load();
        }
        template <class TUPLE>
        void set(const V & v, TUPLE & tuple)
        {
            set(tuple);
            assign(v);
        }

        tuple_entry & operator=(tuple_entry const & e)
        {
            ptr_type tmp = new load_type(*e.p);
            delete p;
            p = tmp;
            return *this;
        }
        tuple_entry(tuple_entry const & e)
            : p(0)
        {
            if (e.p) // properly care for empty original
                p = new load_type(*e.p);
        }
        ~tuple_entry()
        {
            delete p;
        }
        friend std::ostream & operator<<(std::ostream & os,
                                                          const tuple_entry & x)
        {
            if (x.p)
                os << x.val(x);
            else
                os << "<nil>";
            return os;
        }
    };

    // pos is the position of type V in the type list of the tuple
    template <class ACX, class T, class Z, class V, unsigned pos>
    struct cond_tuple_entry : public tuple_entry<ACX, T, Z, V>
    {
        template <class TUPLE> // not quite identical to ACX
        void set(TUPLE & tuple)
        {
            this->make_load();
            tuple.template add<V>(this->p, pos);
        }
        template <class TUPLE>
        void reassign(TUPLE & tuple)
        {
            if (this->p)
                tuple.reassign(this->p, pos);
        }
        template <class TUPLE>
        void set(const V & v, TUPLE & tuple)
        {
            set(tuple);
            this->assign(v);
        }
    };
   
    // helper classes for tuples
    
    template <unsigned pos, class X>
    struct at_finder
    {
        typedef at_finder<pos - 1, typename X::rest_type> next_type;
        typedef typename next_type::type type;
        static
        type & at(X & x)
        {
            return next_type::at(x.rest);
        }
    };
    template <class X>
    struct at_finder<0, X>
    {
        typedef typename X::first_type type;
        static
        type & at(X & x) {
            return x.first;
        }
    };

    template <class T, class X>
    struct sub_finder
    {
        typedef typename X::rest_type rest_type;
        typedef sub_finder<T, rest_type>
                next_type;
        typedef typename next_type::type type;
        static type & object(X & x)
        {
            return next_type::object(x.rest);
        }
        static const type & const_object(const X & x)
        {
            return next_type::const_object(x.rest);
        }
    };
    template <class X>
    struct sub_finder<typename X::finder_type, X>
    {
        typedef X type;
        static type & object(X & x)
        {
            return x;
        }
        static const type & const_object(const X & x)
        {
            return x;
        }
    };

    template <class T, class X>
    struct ref_finder
    {
        typedef          sub_finder<T, X>         finder;
        typedef typename finder::type::first_type type;
        static
        type & ref(X & x)
        {
            return finder::object(x).first;
        }
        static
        const type & ref(const X & x)
        {
            return finder::const_object(x).first;
        }
    };

    struct binder_0
    {
        template <class F>
        void operator()(F & first)
        {
            first();
        }
        template <class F>
        void call(F & first)
        {
            first.call();
        }
    };
    template <class A>
    struct binder_1
    {
        A v;
        binder_1(A v_) : v(v_) {}
        template <class F>
        void operator()(F & first)
        {
            first(v);
        }
        template <class F>
        void call(F & first)
        {
            first.call(v);
        }
    };
    template <class A, class B>
    struct binder_2
    {
        A v;
        B w;
        binder_2(A v_, B w_) : v(v_), w(w_) {}
        template <class F>
        void operator()(F & first)
        {
            first(v, w);
        }
        template <class F>
        void call(F & first)
        {
            first.call(v, w);
        }
    };
    template <class A, class B, class C>
    struct binder_3
    {
        A v;
        B w;
        C x;
        binder_3(A v_, B w_, C x_) : v(v_), w(w_), x(x_) {}
        template <class F>
        void operator()(F & first)
        {
            first(v, w, x);
        }
        template <class F>
        void call(F & first)
        {
            first.call(v, w, x);
        }
    };

    // mechanism for iterative application of operator() to a tuple
    template <template <class> class TEST>
    struct exec_op_plain
    {
        template <class TUPLE, class B, class TBASE>
        static void exec(TUPLE & tuple, B & binder, const TBASE & z)
        {
            binder(tuple.first);
            tuple.rest.exec_bound_op(binder, z);
        }
        template <class TUPLE, class B, class TBASE>
        static void call(TUPLE & tuple, B & binder, const TBASE & z)
        {
            typedef typename TUPLE::ref_finder_type ref_finder_type;
            if (TEST<ref_finder_type>::value)
                binder.call(static_cast<ref_finder_type &> (tuple.first));
            tuple.rest.call_bound_op(binder, z);
        }
    };

    // policy classes for tuples
    struct plain_global_data
    {
        void reassign() {}
    };
    
    struct plain_chooser // this policy does effectively nothing.
    {
        template <class V, unsigned pos = 0>
        struct use
        {
             typedef V type;
        };

        template <class, template <class> class TEST>
        struct exec_op : public exec_op_plain<TEST> {};

        // "M" & "S" -> bug in cl.exe's parser.
        template <template<class, class, template<class> class M, unsigned>
                  class, class, template<class> class S, unsigned>
        struct global_data : public plain_global_data
        {
        	typedef global_data global_data_type;
        };
        template <class QV, class TUPLE>
        static bool is_set(const QV &, const TUPLE &) { return true; }
        template <class QV, class TUPLE>
        static void set(QV &, TUPLE &) {}
    };

                        // this policy uses the cond_val template to annotate
                        // each tuple member with a bool that steers conditional
                        // execution of each member's operator(), if called via
                        // the tuple's operator().
    struct cond_chooser_plain : public plain_chooser
    {
        template <class V, unsigned pos = 0>
        struct use
        {
            typedef cond_val<V, pos> type;
        };

        template <class, template <class> class TEST>
        struct exec_op
        {
            template <class TUPLE, class B, class TBASE>
            static void exec(TUPLE & tuple, B & binder, const TBASE & z)
            {
                typedef typename TUPLE::ref_finder_type ref_finder_type;
                if (tuple.first.is_set(z))
                    binder(static_cast<ref_finder_type &>(tuple.first));
                tuple.rest.exec_bound_op(binder, z);
            }
            template <class TUPLE, class B, class TBASE>
            static void call(TUPLE & tuple, B & binder, const TBASE & z)
            {
                typedef typename TUPLE::ref_finder_type ref_finder_type;
                if (TEST<ref_finder_type>::value)
                    if (tuple.first.is_set(z))
                        binder.call(static_cast<ref_finder_type &>
                                                                 (tuple.first));
                tuple.rest.call_bound_op(binder, z);
            }
        };

        template <class QV, class TUPLE>
        static bool is_set(const QV & qv, const TUPLE & t)
        {
            return qv.is_set(t);
        }
        template <class QV, class TUPLE>
        static void set(QV & qv, TUPLE & t)
        {
            qv.set(t);
        }
    };
    
    // start the machinery for cond_chooser that produces nested 'if's

    template <class X, class T, class L = typename get_pure_depends_on<T>::type>
    struct depends_on_deep
    {
        static const bool value =
            depends_on_deep<X, T, typename L::rest>::value   // iterate list
            || depends_on_deep<X, typename L::first>::value; // indirect dep.
    };
    template <class T, class L>
    struct depends_on_deep<typename L::first, T, L>
    {
        static const bool value = true;
    };
    template <class X, class T>
    struct depends_on_deep<X, T, nil>
    {
        static const bool value = false;
    };

    template <class T, class R>
    struct first_depends_on
    {
        static const bool value
            =  depends_on_deep<typename R::first, T>::value;
    };
    template <class T>
    struct first_depends_on<T, nil>
    {
        static const bool value = false;
    };

    template <class RRL, class R>
    struct first_depends_on_all_of
    {
        static const bool value
            = ChooseBool<
                  first_depends_on<typename
                      RRL::first,
                      R
                  >,
                  first_depends_on_all_of<typename RRL::rest, R>,
                  VigraFalseType
              >::value;
    };
    template <class R> // end of list RRL: 'success'
    struct first_depends_on_all_of<nil, R>
    {
        static const bool value = true;
    };
    template <class RRL> // 'invalid' input (e.g., at end of cond_op recursion)
    struct first_depends_on_all_of<RRL, nil>
    {
        static const bool value = false;
    };
    template <> // 'invalid' input (e.g., at end of cond_op recursion)
    struct first_depends_on_all_of<nil, nil>
    {
        static const bool value = false;
    };

    // helper structs for cond_op:
    struct null_exec
    {
        template <class TUPLE, class B, class TBASE>
        static void exec(TUPLE &, B &, const TBASE &) {}
        template <class TUPLE, class B, class TBASE>
        static void call(TUPLE &, B &, const TBASE &) {}
        typedef nil iter_leftover_type;
    };
    template <bool cond, class EX>
    struct if_then
    {
        template <class TUPLE, class B, class TBASE>
        static void exec(TUPLE & t, B & b, const TBASE & z)
        {
            IfBool<cond, EX, null_exec>::type::exec(t, b, z);
        }
        template <class TUPLE, class B, class TBASE>
        static void call(TUPLE & t, B & b, const TBASE & z)
        {
            IfBool<cond, EX, null_exec>::type::call(t, b, z);
        }
    };
    template <class ZL, template <class> class TEST, class RRL>
    struct cond_op_inner;
    
    template <class ZL, template <class> class TEST, class RRL = nil>
    struct cond_op
    {
        typedef typename ZL::first             first_type;
        typedef typename ZL::rest              rest_type;
        typedef          cons<first_type, RRL> next_rr_list;

        static const bool recurse_deep
            = first_depends_on<first_type, rest_type>::value;
        typedef cond_op<rest_type, TEST, next_rr_list>
            deep_type;

        typedef typename IfBool<recurse_deep, typename
                                deep_type::iter_leftover_type,
                                rest_type
                               >::type
            deep_leftover_type;

        static const bool iterate
            = first_depends_on_all_of<RRL, deep_leftover_type>::value;

        typedef cond_op_inner<deep_leftover_type, TEST, RRL>
            iter_type;

        // the type list left over from the deep first recursion of exec()
        // and the iteration step: the recursion patterns of both the type
        // 'iter_leftover_type' and the function 'exec()' must match.
        typedef typename IfBool<iterate, typename
                                iter_type::iter_leftover_type,
                                deep_leftover_type
                               >::type
            iter_leftover_type;

        // the code generation templates
        template <class TUPLE, class B, class TBASE>
        static void exec(TUPLE & tuple, B & binder, const TBASE & z)
        {
            if (tuple.first.is_set(z))
            {
                binder(tuple.first);
                if_then<recurse_deep, deep_type>::exec(tuple.rest, binder, z);
            }
            if_then<iterate, iter_type>::exec(tuple, binder, z);
        }
        template <class TUPLE, class B, class TBASE>
        static void call(TUPLE & tuple, B & binder, const TBASE & z)
        {
            if (tuple.first.is_set(z))
            {
                typedef typename TUPLE::ref_finder_type ref_finder_type;
                if (TEST<ref_finder_type>::value)
                    binder.call(static_cast<ref_finder_type &> (tuple.first));
                    
                if_then<recurse_deep, deep_type>::call(tuple.rest, binder, z);
            }
            if_then<iterate, iter_type>::call(tuple, binder, z);
        }
    };
    template <template <class> class TEST, class RRL> // end of type list ZL
    struct cond_op<nil, TEST, RRL> : public null_exec {};

    template <template <class> class TEST, class RRL> // end of type list ZL
    struct cond_op_inner<nil, TEST, RRL> : public null_exec {};

    template <class ZL, template <class> class TEST, class RRL>
    struct cond_op_inner
    {
        typedef cond_op<ZL, TEST, RRL> exec_type;
        typedef typename exec_type::iter_leftover_type iter_leftover_type;

        template <class TUPLE, class B, class TBASE>
        static void exec(TUPLE & tuple, B & binder, const TBASE & z)
        {
            exec_type::
                exec(sub_finder<typename ZL::first, TUPLE>::object(tuple),
                     binder,
                     z);
        }
        template <class TUPLE, class B, class TBASE>
        static void call(TUPLE & tuple, B & binder, const TBASE & z)
        {
            exec_type::
                call(sub_finder<typename ZL::first, TUPLE>::object(tuple),
                     binder,
                     z);
        }
    };

    struct cond_chooser : public cond_chooser_plain
    {
        template <class ZL, template <class> class TEST>
        struct exec_op : public cond_op<ZL, TEST> {};
    };

    struct bit_cond_chooser : public cond_chooser
    {
        template <class V, unsigned pos>
        struct use
        {
            typedef bit_cond_val<V, pos> type;
        };

        // cl.exe wants this -- maybe it is right.
        template <template<class, class, template<class> class M, unsigned>
                  class, class, template<class> class S, unsigned>
        struct global_data : public plain_global_data
        {
        	typedef global_data global_data_type;
        };
        template <template<class, class, template<class> class M, unsigned>
                  class TBASE, class ITL, template<class> class TEST>
        struct global_data<TBASE, ITL, TEST, 0> : public plain_global_data
        {
            // typedef to catch our copy constructor and assignment operator:
            typedef global_data global_data_type;

            BitArray<size<ITL>::of> bit_set;

            void clear() { bit_set.clear(); }
            template <unsigned pos>
            void set_bit() { bit_set.template set<pos>(); }
            template <unsigned pos>
            bool is_bit_set() const { return bit_set.template test<pos>(); }
        };
    };

    template <class ACX, class T, class Z>
    struct virtual_chooser: public cond_chooser_plain
    {
        template <class V, unsigned pos = 0>
        struct use
        {
            typedef tuple_entry<ACX, T, Z, V> type;
        };
    };

    template <class T> struct set_exec
    {
        template <class ACX>
        static void exec(ACX & x)
        {
            x.template set<T>();
        }
    };

    template <class ACX, class T, class Z>
    struct cond_virtual_chooser: public cond_chooser_plain
    {
        template <class V, unsigned pos = 0>
        struct use
        {
            typedef cond_tuple_entry<ACX, T, Z, V, pos> type;
        };

        template <class, template <class> class TEST>
        struct exec_op
        {
            template <class TUPLE, class B, class TBASE>
            static void exec(TUPLE & tuple, B & binder, const TBASE &)
            {
                for (unsigned i = 0; i != tuple.execs.size; ++i)
                    binder(*tuple.execs.pointers[i]);
            }
            template <class TUPLE, class B, class TBASE>
            static void call(TUPLE & tuple, B & binder, const TBASE &)
            {
                for (unsigned i = 0; i != tuple.callers.size; ++i)
                    binder.call(*tuple.callers.pointers[i]);
            }
        };
        // cl.exe wants this -- maybe it is right.
        template <template<class, class, template<class> class M, unsigned>
                  class, class, template<class> class S, unsigned>
        struct global_data : public plain_global_data
        {
        	typedef global_data global_data_type;
        };

        template <class ITL>
        struct global_data_pointers // weak pointers
        {
            typedef typename Z::template member_base_type<ACX, T>* pointer_type;
            static const unsigned array_len = size<ITL>::of;

            unsigned     orders  [array_len]; // consciously use two arrays
            unsigned     size;
            pointer_type pointers[array_len];

            void clear()
            {
                size = 0;
            }
            global_data_pointers()
            {
                clear();
            }
            unsigned* end(unsigned* = 0)
            {
                return orders + size;
            }
            pointer_type* end(pointer_type*)
            {
                return pointers + size;
            }
            template <class E>
            void make_room(E* i)
            {
                std::copy_backward(i, end(i), end(i) + 1);
            }
            typedef std::pair<unsigned*, pointer_type*> finder;
            bool find(unsigned pos, finder & found)
            {
                found.first  = std::lower_bound(orders, end(), pos);
                found.second = pointers + (found.first - orders);
                return found.first != end() && *found.first == pos;
            }
            void add(pointer_type p, unsigned pos)
            {
                // merge pointer according to its topological sort order
                finder found;
                if (find(pos, found))
                    return;
                make_room(found.first);
                make_room(found.second);
                ++size;
                *found.first  = pos;
                *found.second = p;
            }
            void reassign(pointer_type p, unsigned pos)
            {
                // replace pointers -- the present values still belong to
                // the old tuple object that was copied from
                finder found;
                if (find(pos, found))
                    *found.second = p;
            }
        };
        template <template<class, class, template<class> class M, unsigned>
                  class TBASE, class ITL, template<class> class TEST>
        struct global_data<TBASE, ITL, TEST, 0>
        {
            // typedef to catch our copy constructor and assignment operator:
            typedef global_data global_data_type;
            // the derived class:
            typedef TBASE<ITL, cond_virtual_chooser, TEST, 0> crtp_type;

            typedef global_data_pointers<ITL>        ptrs_type;
            typedef typename ptrs_type::pointer_type pointer_type;
            ptrs_type callers;
            ptrs_type execs;

            template <class V>
            void add(pointer_type p, unsigned pos)
            {
                execs.add(p, pos);
                if (TEST<V>::value)
                    callers.add(p, pos);
            }
            void reassign(pointer_type p, unsigned pos)
            {
                execs.  reassign(p, pos);
                callers.reassign(p, pos);
            }

            template <class V>
            struct reassign_op
            {
                static void exec(crtp_type & tuple)
                {
                    tuple.template ref<V>().reassign(tuple);
                }
            };

            void reassign()
            {
                for_all<ITL, reassign_op>(static_cast<crtp_type &>(*this));
            }
        };
    };


    template <class ITL, class Q = plain_chooser,
              template<class> class TEST = true_test, unsigned index = 0>
    struct tuple_base
        : public Q::template global_data<tuple_base, ITL, TEST, index>
    {
        typedef typename tuple_base::global_data_type global_data_base_type;
        typedef nil derived_type; // dummy declaration for static_cast_2<>
        typedef tuple_base tuple_type;
        typedef ITL        list_type;

        // use the original types from the list to find tuple members via ref():
        typedef typename ITL::first                            finder_type;
        typedef typename ITL::rest                             rest_list_type;

        typedef tuple_base<rest_list_type, Q, TEST, index + 1> rest_type;

        // use policy class Q to annotate the types of the type list ITL
        // for use as members of the tuple:
        // -- the actually stored type
        typedef typename ITL::first                            ref_finder_type;
        typedef typename Q::template use<ref_finder_type,
                                         index>::type          first_type;
        first_type first;
        rest_type  rest;

        template <class T>
        struct has_element
        {
            static const bool value = contains<T, ITL>::value;
        };

        template <unsigned pos>
        typename at_finder<pos, tuple_base>::type & at()
        {
            return at_finder<pos, tuple_base>::at(*this);
        }

        template <class T>
        struct ref_returns
        {
            typedef ref_finder<T, tuple_base>     finder;
            typedef typename finder::type         ref_finder_type;
            typedef       ref_finder_type &       type;
            typedef const ref_finder_type & const_type;
        };
        template <class T>
        typename ref_returns<T>::type
        ref()
        {
            return ref_returns<T>::finder::ref(*this);
        }
        template <class T>
        typename ref_returns<T>::const_type
        ref() const
        {
            return ref_returns<T>::finder::ref(*this);
        }

        template <class V>
        bool is_set() const
        {
            return Q::template is_set(this->template ref<V>(), *this);
        }
        template <class RV>
        void set_if_not(RV & rv)
        {
            if (! Q::template is_set(rv, *this))
                Q::template set(rv, *this);
        }
        // recursively set all the cond_val::is_set bits of the depended-on
        // types, or, respectively, take equivalent action.
        template <class V>
        void set()
        {
            set_if_not(this->template ref<V>());
            for_all_used<V, set_exec>(*this);
        }
        // transfer the set bits of *this to another tuple t2:
        // (ITL must be a subset of ITL2.)
        template <class ITL2, class Q2, template<class> class TEST2>
        void transfer_set_to(tuple_base<ITL2, Q2, TEST2> & t2) const
        {
            if (is_set<ref_finder_type>())
                t2.template set<ref_finder_type>();
            rest.transfer_set_to(t2);
        }

        // policy-based application of operator()
        template <class B, class TBASE>
        void exec_bound_op(B & binder, const TBASE & z)
        {
            Q::template exec_op<list_type, true_test>::exec(*this, binder, z);
        }
        template <class B>
        void exec_bound(B & binder)
        {
            exec_bound_op(binder, *this);
        }

        void operator()()
        {
            binder_0 binder;
            exec_bound(binder);
        }

        template <class A>
        void operator()(const A & v)
        {
            binder_1<const A &> binder(v);
            exec_bound(binder);
        }
        template <class A>
        void operator()(A & v)
        {
            binder_1<A &> binder(v);
            exec_bound(binder);
        }
        template <class A>
        void operator()(A & v) const
        {
            binder_1<A &> binder(v);
            exec_bound(binder);
        }

        template <class A, class B>
        void operator()(const A & v, const B & w)
        {
            binder_2<const A &, const B &> binder(v, w);
            exec_bound(binder);
        }
        template <class A, class B>
        void operator()(A & v, const B & w)
        {
            binder_2<A &, const B &> binder(v, w);
            exec_bound(binder);
        }
        template <class A, class B>
        void operator()(A & v, const B & w) const
        {
            binder_2<A &, const B &> binder(v, w);
            exec_bound(binder);
        }

        template <class A, class B, class C>
        void operator()(const A & v, const B & w, const C & x)
        {
            binder_3<const A &, const B &, const C &> binder(v, w, x);
            exec_bound(binder);
        }
        template <class A, class B, class C>
        void operator()(A & v, const B & w, const C & x)
        {
            binder_3<A &, const B &, const C &> binder(v, w, x);
            exec_bound(binder);
        }
        template <class A, class B, class C>
        void operator()(A & v, const B & w, const C & x) const
        {
            binder_3<A &, const B &, const C &> binder(v, w, x);
            exec_bound(binder);
        }

        // policy-based application of call()
        template <class B, class TBASE>
        void call_bound_op(B & binder, const TBASE & z)
        {
            Q::template exec_op<list_type, TEST>::call(*this, binder, z);
        }
        template <class B>
        void call_bound(B & binder)
        {
            call_bound_op(binder, *this);
        }

        template <class A, class B>
        void call(const A & v, const B & w)
        {
            binder_2<const A &, const B &> binder(v, w);
            call_bound(binder);
        }
        template <class A, class B>
        void call(A & v, const B & w)
        {
            binder_2<A &, const B &> binder(v, w);
            call_bound(binder);
        }
        template <class A, class B>
        void call(A & v, const B & w) const
        {
            binder_2<A &, const B &> binder(v, w);
            call_bound(binder);
        }

        template <class A, class B, class C>
        void call(const A & v, const B & w, const C & x)
        {
            binder_3<const A &, const B &, const C &> binder(v, w, x);
            call_bound(binder);
        }
        template <class A, class B, class C>
        void call(A & v, const B & w, const C & x)
        {
            binder_3<A &, const B &, const C &> binder(v, w, x);
            call_bound(binder);
        }
        template <class A, class B, class C>
        void call(A & v, const B & w, const C & x) const
        {
            binder_3<A &, const B &, const C &> binder(v, w, x);
            call_bound(binder);
        }

        // the possible reassign() of global data requires handmade copies here
        tuple_base() {}
        tuple_base(const tuple_base & x)
            : global_data_base_type(x), first(x.first), rest(x.rest)
        {
            this->reassign();
        }
        tuple_base & operator=(const tuple_base & x)
        {
            global_data_base_type::operator=(x);
            first = x.first;
            rest  = x.rest;
            this->reassign();
            return *this;
        }
    };
    template <class Q, template <class> class TEST, unsigned index>
    struct tuple_base<nil, Q, TEST, index>
    {
        template <class>
        struct has_element
        {
            static const bool value = false;
        };
        template <class B, class TBASE>
        void exec_bound_op(B &, const TBASE &) {}
        template <class B, class TBASE>
        void call_bound_op(B &, const TBASE &) {}
        template <class ITL2, class Q2, template<class> class TEST2>
        void transfer_set_to(tuple_base<ITL2, Q2, TEST2> &) const {}
    };

    template <class V, class L, class A, class B>
    struct choose_tuple
    {
        static const bool value = contains<V, L>::value;
        typedef typename IfBool<value, A, B>::type type;
        static type & at(A & a, B & b)
        {
            return choose_type<value>::at(a, b);
        }
        static const type & at(const A & a, const B & b)
        {
            return choose_type<value>::at(a, b);
        }

        typedef typename type::template ref_returns<V>::      type
            ref_type;
        typedef typename type::template ref_returns<V>::const_type
            const_ref_type;

        ref_type
        static ref(A & a, B & b)
        {
            return at(a, b).template ref<V>();
        }
        const_ref_type
        static ref(const A & a, const B & b)
        {
            return at(a, b).template ref<V>();
        }
    };

    template <class ITL, template<class> class TEST = true_test>
    struct tuple : public tuple_base<ITL, plain_chooser, TEST> {};

    template <class ITL, template<class> class TEST = true_test>
    struct cond_tuple_plain
        : public tuple_base<ITL, cond_chooser_plain, TEST> {};

    template <class ITL, template<class> class TEST = true_test>
    struct cond_tuple : public tuple_base<ITL, cond_chooser, TEST> {};

    template <class ITL, template<class> class TEST = true_test>
    struct bit_cond_tuple : public tuple_base<ITL, bit_cond_chooser, TEST> {};

    template <class ITL, class ACX, class T, class Z,
              template<class> class TEST = true_test>
    struct virtual_tuple_base
        : public tuple_base<ITL, virtual_chooser<ACX, T, Z>, TEST> {};

    template <template <class, class> class ACXTT, class T,
              class Z = simple_member_dispatch>
    struct virtual_tuple
    {
        template <class ITL, template<class> class TEST = true_test>
        struct type
            : public virtual_tuple_base<ITL, ACXTT<T, ITL>, T, Z, TEST> {};
    };

    template <class ITL, class ACX, class T, class Z,
              template<class> class TEST = true_test>
    struct cond_virtual_tuple_base
        : public tuple_base<ITL, cond_virtual_chooser<ACX, T, Z>, TEST>
    {
        typedef ACX derived_type;
    };

    template <template <class, class> class ACXTT, class T,
              class Z = simple_member_dispatch>
    struct cond_virtual_tuple
    {
        template <class ITL, template<class> class TEST = true_test>
        struct type
            : public cond_virtual_tuple_base<ITL, ACXTT<T, ITL>, T, Z, TEST> {};
    };

} // namespace type_lists

// mask cl.exe shortcomings [end]
#if defined(_MSC_VER)
#pragma warning( pop )
#endif

} // namespace vigra

#endif // VIGRA_TYPE_LISTS_HXX
