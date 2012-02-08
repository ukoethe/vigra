#ifndef VIGRA_INSPECTOR_PASSES_HXX
#define VIGRA_INSPECTOR_PASSES_HXX

#include "metaprogramming.hxx"

namespace vigra {

// test and accomodate for functors that require extra passes over arrays / etc.

namespace detail {

template <bool>
struct extra_passes_selector
{
    template <class Inspector, class Functor>
    static void
    call(Inspector g, Functor & f) {}
};
template <>
struct extra_passes_selector<true>
{
    template <class Inspector, class Functor_n>
    static void
    call_n(Inspector g, Functor_n f_n)
    {
        g(f_n);
    }
    template <class Inspector, class Functor>
    static void
    call(Inspector g, Functor & f)
    {
        for (unsigned n = 2; n <= Functor::max_passes; ++n)
        {
            f.calc_sync();
            call_n(g, f.pass_n(n));
        }
    }
};

template <class T>
struct has_extra_passes : public sfinae_test<T, has_extra_passes>
{
	template <class U> has_extra_passes(U*, typename U::extra_passes* = 0);
};

template <class Functor, bool extra = has_extra_passes<Functor>::value>
struct get_extra_passes
    : public VigraFalseType
{
    void sync(Functor &) {}
};

template <class Functor>
struct get_extra_passes<Functor, true>
{
    typedef get_extra_passes extra_passes;
    static const unsigned max_passes = Functor::max_passes;
    static const bool value = Functor::max_passes >= 2;

    void sync(Functor & f)
    {
        f.calc_sync();
    }
};

template <class Inspector, class Functor>
inline
void
extra_passes_select(Inspector g, Functor & f)
{
    g(f);
    extra_passes_selector<get_extra_passes<Functor>::value>::call(g, f);
}

} // namespace detail

} // namespace vigra

#endif // VIGRA_INSPECTOR_PASSES_HXX
