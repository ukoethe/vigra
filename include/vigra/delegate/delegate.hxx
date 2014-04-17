/*
	(c) Sergey Ryazanov (http://home.onego.ru/~ryazanov)

	Fast delegate compatible with C++ Standard.
*/
#ifndef VIGRA_DELEGATE_INCLUDED
#define VIGRA_DELEGATE_INCLUDED

namespace vigra
{
#ifdef VIGRA_DELEGATE_PREFERRED_SYNTAX
	template <typename TSignature> class delegate;
	template <typename TSignature> class delegate_invoker;
#endif
}

#ifdef _MSC_VER
#define VIGRA_DELEGATE_CALLTYPE __fastcall
#else
#define VIGRA_DELEGATE_CALLTYPE
#endif

#include "detail/delegate_list.hxx"

#undef VIGRA_DELEGATE_CALLTYPE

#endif//VIGRA_DELEGATE_INCLUDED
