/*
	(c) Sergey Ryazanov (http://home.onego.ru/~ryazanov)

	Template file. May be included many times with different predefined macros.
*/
#if VIGRA_DELEGATE_PARAM_COUNT > 0
#define VIGRA_DELEGATE_SEPARATOR ,
#else
#define VIGRA_DELEGATE_SEPARATOR
#endif

// see BOOST_JOIN for explanation
#define VIGRA_DELEGATE_JOIN_MACRO( X, Y) VIGRA_DELEGATE_DO_JOIN( X, Y )
#define VIGRA_DELEGATE_DO_JOIN( X, Y ) VIGRA_DELEGATE_DO_JOIN2(X,Y)
#define VIGRA_DELEGATE_DO_JOIN2( X, Y ) X##Y

namespace vigra
{
#ifdef VIGRA_DELEGATE_PREFERRED_SYNTAX
#define VIGRA_DELEGATE_CLASS_NAME delegate
#define VIGRA_DELEGATE_INVOKER_CLASS_NAME delegate_invoker
#else
#define VIGRA_DELEGATE_CLASS_NAME VIGRA_DELEGATE_JOIN_MACRO(delegate,VIGRA_DELEGATE_PARAM_COUNT)
#define VIGRA_DELEGATE_INVOKER_CLASS_NAME VIGRA_DELEGATE_JOIN_MACRO(delegate_invoker,VIGRA_DELEGATE_PARAM_COUNT)
	template <typename R VIGRA_DELEGATE_SEPARATOR VIGRA_DELEGATE_TEMPLATE_PARAMS>
	class VIGRA_DELEGATE_INVOKER_CLASS_NAME;
#endif

	template <typename R VIGRA_DELEGATE_SEPARATOR VIGRA_DELEGATE_TEMPLATE_PARAMS>
#ifdef VIGRA_DELEGATE_PREFERRED_SYNTAX
	class VIGRA_DELEGATE_CLASS_NAME<R (VIGRA_DELEGATE_TEMPLATE_ARGS)>
#else
	class VIGRA_DELEGATE_CLASS_NAME
#endif
	{
	public:
		typedef R return_type;
#ifdef VIGRA_DELEGATE_PREFERRED_SYNTAX
		typedef return_type (VIGRA_DELEGATE_CALLTYPE *signature_type)(VIGRA_DELEGATE_TEMPLATE_ARGS);
		typedef VIGRA_DELEGATE_INVOKER_CLASS_NAME<signature_type> invoker_type;
#else
		typedef VIGRA_DELEGATE_INVOKER_CLASS_NAME<R VIGRA_DELEGATE_SEPARATOR VIGRA_DELEGATE_TEMPLATE_ARGS> invoker_type;
#endif

		VIGRA_DELEGATE_CLASS_NAME()
			: object_ptr(0)
			, stub_ptr(0)
		{}

		template <return_type (*TMethod)(VIGRA_DELEGATE_TEMPLATE_ARGS)>
		static VIGRA_DELEGATE_CLASS_NAME from_function()
		{
			return from_stub(0, &function_stub<TMethod>);
		}

		template <class T, return_type (T::*TMethod)(VIGRA_DELEGATE_TEMPLATE_ARGS)>
		static VIGRA_DELEGATE_CLASS_NAME from_method(T* object_ptr)
		{
			return from_stub(object_ptr, &method_stub<T, TMethod>);
		}

		template <class T, return_type (T::*TMethod)(VIGRA_DELEGATE_TEMPLATE_ARGS) const>
		static VIGRA_DELEGATE_CLASS_NAME from_const_method(T const* object_ptr)
		{
			return from_stub(const_cast<T*>(object_ptr), &const_method_stub<T, TMethod>);
		}

		return_type operator()(VIGRA_DELEGATE_PARAMS) const
		{
			return (*stub_ptr)(object_ptr VIGRA_DELEGATE_SEPARATOR VIGRA_DELEGATE_ARGS);
		}

		operator bool () const
		{
			return stub_ptr != 0;
		}

		bool operator!() const
		{
			return !(operator bool());
		}

	private:
		
		typedef return_type (VIGRA_DELEGATE_CALLTYPE *stub_type)(void* object_ptr VIGRA_DELEGATE_SEPARATOR VIGRA_DELEGATE_PARAMS);

		void* object_ptr;
		stub_type stub_ptr;

		static VIGRA_DELEGATE_CLASS_NAME from_stub(void* object_ptr, stub_type stub_ptr)
		{
			VIGRA_DELEGATE_CLASS_NAME d;
			d.object_ptr = object_ptr;
			d.stub_ptr = stub_ptr;
			return d;
		}

		template <return_type (*TMethod)(VIGRA_DELEGATE_TEMPLATE_ARGS)>
		static return_type VIGRA_DELEGATE_CALLTYPE function_stub(void* VIGRA_DELEGATE_SEPARATOR VIGRA_DELEGATE_PARAMS)
		{
			return (TMethod)(VIGRA_DELEGATE_ARGS);
		}

		template <class T, return_type (T::*TMethod)(VIGRA_DELEGATE_TEMPLATE_ARGS)>
		static return_type VIGRA_DELEGATE_CALLTYPE method_stub(void* object_ptr VIGRA_DELEGATE_SEPARATOR VIGRA_DELEGATE_PARAMS)
		{
			T* p = static_cast<T*>(object_ptr);
			return (p->*TMethod)(VIGRA_DELEGATE_ARGS);
		}

		template <class T, return_type (T::*TMethod)(VIGRA_DELEGATE_TEMPLATE_ARGS) const>
		static return_type VIGRA_DELEGATE_CALLTYPE const_method_stub(void* object_ptr VIGRA_DELEGATE_SEPARATOR VIGRA_DELEGATE_PARAMS)
		{
			T const* p = static_cast<T*>(object_ptr);
			return (p->*TMethod)(VIGRA_DELEGATE_ARGS);
		}
	};

	template <typename R VIGRA_DELEGATE_SEPARATOR VIGRA_DELEGATE_TEMPLATE_PARAMS>
#ifdef VIGRA_DELEGATE_PREFERRED_SYNTAX
	class VIGRA_DELEGATE_INVOKER_CLASS_NAME<R (VIGRA_DELEGATE_TEMPLATE_ARGS)>
#else
	class VIGRA_DELEGATE_INVOKER_CLASS_NAME
#endif
	{
		VIGRA_DELEGATE_INVOKER_DATA

	public:
		VIGRA_DELEGATE_INVOKER_CLASS_NAME(VIGRA_DELEGATE_PARAMS)
#if VIGRA_DELEGATE_PARAM_COUNT > 0
			:
#endif
			VIGRA_DELEGATE_INVOKER_INITIALIZATION_LIST
		{
		}

		template <class TDelegate>
		R operator()(TDelegate d) const
		{
			return d(VIGRA_DELEGATE_ARGS);
		}
	};
}

#undef VIGRA_DELEGATE_CLASS_NAME
#undef VIGRA_DELEGATE_SEPARATOR
#undef VIGRA_DELEGATE_JOIN_MACRO
#undef VIGRA_DELEGATE_DO_JOIN
#undef VIGRA_DELEGATE_DO_JOIN2
