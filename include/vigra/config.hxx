/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
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
 
 
#ifndef VIGRA_CONFIG_HXX 
#define VIGRA_CONFIG_HXX

#include <stdexcept>

///////////////////////////////////////////////////////////
//                                                       //
//                   VisualC++ 5.0                       //
//                                                       //
///////////////////////////////////////////////////////////

#ifdef _MSC_VER
    #if(_MSC_VER < 1100)    // before VisualC++ 5.0
        #error "Need VisualC++ 5.0, Service Pack 2, or later"
    #endif // _MSC_VER < 1100
    
    #pragma warning( disable : 4786 4250 4244 4305)
    #define NO_TYPENAME         // no 'typename' keyword
    #define TEMPLATE_COPY_CONSTRUCTOR_BUG
    #define NO_PARTIAL_TEMPLATE_SPECIALIZATION
    #define NO_STL_MEMBER_TEMPLATES
    #define NO_INLINE_STATIC_CONST_DEFINITION
    #define CMATH_NOT_IN_STD
    #define NO_COVARIANT_RETURN_TYPES
    
    #ifdef VIGRA_NO_STD_MINMAX  // activate if necessary
    namespace std {
    
	template<class T>
	const T& min(const T& x, const T& y)
	{
		return (y < x) 
			? y
			: x;
	}

	template<class T>
	const T& max(const T& x, const T& y)
	{
		return (x < y) 
			? y
			: x;
	}
    }

    #endif // VIGRA_NO_STD_MINMAX
#endif // _MSC_VER

///////////////////////////////////////////////////////////
//                                                       //
//                        egcs 1.1                       //
//                                                       //
///////////////////////////////////////////////////////////

#if defined(__GNUC__) 
    #if  __GNUC__ < 2 || ((__GNUC__ == 2) && (__GNUC_MINOR__ <= 8))
        #error "Need at least egcs 1.1 or g++ 2.95"
    #endif 
    #define HAS_HASH_CONTAINERS
#endif  // __GNUC__

///////////////////////////////////////////////////////////
//                                                       //
//                      SGI C++ 7.2                      //
//                                                       //
///////////////////////////////////////////////////////////

#if defined(__sgi) && !defined(__GNUC__) 
    #if _COMPILER_VERSION < 720
        #error "Need SGI C++ 7.2 or later"
    #endif
    #if (_COMPILER_VERSION  == 720) || (_COMPILER_VERSION  == 721)
        #define SPECIAL_STDEXCEPTION_DEFINITION_NEEDED

        namespace vigra {
            typedef std::exception StdException; // must be above next #define !!
        }
        #define std
        #define NO_NAMESPACE_STD
    #endif // _COMPILER_VERSION
    #define HAS_HASH_CONTAINERS
#endif // __sgi

///////////////////////////////////////////////////////////
//                                                       //
//                        general                        //
//                                                       //
///////////////////////////////////////////////////////////

#ifdef CMATH_NOT_IN_STD
    #define VIGRA_CSTD
#else
    #define VIGRA_CSTD std
#endif

#ifdef NO_TYPENAME
    #define typename
#endif

#ifdef NO_EXPLICIT
    #define explicit
#endif

namespace vigra {

#ifndef SPECIAL_STDEXCEPTION_DEFINITION_NEEDED
     typedef std::exception StdException;
#endif

} // namespace vigra

#endif // VIGRA_CONFIG_HXX
