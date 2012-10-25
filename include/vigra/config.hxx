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


#ifndef VIGRA_CONFIG_HXX
#define VIGRA_CONFIG_HXX

#include <vigra/configVersion.hxx>
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

    #if (_MSC_VER < 1300)
        #define NO_TYPENAME         // no 'typename' keyword
        #define TEMPLATE_COPY_CONSTRUCTOR_BUG
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
    #endif // (_MSC_VER < 1300)

    #if _MSC_VER < 1310
        #pragma warning( disable : 4786 4250 4244 4305)

        #define NO_PARTIAL_TEMPLATE_SPECIALIZATION
        #define NO_OUT_OF_LINE_MEMBER_TEMPLATES
        #include <cmath>

        #ifdef _MSC_EXTENSIONS
        #ifndef CMATH_NOT_IN_STD
                namespace std {
        #endif // CMATH_NOT_IN_STD
                inline double abs(double v) { return fabs(v); }
                inline float  abs(float v)  { return fabs(v); }
        #ifndef CMATH_NOT_IN_STD
                }
        #endif // CMATH_NOT_IN_STD
        #endif // _MSC_EXTENSIONS
    #endif // _MSC_VER < 1310

    #if _MSC_VER < 1400
        #define VIGRA_NO_WORKING_STRINGSTREAM
    #endif
    
    #if _MSC_VER < 1600
        #define VIGRA_NO_UNIQUE_PTR
    #endif
    
    #define VIGRA_NEED_BIN_STREAMS
    
    #ifndef VIGRA_ENABLE_ANNOYING_WARNINGS
        #pragma warning ( disable: 4244 4267) // implicit integer conversion warnings
    #endif

    #ifdef VIGRA_DLL
        #define VIGRA_EXPORT __declspec(dllexport)
    #elif defined(VIGRA_STATIC_LIB)
        #define VIGRA_EXPORT
    #else
        #define VIGRA_EXPORT __declspec(dllimport)
    #endif
#endif // _MSC_VER

///////////////////////////////////////////////////////////
//                                                       //
//                           gcc                         //
//                                                       //
///////////////////////////////////////////////////////////

#if defined(__GNUC__)
    #if  __GNUC__ < 2 || ((__GNUC__ == 2) && (__GNUC_MINOR__ <= 8))
        #error "Need at least g++ 2.95"
    #endif
    #if __GNUC__ < 3
        #define VIGRA_NO_WORKING_STRINGSTREAM
    #endif
    #define HAS_HASH_CONTAINERS
    
    // these warnings produce too many false positives to be useful
    #pragma GCC diagnostic ignored "-Wstrict-aliasing"  
    #pragma GCC diagnostic ignored "-Wshadow"  
    
    #if !defined(__GXX_EXPERIMENTAL_CXX0X) && __cplusplus < 201103L
        #define VIGRA_NO_UNIQUE_PTR
    #endif

#endif  // __GNUC__

///////////////////////////////////////////////////////////
//                                                       //
//                         MingW                         //
//                                                       //
///////////////////////////////////////////////////////////

#if defined(__MINGW32__)
    #define VIGRA_NEED_BIN_STREAMS

    #ifdef VIGRA_DLL
        #define VIGRA_EXPORT __declspec(dllexport)
    #elif defined(VIGRA_STATIC_LIB)
        #define VIGRA_EXPORT
    #else
        #define VIGRA_EXPORT __declspec(dllimport)
    #endif
#endif  // __MINGW32__

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
//                      Sun C++ ???                      //
//                                                       //
///////////////////////////////////////////////////////////

#if defined(__sun) && !defined(__GNUC__)
    #define VIGRA_HAS_ERF
#endif // __sun

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

#ifndef VIGRA_EXPORT
    #define VIGRA_EXPORT
#endif

#ifdef VIGRA_NO_UNIQUE_PTR
#  define VIGRA_UNIQUE_PTR  std::auto_ptr
#else
#  define VIGRA_UNIQUE_PTR  std::unique_ptr
#endif

namespace vigra {

#ifndef SPECIAL_STDEXCEPTION_DEFINITION_NEEDED
     typedef std::exception StdException;
#endif

} // namespace vigra

#ifdef DOXYGEN
#  define doxygen_overloaded_function(fun) fun(...);
#else
#  define doxygen_overloaded_function(fun)
#endif


#endif // VIGRA_CONFIG_HXX
