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


#ifndef VIGRA_BASICS_HXX
#define VIGRA_BASICS_HXX

#include "config.hxx"
#include "error.hxx"
#include "metaprogramming.hxx"
#include "tuple.hxx"
#include "diff2d.hxx"
#include "mathutil.hxx"
#include <string>
#include <sstream>
#include <cctype>

namespace vigra {

/** Convert a value to a string. Available for integral and floating point types
    and void *.
*/
doxygen_overloaded_function(template <class T> std::string asString(T t))

#define VIGRA_AS_STRING(T) \
inline std::string asString(T t) \
{ \
    std::stringstream s; \
    s << t; \
    return s.str(); \
}

VIGRA_AS_STRING(bool)
VIGRA_AS_STRING(signed char)
VIGRA_AS_STRING(unsigned char)
VIGRA_AS_STRING(signed short)
VIGRA_AS_STRING(unsigned short)
VIGRA_AS_STRING(signed long)
VIGRA_AS_STRING(unsigned long)
VIGRA_AS_STRING(signed long long)
VIGRA_AS_STRING(unsigned long long)
VIGRA_AS_STRING(signed int)
VIGRA_AS_STRING(unsigned int)
VIGRA_AS_STRING(float)
VIGRA_AS_STRING(double)
VIGRA_AS_STRING(long double)
VIGRA_AS_STRING(void *)

#undef VIGRA_AS_STRING

template <class T>
std::string operator<<(std::string const & s, T const & t)
{
    std::stringstream ss;
    ss << t; 
    return s + ss.str();
}

    /** Convert string to lower case.
    */
inline std::string tolower(std::string s)
{
    for(unsigned int k=0; k<s.size(); ++k)
        s[k] = (std::string::value_type)std::tolower(s[k]);
    return s;
}

inline std::string tolower(const char * s)
{
    return tolower(std::string(s));
}

    /** Convert string to lower case and remove any white space characters.
    */
inline std::string normalizeString(std::string const & s)
{
    std::string res;
    
    for(unsigned int k=0; k<s.size(); ++k)
    {
        if(std::isspace(s[k]))
            continue;
        res += (std::string::value_type)std::tolower(s[k]);
    }
    return res;
}

inline std::string normalizeString(const char * s)
{
    return normalizeString(std::string(s));
}

namespace detail {

template <class T>
struct FinallyImpl
{
    T & destructor_;
    
    FinallyImpl(T & destructor)
    : destructor_(destructor)
    {}
    
    ~FinallyImpl()
    {
        destructor_();
    }
};

} // namespace detail

} // namespace vigra

#define VIGRA_TOKEN_PASTE_IMPL(x, y) x##y
#define VIGRA_TOKEN_PASTE(x, y) VIGRA_TOKEN_PASTE_IMPL(x, y)

#define VIGRA_FINALLY_IMPL(destructor, counter) \
    auto VIGRA_TOKEN_PASTE(_vigra_finally_impl_, counter) = [&]() { destructor; }; \
    ::vigra::detail::FinallyImpl<decltype(VIGRA_TOKEN_PASTE(_vigra_finally_impl_, counter))> \
        VIGRA_TOKEN_PASTE(_vigra_finally_, counter)(VIGRA_TOKEN_PASTE(_vigra_finally_impl_, counter))

    /** Emulate the 'finally' keyword as known from Python and other languages.
    
        This macro improves upon the famous 
        <a href="http://en.wikipedia.org/wiki/Resource_Acquisition_Is_Initialization">Resource Acquisition Is Initialization</a> idiom, where a resource (e.g. heap memory or a mutex) is automatically free'ed when program execution leaves the current scope. Normally, this is implemented by calling a suitable function in the destructor of a dedicated helper class (e.g. <tt>std::unique_ptr</tt> or <tt>std::lock_guard<std::mutex></tt>). 
        
        Traditionally, a separate helper class has to be implemented for each kind of resource to be handled. In contrast, the macro <tt>VIGRA_FINALLY</tt> creates such a class on the fly by means of an embedded lambda expression.
        
        <b>Usage:</b>
        
        <b>\#include</b> \<vigra/utilities.hxx\><br/>

        \code
        std::mutex my_mutex;
        ...
        {
            // the following two lines are equivalent to 
            //     std::unique_ptr<std::string> my_string = new std::string("foo");
            std::string * my_string = new std::string("foo");
            VIGRA_FINALLY(delete my_string);
      
            // the following two lines are equivalent to 
            //     std::lock_guard<std::mutex> lock(my_mutex);
            my_mutex.lock();
            VIGRA_FINALLY(my_mutex.unlock());
      
            ...
        }
        // the string has been deallocated and the mutex is unlocked
        \endcode
        
        You can pass any code to this macro. Multiple statements must be enclosed in braces as usual. Arbitrary many calls to <tt>VIGRA_FINALLY</tt> can be placed in the same scope. Their actions will be executed in the reversed order of declaration:

        \code
        int i = 0;
        ...
        {
            VIGRA_FINALLY({           // execute multiple statements
                i = i*i;
                ++i;
            });
      
            VIGRA_FINALLY( i += 2 );  // this executes first
            
            assert(i == 0);           // as yet, nothing happend
        }
        assert(i == 5);               // 'finally' code was executed in reversed order at end-of-scope
        \endcode
        
        This idea was popularized by Marko Tintor in "<a href="http://blog.memsql.com/c-error-handling-with-auto/">The Auto Macro: A Clean Approach to C++ Error Handling</a>".
    */
#define VIGRA_FINALLY(destructor) \
    VIGRA_FINALLY_IMPL(destructor, __COUNTER__)
    
namespace std {

template <class T1, class T2>
ostream & operator<<(ostream & s, std::pair<T1, T2> const & p)
{
    s << "(" << p.first << ", " << p.second << ")";
    return s;
}

}

/** \page Utilities Utilities
    Basic helper functionality needed throughout.

    <UL style="list-style-image:url(documents/bullet.gif)">
    <LI> \ref vigra::ArrayVector
         <BR>&nbsp;&nbsp;&nbsp;<em>replacement for std::vector (always uses consecutive memory)</em>
    <LI> \ref vigra::BucketQueue and \ref vigra::MappedBucketQueue
         <BR>&nbsp;&nbsp;&nbsp;<em>efficient priority queues for integer priorities</em>
    <LI> \ref RangesAndPoints
         <BR>&nbsp;&nbsp;&nbsp;<em>2-D and N-D positions, extents, and boxes</em>
    <LI> \ref PixelNeighborhood
         <BR>&nbsp;&nbsp;&nbsp;<em>4- and 8-neighborhood definitions and circulators</em>
    <LI> \ref VoxelNeighborhood
         <BR>&nbsp;&nbsp;&nbsp;<em>6- and 26-neighborhood definitions and circulators</em>
    <LI> \ref vigra::IteratorAdaptor
         <BR>&nbsp;&nbsp;&nbsp;<em>Quickly create STL-compatible 1D iterator adaptors</em>
    <LI> \ref TupleTypes
         <BR>&nbsp;&nbsp;&nbsp;<em>pair, triple, tuple4, tuple5</em>
    <LI> \ref MathConstants
         <BR>&nbsp;&nbsp;&nbsp;<em>M_PI, M_SQRT2</em>
    <LI> \ref TimingMacros
         <BR>&nbsp;&nbsp;&nbsp;<em>Macros for taking execution speed measurements</em>
    <LI> \ref VIGRA_FINALLY
         <BR>&nbsp;&nbsp;&nbsp;<em>Emulation of the 'finally' keyword from Python</em>
    </UL>
*/

#endif // VIGRA_BASICS_HXX
