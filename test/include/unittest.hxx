/************************************************************************/
/*                                                                      */
/* a simple unit test framework, similar to Kent Beck's JUnit           */
/*                                                                      */
/*               Copyright 2002-2004 by Ullrich Koethe                  */
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

#ifndef VIGRA_UNIT_TEST_HPP
#define VIGRA_UNIT_TEST_HPP

#include <vector>
#include <string>
#include <new>                // for bad_alloc
#include <typeinfo>           // for bad_cast, bad_typeid
#include <exception>          // for exception, bad_exception
#include <stdexcept>

#if __GNUC__ >= 3  // why does this not work with MSVC 7.1 ???
#include <sstream>
#define VIGRA_SSTREAM std::basic_stringstream<char>
#define VIGRA_SSTREAM_STR(s) s.str().c_str()
#else
#include <strstream>
#define VIGRA_SSTREAM std::strstream
#define VIGRA_SSTREAM_STR(s) ((s << char()), s.str())
#endif

#include <iostream>
#include <limits.h>
#include <cfloat>
#include <cmath>
#include "vigra/config.hxx"
#include "vigra/error.hxx"

#ifdef _MSC_VER

#include <wtypes.h>
#include <winbase.h>
#include <excpt.h>

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif
#ifdef DIFFERENCE
#undef DIFFERENCE
#endif
#ifdef RGB
#undef RGB
#endif

#elif defined(__CYGWIN__)

#define VIGRA_CANT_CATCH_SIGNALS

#elif defined(__unix) || defined(unix)

#include <unistd.h>
#include <signal.h>
#include <setjmp.h>

#else

#define VIGRA_CANT_CATCH_SIGNALS

#endif

#define VIGRA_TEST_CASE(function)  vigra::create_test_case(function, #function "()")

#define testCase VIGRA_TEST_CASE

#define VIGRA_TEST_SUITE(testsuite)  ( new testsuite )

#define VIGRA_CHECKPOINT(message) \
    vigra::detail::checkpoint_impl(message, __FILE__, __LINE__)

#define VIGRA_ASSERT(predicate) \
    vigra::detail::should_impl((predicate), #predicate, __FILE__, __LINE__)

#define should VIGRA_ASSERT

#define VIGRA_ASSERT_MESSAGE(predicate, message) \
    vigra::detail::should_impl((predicate), message, __FILE__, __LINE__)

#define shouldMsg VIGRA_ASSERT_MESSAGE

#define shouldEqual(left, right) \
    vigra::detail::equal_impl(left, right, #left " == " #right, __FILE__, __LINE__)

#define shouldEqualTolerance(left, right, eps) \
    vigra::detail::tolerance_equal_impl(left, right, eps, #left " == " #right, __FILE__, __LINE__)

#define shouldEqualSequence(begin1, end1, begin2) \
    vigra::detail::sequence_equal_impl(begin1, end1, begin2, __FILE__, __LINE__)

#define shouldEqualSequenceTolerance(begin1, end1, begin2, eps) \
    vigra::detail::sequence_equal_tolerance_impl(begin1, end1, begin2, eps, __FILE__, __LINE__)

#define VIGRA_ERROR(message) \
    vigra::detail::should_impl(false, message, __FILE__, __LINE__)

#define failTest VIGRA_ERROR

namespace vigra {

class test_suite;

namespace detail {

struct errstream
{
    VIGRA_SSTREAM buf;
    char const * str() { return VIGRA_SSTREAM_STR(buf); }
    template <class T>
    errstream & operator<<(T t) { buf << t;  return *this; }
};

inline std::string & exception_checkpoint()
{
    static std::string test_checkpoint_;
    return test_checkpoint_;
}

//  A separate reporting function was requested during formal review.
inline void report_exception( detail::errstream & os,
                       const char * name, const char * info )
{
    os << "Unexpected " << name << " " << info << "\n";
    if(exception_checkpoint().size() > 0)
    {
        os << "Last checkpoint: " << exception_checkpoint().c_str() << "\n";
    }
}

enum {
    unexpected_exception = -1,
    os_exception = -2,
    memory_access_violation = -3,
    destructor_failure = -4
};

inline bool critical_error(int i)
{ return i <= memory_access_violation; }

inline bool unexpected_error(int i)
{ return i < 0; }

#ifndef VIGRA_CANT_CATCH_SIGNALS

#ifdef _MSC_VER

inline long handle_signal_here(long code)
{
    switch (code)
    {
        case EXCEPTION_ACCESS_VIOLATION:
        case EXCEPTION_INT_DIVIDE_BY_ZERO:
            return EXCEPTION_EXECUTE_HANDLER;
        default:
            return EXCEPTION_CONTINUE_SEARCH;
    }
}

template< class Generator >  // Generator is function object returning int
int catch_signals( Generator function_object, detail::errstream & err, int timeout )
{
    int result = 0;
    int code;
    __try
    {
        result = function_object();
    }
    __except (handle_signal_here(code = GetExceptionCode()))
    {
        switch (code)
        {
            case EXCEPTION_ACCESS_VIOLATION:
                report_exception(err, "operating system exception:", "memory access violation");
                result = memory_access_violation;
                break;
            case EXCEPTION_INT_DIVIDE_BY_ZERO:
                report_exception(err, "operating system exception:", "integer divide by zero");
                result = os_exception;
                break;
            default:
                report_exception(err, "operating system exception:", "unrecognized exception or signal");
                result = os_exception;
        }
    }
    return result;

}


#elif defined(__unix)

inline jmp_buf & unit_test_jump_buffer()
{
    static jmp_buf unit_test_jump_buffer_;
    return unit_test_jump_buffer_;
}

static void unit_test_signal_handler(int sig)
{
    longjmp(unit_test_jump_buffer(), sig);
}

template< class Generator >  // Generator is function object returning int
int catch_signals( Generator function_object, detail::errstream & err, int timeout)
{
    volatile int sigtype;
    int result;

#if defined(linux) || defined(__linux)
    signal(SIGFPE, &unit_test_signal_handler);
    signal(SIGTRAP, &unit_test_signal_handler);
    signal(SIGSEGV, &unit_test_signal_handler);
    signal(SIGBUS, &unit_test_signal_handler);
#else
    sigset(SIGFPE, &unit_test_signal_handler);
    sigset(SIGTRAP, &unit_test_signal_handler);
    sigset(SIGSEGV, &unit_test_signal_handler);
    sigset(SIGBUS, &unit_test_signal_handler);
#endif

    if(timeout)
    {
#if defined(linux) || defined(__linux)
        signal(SIGALRM, &unit_test_signal_handler);
#else
        sigset(SIGALRM, &unit_test_signal_handler);
#endif
        alarm(timeout);
    }

    sigtype = setjmp(unit_test_jump_buffer());
    if(sigtype == 0)
    {
       result = function_object();
    }
    else
    {
        switch(sigtype)
        {
            case SIGALRM:
                report_exception(err, "signal:", "SIGALRM (timeout while executing function)");
                result = os_exception;
                break;
            case SIGTRAP:
                report_exception(err, "signal:", "SIGTRAP (perhaps integer divide by zero)");
                result = os_exception;
                break;
            case SIGFPE:
                report_exception(err, "signal:", "SIGFPE (arithmetic exception)");
                result = os_exception;
                break;
            case SIGSEGV:
            case SIGBUS:
                report_exception(err, "signal:", "memory access violation");
                result = memory_access_violation;
                break;
            default:
                report_exception(err, "signal:", "unrecognized signal");
                result = os_exception;
        }
    }

    if(timeout)
    {
        alarm(0);
#if defined(linux) || defined(__linux)
#else
        sigrelse(SIGALRM);
#endif
    }

#if defined(linux) || defined(__linux)
#else
    sigrelse(SIGFPE);
    sigrelse(SIGTRAP);
    sigrelse(SIGSEGV);
    sigrelse(SIGBUS);
#endif

    return result;
}

#endif  /* _MSC_VER || __unix */

#else  /* VIGRA_CANT_CATCH_SIGNALS */

template< class Generator >  // Generator is function object returning int
int catch_signals( Generator function_object, detail::errstream & err , int)
{
     return function_object();
}

#endif /* VIGRA_CANT_CATCH_SIGNALS */

} // namespace detail

template< class Generator >  // Generator is function object returning int
int catch_exceptions( Generator function_object, detail::errstream & err, int timeout )
{
    int result = detail::unexpected_exception;

    try
    {
        result = detail::catch_signals(function_object, err, timeout);
    }

    //  As a result of hard experience with strangely interleaved output
    //  under some compilers, there is a lot of use of endl in the code below
    //  where a simple '\n' might appear to do.

    //  The rules for catch & arguments are a bit different from function
    //  arguments (ISO 15.3 paragraphs 18 & 19). Apparently const isn't
    //  required, but it doesn't hurt and some programmers ask for it.

    catch ( vigra::ContractViolation & ex )
      { detail::report_exception( err, "Contract exception: ", ex.what() ); }
    catch ( const char * ex )
      { detail::report_exception( err, "string exception: ", ex ); }
    catch ( const std::string & ex )
      { detail::report_exception( err, "string exception: ", ex.c_str() ); }

    //  std:: exceptions
    catch ( const std::bad_alloc & ex )
      { detail::report_exception( err, "exception: std::bad_alloc:", ex.what() ); }

# if !defined(__BORLANDC__) || __BORLANDC__ > 0x0551
    catch ( const std::bad_cast & ex )
      { detail::report_exception( err, "exception: std::bad_cast:", ex.what() ); }
    catch ( const std::bad_typeid & ex )
      { detail::report_exception( err, "exception: std::bad_typeid:", ex.what() ); }
# else
    catch ( const std::bad_cast & ex )
      { detail::report_exception( err, "exception: std::bad_cast", "" ); }
    catch ( const std::bad_typeid & ex )
      { detail::report_exception( err, "exception: std::bad_typeid", "" ); }
# endif

    catch ( const std::bad_exception & ex )
      { detail::report_exception( err, "exception: std::bad_exception:", ex.what() ); }
    catch ( const std::domain_error & ex )
      { detail::report_exception( err, "exception: std::domain_error:", ex.what() ); }
    catch ( const std::invalid_argument & ex )
      { detail::report_exception( err, "exception: std::invalid_argument:", ex.what() ); }
    catch ( const std::length_error & ex )
      { detail::report_exception( err, "exception: std::length_error:", ex.what() ); }
    catch ( const std::out_of_range & ex )
      { detail::report_exception( err, "exception: std::out_of_range:", ex.what() ); }
    catch ( const std::range_error & ex )
      { detail::report_exception( err, "exception: std::range_error:", ex.what() ); }
    catch ( const std::overflow_error & ex )
      { detail::report_exception( err, "exception: std::overflow_error:", ex.what() ); }
    catch ( const std::underflow_error & ex )
      { detail::report_exception( err, "exception: std::underflow_error:", ex.what() ); }
    catch ( const std::logic_error & ex )
      { detail::report_exception( err, "exception: std::logic_error:", ex.what() ); }
    catch ( const std::runtime_error & ex )
      { detail::report_exception( err, "exception: std::runtime_error:", ex.what() ); }
   catch ( const std::exception & ex )
      { detail::report_exception( err, "exception: std::exception:", ex.what() ); }

    catch ( ... )
      { throw;detail::report_exception( err, "unknown exception", "" ); }

    return result;
} // catch_exceptions

template< class Generator >  // Generator is function object returning int
inline
int catch_exceptions( Generator function_object, detail::errstream & err)
{
    return catch_exceptions(function_object, err, 0);
}

namespace detail {

struct unit_test_failed
: public std::exception
{
    unit_test_failed(char const * message)
    : what_(message)
    {}

    virtual ~unit_test_failed() throw()
    {
    }

    virtual const char * what() const throw()
    {
        return what_.c_str();
    }

    std::string what_;
};

inline void
checkpoint_impl(const char * message, const char * file, int line)
{
    detail::errstream buf;
    buf << message << " (" << file <<":" << line << ")";
    exception_checkpoint() = buf.str();
}

inline void
should_impl(bool predicate, const char * message, const char * file, int line)
{
    checkpoint_impl(message, file, line);
    if(!predicate)
    {
        detail::errstream buf;
        buf << message << " (" << file <<":" << line << ")";
        throw unit_test_failed(buf.str()); 
    } 
}

template <class Iter1, class Iter2>
void 
sequence_equal_impl(Iter1 i1, Iter1 end1, Iter2 i2, const char * file, int line)
{
    for(int counter = 0; i1 != end1; ++i1, ++i2, ++counter)
    {
        if(*i1 != *i2)
        {
            detail::errstream buf;
            buf << "Sequence items differ at index " << counter <<
                   " ["<< *i1 << " != " << *i2 << "]";
            should_impl(false, buf.str(), file, line); 
        }
    }
}

/******************Floating point comparison********************************/
/**
* See Knuth "The art of computer programming" (Vol II, Ch.4.2)
*/
struct ScalarType {};
struct VectorType {};

template<class T>
struct FloatTraits
{
    typedef VectorType ScalarOrVector;
};

template<>
struct FloatTraits<float>
{
    typedef ScalarType ScalarOrVector;
    static float epsilon() { return FLT_EPSILON; }
    static float smallestPositive() { return FLT_MIN; }
    static float min() { return -FLT_MAX; }
    static float max() { return FLT_MAX; }
};

template<>
struct FloatTraits<double>
{
    typedef ScalarType ScalarOrVector;
    static double epsilon() { return DBL_EPSILON; }
    static double smallestPositive() { return DBL_MIN; }
    static double min() { return -DBL_MAX; }
    static double max() { return DBL_MAX; }
};

template<>
struct FloatTraits<long double>
{
    typedef ScalarType ScalarOrVector;
    static long double epsilon() { return LDBL_EPSILON; }
    static long double smallestPositive() { return LDBL_MIN; }
    static long double min() { return -LDBL_MAX; }
    static long double max() { return LDBL_MAX; }
};

template<class FPT>
inline
FPT fpt_abs( FPT arg )
{
    return arg < 0 ? -arg : arg;
}


/***********************************************************************/

// both f1 and f2 are unsigned here
template<class FPT>
inline
FPT safe_fpt_division( FPT f1, FPT f2 )
{
        /*  ist f1 das absolute minimum (in diesem Fall einfach nur sehr kleine Zahl)
        *   aber nicht null (1.65242e-28) und f2 = 0,
        *   dann tritt die erste Bedingung in Kraft 0<1 && 1.65242e-28 > 0*1.79769e+308 (max)
        *   deshalb schlaegt es fehl sogar wenn min closed at tolarance zu 0 ist ???
        *   Der Vergleich aller Zahlen closed at tolarance zu 0 wuerden fehlschlagen;
        *   Sie umzudrehen bringt nichts, denn diese Funktion wird symetrisch fuer beide
        *   angewendet wird.
        *   0 mit 0 zu Vergleichen bereitet keine Probleme.
        *   Ausweg: evl. eine extra Behandlung der F = 0 ???
        */
    return  ((f2 < 1) && (f1 > (f2 *    FloatTraits<FPT>::max()))) ? 
            FloatTraits<FPT>::max() :
            (((f2 > 1) && (f1 < (f2 * FloatTraits<FPT>::smallestPositive())) || (f1 == 0)) ? 0 : f1/f2 );
        /*  Die Multiplikation mit max in 1.ten Bedingung und mit min in der 2.ten ist eine Absicherung gegen
        *   die Owerflow bzw Underflow ???
        */
}

/***********************************************************************/

template<class FPT>
class close_at_tolerance {
public:
    explicit    close_at_tolerance( FPT tolerance, bool strong_test = true )
        : m_strong_test( strong_test ),
          m_tolerance( tolerance ) {}

    explicit    close_at_tolerance( int number_of_rounding_errors, bool strong_test = true )
        : m_strong_test( strong_test ),
          m_tolerance( FloatTraits<FPT>::epsilon() * number_of_rounding_errors / 2.0 ) {}

    bool        operator()( FPT left, FPT right ) const
    {
        if (left == 0 && right != 0)
        {
                return (fpt_abs(right) <= m_tolerance);
        }
        if (right == 0 && left != 0)
        {
                return (fpt_abs(left) <= m_tolerance);
        }
        FPT diff = fpt_abs( left - right );
        FPT d1   = safe_fpt_division( diff, fpt_abs( right ) );
        FPT d2   = safe_fpt_division( diff, fpt_abs( left ) );

        return m_strong_test ? (d1 <= m_tolerance && d2 <= m_tolerance)
                                : (d1 <= m_tolerance || d2 <= m_tolerance);
    }

private:
    bool        m_strong_test;
    FPT         m_tolerance;
};

/*****************end of float comparison***********************************/

template <class T1, class T2, class T3>
void
tolerance_equal_impl(T1 left, T2 right, T3 epsilon, 
                     const char * message, const char * file, int line, ScalarType)
{
    detail::errstream buf;
    buf << message << " [" << left << " != " << right << "]";

    close_at_tolerance<T3> fcomparator( epsilon );
    bool compare = fcomparator ( left , right );
    should_impl(compare, buf.str(), file, line);

}

template <class T1, class T2, class T3>
void
tolerance_equal_impl(T1 left, T2 right, T3 epsilon, 
                     const char * message, const char * file, int line, VectorType)
{
    detail::errstream buf;
    buf << message << " [" << left << " != " << right << "]";

    bool compare = true;
    for(unsigned int i=0; i<epsilon.size(); ++i)
    {
        close_at_tolerance<typename T3::value_type> fcomparator( epsilon[i] );
        compare = compare && fcomparator ( left[i] , right[i] );
    }
    should_impl(compare, buf.str(), file, line);
}

template <class T1, class T2, class T3>
void
tolerance_equal_impl(T1 left, T2 right, T3 epsilon, const char * message, const char * file, int line)
{
    tolerance_equal_impl(left, right, epsilon, 
                         message, file, line, typename FloatTraits<T3>::ScalarOrVector());
}

template <class Iter1, class Iter2, class T>
void 
sequence_equal_tolerance_impl(Iter1 i1, Iter1 end1, Iter2 i2, T epsilon, const char * file, int line)
{
    for(int counter = 0; i1 != end1; ++i1, ++i2, ++counter)
    {
        detail::errstream buf;
        buf << "Sequence items differ at index " << counter;
        tolerance_equal_impl(*i1, *i2, epsilon, buf.str(), file, line); 
    }
}

template <class Left, class Right>
void
equal_impl(Left left, Right right, const char * message, const char * file, int line)
{
    detail::errstream buf;
    buf << message << " [" << left << " != " << right << "]";
    should_impl(left == right, buf.str(), file, line);
}

inline void
equal_impl(double left, double right, const char * message, const char * file, int line)
{
    tolerance_equal_impl(left, right, 1.0e-16, message, file, line);
}

inline void
equal_impl(float left, float right, const char * message, const char * file, int line)
{
    tolerance_equal_impl(left, right, 1.0e-6f, message, file, line);
}

inline void
equal_impl(float left, double right, const char * message, const char * file, int line)
{
    tolerance_equal_impl(left, right, 1.0e-6f, message, file, line);
}

inline void
equal_impl(double left, float right, const char * message, const char * file, int line)
{
    tolerance_equal_impl(left, right, 1.0e-6f, message, file, line);
}

class test_case
{
  public:

    test_case(char const * name = "Unnamed")
    : name_(name), timeout(0)
    {}

    virtual ~test_case() {}

    virtual int run() = 0;
    virtual void do_init() {}
    virtual void do_run() {}
    virtual void do_destroy() {}

    virtual char const * name() { return name_.c_str(); }
    virtual int size() const { return 1; }

    std::string name_;
    std::string report_;
    int timeout;
};


} // namespace detail

class test_suite
: public detail::test_case
{
  public:
    test_suite(char const * name = "TopLevel")
    : detail::test_case(name),
      size_(0)
    {}

    virtual ~test_suite()
    {
        for(unsigned int i=0; i != testcases_.size(); ++i)
            delete testcases_[i];
    }

    virtual void add(detail::test_case * t, int timeout = 0)
    {
        t->timeout = timeout;
        testcases_.push_back(t);
        size_ += t->size();
    }

    virtual int run()
    {
        int failed = 0;
        report_ = std::string("Entering test suite ") + name() + "\n";

        for(unsigned int i=0; i != testcases_.size(); ++i)
        {
            int result = testcases_[i]->run();
            report_ += testcases_[i]->report_;

            if(detail::critical_error(result))
            {
                report_ += std::string("\nFatal error - aborting test suite ") + name() + ".\n";
                return result;
            }
            else if(detail::unexpected_error(result))
                failed++;
            else
                failed += result;
        }

        if(failed)
        {
            detail::errstream buf;
            buf << "\n" << failed << " of " << size() <<
                " tests failed in test suite " << name() << "\n";
            report_ += buf.str();
        }
        else
        {
            detail::errstream buf;
            buf << "All (" << size() <<
               ") tests passed in test suite " << name() << "\n";
            report_ += buf.str();
        }

        report_ += std::string("Leaving test suite ") + name() + "\n";

        return failed;
    }

    virtual int size() const { return size_; }
    virtual std::string report() { return report_; }

    std::vector<detail::test_case *> testcases_;
    int size_;
};

namespace detail {

struct test_case_init_functor
{
    detail::errstream & buf_;
    test_case * test_case_;

    test_case_init_functor(detail::errstream & b, test_case * tc)
    : buf_(b), test_case_(tc)
    {}

    int operator()()
    {
        try
        {
            test_case_->do_init();
            return 0;
        }
        catch(unit_test_failed & e)
        {
            buf_ << "Assertion failed: " << e.what() << "\n";
            return 1;
        }
    }
};

struct test_case_run_functor
{
    detail::errstream & buf_;
    test_case * test_case_;

    test_case_run_functor(detail::errstream & b, test_case * tc)
    : buf_(b), test_case_(tc)
    {}

    int operator()()
    {
        try
        {
            test_case_->do_run();
            return 0;
        }
        catch(unit_test_failed & e)
        {
            buf_ << "Assertion failed: " << e.what() << "\n";
            return 1;
        }
    }
};

struct test_case_destroy_functor
{
    detail::errstream & buf_;
    test_case *  test_case_;

    test_case_destroy_functor(detail::errstream & b, test_case * tc)
    : buf_(b), test_case_(tc)
    {}

    int operator()()
    {
        try
        {
            test_case_->do_destroy();
            return 0;
        }
        catch(unit_test_failed & e)
        {
            buf_ << "Assertion failed: " << e.what() << "\n";
            return 1;
        }
    }
};

template <class TESTCASE>
class class_test_case
: public test_case
{
  public:

    class_test_case(void (TESTCASE::*fct)(), char const * name)
    : test_case(name),
      fct_(fct),
      testcase_(0)
    {}

    virtual ~class_test_case()
    {
        delete testcase_;
    }

    virtual void do_init()
    {
        testcase_ = new TESTCASE;
    }

    int init()
    {
        exception_checkpoint() = "";
        report_ = "";
        int failed = 0;

        detail::errstream buf;
        buf << "\nFailure in initialization of " << name() << "\n";
        if(testcase_ != 0)
        {
            buf << "Test case failed to clean up after previous run.\n";
            failed = 1;
        }
        else
        {
            failed = catch_exceptions(
                detail::test_case_init_functor(buf, this), buf, timeout);
        }

        if(failed)
        {
            report_ += buf.str();
        }

        return failed;
    }

    virtual void do_run()
    {
        if(testcase_ != 0)
            (testcase_->*fct_)();
    }

    virtual int run()
    {
        int failed = init();

        if(failed)
            return failed;

        detail::errstream buf;
        buf << "\nFailure in " << name() << "\n";

        failed = catch_exceptions(
            detail::test_case_run_functor(buf, this), buf, timeout);
        if(failed)
            report_ += buf.str();

        if(critical_error(failed))
            return failed;

        int destruction_failed = destroy();

        return destruction_failed ?
                destruction_failed :
                failed;
    }

    virtual void do_destroy()
    {
        delete testcase_;
        testcase_ = 0;
    }

    int destroy()
    {
        detail::errstream buf;
        buf << "\nFailure in destruction of " << "\n";

        int failed = catch_exceptions(
            detail::test_case_destroy_functor(buf, this), buf, timeout);
        if(failed)
        {
            report_ += buf.str();
            return destructor_failure;
        }
        else
        {
            return 0;
        }
    }

    void (TESTCASE::*fct_)();
    TESTCASE * testcase_;
};

class function_test_case
: public test_case
{
  public:

    function_test_case(void (*fct)(), char const * name)
    : test_case(name),
      fct_(fct)
    {}

    virtual void do_run()
    {
        (*fct_)();
    }

    virtual int run()
    {
        report_ = "";
        exception_checkpoint() = "";

        detail::errstream buf;
        buf << "\nFailure in " << name() << "\n";

        int failed = catch_exceptions(
            detail::test_case_run_functor(buf, this), buf, timeout);
        if(failed)
        {
            report_ += buf.str();
        }

        return failed;
    }

    void (*fct_)();
};

template <class FCT>
struct test_functor
{
    virtual ~test_functor() {}
    virtual void operator()() = 0;
    
    FCT clone() const 
        { return FCT(static_cast<FCT const &>(*this)); }
};
    
template <class FCT>
class functor_test_case
: public test_case
{
  public:
    
    functor_test_case(FCT const & fct, char const * name)
    : test_case(name),
      fct_(fct)
    {}
    
    virtual void do_run()
    {
        fct_();
    }
    
    virtual int run()
    {
        report_ = "";
        exception_checkpoint() = "";
        
        detail::errstream buf;
        buf << "\nFailure in " << name() << "\n";
        
        int failed = catch_exceptions(
            detail::test_case_run_functor(buf, this), buf, timeout);
        if(failed)
        {
            report_ += buf.str();
        }

        return failed;
    }
    
    FCT fct_;
};
    
} // namespace detail

template <class TESTCASE>
inline
detail::test_case *
create_test_case(void (TESTCASE::*fct)(), char const * name)
{
    if(*name == '&') ++name;
    return new detail::class_test_case<TESTCASE>(fct, name);
}

inline
detail::test_case *
create_test_case(void (*fct)(), char const * name)
{
    if(*name == '&') ++name;
    return new detail::function_test_case(fct, name);
}

template <class FCT>
inline 
detail::test_case * 
create_test_case(detail::test_functor<FCT> const & fct, char const * name)
{
    if(*name == '&') ++name;
    return new detail::functor_test_case<FCT>(fct.clone(), name);
}

} // namespace vigra


#if !defined(__GNUC__) || __GNUC__ >= 3

// provide more convenient output functions, used like:
// std::cerr << 1, 2, 3, 4, "\n";
template <class E, class T, class V>
inline
std::basic_ostream<E,T> & operator,(std::basic_ostream<E,T> & o, V const & t)
{
    return (o << ' ' << t);
}

template <class E, class T>
inline
std::basic_ostream<E,T> & operator,(std::basic_ostream<E,T> & o, 
              std::basic_ostream<E,T> & (*t)(std::basic_ostream<E,T> &))
{
    return (o << t);
}

#else

template <class V>
inline
std::ostream & operator,(std::ostream & o, V const & t)
{
    return (o << ' ' << t);
}

inline
std::ostream & operator,(std::ostream & o, 
              std::ostream & (*t)(std::ostream &))
{
    return (o << t);
}

#endif


#endif /* VIGRA_UNIT_TEST_HPP */
