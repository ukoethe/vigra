/**********************************************************************/
/*                                                                    */
/* a simple unit test framework, similar to Kent Beck's JUnit         */
/*                                                                    */
/*              (C) Copyright Ullrich Koethe 2001                     */
/* Permission to copy, use, modify, sell and distribute this software */
/* is granted provided this copyright notice appears in all copies.   */
/* This software is provided "as is" without express or implied       */
/* warranty, and with no claim as to its suitability for any purpose. */
/*                                                                    */
/**********************************************************************/
 
#ifndef VIGRA_UNIT_TEST_HPP
#define VIGRA_UNIT_TEST_HPP

#include <vector>
#include <string>
#include <new>                // for bad_alloc
#include <typeinfo>           // for bad_cast, bad_typeid
#include <exception>          // for exception, bad_exception
#include <stdexcept>
#include <strstream>
#include <iostream>
#include <cmath>
#include "vigra/config.hxx"

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

#elif defined(__unix)

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

#define shouldEqualEps(left, right, eps) \
    vigra::detail::eps_equal_impl(left, right, eps, #left " == " #right, __FILE__, __LINE__)

#define VIGRA_ERROR(message) \
    vigra::detail::should_impl(false, message, __FILE__, __LINE__)

#define failTest VIGRA_ERROR

namespace vigra {

class test_suite;

namespace detail {

struct errstream
{
    std::strstream buf;
    char const * str() { buf << char(); return buf.str(); }
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

template< class Generator >  // Generator is function object returning int
int catch_signals( Generator function_object, detail::errstream & err, int timeout )
{
    int result = 0;
    __try 
    {
        result = function_object();
    }
    __except (true) 
    {
        switch (GetExceptionCode())
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

    sigset(SIGFPE, &unit_test_signal_handler); 
    sigset(SIGTRAP, &unit_test_signal_handler); 
    sigset(SIGSEGV, &unit_test_signal_handler); 
    sigset(SIGBUS, &unit_test_signal_handler); 
    if(timeout)
    {
        sigset(SIGALRM, &unit_test_signal_handler); 
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
        sigrelse(SIGALRM); 
    }
    sigrelse(SIGFPE); 
    sigrelse(SIGTRAP); 
    sigrelse(SIGSEGV); 
    sigrelse(SIGBUS); 

    return result;
}

#endif  /* _MSC_VER || __unix */

#else  /* VIGRA_CANT_CATCH_SIGNALS */

template< class Generator >  // Generator is function object returning int
int catch_signals( Generator function_object, detail::errstream & err )
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
      { detail::report_exception( err, "unknown exception", "" ); }

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
    
    virtual const char * what() const
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

inline void 
eps_equal_impl(double left, double right, double epsilon, const char * message, const char * file, int line)
{
    detail::errstream buf;
    buf << message << " [" << left << " != " << right << "]";
    should_impl(VIGRA_CSTD::fabs(left - right) < epsilon, buf.str(), file, line); 
}

template <class Left, class Right>
inline void 
equal_impl(Left left, Right right, const char * message, const char * file, int line)
{
    detail::errstream buf;
    buf << message << " [" << left << " != " << right << "]";
    should_impl(left == right, buf.str(), file, line); 
}

inline void 
equal_impl(double left, double right, const char * message, const char * file, int line)
{
    eps_equal_impl(left, right, 1.0e-16, message, file, line); 
}

inline void 
equal_impl(float left, float right, const char * message, const char * file, int line)
{
    eps_equal_impl(left, right, 1.0e-6, message, file, line); 
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
    
} // namespace detail

template <class TESTCASE>
inline 
detail::test_case * 
create_test_case(void (TESTCASE::*fct)(), char const * name)
{
    if(*name == '&') ++name;
    return new detail::class_test_case<TESTCASE>(fct, name);
};

inline 
detail::test_case * 
create_test_case(void (*fct)(), char const * name)
{
    if(*name == '&') ++name;
    return new detail::function_test_case(fct, name);
};

} // namespace vigra

#endif /* VIGRA_UNIT_TEST_HPP */
