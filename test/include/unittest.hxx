/**************************************************************/
/*                                                            */
/*                    Copyright U. Koethe                     */
/*    Fraunhoferinstitut fuer Graphische Datenverarbeitung    */
/*                     Rostock, Germany                       */
/*               Wed Jun 24 15:11:45 MET 1998                 */
/*                                                            */
/**************************************************************/
 
#ifndef UNITTEST_HXX
#define UNITTEST_HXX

#include <stdio.h>
#include <list>
#include <string>
#include <stdexcept>

#ifdef _MSC_VER

#include <wtypes.h>
#include <winbase.h>
#include <excpt.h>
#define snprintf _snprintf

#elif defined(__unix)

#include <signal.h>
#include <setjmp.h>

#else

#define CANT_CATCH_SIGNALS

#endif

#define testCase(function)  createTestCase(function, #function "()")
#define testSuite(testsuite)  ( new testsuite )

#define checkpoint(message) \
    checkpointImpl(message, __FILE__, __LINE__)

#define should(predicate) \
    shouldImpl((predicate), #predicate, __FILE__, __LINE__)

#define shouldMsg(predicate, message) \
    shouldImpl((predicate), message, __FILE__, __LINE__)
    
#define failTest(message) \
    shouldImpl(false, message, __FILE__, __LINE__)

struct UnitTestFailed 
{
    UnitTestFailed(char * message)
    : what_(message)
    {}
    
    virtual const char * what() const
    {
        return what_.c_str();
    }
    
    std::string what_;
};

inline std::string & testCheckpoint()
{
    static std::string testCheckpoint_;
    return testCheckpoint_;
}

inline void 
checkpointImpl(const char * message, const char * file, int line)
{
    char buf[1000]; 
    snprintf(buf, 1000, "\"%s\" (%s:%d)", message, file, line); 
    testCheckpoint() = buf;
}

inline void 
shouldImpl(bool predicate, const char * message, const char * file, int line)
{
    checkpointImpl(message, file, line);
    if(!predicate)
    { 
        char buf[1000]; 
        snprintf(buf, 1000, "\"%s\" (%s:%d)", message, file, line); 
        throw UnitTestFailed(buf); 
    } 
}

class TestCase
{
  public:

    TestCase(char const * name = "Unnamed")
    : name_(name)
    {}

    virtual ~TestCase() {}
    
    virtual int run() = 0;
    virtual void init() {}
    virtual void destroy() {}
    virtual bool hasCheckpoint() const 
    {
        return testCheckpoint().size() > 0;
    }
    
    virtual std::string getCheckpointMessage() const 
    {
        return testCheckpoint();
    }

    virtual char const * name() { return name_.c_str(); } 
    virtual void setName(char const * name) { name_ = name; }
    virtual int size() const { return 1; }
    virtual std::string report() { return std::string(""); }

    std::string name_;
    
};

class TestSuite;

static int doRun(TestSuite *, TestCase *);

class TestSuite
: public TestCase
{
  public:
    TestSuite(char const * name = "TopLevel")
    : TestCase(name),
      size_(0)
    {}
    
    virtual ~TestSuite()
    {
        std::list<TestCase *>::iterator i = testcases_.begin();
        for(; i != testcases_.end(); ++i) delete (*i);        
    }
    
    virtual void add(TestCase * t)
    {
        testcases_.push_back(t);
        size_ += t->size();
    }
   
    virtual int run()
    {
        int failed = 0;
        
        std::list<TestCase *>::iterator i = testcases_.begin();
        for(; i != testcases_.end(); ++i) 
        {
            failed += doRun(this, *i);
            report_ += (*i)->report();
        }
        
        char buf[100];
        if(failed) 
        {
            sprintf(buf, "%d of %d", failed, size());
            report_ += std::string("\n") + buf +
                       " tests failed in TestSuite " + name() + "\n";
        }
        else 
        {
            sprintf(buf, "%d", size());
            report_ += std::string("\nAll (") + buf +
                       ") tests passed in TestSuite " + name() + "\n";
        }

        return failed;
    }
    
    virtual void handleError(const char * where, TestCase * t)
    {
        try
        {
            throw;
        }
        catch(UnitTestFailed & e)
        {
            report_ += std::string("\nFailure in ") + where + t->name() + 
                         " - assertion failed: " + e.what() + "\n";
        }
        catch(std::exception & e)
        {
            report_ += std::string("\nFailure in ") + where + t->name() + 
                         " - unexpected exception: " + e.what() + "\n";

            if(t->hasCheckpoint())
            {
                report_ += "Last checkpoint: " + 
                                t->getCheckpointMessage() + "\n";
            }
        }
    }
    
    virtual int tryInit(TestCase * t)
    {
        try
        {
            t->init();
        }
        catch(...)
        {
            handleError("initialization of ", t);
            return 1;
        }

        return 0;
    }
    
    virtual int tryRun(TestCase * t)
    {
        int res = 0;
        
        try
        {
            res = t->run();
        }
        catch(...)
        {
            handleError("", t);
            return 1;
        }

        return res;
    }
    
    virtual int tryDestroy(TestCase * t)
    {
        try
        {
            t->destroy();
        }
        catch(...)
        {
            handleError("destruction of ", t);
            return 1;
        }

        return 0;
    }
    
    virtual int size() const { return size_; }
    virtual std::string report() { return report_; }
  
    std::list<TestCase *> testcases_;
    int size_;
    std::string report_;
};

static void handleUnrecognizedException(const char * where, TestSuite * ts, TestCase * tc)
{
    ts->report_ += std::string("\nFailure in ") + where + tc->name() + 
                               " - unrecognized exception\n";

    if(tc->hasCheckpoint())
    {
        ts->report_ += "Last checkpoint: " + 
                        tc->getCheckpointMessage() + "\n";
    }
}

#ifndef CANT_CATCH_SIGNALS

#ifdef _MSC_VER

static void handleSignal(const char * where, TestSuite * ts, TestCase * tc)
{
    ts->report_ += std::string("\nFailure in ") + where + tc->name();

    switch (GetExceptionCode())
    {
        case EXCEPTION_ACCESS_VIOLATION:
            ts->report_ += " - memory access violation\n";
            break;
        case EXCEPTION_INT_DIVIDE_BY_ZERO:
            ts->report_ += " - integer divide by zero\n";
            break;
        default:
            ts->report_ += " - unrecognized exception or signal\n";
    }

    if(tc->hasCheckpoint())
    {
        ts->report_ += "Last checkpoint: " + 
                        tc->getCheckpointMessage() + "\n";
    }
}
   
static int doRun(TestSuite * ts, TestCase * tc)
{
    int current_failed = 0;
    
    __try 
    {
        current_failed += ts->tryInit(tc);
    }
    __except (true) 
    {
        handleSignal("initialization of ", ts, tc);
        current_failed++;
    }
   
    __try 
    {
        if(current_failed == 0) current_failed += ts->tryRun(tc);
    }
    __except (true) 
    {
        handleSignal("", ts, tc);
        current_failed++;
    }
    
    __try 
    {
        if(current_failed == 0) current_failed += ts->tryDestroy(tc);
    }
    __except (true) 
    {
        handleSignal("destruction of ", ts, tc);
        current_failed++;
    }

    return current_failed;
}

#elif defined(__unix)

inline jmp_buf & unitTestJumpBuffer()
{
    static jmp_buf unitTestJumpBuffer_;
    return unitTestJumpBuffer_;
}

static void unitTestSignalHandler(int sig, int code, struct sigcontext * c)
{
    longjmp(unitTestJumpBuffer(), sig);
}

static void handleSignal(int sigtype, const char * where, TestSuite * ts, TestCase * tc)
{
    ts->report_ += std::string("\nFailure in ") + where + tc->name();

    switch(sigtype)
    {
        case SIGTRAP:
            ts->report_ += " - SIGTRAP (perhaps integer divide by zero)\n";
            break;
        case SIGFPE:
            ts->report_ += " - SIGFPE (arithmetic exception)\n";
            break;
        case SIGSEGV:
        case SIGBUS:
            ts->report_ += " - memory access violation\n";
            break;
        default:
            ts->report_ += " - unrecognized signal\n";
    }

    if(tc->hasCheckpoint())
    {
        ts->report_ += "Last checkpoint: " + 
                        tc->getCheckpointMessage() + "\n";
    }
}
   
static int doRun(TestSuite * ts, TestCase * tc)
{
    int current_failed = 0;
    int sigtype;


    sigset(SIGFPE, (SIG_PF)&unitTestSignalHandler); 
    sigset(SIGTRAP, (SIG_PF)&unitTestSignalHandler); 
    sigset(SIGSEGV, (SIG_PF)&unitTestSignalHandler); 
    sigset(SIGBUS, (SIG_PF)&unitTestSignalHandler); 

    try 
    {
        sigtype = setjmp(unitTestJumpBuffer());
        if(sigtype == 0)
        {
           current_failed += ts->tryInit(tc);
        }
        else
        {
            handleSignal(sigtype, "initialization of ", ts, tc);
            ++current_failed;
        }
    }
    catch(...) 
    {
        handleUnrecognizedException("initialization of ", ts, tc);
        ++current_failed;
    }
    
    if(current_failed) return current_failed;
   
    try 
    {
        sigtype = setjmp(unitTestJumpBuffer());
        if(sigtype == 0)
        {
           current_failed += ts->tryRun(tc);
        }
        else
        {
            handleSignal(sigtype, "", ts, tc);
            ++current_failed;
        }
    }
    catch(...) 
    {
        handleUnrecognizedException("", ts, tc);
        ++current_failed;
    }
    
    if(current_failed) return current_failed;
   
    try 
    {
        sigtype = setjmp(unitTestJumpBuffer());
        if(sigtype == 0)
        {
           current_failed += ts->tryDestroy(tc);
        }
        else
        {
            handleSignal(sigtype, "destruction of ", ts, tc);
            ++current_failed;
        }
    }
    catch(...) 
    {
        handleUnrecognizedException("destruction of ", ts, tc);
        ++current_failed;
    }

    return current_failed;
}

#endif  /* _MSC_VER || __unix */

#else  /* CANT_CATCH_SIGNALS */

static int doRun(TestSuite * ts, TestCase * tc)
{
    int current_failed = 0;

    try 
    {
        current_failed += ts->tryInit(tc);
    }
    catch(...) 
    {
        handleUnrecognizedException("initialization of ", ts, tc);
        ++current_failed;
    }
    
    if(current_failed) return current_failed;
   
    try 
    {
        current_failed += ts->tryRun(tc);
    }
    catch(...) 
    {
        handleUnrecognizedException("", ts, tc);
        ++current_failed;
    }
    
    if(current_failed) return current_failed;
   
    try 
    {
        current_failed += ts->tryDestroy(tc);
    }
    catch(...) 
    {
        handleUnrecognizedException("destruction of ", ts, tc);
        ++current_failed;
    }

    return current_failed;
}

#endif /* CANT_CATCH_SIGNALS */

template <class TESTCASE>
class ClassTestCase
: public TestCase
{
  public:
    
    ClassTestCase(void (TESTCASE::*fct)(), char const * name)
    : TestCase(name),
      fct_(fct),
      testcase_(0)
    {}
    
    virtual ~ClassTestCase()
    {
        delete testcase_;
    }
    
    virtual void init()
    {
        if(testcase_ != 0)
        {
            shouldMsg(0, "Attempt to run test case which failed "
                         "to clean up after previous run.");
        }
        
        testCheckpoint() = "";
        testcase_ = new TESTCASE;
    }
    
    virtual int run()
    {
        testCheckpoint() = "";
        if(testcase_ != 0) (testcase_->*fct_)();
        
        return 0;
    }
    
    virtual void destroy()
    {
        testCheckpoint() = "";
        TESTCASE * toDelete = testcase_;
        testcase_ = 0;
        delete toDelete; 
    }
    
    void (TESTCASE::*fct_)();
    TESTCASE * testcase_;
};
    
class FunctionTestCase
: public TestCase
{
  public:
    
    FunctionTestCase(void (*fct)(), char const * name)
    : TestCase(name),
      fct_(fct)
    {}
    
    virtual int run()
    {
        testCheckpoint() = "";
        (*fct_)();
        
        return 0;
    }
    
    void (*fct_)();
};
    
template <class TESTCASE>
inline 
TestCase * 
createTestCase(void (TESTCASE::*fct)(), char const * name)
{
    if(*name == '&') ++name;
    return new ClassTestCase<TESTCASE>(fct, name);
};

inline 
TestCase * 
createTestCase(void (*fct)(), char const * name)
{
    if(*name == '&') ++name;
    return new FunctionTestCase(fct, name);
};


#endif /* UNITTEST_HXX */
