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

#elif defined(__unix)

#include <signal.h>
#include <setjmp.h>

#else

#error unsupported system, sorry

#endif

#include <vigra/config.hxx>


#define testCase(function)  createTestCase(function, #function "()")
#define testSuite(testsuite)  ( new testsuite )

#define checkpoint(message) \
    { \
        char buf[1000]; \
        sprintf(buf, "%s (" __FILE__ ":%d)", (message), __LINE__); \
        testCheckpoint_ = buf; \
    } 

#define should(p) \
    {\
        checkpoint(#p) \
        if(p); else\
        { \
            char buf[1000]; \
            sprintf(buf, "\"" #p "\" (" __FILE__ ":%d)", __LINE__); \
            throw UnitTestFailed(buf); \
        } \
    }

#define shouldMsg(p, message) \
    {\
        checkpoint(message) \
        if(p); else \
        { \
            char buf[1000]; \
            sprintf(buf, "%s (" __FILE__ ":%d)", (message), __LINE__); \
            throw UnitTestFailed(buf); \
        } \
    } 

#define failTest(message) shouldMsg(false, message)

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

static std::string testCheckpoint_;

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
        return testCheckpoint_.size() > 0;
    }
    
    virtual std::string getCheckpointMessage() const 
    {
        return testCheckpoint_;
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
    
    virtual int tryInit(TestCase * t)
    {
        try
        {
            t->init();
        }
        catch(UnitTestFailed & e)
        {
            report_ += std::string("\nFailure in initialization of ") + t->name() + 
                         " - assertion failed: " + e.what() + "\n";
            return 1;
        }
        catch(VigraStdException & e)
        {
            report_ += std::string("\nFailure in initialization of ") + t->name() + 
                         " - unexpected exception: " + e.what() + "\n";

            if(t->hasCheckpoint())
            {
                report_ += "Last checkpoint: " + 
                                t->getCheckpointMessage() + "\n";
            }
        
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
        catch(UnitTestFailed & e)
        {
            report_ += std::string("\nFailure in ") + t->name() + 
                         " - assertion failed: " + e.what() + "\n";
            return 1;
        }
        catch(VigraStdException & e)
        {
            report_ += std::string("\nFailure in ") + t->name() + 
                         " - unexpected exception: " + e.what() + "\n";
            if(t->hasCheckpoint())
            {
                report_ += "Last checkpoint: " + 
                                t->getCheckpointMessage() + "\n";
            }
        
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
        catch(UnitTestFailed & e)
        {
            report_ += std::string("\nFailure in destruction of ") + t->name() + 
                         " - assertion failed: " + e.what() + "\n";
            return 1;
        }
        catch(VigraStdException & e)
        {
            report_ += std::string("\nFailure in destruction of ") + t->name() + 
                         " - unexpected exception: " + e.what() + "\n";
            if(t->hasCheckpoint())
            {
                report_ += "Last checkpoint: " + 
                                t->getCheckpointMessage() + "\n";
            }
        
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

#ifdef _MSC_VER
   
static int doRun(TestSuite * ts, TestCase * tc)
{
    int current_failed = 0;
    
    __try 
    {
        current_failed += ts->tryInit(tc);
    }
    __except (true) 
    {
        ts->report_ += "\nFailure in initialization of ";
        ts->report_ += tc->name();

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
        current_failed++;
    }
   
    __try 
    {
        if(current_failed == 0) current_failed += ts->tryRun(tc);
    }
    __except (true) 
    {
        ts->report_ += "\nFailure in ";
        ts->report_ += tc->name();
        
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
        current_failed++;
    }
    
    __try 
    {
        if(current_failed == 0) current_failed += ts->tryDestroy(tc);
    }
    __except (true) 
    {
        ts->report_ += "\nFailure in destruction of ";
        ts->report_ += tc->name();
        
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
        current_failed++;
    }

    return current_failed;
}

#elif defined(__unix)

static jmp_buf unitTestJumpBuffer;

static void unitTestSignalHandler(int sig, int code, struct sigcontext * c)
{
    longjmp(unitTestJumpBuffer, sig);
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
        sigtype = setjmp(unitTestJumpBuffer);
        if(sigtype == 0)
        {
           current_failed += ts->tryInit(tc);
        }
        else
        {
            ts->report_ += "\nFailure in initialization of ";
            ts->report_ += tc->name();
            
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
            }
            
            if(tc->hasCheckpoint())
            {
                ts->report_ += "Last checkpoint: " + 
                                tc->getCheckpointMessage() + "\n";
            }
        
            ++current_failed;
        }
    }
    catch(...) 
    {
        ts->report_ += "\nFailure in initialization of ";
        ts->report_ += tc->name();
        ts->report_ += " - unrecognized exception\n";
        
        if(tc->hasCheckpoint())
        {
            ts->report_ += "Last checkpoint: " + 
                            tc->getCheckpointMessage() + "\n";
        }
        
        ++current_failed;
    }
    
    if(current_failed) return current_failed;
   
    try 
    {
        sigtype = setjmp(unitTestJumpBuffer);
        if(sigtype == 0)
        {
           current_failed += ts->tryRun(tc);
        }
        else
        {
            ts->report_ += "\nFailure in ";
            ts->report_ += tc->name();
            
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
            }
        
            if(tc->hasCheckpoint())
            {
                ts->report_ += "Last checkpoint: " + 
                                tc->getCheckpointMessage() + "\n";
            }
        
            ++current_failed;
        }
    }
    catch(...) 
    {
        ts->report_ += "\nFailure in ";
        ts->report_ += tc->name();
        ts->report_ += " - unrecognized exception\n";
        
        if(tc->hasCheckpoint())
        {
            ts->report_ += "Last checkpoint: " + 
                            tc->getCheckpointMessage() + "\n";
        }
        
        ++current_failed;
    }
    
    if(current_failed) return current_failed;
   
    try 
    {
        sigtype = setjmp(unitTestJumpBuffer);
        if(sigtype == 0)
        {
           current_failed += ts->tryDestroy(tc);
        }
        else
        {
            ts->report_ += "\nFailure in destruction of ";
            ts->report_ += tc->name();
            
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
            }
        
            if(tc->hasCheckpoint())
            {
                ts->report_ += "Last checkpoint: " + 
                                tc->getCheckpointMessage() + "\n";
            }
        
            ++current_failed;
        }
    }
    catch(...) 
    {
        ts->report_ += "\nFailure in destruction of ";
        ts->report_ += tc->name();
        ts->report_ += " - unrecognized exception\n";
        
        if(tc->hasCheckpoint())
        {
            ts->report_ += "Last checkpoint: " + 
                            tc->getCheckpointMessage() + "\n";
        }
        
        ++current_failed;
    }

    return current_failed;
}

#endif

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
        
        testCheckpoint_ = "";
        testcase_ = new TESTCASE;
    }
    
    virtual int run()
    {
        testCheckpoint_ = "";
        if(testcase_ != 0) (testcase_->*fct_)();
        
        return 0;
    }
    
    virtual void destroy()
    {
        testCheckpoint_ = "";
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
        testCheckpoint_ = "";
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
