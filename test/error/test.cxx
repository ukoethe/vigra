#include <iostream>

#ifdef NDEBUG
#undef NDEBUG
#endif

#include "unittest.hxx"
#include "vigra/error.hxx"

struct ErrorTest
{
    void testPrecondition()
    {
        try
        {
            vigra_precondition(0, "Intentional error");
            failTest("no exception thrown");
        }
        catch(vigra::ContractViolation & c)
        {
            std::string expected("\nPrecondition violation!\nIntentional error");
            std::string message(c.what());
            should(0 == expected.compare(message.substr(0,expected.size())));
        }
        try
        {
            vigra_assert(0, "Intentional error");
            failTest("no exception thrown");
        }
        catch(vigra::ContractViolation & c)
        {
            std::string expected("\nPrecondition violation!\nIntentional error");
            std::string message(c.what());
            should(0 == expected.compare(message.substr(0,expected.size())));
        }
    }

    void testPostcondition()
    {
        try
        {
            vigra_postcondition(0, "Intentional error");
            failTest("no exception thrown");
        }
        catch(vigra::ContractViolation & c)
        {
            std::string expected("\nPostcondition violation!\nIntentional error");
            std::string message(c.what());
            should(0 == expected.compare(message.substr(0,expected.size())));
        }
    }

    void testInvariant()
    {
        try
        {
            vigra_invariant(0, "Intentional error");
            failTest("no exception thrown");
        }
        catch(vigra::ContractViolation & c)
        {
            std::string expected("\nInvariant violation!\nIntentional error");
            std::string message(c.what());
            should(0 == expected.compare(message.substr(0,expected.size())));
        }
    }
};

struct ErrorTestSuite
: public vigra::test_suite
{
    ErrorTestSuite()
    : vigra::test_suite("ErrorTest")
    {
        add( testCase(&ErrorTest::testPrecondition));
        add( testCase(&ErrorTest::testPostcondition));
        add( testCase(&ErrorTest::testInvariant));
    }
};

int main()
{
    ErrorTestSuite test;

    int failed = test.run();

    std::cout << test.report() << std::endl;

    return (failed != 0);
}
