#include <iostream>
#include "unittest.h"
#include "vigra/error.hxx"

struct ErrorTest
{
	void testPrecondition()
	{
		bool caughtIt = false;
		
		try
		{
			precondition(0, "Intentional error");
		}
		catch(ContractViolation & c)
		{
			caughtIt = true;
			std::cout << "Testing error reporting: " << c.what() << std::endl;
		}

		should(caughtIt == true);
	}

	void testPostcondition()
	{
		bool caughtIt = false;
		
		try
		{
			postcondition(0, "Intentional error");
		}
		catch(ContractViolation & c)
		{
			caughtIt = true;
			std::cout << "Testing error reporting: " << c.what() << std::endl;
		}

		should(caughtIt == true);
	}

	void testInvariant()
	{
		bool caughtIt = false;
		
		try
		{
			invariant(0, "Intentional error");
		}
		catch(ContractViolation & c)
		{
			caughtIt = true;
			std::cout << "Testing error reporting: " << c.what() << std::endl;
		}

		should(caughtIt == true);
	}
};

struct ErrorTestSuite
: public TestSuite
{
    ErrorTestSuite()
	: TestSuite("")
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
