/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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
 
 
#ifndef VIGRA_ERROR_HXX
#define VIGRA_ERROR_HXX

#include <stdexcept>
#include <stdio.h>
#include "vigra/config.hxx"
          
/*! \page ErrorReporting Error Reporting
    Exceptions and assertions provided by VIGRA

    <b>\#include</b> "<a href="error_8hxx-source.html">vigra/error.hxx</a>"
    
    VIGRA defines the following exception classes:
    
    \code
    namespace vigra {
        class ContractViolation : public std::exception;
        class PreconditionViolation : public ContractViolation;
        class PostconditionViolation : public ContractViolation;
        class InvariantViolation : public ContractViolation;
    }
    \endcode
    
    The following associated macros throw the corresponding exception if 
    their PREDICATE evaluates to '<TT>false</TT>':
    
    \code
    vigra_precondition(PREDICATE, MESSAGE);
    vigra_postcondition(PREDICATE, MESSAGE);
    vigra_invariant(PREDICATE, MESSAGE);
    \endcode
    
    The MESSAGE is passed to the exception and can be retrieved via
    the overloaded member function '<TT>exception.what()</TT>'. If the compiler
    flag '<TT>NDEBUG</TT>' is <em>not</em> defined, the file name and line number of 
    the error are automatically included in the message.
    
    The following macro
    
    \code
    vigra_fail(MESSAGE);
    \endcode
    
    unconditionally throws a '<TT>std::runtime_error</TT>' constructed from the message 
    (along with file name and line number, if NDEBUG is not set).
    
    <b> Usage:</b>
    
    Include-File:
    "<a href="error_8hxx-source.html">vigra/error.hxx</a>"
    <p>
    Namespace: vigra (except for the macros, of course)
    
    \code
    int main(int argc, char ** argv)
    {
        try
        {
            const char* input_file_name = argv[1];

            // read input image
            vigra::ImageImportInfo info(input_file_name);

            // fail if input image is not grayscale
            vigra_precondition(info.isGrayscale(), "Input image must be grayscale");

            ...// process image
        }
        catch (std::exception & e)
        {
            std::cerr << e.what() << std::endl; // print message
            return 1;
        }

        return 0;
    }
    \endcode
**/

namespace vigra {

class ContractViolation : public StdException
{
  public:
    ContractViolation(char const * prefix, char const * message, 
                      char const * file, int line)
    {
        sprintf(what_, "\n%.30s\n%.900s\n(%.100s:%d)\n", prefix, message, file, line);
    }
    
    ContractViolation(char const * prefix, char const * message)
    {
        sprintf(what_, "\n%.30s\n%.900s\n", prefix, message);
    }
    
    virtual const char * what() const throw()
    {
        return what_;
    }
  
  private:
    enum { bufsize_ = 1100 };
    char what_[bufsize_];
};

class PreconditionViolation : public ContractViolation
{
  public:
    PreconditionViolation(char const * message, const char * file, int line)
    : ContractViolation("Precondition violation!", message, file, line)
    {}
    
    PreconditionViolation(char const * message)
    : ContractViolation("Precondition violation!", message)
    {}
};

class PostconditionViolation : public ContractViolation
{
  public:
    PostconditionViolation(char const * message, const char * file, int line)
    : ContractViolation("Postcondition violation!", message, file, line)
    {}
    
    PostconditionViolation(char const * message)
    : ContractViolation("Postcondition violation!", message)
    {}
};

class InvariantViolation : public ContractViolation
{
  public:
    InvariantViolation(char const * message, const char * file, int line)
    : ContractViolation("Invariant violation!", message, file, line)
    {}
    
    InvariantViolation(char const * message)
    : ContractViolation("Invariant violation!", message)
    {}
};

#ifndef NDEBUG

inline
void throw_invariant_error(bool predicate, char const * message, char const * file, int line)
{
    if(!predicate)
	   throw vigra::InvariantViolation(message, file, line); 
}

inline
void throw_precondition_error(bool predicate, char const * message, char const * file, int line)
{
    if(!predicate)
	   throw vigra::PreconditionViolation(message, file, line); 
}

inline
void throw_postcondition_error(bool predicate, char const * message, char const * file, int line)
{
    if(!predicate)
	   throw vigra::PostconditionViolation(message, file, line); 
}

inline
void throw_runtime_error(char const * message, char const * file, int line)
{
    char what_[1100];
    sprintf(what_, "\n%.900s\n(%.100s:%d)\n", message, file, line);
    throw std::runtime_error(what_); 
}

#define vigra_precondition(PREDICATE, MESSAGE) vigra::throw_precondition_error((PREDICATE), MESSAGE, __FILE__, __LINE__)

#define vigra_postcondition(PREDICATE, MESSAGE) vigra::throw_postcondition_error((PREDICATE), MESSAGE, __FILE__, __LINE__)

#define vigra_invariant(PREDICATE, MESSAGE) vigra::throw_invariant_error((PREDICATE), MESSAGE, __FILE__, __LINE__)
            
#define vigra_fail(MESSAGE) vigra::throw_runtime_error(MESSAGE, __FILE__, __LINE__)

#else // NDEBUG

inline
void throw_invariant_error(bool predicate, char const * message)
{
    if(!predicate)
	   throw vigra::InvariantViolation(message); 
}

inline
void throw_precondition_error(bool predicate, char const * message)
{
    if(!predicate)
	   throw vigra::PreconditionViolation(message); 
}

inline
void throw_postcondition_error(bool predicate, char const * message)
{
    if(!predicate)
	   throw vigra::PostconditionViolation(message); 
}

#define vigra_precondition(PREDICATE, MESSAGE) vigra::throw_precondition_error((PREDICATE), MESSAGE)

#define vigra_postcondition(PREDICATE, MESSAGE) vigra::throw_postcondition_error((PREDICATE), MESSAGE)

#define vigra_invariant(PREDICATE, MESSAGE) vigra::throw_invariant_error((PREDICATE), MESSAGE)
            
#define vigra_fail(MESSAGE) throw std::runtime_error(MESSAGE)

#endif // NDEBUG

} // namespace vigra

#endif // VIGRA_ERROR_HXX
