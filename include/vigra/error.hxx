/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2000 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    ( Version 1.1.1, Feb 12 2001 )                                    */
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

#include <stdio.h>
#include "vigra/config.hxx"
          
#ifdef _MSC_VER
#  define snprintf _snprintf
#endif

/*! \page ErrorReporting Error Reporting

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
    their PREDICATE evaluates to '<TT>false<TT>':
    
    \code
    vigra_precondition(PREDICATE, MESSAGE);
    vigra_postcondition(PREDICATE, MESSAGE);
    vigra_invariant(PREDICATE, MESSAGE);
    \endcode
    
    The MESSAGE is passed to the exception and can be retrieved via
    the overloaded member function '<TT>exception.what()<TT>'. If the compiler
    flag '<TT>NDEBUG<TT>' is <em> not</em> defined, the file name and line number of 
    the error are automatically included in the message.
    
    The following macro
    
    \code
    vigra_fail(MESSAGE);
    \endcode
    
    unconditionally throws a '<TT>std::runtime_error<TT>' constructed from the message 
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
    
    \brief Exceptions and assertions provided by VIGRA
**/

namespace vigra {


class ContractViolation : public StdException
{
  protected:
#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { bufsize_ = 1100 };
#else
    static const int bufsize_ = 1100;
#endif
    
  public:
    ContractViolation(char const * message)
    {
        snprintf(what_, bufsize_, "%s", message);
    }
    
    ContractViolation()
    {
        what_[0] = 0;
    }
    
    ContractViolation(ContractViolation const & o)
    : StdException(o)
    {
        snprintf(what_, bufsize_, "%s", o.what_);
    }

    virtual const char * what() const throw()
    {
        return what_;
    }
    
    char what_[bufsize_];
};

class PreconditionViolation : public ContractViolation
{
  public:
    PreconditionViolation(char const * message, const char * file, int line)
    {
        snprintf( what_, bufsize_, "\nPrecondition violation!\n%s\n(%s:%d)\n", message, file, line);
    }
    
    PreconditionViolation(char const * message)
    {
        snprintf( what_, bufsize_, "\nPrecondition violation!\n%s\n", message);
    }
};

class PostconditionViolation : public ContractViolation
{
  public:
    PostconditionViolation(char const * message, const char * file, int line)
    {
        snprintf( what_, bufsize_, "\nPostcondition violation!\n%s\n(%s:%d)\n", message, file, line);
    }
    
    PostconditionViolation(char const * message)
    {
        snprintf( what_, bufsize_, "\nPostcondition violation!\n%s\n", message);
    }
};

class InvariantViolation : public ContractViolation
{
  public:
    InvariantViolation(char const * message, const char * file, int line)
    {
        snprintf( what_, bufsize_, "\nInvariant violation!\n%s\n(%s:%d)\n", message, file, line);
    }
    
    InvariantViolation(char const * message)
    {
        snprintf( what_, bufsize_, "\nInvariant violation!\n%s\n", message);
    }
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
    char buf[1000]; 
    snprintf(buf, 1000, "%s (%s:%d)", message, file, line); 
    throw std::runtime_error(buf); 
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
