/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2000 by Ullrich Koethe                  */
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

#include <stdio.h>
#include "vigra/config.hxx"
          
/** @heading Error Reporting

    VIGRA defines the following exception classes:
    
    \begin{verbatim}
    namespace vigra {
        class ContractViolation : public std::exception;
        class PreconditionViolation : public ContractViolation;
        class PostconditionViolation : public ContractViolation;
        class InvariantViolation : public ContractViolation;
    }
    \end{verbatim}
    
    The following associated macros throw the corresponding exception if 
    their PREDICATE evaluates to '#false#':
    
    \begin{verbatim}
    vigra_precondition(PREDICATE, MESSAGE);
    vigra_postcondition(PREDICATE, MESSAGE);
    vigra_invariant(PREDICATE, MESSAGE);
    \end{verbatim}
    
    The MESSAGE is passed to the exception and can be retrieved via
    the overloaded member function '#exception.what()#'. If the compiler
    flag '#NDEBUD#' is {\em not} defined, the file name and line number of 
    the error are automatically included in the message.
    
    The following macro
    
    \begin{verbatim}
    vigra_fail(MESSAGE);
    \end{verbatim}
    
    unconditionally throws a '#std::runtime_error#' constructed from the message 
    (along with file name and line number, if NDEBUG is not set).
    
    {\bf Usage:}
    
    Include-File:
    \URL[vigra/error.hxx]{../include/vigra/error.hxx}\\
    Namespace: vigra (except for the macros, of course)
    
    \begin{verbatim}
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
    \end{verbatim}
    
    @memo Exceptions and assertions provided by VIGRA
**/

namespace vigra {


class ContractViolation : public StdException
{
  public:
    ContractViolation(char const * message)
    {
        sprintf(what_, "%.1099s", message);
    }
    
    ContractViolation()
    {
        what_[0] = 0;
    }
    
    ContractViolation(ContractViolation const & o)
    {
        sprintf(what_, "%.1099s", o.what_);
    }

    virtual const char * what() const throw()
    {
        return what_;
    }
    
    char what_[1100];
};

class PreconditionViolation : public ContractViolation
{
  public:
    PreconditionViolation(char const * message, const char * file, int line)
    {
        sprintf( what_, "\nPrecondition Violation!\n%.900s\n(%.150s:%d)\n", message, file, line);
    }
    
    PreconditionViolation(char const * message)
    {
        sprintf( what_, "\nPrecondition Violation!\n%.900s\n", message);
    }
};

class PostconditionViolation : public ContractViolation
{
  public:
    PostconditionViolation(char const * message, const char * file, int line)
    {
        sprintf( what_, "\nPostcondition Violation!\n%.900s\n(%.150s:%d)\n", message, file, line);
    }
    
    PostconditionViolation(char const * message)
    {
        sprintf( what_, "\nPostcondition Violation!\n%.900s\n", message);
    }
};

class InvariantViolation : public ContractViolation
{
  public:
    InvariantViolation(char const * message, const char * file, int line)
    {
        sprintf( what_, "\nInvariant Violation!\n%.900s\n(%.150s:%d)\n", message, file, line);
    }
    
    InvariantViolation(char const * message)
    {
        sprintf( what_, "\nInvariant Violation!\n%.900s\n", message);
    }
};

#ifndef NDEBUG

#define vigra_precondition(PREDICATE, MESSAGE) \
        if(PREDICATE); else throw ::vigra::PreconditionViolation(MESSAGE, __FILE__, __LINE__)

#define vigra_postcondition(PREDICATE, MESSAGE) \
        if(PREDICATE); else throw ::vigra::PostconditionViolation(MESSAGE, __FILE__, __LINE__)

#define vigra_invariant(PREDICATE, MESSAGE) \
        if(PREDICATE); else throw ::vigra::InvariantViolation(MESSAGE, __FILE__, __LINE__)
            
#define vigra_fail(MESSAGE) \
        { \
            char buf[1000]; \
            sprintf(buf, "%.900s (" __FILE__ ":%d)", (MESSAGE), __LINE__); \
            throw std::runtime_error(buf); \
        } 

#else // NDEBUG

#define vigra_precondition(PREDICATE, MESSAGE) \
        if(PREDICATE); else throw ::vigra::PreconditionViolation(MESSAGE)

#define vigra_postcondition(PREDICATE, MESSAGE) \
        if(PREDICATE); else throw ::vigra::PostconditionViolation(MESSAGE)

#define vigra_invariant(PREDICATE, MESSAGE) \
        if(PREDICATE); else throw ::vigra::InvariantViolation(MESSAGE)
            
#define vigra_fail(MESSAGE) \
        { \
            throw std::runtime_error(MESSAGE); \
        } 

#endif // NDEBUG

} // namespace vigra

#endif // VIGRA_ERROR_HXX
