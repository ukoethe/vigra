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
 
 
#ifndef VIGRA_ERROR_HXX
#define VIGRA_ERROR_HXX

#include <stdexcept>
#include <stdio.h>
#include <string>
#include "config.hxx"
          
/*! \page ErrorReporting Error Reporting
    Exceptions and assertions provided by VIGRA

    <b>\#include</b> \<vigra/error.hxx\>
    
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
    the error are automatically included in the message. The macro
    
    \code
    vigra_assert(PREDICATE, MESSAGE);
    \endcode
    
    is identical to <tt>vigra_precondition()</tt> except that it is completely removed
    when '<TT>NDEBUG</TT>' is defined. This is useful for test that are only needed during 
    debugging, such as array index bound checking. The following macro
    
    \code
    vigra_fail(MESSAGE);
    \endcode
    
    unconditionally throws a '<TT>std::runtime_error</TT>' constructed from the message 
    (along with file name and line number, if NDEBUG is not set).
    
    <b> Usage:</b>
    
    Include-File:
    \<vigra/error.hxx\>
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
void throw_invariant_error(bool predicate, std::string message, char const * file, int line)
{
    if(!predicate)
	   throw vigra::InvariantViolation(message.c_str(), file, line); 
}

inline
void throw_precondition_error(bool predicate, char const * message, char const * file, int line)
{
    if(!predicate)
	   throw vigra::PreconditionViolation(message, file, line); 
}

inline
void throw_precondition_error(bool predicate, std::string message, char const * file, int line)
{
    if(!predicate)
	   throw vigra::PreconditionViolation(message.c_str(), file, line); 
}

inline
void throw_postcondition_error(bool predicate, char const * message, char const * file, int line)
{
    if(!predicate)
	   throw vigra::PostconditionViolation(message, file, line); 
}

inline
void throw_postcondition_error(bool predicate, std::string message, char const * file, int line)
{
    if(!predicate)
	   throw vigra::PostconditionViolation(message.c_str(), file, line); 
}

inline
void throw_runtime_error(char const * message, char const * file, int line)
{
    char what_[1100];
    sprintf(what_, "\n%.900s\n(%.100s:%d)\n", message, file, line);
    throw std::runtime_error(what_); 
}

inline
void throw_runtime_error(std::string message, char const * file, int line)
{
    char what_[1100];
    sprintf(what_, "\n%.900s\n(%.100s:%d)\n", message.c_str(), file, line);
    throw std::runtime_error(what_); 
}

#define vigra_precondition(PREDICATE, MESSAGE) vigra::throw_precondition_error((PREDICATE), MESSAGE, __FILE__, __LINE__)

#define vigra_assert(PREDICATE, MESSAGE) vigra_precondition(PREDICATE, MESSAGE)

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

inline
void throw_invariant_error(bool predicate, std::string message)
{
    if(!predicate)
	   throw vigra::InvariantViolation(message.c_str()); 
}

inline
void throw_precondition_error(bool predicate, std::string message)
{
    if(!predicate)
	   throw vigra::PreconditionViolation(message.c_str()); 
}

inline
void throw_postcondition_error(bool predicate, std::string message)
{
    if(!predicate)
	   throw vigra::PostconditionViolation(message.c_str()); 
}

#define vigra_precondition(PREDICATE, MESSAGE) vigra::throw_precondition_error((PREDICATE), MESSAGE)

#define vigra_assert(PREDICATE, MESSAGE)

#define vigra_postcondition(PREDICATE, MESSAGE) vigra::throw_postcondition_error((PREDICATE), MESSAGE)

#define vigra_invariant(PREDICATE, MESSAGE) vigra::throw_invariant_error((PREDICATE), MESSAGE)
            
#define vigra_fail(MESSAGE) throw std::runtime_error(MESSAGE)

#endif // NDEBUG

} // namespace vigra

#endif // VIGRA_ERROR_HXX
