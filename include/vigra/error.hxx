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
          
class ContractViolation : public VigraStdException
{
  public:
    ContractViolation(char const * message)
    {
        sprintf(what_, "%.1099s", message);
    }
    
    ContractViolation()
    {}
    
    ContractViolation(ContractViolation const & o)
    {
        sprintf(what_, "%.1099s", o.what_);
    }

//    virtual ~ContractViolation() {}
   
    virtual const char * what() const throw()
    {
        return what_;
    }
    
    char what_[1100];
};

class PreconditionViolation : public ContractViolation
{
  public:
//    virtual ~PreconditionViolation() {}
   
    PreconditionViolation(char const * message, const char * file, int line)
    {
        sprintf( what_, "\nPrecondition Violation!\n%.900s\n(%.150s:%d)\n", message, file, line);
    }
};

class PostconditionViolation : public ContractViolation
{
  public:
//    virtual ~PostconditionViolation() {}
   
    PostconditionViolation(char const * message, const char * file, int line)
    {
        sprintf( what_, "\nPostcondition Violation!\n%.900s\n(%.150s:%d)\n", message, file, line);
    }
};

class InvariantViolation : public ContractViolation
{
  public:
//    virtual ~InvariantViolation() {}
   
    InvariantViolation(char const * message, const char * file, int line)
    {
        sprintf( what_, "\nInvariant Violation!\n%.900s\n(%.150s:%d)\n", message, file, line);
    }
};


#define precondition(PREDICATE, MESSAGE) \
        if(PREDICATE); else throw PreconditionViolation(MESSAGE, __FILE__, __LINE__)

#define postcondition(PREDICATE, MESSAGE) \
        if(PREDICATE); else throw PostconditionViolation(MESSAGE, __FILE__, __LINE__)

#define invariant(PREDICATE, MESSAGE) \
        if(PREDICATE); else throw InvariantViolation(MESSAGE, __FILE__, __LINE__)
            
#define fail(MESSAGE) \
        { \
            char buf[1000]; \
            sprintf(buf, "%.900s (" __FILE__ ":%d)", (MESSAGE), __LINE__); \
            throw std::runtime_error(buf); \
        } 

#endif // VIGRA_ERROR_HXX
