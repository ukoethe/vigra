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
 
 
#ifndef VIGRA_NUMERICTRAITS_HXX
#define VIGRA_NUMERICTRAITS_HXX

#include <limits.h>
#include <vigra/utilities.hxx>

/********************************************************/
/*                                                      */
/*                      NumericTraits                   */
/*                                                      */
/********************************************************/


/**@name Numeric and Promotion Traits

    These traits classes contain information that is used by generic
    algorithms and data structures to determine intermediate and result
    types of numerical calculations, to convert between different 
    representations of arithmetic types, and to create certain important
    constants of each type. Thus, algorithms and data structures
    operating that need arithmetic operations can be made more
    independent from the actual data representation.
    
    NumericTraits are implemented as template specializations of one
    arithmetic type, while PromoteTraits are specialized for a pair of
    arithmetic types that shall be combined in one operation.
    
* @memo   Meta-information about arithmetic types
*/
//@{
/**@name template<> struct NumericTraits<ArithmeticType>

    This traits class is used derive important properties of
    an arithmetic type. Consider the following algorithm:
    
    \begin{verbatim}
    // calculate the sum of a sequence of bytes
    int sumBytes(unsigned char * begin, unsigned char * end)
    {
        int result = 0;
    for(; begin != end; ++begin)  result += *begin;
    return result;
    }
    \end{verbatim} 
    
    The return type of this function can not be 'unsigned char' because
    the summation would very likely overflow. Since we know the source
    type, we can easily choose 'int' as an appropriate return type.
    Likewise, we would have choosen 'float' if we had to sum a 
    sequence of floats. If we want to make this 
    algorithm generic, we would like to derive the appropriate return 
    type automatically. This can be done with NumericTraits. 
    The code would look like this (we use \Ref{Data Accessors} to 
    read the data from the sequence):
    
    \begin{verbatim}
    // calculate the sum of any sequence
    template <class Iterator, class Accessor>
    NumericTraits<typename Accessor::value_type>::Promote
    sumSequence(Iterator begin, Iterator end, Accessor a)
    {
        // an abbraviation
    typedef NumericTraits<typename Accessor::value_type>  SrcTraits;
        
    // find out result type
    typedef typename SrcTraits::Promote ResultType;
      
    // init result to zero
    ResultType result = NumericTraits<ResultType>::zero();
    
    for(; begin != end; ++begin)
    {  
        // cast current item to ResultType and add
        result += SrcTraits::toPromote(a(begin));
    }
        
    return result;
    }
    \end{verbatim}
    
    In this example NumericTraits is not only used to deduce the 
    ReturnType of the operation, but also to initialize it with the
    constant 'zero'. This is necessary since we do not know in general,
    which expression must be used to obtain a zero of some arbitrary
    type - '#ResultType result = 0;#' would only work if the 
    ResultType had an constructor taking an '#int#' argument, and we 
    would not even have any guarantee as to what the semantics of this
    constructor are. In addition, the traits are used to cast the 
    source type into the promote type.
    
    Similarly, an algorithm that needs multiplication would use the 
    return type #RealPromote# and the functions #one()# and
    #toRealPromote()#. The following members are defined in 
    {\bf #NumericTraits<ArithmeticType>#}:
    
    \begin{tabular}{ll}
    {\bf #typedef ... Promote;#} & 
    
            promote type for addition and subtraction 
        
        \\

    {\bf #typedef ... RealPromote;#} &
    
        promote type for multiplication and division
    
    (only defined if #ArithmeticType# supports these operations) 
    
    \\

    {\bf #static Promote toPromote(ArithmeticType v);#} &
    
        convert to #Promote# type 
    
    \\
    
    {\bf #static RealPromote toRealPromote(ArithmeticType v);#} &
    
        convert to #RealPromote# type 

    (only defined if #ArithmeticType# supports multiplication) 
    
    \\
    
    {\bf #static ArithmeticType fromPromote(Promote v);#} &
    
        convert from #Promote# type
    
    if #v# is outside the range of #ArithmeticType# it is clipped;

    \\
    
    {\bf #static ArithmeticType fromRealPromote(RealPromote v);#} &
    
        convert from #RealPromote# type 
    
    (only defined if 
    #ArithmeticType# supports multiplication)
    
    if #ArithmeticType# is an integral type, the result is rounded 
    
    if #v# is outside the range of #ArithmeticType# it is clipped
    
    \\

    {\bf #static ArithmeticType zero();#} & 
    
    create neutral element of addition
    
    i.e. #(ArithmeticType a = ...,# 
    #  a + NumericTraits<ArithmeticType>::zero() == a)# 
    must always yield #true# 
    
    \\

    {\bf #static ArithmeticType nonZero();#} & 
    
    create a non-zero element (if multiplication is defined, this yields one())
    
    i.e. #(ArithmeticType a = ...,# 
    #  a + NumericTraits<ArithmeticType>::nonZero() == a)# 
    must always yield #false# 
    
    \\
    
    {\bf #static ArithmeticType one();#} & 
    
    create neutral element of multiplication 
    
    (only defined if #ArithmeticType# supports multiplication)
    
    i.e. #(ArithmeticType a = ...,# 
    #  a * NumericTraits<ArithmeticType>::one() == a)# 
    must always yield #true# 
    
    \\
    
    {\bf #static const bool is_integral;#} & 
    
        true if #ArithmeticType# is an integral type, false otherwise 
    
    \\
    
    {\bf #static const bool is_scalar;#} &
    
        true if #ArithmeticType# is a scalar type, false otherwise 
    
    \\
    
    \end{tabular}
    
    NumericTraits for the built-in types are defined in Include-File: 
    \URL[numerictraits.hxx]{../include/numerictraits.hxx}
    
* @memo Unary traits for promotion, conversion, creation of arithmetic objects
*/

/**@name template<> struct PromoteTraits<ArithmeticType1, ArithmeticType2>

    This traits class is used to determine the appropriate result type
    of arithmetic expressions which depend of two arguments. Consider
    the following function:
    
    \begin{verbatim}
    template <class T>
    T min(T t1, T t2)
    {
        return (t1 < t2) ? t1 : t2;
    }
    \end{verbatim}
    
    This template is only applicable if both arguments have the same
    type. However, sometimes we may want to use the function in cases
    where the argument types differ. The we can deduce the approrpiate
    return type by using #PromoteTraits#:
    
    \begin{verbatim}
    template <class T1, class T2>
    PromoteTraits<T1, T2>::Promote
    min(T1 t1, T2 t2)
    {
        return (t1 < t2) ? PromoteTraits<T1, T2>::toPromote(t1) : 
                       PromoteTraits<T1, T2>::toPromote(t2);
    }    
    \end{verbatim}
    
    In addition, the traits class provide static functions to cast the
    arguments to the promote type. For example, if #T1# were #int# and 
    #T2# were #float#, the #Promote# type would be #float#. 
    The following members are defined in 
    {\bf #PromoteTraits<ArithmeticType1, ArithmeticType2>#}:
    
    \begin{tabular}{ll}
    {\bf #typedef ... Promote;#} & 
    
            promote type 
        
        \\

    {\bf #static Promote toPromote(ArithmeticType1 v);#} 
    
    {\bf #static Promote toPromote(ArithmeticType2 v);#}  &
    
        convert to #Promote# type 
    
    \\
    
    \end{tabular}
    
    PromoteTraits for the built-in types are defined in Include-File: 
    \URL[numerictraits.hxx]{../include/numerictraits.hxx}
    
* @memo   Binary traits for promotion of arithmetic objects
*/
//@}

struct Error_NumericTraits_not_specialized_for_this_case { };

template<class A>
struct NumericTraits
{
    typedef Error_NumericTraits_not_specialized_for_this_case Promote;
    typedef Error_NumericTraits_not_specialized_for_this_case RealPromote;
    typedef Error_NumericTraits_not_specialized_for_this_case isScalar;
    typedef Error_NumericTraits_not_specialized_for_this_case isIntegral;
    typedef Error_NumericTraits_not_specialized_for_this_case isOrdered;
};

#ifndef NO_BOOL
template<>
struct NumericTraits<bool>
{
    typedef bool Type;
    typedef int Promote;
    typedef double RealPromote;
    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isOrdered;
    
    static bool zero() { return false; }
    static bool one() { return true; }
    static bool nonZero() { return true; }
    static bool min() { return false; }
    static bool max() { return true; }
    
    static const bool minConst = false;
    static const bool maxConst = true;
    
    static Promote toPromote(bool v) { return v ? 1 : 0; }
    static RealPromote toRealPromote(bool v) { return v ? 1.0 : 0.0; }
    static bool fromPromote(Promote v) { 
        return (v == 0) ? false : true; 
    }
    static bool fromRealPromote(RealPromote v) {
        return (v == 0.0) ? false : true; 
    }
};
#endif

template<>
struct NumericTraits<signed char>
{
    typedef signed char Type;
    typedef int Promote;
    typedef double RealPromote;
    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isOrdered;
    
    static signed char zero() { return 0; }
    static signed char one() { return 1; }
    static signed char nonZero() { return 1; }
    static signed char min() { return SCHAR_MIN; }
    static signed char max() { return SCHAR_MAX; }
    
    static const signed char minConst = SCHAR_MIN;
    static const signed char maxConst = SCHAR_MIN;
    
    static Promote toPromote(signed char v) { return v; }
    static RealPromote toRealPromote(signed char v) { return v; }
    static signed char fromPromote(Promote v) { 
        return ((v < SCHAR_MIN) ? SCHAR_MIN : (v > SCHAR_MAX) ? SCHAR_MAX : v); 
    }
    static signed char fromRealPromote(RealPromote v) {
        return ((v < 0.0) ? ((v < (float)SCHAR_MIN) ? SCHAR_MIN : static_cast<signed char>(v - 0.5)) : 
                (v > SCHAR_MAX) ? SCHAR_MAX : static_cast<signed char>(v + 0.5)); 
    }
};

template<>
struct NumericTraits<unsigned char>
{
    typedef unsigned char Type;
    typedef int Promote;
    typedef double RealPromote;
    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isOrdered;
    
    static unsigned char zero() { return 0; }
    static unsigned char one() { return 1; }
    static unsigned char nonZero() { return 1; }
    static unsigned char min() { return 0; }
    static unsigned char max() { return UCHAR_MAX; }
    
    static const unsigned char minConst = 0;
    static const unsigned char maxConst = UCHAR_MAX;
    
    static Promote toPromote(unsigned char v) { return v; }
    static RealPromote toRealPromote(unsigned char v) { return v; }
    static unsigned char fromPromote(Promote const & v) { 
        return ((v < 0) ? 0 : (v > UCHAR_MAX) ? UCHAR_MAX : v); 
    }
    static unsigned char fromRealPromote(RealPromote const & v) {
            return ((v < 0.0) ? 0 : ((v > (float)UCHAR_MAX) ? UCHAR_MAX : static_cast<unsigned char>(v + 0.5)));
    }
};

template<>
struct NumericTraits<short int>
{
    typedef short int Type;
    typedef int Promote;
    typedef double RealPromote;
    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isOrdered;
    
    static short int zero() { return 0; }
    static short int one() { return 1; }
    static short int nonZero() { return 1; }
    static short int min() { return SHRT_MIN; }
    static short int max() { return SHRT_MAX; }
    
    static const short int minConst = SHRT_MIN;
    static const short int maxConst = SHRT_MAX;
    
    static Promote toPromote(short int v) { return v; }
    static RealPromote toRealPromote(short int v) { return v; }
    static short int fromPromote(Promote v) { 
        return ((v < SHRT_MIN) ? SHRT_MIN : 
                (v > SHRT_MAX) ? SHRT_MAX : v); 
    }
    static short int fromRealPromote(RealPromote v) {
        return ((v < 0.0) ? 
                ((v < (float)SHRT_MIN) ? SHRT_MIN : static_cast<short int>(v - 0.5)) : 
                ((v > (float)SHRT_MAX) ? SHRT_MAX : static_cast<short int>(v + 0.5))); 
    }
};

template<>
struct NumericTraits<short unsigned int>
{
    typedef short unsigned int Type;
    typedef int Promote;
    typedef double RealPromote;

    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isOrdered;

    static short unsigned int zero() { return 0; }
    static short unsigned int one() { return 1; }
    static short unsigned int nonZero() { return 1; }
    static short unsigned int min() { return 0; }
    static short unsigned int max() { return USHRT_MAX; }
    
    static const short unsigned int minConst = 0;
    static const short unsigned int maxConst = USHRT_MAX;

    static Promote toPromote(short unsigned int v) { return v; }
    static RealPromote toRealPromote(short unsigned int v) { return v; }
    static short unsigned int fromPromote(Promote v) { 
        return ((v < 0) ? 0 : (v > USHRT_MAX) ? USHRT_MAX : v); 
    }
    static short unsigned int fromRealPromote(RealPromote v) {
            return ((v < 0.0) ? 
              0 : ((v > (float)USHRT_MAX) ? USHRT_MAX : static_cast<short unsigned int>(v + 0.5)));
    }
};

template<>
struct NumericTraits<int>
{
    typedef int Type;
    typedef int Promote;
    typedef double RealPromote;
    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isOrdered;

    static int zero() { return 0; }
    static int one() { return 1; }
    static int nonZero() { return 1; }
    static int min() { return INT_MIN; }
    static int max() { return INT_MAX; }
    
    static const int minConst = INT_MIN;
    static const int maxConst = INT_MAX;

    static Promote toPromote(int v) { return v; }
    static RealPromote toRealPromote(int v) { return v; }
    static int fromPromote(Promote v) { return v; }
    static int fromRealPromote(RealPromote v) {
        return ((v < 0.0) ? 
                ((v < (float)INT_MIN) ? INT_MIN : static_cast<int>(v - 0.5)) : 
                ((v > (float)INT_MAX) ? INT_MAX : static_cast<int>(v + 0.5))); 
    }
};

template<>
struct NumericTraits<unsigned int>
{
    typedef unsigned int Type;
    typedef unsigned int Promote;
    typedef double RealPromote;
    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isOrdered;
    
    static unsigned int zero() { return 0; }
    static unsigned int one() { return 1; }
    static unsigned int nonZero() { return 1; }
    static unsigned int min() { return 0; }
    static unsigned int max() { return UINT_MAX; }
    
    static const unsigned int minConst = 0;
    static const unsigned int maxConst = UINT_MAX;

    static Promote toPromote(unsigned int v) { return v; }
    static RealPromote toRealPromote(unsigned int v) { return v; }
    static unsigned int fromPromote(Promote v) { return v; }
    static unsigned int fromRealPromote(RealPromote v) {
            return ((v < 0.0) ? 0 : 
               ((v > (float)UINT_MAX) ? 
                         UINT_MAX : static_cast<unsigned int>(v + 0.5)));
    }
};

template<>
struct NumericTraits<long>
{
    typedef long Type;
    typedef long Promote;
    typedef double RealPromote;
    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isOrdered;
    
    static long zero() { return 0; }
    static long one() { return 1; }
    static long nonZero() { return 1; }
    static long min() { return LONG_MIN; }
    static long max() { return LONG_MAX; }
    
    static const long minConst = LONG_MIN;
    static const long maxConst = LONG_MAX;

    static Promote toPromote(long v) { return v; }
    static RealPromote toRealPromote(long v) { return v; }
    static long fromPromote(Promote v) { return v; }
    static long fromRealPromote(RealPromote v) {
        return ((v < 0.0) ? 
                ((v < (float)LONG_MIN) ? LONG_MIN : static_cast<long>(v - 0.5)) : 
                ((v > (float)LONG_MAX) ? LONG_MAX : static_cast<long>(v + 0.5))); 
    }
};

template<>
struct NumericTraits<unsigned long>
{
    typedef unsigned long Type;
    typedef unsigned long Promote;
    typedef double RealPromote;
    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isOrdered;
    
    static unsigned long zero() { return 0; }
    static unsigned long one() { return 1; }
    static unsigned long nonZero() { return 1; }
    static unsigned long min() { return 0; }
    static unsigned long max() { return ULONG_MAX; }
    
    static const unsigned long minConst = 0;
    static const unsigned long maxConst = ULONG_MAX;

    static Promote toPromote(unsigned long v) { return v; }
    static RealPromote toRealPromote(unsigned long v) { return v; }
    static unsigned long fromPromote(Promote v) { return v; }
    static unsigned long fromRealPromote(RealPromote v) {
            return ((v < 0.0) ? 0 : 
               ((v > (float)ULONG_MAX) ? 
                            ULONG_MAX : static_cast<unsigned long>(v + 0.5)));
    }
};

template<>
struct NumericTraits<float>
{
    typedef float Type;
    typedef float Promote;
    typedef double RealPromote;
    typedef VigraFalseType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isOrdered;
    
    static float zero() { return 0.0; }
    static float one() { return 1.0; }
    static float nonZero() { return 1.0; }
    static float min() { return -FLT_MAX; }
    static float max() { return FLT_MAX; }
    
    static Promote toPromote(float v) { return v; }
    static RealPromote toRealPromote(float v) { return v; }
    static float fromPromote(Promote v) { return v; }
    static float fromRealPromote(RealPromote v) { return v; }
};

template<>
struct NumericTraits<double>
{
    typedef double Type;
    typedef double Promote;
    typedef double RealPromote;
    typedef VigraFalseType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isOrdered;
    
    static double zero() { return 0.0; }
    static double one() { return 1.0; }
    static double nonZero() { return 1.0; }
    static double min() { return -DBL_MAX; }
    static double max() { return DBL_MAX; }

    static Promote toPromote(double v) { return v; }
    static RealPromote toRealPromote(double v) { return v; }
    static double fromPromote(Promote v) { return v; }
    static double fromRealPromote(RealPromote v) { return v; }
};

template<>
struct NumericTraits<long double>
{
    typedef long double Type;
    typedef long double Promote;
    typedef long double RealPromote;
    typedef VigraFalseType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isOrdered;
    
    static long double zero() { return 0.0; }
    static long double one() { return 1.0; }
    static long double nonZero() { return 1.0; }

    static Promote toPromote(long double v) { return v; }
    static RealPromote toRealPromote(long double v) { return v; }
    static long double fromPromote(Promote v) { return v; }
    static long double fromRealPromote(RealPromote v) { return v; }
};

/********************************************************/
/*                                                      */
/*                      PromoteTraits                  */
/*                                                      */
/********************************************************/

struct Error_PromoteTraits_not_specialized_for_this_case { };

template<class A, class B>
struct PromoteTraits
{
        typedef Error_PromoteTraits_not_specialized_for_this_case Promote;
};

template<>
struct PromoteTraits<char, char>
{
    typedef int Promote;
    static Promote toPromote(char v) { return v; }
};

template<>
struct PromoteTraits<char, unsigned char>
{
    typedef int Promote;
    static Promote toPromote(char v) { return v; }
    static Promote toPromote(unsigned char v) { return v; }
};

template<>
struct PromoteTraits<char, short int>
{
    typedef int Promote;
    static Promote toPromote(char v) { return v; }
    static Promote toPromote(short int v) { return v; }
};

template<>
struct PromoteTraits<char, short unsigned int>
{
    typedef unsigned int Promote;
    static Promote toPromote(char v) { return v; }
    static Promote toPromote(short unsigned int v) { return v; }
};

template<>
struct PromoteTraits<char, int>
{
    typedef int Promote;
    static Promote toPromote(char v) { return v; }
    static Promote toPromote(int v) { return v; }
};

template<>
struct PromoteTraits<char, unsigned int>
{
    typedef unsigned int Promote;
    static Promote toPromote(char v) { return v; }
    static Promote toPromote(unsigned int v) { return v; }
};

template<>
struct PromoteTraits<char, long>
{
    typedef long Promote;
    static Promote toPromote(char v) { return v; }
    static Promote toPromote(long v) { return v; }
};

template<>
struct PromoteTraits<char, unsigned long>
{
    typedef unsigned long Promote;
    static Promote toPromote(char v) { return v; }
    static Promote toPromote(unsigned long v) { return v; }
};

template<>
struct PromoteTraits<char, float>
{
    typedef double Promote;
    static Promote toPromote(char v) { return v; }
    static Promote toPromote(float v) { return v; }
};

template<>
struct PromoteTraits<char, double>
{
    typedef double Promote;
    static Promote toPromote(char v) { return v; }
    static Promote toPromote(double v) { return v; }
};

template<>
struct PromoteTraits<char, long double>
{
    typedef long double Promote;
    static Promote toPromote(char v) { return v; }
    static Promote toPromote(long double v) { return v; }
};

template<>
struct PromoteTraits<unsigned char, char>
{
    typedef int Promote;
    static Promote toPromote(unsigned char v) { return v; }
    static Promote toPromote(char v) { return v; }
};

template<>
struct PromoteTraits<unsigned char, unsigned char>
{
    typedef int Promote;
    static Promote toPromote(unsigned char v) { return v; }
};

template<>
struct PromoteTraits<unsigned char, short int>
{
    typedef int Promote;
    static Promote toPromote(unsigned char v) { return v; }
    static Promote toPromote(short int v) { return v; }
};

template<>
struct PromoteTraits<unsigned char, short unsigned int>
{
    typedef unsigned int Promote;
    static Promote toPromote(unsigned char v) { return v; }
    static Promote toPromote(short unsigned int v) { return v; }
};

template<>
struct PromoteTraits<unsigned char, int>
{
    typedef int Promote;
    static Promote toPromote(unsigned char v) { return v; }
    static Promote toPromote(int v) { return v; }
};

template<>
struct PromoteTraits<unsigned char, unsigned int>
{
    typedef unsigned int Promote;
    static Promote toPromote(unsigned char v) { return v; }
    static Promote toPromote(unsigned int v) { return v; }
};

template<>
struct PromoteTraits<unsigned char, long>
{
    typedef long Promote;
    static Promote toPromote(unsigned char v) { return v; }
    static Promote toPromote(long v) { return v; }
};

template<>
struct PromoteTraits<unsigned char, unsigned long>
{
    typedef unsigned long Promote;
    static Promote toPromote(unsigned char v) { return v; }
    static Promote toPromote(unsigned long v) { return v; }
};

template<>
struct PromoteTraits<unsigned char, float>
{
    typedef double Promote;
    static Promote toPromote(unsigned char v) { return v; }
    static Promote toPromote(float v) { return v; }
};

template<>
struct PromoteTraits<unsigned char, double>
{
    typedef double Promote;
    static Promote toPromote(unsigned char v) { return v; }
    static Promote toPromote(double v) { return v; }
};

template<>
struct PromoteTraits<unsigned char, long double>
{
    typedef long double Promote;
    static Promote toPromote(unsigned char v) { return v; }
    static Promote toPromote(long double v) { return v; }
};

template<>
struct PromoteTraits<short int, char>
{
    typedef int Promote;
    static Promote toPromote(short int v) { return v; }
    static Promote toPromote(char v) { return v; }
};

template<>
struct PromoteTraits<short int, unsigned char>
{
    typedef int Promote;
    static Promote toPromote(short int v) { return v; }
    static Promote toPromote(unsigned char v) { return v; }
};

template<>
struct PromoteTraits<short int, short int>
{
    typedef int Promote;
    static Promote toPromote(short int v) { return v; }
};

template<>
struct PromoteTraits<short int, short unsigned int>
{
    typedef unsigned int Promote;
    static Promote toPromote(short int v) { return v; }
    static Promote toPromote(short unsigned int v) { return v; }
};

template<>
struct PromoteTraits<short int, int>
{
    typedef int Promote;
    static Promote toPromote(short int v) { return v; }
    static Promote toPromote(int v) { return v; }
};

template<>
struct PromoteTraits<short int, unsigned int>
{
    typedef unsigned int Promote;
    static Promote toPromote(short int v) { return v; }
    static Promote toPromote(unsigned int v) { return v; }
};

template<>
struct PromoteTraits<short int, long>
{
    typedef long Promote;
    static Promote toPromote(short int v) { return v; }
    static Promote toPromote(long v) { return v; }
};

template<>
struct PromoteTraits<short int, unsigned long>
{
    typedef unsigned long Promote;
    static Promote toPromote(short int v) { return v; }
    static Promote toPromote(unsigned long v) { return v; }
};

template<>
struct PromoteTraits<short int, float>
{
    typedef double Promote;
    static Promote toPromote(short int v) { return v; }
    static Promote toPromote(float v) { return v; }
};

template<>
struct PromoteTraits<short int, double>
{
    typedef double Promote;
    static Promote toPromote(short int v) { return v; }
    static Promote toPromote(double v) { return v; }
};

template<>
struct PromoteTraits<short int, long double>
{
    typedef long double Promote;
    static Promote toPromote(short int v) { return v; }
    static Promote toPromote(long double v) { return v; }
};

template<>
struct PromoteTraits<short unsigned int, char>
{
    typedef unsigned int Promote;
    static Promote toPromote(short unsigned int v) { return v; }
    static Promote toPromote(char v) { return v; }
};

template<>
struct PromoteTraits<short unsigned int, unsigned char>
{
    typedef unsigned int Promote;
    static Promote toPromote(short unsigned int v) { return v; }
    static Promote toPromote(unsigned char v) { return v; }
};

template<>
struct PromoteTraits<short unsigned int, short int>
{
    typedef unsigned int Promote;
    static Promote toPromote(short unsigned int v) { return v; }
    static Promote toPromote(short int v) { return v; }
};

template<>
struct PromoteTraits<short unsigned int, short unsigned int>
{
    typedef unsigned int Promote;
    static Promote toPromote(short unsigned int v) { return v; }
};

template<>
struct PromoteTraits<short unsigned int, int>
{
    typedef unsigned int Promote;
    static Promote toPromote(short unsigned int v) { return v; }
    static Promote toPromote(int v) { return v; }
};

template<>
struct PromoteTraits<short unsigned int, unsigned int>
{
    typedef unsigned int Promote;
    static Promote toPromote(short unsigned int v) { return v; }
    static Promote toPromote(unsigned int v) { return v; }
};

template<>
struct PromoteTraits<short unsigned int, long>
{
    typedef long Promote;
    static Promote toPromote(short unsigned int v) { return v; }
    static Promote toPromote(long v) { return v; }
};

template<>
struct PromoteTraits<short unsigned int, unsigned long>
{
    typedef unsigned long Promote;
    static Promote toPromote(short unsigned int v) { return v; }
    static Promote toPromote(unsigned long v) { return v; }
};

template<>
struct PromoteTraits<short unsigned int, float>
{
    typedef double Promote;
    static Promote toPromote(short unsigned int v) { return v; }
    static Promote toPromote(float v) { return v; }
};

template<>
struct PromoteTraits<short unsigned int, double>
{
    typedef double Promote;
    static Promote toPromote(short unsigned int v) { return v; }
    static Promote toPromote(double v) { return v; }
};

template<>
struct PromoteTraits<short unsigned int, long double>
{
    typedef long double Promote;
    static Promote toPromote(short unsigned int v) { return v; }
    static Promote toPromote(long double v) { return v; }
};

template<>
struct PromoteTraits<int, char>
{
    typedef int Promote;
    static Promote toPromote(int v) { return v; }
    static Promote toPromote(char v) { return v; }
};

template<>
struct PromoteTraits<int, unsigned char>
{
    typedef int Promote;
    static Promote toPromote(int v) { return v; }
    static Promote toPromote(unsigned char v) { return v; }
};

template<>
struct PromoteTraits<int, short int>
{
    typedef int Promote;
    static Promote toPromote(int v) { return v; }
    static Promote toPromote(short int v) { return v; }
};

template<>
struct PromoteTraits<int, short unsigned int>
{
    typedef unsigned int Promote;
    static Promote toPromote(int v) { return v; }
    static Promote toPromote(short unsigned int v) { return v; }
};

template<>
struct PromoteTraits<int, int>
{
    typedef int Promote;
    static Promote toPromote(int v) { return v; }
};

template<>
struct PromoteTraits<int, unsigned int>
{
    typedef unsigned int Promote;
    static Promote toPromote(int v) { return v; }
    static Promote toPromote(unsigned int v) { return v; }
};

template<>
struct PromoteTraits<int, long>
{
    typedef long Promote;
    static Promote toPromote(int v) { return v; }
    static Promote toPromote(long v) { return v; }
};

template<>
struct PromoteTraits<int, unsigned long>
{
    typedef unsigned long Promote;
    static Promote toPromote(int v) { return v; }
    static Promote toPromote(unsigned long v) { return v; }
};

template<>
struct PromoteTraits<int, float>
{
    typedef double Promote;
    static Promote toPromote(int v) { return v; }
    static Promote toPromote(float v) { return v; }
};

template<>
struct PromoteTraits<int, double>
{
    typedef double Promote;
    static Promote toPromote(int v) { return v; }
    static Promote toPromote(double v) { return v; }
};

template<>
struct PromoteTraits<int, long double>
{
    typedef long double Promote;
    static Promote toPromote(int v) { return v; }
    static Promote toPromote(long double v) { return v; }
};

template<>
struct PromoteTraits<unsigned int, char>
{
    typedef unsigned int Promote;
    static Promote toPromote(unsigned int v) { return v; }
    static Promote toPromote(char v) { return v; }
};

template<>
struct PromoteTraits<unsigned int, unsigned char>
{
    typedef unsigned int Promote;
    static Promote toPromote(unsigned int v) { return v; }
    static Promote toPromote(unsigned char v) { return v; }
};

template<>
struct PromoteTraits<unsigned int, short int>
{
    typedef unsigned int Promote;
    static Promote toPromote(unsigned int v) { return v; }
    static Promote toPromote(short int v) { return v; }
};

template<>
struct PromoteTraits<unsigned int, short unsigned int>
{
    typedef unsigned int Promote;
    static Promote toPromote(unsigned int v) { return v; }
    static Promote toPromote(short unsigned int v) { return v; }
};

template<>
struct PromoteTraits<unsigned int, int>
{
    typedef unsigned int Promote;
    static Promote toPromote(unsigned int v) { return v; }
    static Promote toPromote(int v) { return v; }
};

template<>
struct PromoteTraits<unsigned int, unsigned int>
{
    typedef unsigned int Promote;
    static Promote toPromote(unsigned int v) { return v; }
};

template<>
struct PromoteTraits<unsigned int, long>
{
    typedef long Promote;
    static Promote toPromote(unsigned int v) { return v; }
    static Promote toPromote(long v) { return v; }
};

template<>
struct PromoteTraits<unsigned int, unsigned long>
{
    typedef unsigned long Promote;
    static Promote toPromote(unsigned int v) { return v; }
    static Promote toPromote(unsigned long v) { return v; }
};

template<>
struct PromoteTraits<unsigned int, float>
{
    typedef double Promote;
    static Promote toPromote(unsigned int v) { return v; }
    static Promote toPromote(float v) { return v; }
};

template<>
struct PromoteTraits<unsigned int, double>
{
    typedef double Promote;
    static Promote toPromote(unsigned int v) { return v; }
    static Promote toPromote(double v) { return v; }
};

template<>
struct PromoteTraits<unsigned int, long double>
{
    typedef long double Promote;
    static Promote toPromote(unsigned int v) { return v; }
    static Promote toPromote(long double v) { return v; }
};

template<>
struct PromoteTraits<long, char>
{
    typedef long Promote;
    static Promote toPromote(long v) { return v; }
    static Promote toPromote(char v) { return v; }
};

template<>
struct PromoteTraits<long, unsigned char>
{
    typedef long Promote;
    static Promote toPromote(long v) { return v; }
    static Promote toPromote(unsigned char v) { return v; }
};

template<>
struct PromoteTraits<long, short int>
{
    typedef long Promote;
    static Promote toPromote(long v) { return v; }
    static Promote toPromote(short int v) { return v; }
};

template<>
struct PromoteTraits<long, short unsigned int>
{
    typedef long Promote;
    static Promote toPromote(long v) { return v; }
    static Promote toPromote(short unsigned int v) { return v; }
};

template<>
struct PromoteTraits<long, int>
{
    typedef long Promote;
    static Promote toPromote(long v) { return v; }
    static Promote toPromote(int v) { return v; }
};

template<>
struct PromoteTraits<long, unsigned int>
{
    typedef long Promote;
    static Promote toPromote(long v) { return v; }
    static Promote toPromote(unsigned int v) { return v; }
};

template<>
struct PromoteTraits<long, long>
{
    typedef long Promote;
    static Promote toPromote(long v) { return v; }
};

template<>
struct PromoteTraits<long, unsigned long>
{
    typedef unsigned long Promote;
    static Promote toPromote(long v) { return v; }
    static Promote toPromote(unsigned long v) { return v; }
};

template<>
struct PromoteTraits<long, float>
{
    typedef double Promote;
    static Promote toPromote(long v) { return v; }
    static Promote toPromote(float v) { return v; }
};

template<>
struct PromoteTraits<long, double>
{
    typedef double Promote;
    static Promote toPromote(long v) { return v; }
    static Promote toPromote(double v) { return v; }
};

template<>
struct PromoteTraits<long, long double>
{
    typedef long double Promote;
    static Promote toPromote(long v) { return v; }
    static Promote toPromote(long double v) { return v; }
};

template<>
struct PromoteTraits<unsigned long, char>
{
    typedef unsigned long Promote;
    static Promote toPromote(unsigned long v) { return v; }
    static Promote toPromote(char v) { return v; }
};

template<>
struct PromoteTraits<unsigned long, unsigned char>
{
    typedef unsigned long Promote;
    static Promote toPromote(unsigned long v) { return v; }
    static Promote toPromote(unsigned char v) { return v; }
};

template<>
struct PromoteTraits<unsigned long, short int>
{
    typedef unsigned long Promote;
    static Promote toPromote(unsigned long v) { return v; }
    static Promote toPromote(short int v) { return v; }
};

template<>
struct PromoteTraits<unsigned long, short unsigned int>
{
    typedef unsigned long Promote;
    static Promote toPromote(unsigned long v) { return v; }
    static Promote toPromote(short unsigned int v) { return v; }
};

template<>
struct PromoteTraits<unsigned long, int>
{
    typedef unsigned long Promote;
    static Promote toPromote(unsigned long v) { return v; }
    static Promote toPromote(int v) { return v; }
};

template<>
struct PromoteTraits<unsigned long, unsigned int>
{
    typedef unsigned long Promote;
    static Promote toPromote(unsigned long v) { return v; }
    static Promote toPromote(unsigned int v) { return v; }
};

template<>
struct PromoteTraits<unsigned long, long>
{
    typedef unsigned long Promote;
    static Promote toPromote(unsigned long v) { return v; }
    static Promote toPromote(long v) { return v; }
};

template<>
struct PromoteTraits<unsigned long, unsigned long>
{
    typedef unsigned long Promote;
    static Promote toPromote(unsigned long v) { return v; }
};

template<>
struct PromoteTraits<unsigned long, float>
{
    typedef double Promote;
    static Promote toPromote(unsigned long v) { return v; }
    static Promote toPromote(float v) { return v; }
};

template<>
struct PromoteTraits<unsigned long, double>
{
    typedef double Promote;
    static Promote toPromote(unsigned long v) { return v; }
    static Promote toPromote(double v) { return v; }
};

template<>
struct PromoteTraits<unsigned long, long double>
{
    typedef long double Promote;
    static Promote toPromote(unsigned long v) { return v; }
    static Promote toPromote(long double v) { return v; }
};

template<>
struct PromoteTraits<float, char>
{
    typedef double Promote;
    static Promote toPromote(float v) { return v; }
    static Promote toPromote(char v) { return v; }
};

template<>
struct PromoteTraits<float, unsigned char>
{
    typedef double Promote;
    static Promote toPromote(float v) { return v; }
    static Promote toPromote(unsigned char v) { return v; }
};

template<>
struct PromoteTraits<float, short int>
{
    typedef double Promote;
    static Promote toPromote(float v) { return v; }
    static Promote toPromote(short int v) { return v; }
};

template<>
struct PromoteTraits<float, short unsigned int>
{
    typedef double Promote;
    static Promote toPromote(float v) { return v; }
    static Promote toPromote(short unsigned int v) { return v; }
};

template<>
struct PromoteTraits<float, int>
{
    typedef double Promote;
    static Promote toPromote(float v) { return v; }
    static Promote toPromote(int v) { return v; }
};

template<>
struct PromoteTraits<float, unsigned int>
{
    typedef double Promote;
    static Promote toPromote(float v) { return v; }
    static Promote toPromote(unsigned int v) { return v; }
};

template<>
struct PromoteTraits<float, long>
{
    typedef double Promote;
    static Promote toPromote(float v) { return v; }
    static Promote toPromote(long v) { return v; }
};

template<>
struct PromoteTraits<float, unsigned long>
{
    typedef double Promote;
    static Promote toPromote(float v) { return v; }
    static Promote toPromote(unsigned long v) { return v; }
};

template<>
struct PromoteTraits<float, float>
{
    typedef float Promote;
    static Promote toPromote(float v) { return v; }
};

template<>
struct PromoteTraits<float, double>
{
    typedef double Promote;
    static Promote toPromote(float v) { return v; }
    static Promote toPromote(double v) { return v; }
};

template<>
struct PromoteTraits<float, long double>
{
    typedef long double Promote;
    static Promote toPromote(float v) { return v; }
    static Promote toPromote(long double v) { return v; }
};

template<>
struct PromoteTraits<double, char>
{
    typedef double Promote;
    static Promote toPromote(double v) { return v; }
    static Promote toPromote(char v) { return v; }
};

template<>
struct PromoteTraits<double, unsigned char>
{
    typedef double Promote;
    static Promote toPromote(double v) { return v; }
    static Promote toPromote(unsigned char v) { return v; }
};

template<>
struct PromoteTraits<double, short int>
{
    typedef double Promote;
    static Promote toPromote(double v) { return v; }
    static Promote toPromote(short int v) { return v; }
};

template<>
struct PromoteTraits<double, short unsigned int>
{
    typedef double Promote;
    static Promote toPromote(double v) { return v; }
    static Promote toPromote(short unsigned int v) { return v; }
};

template<>
struct PromoteTraits<double, int>
{
    typedef double Promote;
    static Promote toPromote(double v) { return v; }
    static Promote toPromote(int v) { return v; }
};

template<>
struct PromoteTraits<double, unsigned int>
{
    typedef double Promote;
    static Promote toPromote(double v) { return v; }
    static Promote toPromote(unsigned int v) { return v; }
};

template<>
struct PromoteTraits<double, long>
{
    typedef double Promote;
    static Promote toPromote(double v) { return v; }
    static Promote toPromote(long v) { return v; }
};

template<>
struct PromoteTraits<double, unsigned long>
{
    typedef double Promote;
    static Promote toPromote(double v) { return v; }
    static Promote toPromote(unsigned long v) { return v; }
};

template<>
struct PromoteTraits<double, float>
{
    typedef double Promote;
    static Promote toPromote(double v) { return v; }
    static Promote toPromote(float v) { return v; }
};

template<>
struct PromoteTraits<double, double>
{
    typedef double Promote;
    static Promote toPromote(double v) { return v; }
};

template<>
struct PromoteTraits<double, long double>
{
    typedef long double Promote;
    static Promote toPromote(double v) { return v; }
    static Promote toPromote(long double v) { return v; }
};

template<>
struct PromoteTraits<long double, char>
{
    typedef long double Promote;
    static Promote toPromote(long double v) { return v; }
    static Promote toPromote(char v) { return v; }
};

template<>
struct PromoteTraits<long double, unsigned char>
{
    typedef long double Promote;
    static Promote toPromote(long double v) { return v; }
    static Promote toPromote(unsigned char v) { return v; }
};

template<>
struct PromoteTraits<long double, short int>
{
    typedef long double Promote;
    static Promote toPromote(long double v) { return v; }
    static Promote toPromote(short int v) { return v; }
};

template<>
struct PromoteTraits<long double, short unsigned int>
{
    typedef long double Promote;
    static Promote toPromote(long double v) { return v; }
    static Promote toPromote(short unsigned int v) { return v; }
};

template<>
struct PromoteTraits<long double, int>
{
    typedef long double Promote;
    static Promote toPromote(long double v) { return v; }
    static Promote toPromote(int v) { return v; }
};

template<>
struct PromoteTraits<long double, unsigned int>
{
    typedef long double Promote;
    static Promote toPromote(long double v) { return v; }
    static Promote toPromote(unsigned int v) { return v; }
};

template<>
struct PromoteTraits<long double, long>
{
    typedef long double Promote;
    static Promote toPromote(long double v) { return v; }
    static Promote toPromote(long v) { return v; }
};

template<>
struct PromoteTraits<long double, unsigned long>
{
    typedef long double Promote;
    static Promote toPromote(long double v) { return v; }
    static Promote toPromote(unsigned long v) { return v; }
};

template<>
struct PromoteTraits<long double, float>
{
    typedef long double Promote;
    static Promote toPromote(long double v) { return v; }
    static Promote toPromote(float v) { return v; }
};

template<>
struct PromoteTraits<long double, double>
{
    typedef long double Promote;
    static Promote toPromote(long double v) { return v; }
    static Promote toPromote(double v) { return v; }
};

template<>
struct PromoteTraits<long double, long double>
{
    typedef long double Promote;
    static Promote toPromote(long double v) { return v; }
};



#endif // VIGRA_NUMERICTRAITS_HXX

