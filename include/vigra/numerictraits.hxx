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
 
 
#ifndef VIGRA_NUMERICTRAITS_HXX
#define VIGRA_NUMERICTRAITS_HXX

#include <climits>
#include <limits>
#include <cfloat>
#include <complex>
#include "metaprogramming.hxx"
#include "sized_int.hxx"

/********************************************************/
/*                                                      */
/*                      NumericTraits                   */
/*                                                      */
/********************************************************/


/** \page NumericPromotionTraits Numeric and Promotion Traits

    Meta-information about arithmetic types.
    
    <UL style="list-style-image:url(documents/bullet.gif)">
    <LI> \ref NumericTraits
         <BR>&nbsp;&nbsp;&nbsp;<em>Unary traits for promotion, conversion, creation of arithmetic objects</em>
    <LI> \ref PromoteTraits
         <BR>&nbsp;&nbsp;&nbsp;<em>Binary traits for promotion of arithmetic objects</em>
    <LI> \ref SquareRootTraits
         <BR>&nbsp;&nbsp;&nbsp;<em>Unary traits for the calculation of the square root of arithmetic objects</em>
    <LI> \ref NormTraits
         <BR>&nbsp;&nbsp;&nbsp;<em>Unary traits for the calculation of the norm and squared norm of arithmetic objects</em>
    </UL>
    
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
*/

/** \page NumericTraits template<> struct NumericTraits<ArithmeticType>

    Unary traits for promotion, conversion, creation of arithmetic objects.

    <b>\#include</b> 
    \<vigra/numerictraits.hxx\>

    This traits class is used derive important properties of
    an arithmetic type. Consider the following algorithm:
    
    \code
    // calculate the sum of a sequence of bytes
    int sumBytes(unsigned char * begin, unsigned char * end)
    {
        int result = 0;
        for(; begin != end; ++begin)  result += *begin;
        return result;
    }
    \endcode 
    
    The return type of this function can not be 'unsigned char' because
    the summation would very likely overflow. Since we know the source
    type, we can easily choose 'int' as an appropriate return type.
    Likewise, we would have chosen 'float' if we had to sum a 
    sequence of floats. If we want to make this 
    algorithm generic, we would like to derive the appropriate return 
    type automatically. This can be done with NumericTraits. 
    The code would look like this (we use \ref DataAccessors to 
    read the data from the sequence):
    
    \code
    // calculate the sum of any sequence
    template <class Iterator, class Accessor>
    typename vigra::NumericTraits<typename Accessor::value_type>::Promote
    sumSequence(Iterator begin, Iterator end, Accessor a)
    {
        // an abbreviation
        typedef vigra::NumericTraits<typename Accessor::value_type>  SrcTraits;
        
        // find out result type
        typedef typename SrcTraits::Promote ResultType;
      
        // init result to zero
        ResultType result = vigra::NumericTraits<ResultType>::zero();
    
        for(; begin != end; ++begin)
        {  
            // cast current item to ResultType and add
            result += SrcTraits::toPromote(a(begin));
        }
        
        return result;
    }
    \endcode
    
    In this example NumericTraits is not only used to deduce the 
    ReturnType of the operation, but also to initialize it with the
    constant 'zero'. This is necessary since we do not know in general,
    which expression must be used to obtain a zero of some arbitrary
    type - '<TT>ResultType result = 0;</TT>' would only work if the 
    ResultType had an constructor taking an '<TT>int</TT>' argument, and we 
    would not even have any guarantee as to what the semantics of this
    constructor are. In addition, the traits are used to cast the 
    source type into the promote type.
    
    Similarly, an algorithm that needs multiplication would use the 
    return type <TT>RealPromote</TT> and the functions <TT>one()</TT> and
    <TT>toRealPromote()</TT>. The following members are defined in 
    <b> <TT>NumericTraits<ArithmeticType></TT></b>:
    
    <table>
    <tr><td>
    <b> <TT>typedef ... Type;</TT></b>
    </td><td>
    
            the type itself 
        
    </td></tr>
    <tr><td>
    <b> <TT>typedef ... Promote;</TT></b>
    </td><td>
    
            promote type for addition and subtraction 
        
    </td></tr>
    <tr><td>
    <b> <TT>typedef ... RealPromote;</TT></b>
    </td><td>
            promote type for multiplication and division with a real number
    
    (only defined if <TT>ArithmeticType</TT> supports these operations) 
    
    </td></tr>
    <tr><td>
    <b> <TT>typedef ... ComplexPromote;</TT></b>
    </td><td>
    
            promote type for complex arithmetic 
        
    </td></tr>
    <tr><td>
    <b> <TT>typedef ... ValueType;</TT></b>
    </td><td>
    
            for scalar types: the type itself<br>
            otherwise: typename Type::value_type (if defined)
        
    </td></tr>
    <tr><td>
    <b> <TT>static Promote toPromote(ArithmeticType v);</TT></b>
    </td><td>
        convert to <TT>Promote</TT> type 
    
    </td></tr>
    <tr><td>
    <b> <TT>static RealPromote toRealPromote(ArithmeticType v);</TT></b>
    </td><td>
        convert to <TT>RealPromote</TT> type 

    (only defined if <TT>ArithmeticType</TT> supports multiplication) 
    
    </td></tr>
    <tr><td>
    <b> <TT>static ArithmeticType fromPromote(Promote v);</TT></b> 
    </td><td>
        convert from <TT>Promote</TT> type
    
    if <TT>v</TT> is outside the range of <TT>ArithmeticType</TT> it is clipped;

    </td></tr>
    <tr><td>
    <b> <TT>static ArithmeticType fromRealPromote(RealPromote v);</TT></b>
    </td><td>
        convert from <TT>RealPromote</TT> type 
    
    (only defined if 
    <TT>ArithmeticType</TT> supports multiplication)
    
    if <TT>ArithmeticType</TT> is an integral type, the result is rounded 
    
    if <TT>v</TT> is outside the range of <TT>ArithmeticType</TT> it is clipped
    
    </td></tr>
    <tr><td>
    <b> <TT>static ArithmeticType zero();</TT></b>
    </td><td>
    create neutral element of addition
    
    i.e. <TT>(ArithmeticType a = ...,</TT> 
    <TT>  a + NumericTraits<ArithmeticType>::zero() == a)</TT> 
    must always yield <TT>true</TT> 
    
    </td></tr>
    <tr><td>
    <b> <TT>static ArithmeticType nonZero();</TT></b>
    </td><td>
    create a non-zero element (if multiplication is defined, this yields one())
    
    i.e. <TT>(ArithmeticType a = ...,</TT> 
    <TT>  a + NumericTraits<ArithmeticType>::nonZero() == a)</TT> 
    must always yield <TT>false</TT> 
    
    </td></tr>
    <tr><td>
    <b> <TT>static ArithmeticType min();</TT></b>
    </td><td>
    the smallest number representable in this type.<br>
    Only available if isOrdered is VigraTrueType. For integral types,
    this equals <TT>INT_MIN</TT> etc., for real valued types it is <TT>-FLT_MAX</TT>
    etc. (<b>not</b> <TT>FLT_MIN</TT> -- this is the smallest positive <tt>float</tt>)
    
    </td></tr>
    <tr><td>
    <b> <TT>static ArithmeticType max();</TT></b>
    </td><td>
    the largest number representable in this type.<br>
    Only available if isOrdered is VigraTrueType. For integral types,
    this equals <TT>INT_MAX</TT> etc., for real valued types it is <TT>FLT_MAX</TT>
    etc.
    
    </td></tr>
    <tr><td>
    <b> <TT>static ArithmeticType one();</TT></b>
    </td><td>
    create neutral element of multiplication 
    
    (only defined if <TT>ArithmeticType</TT> supports multiplication)
    
    i.e. <TT>(ArithmeticType a = ...,</TT> 
    <TT>  a * NumericTraits<ArithmeticType>::one() == a)</TT> 
    must always yield <TT>true</TT> 
    
    </td></tr>
    <tr><td>
    <b> <TT>typedef ... isIntegral;</TT></b>
    </td><td>
        VigraTrueType if <TT>ArithmeticType</TT> is an integral type, 
        VigraFalseType otherwise 
    
    </td></tr>
    <tr><td>
    <b> <TT>typedef ... isScalar;</TT></b>
    </td><td>
        VigraTrueType if <TT>ArithmeticType</TT> is a scalar type, 
        VigraFalseType otherwise 
    
    </td></tr>
    <tr><td>
    <tr><td>
    <b> <TT>typedef ... isSigned;</TT></b>
    </td><td>
        VigraTrueType if <TT>ArithmeticType</TT> is a signed type, 
        VigraFalseType otherwise 
    
    </td></tr>
    <tr><td>
    <tr><td>
    <b> <TT>typedef ... isOrdered;</TT></b>
    </td><td>
        VigraTrueType if <TT>ArithmeticType</TT> supports operator<(), 
        VigraFalseType otherwise 
    
    </td></tr>
    <tr><td>
    <b> <TT>typedef ... isComplex;</TT></b>
    </td><td>
        VigraTrueType if <TT>ArithmeticType</TT> is a complex number, 
        VigraFalseType otherwise 
    
    </td></tr>
    <tr><td>
    </table>
    
    NumericTraits for the built-in types are defined in <b>\#include</b> 
    \<vigra/numerictraits.hxx\>
    
    Namespace: vigra
    
*/

/** \page PromoteTraits template<> struct PromoteTraits<ArithmeticType1, ArithmeticType2>

    Binary traits for promotion of arithmetic objects.
    
    <b>\#include</b> 
    \<vigra/numerictraits.hxx\>

    This traits class is used to determine the appropriate result type
    of arithmetic expressions which depend of two arguments. Consider
    the following function:
    
    \code
    template <class T>
    T min(T t1, T t2)
    {
        return (t1 < t2) ? t1 : t2;
    }
    \endcode
    
    This template is only applicable if both arguments have the same
    type. However, sometimes we may want to use the function in cases
    where the argument types differ. Then we can deduce the appropriate
    return type by using <TT>PromoteTraits</TT>:
    
    \code
    template <class T1, class T2>
    typename vigra::PromoteTraits<T1, T2>::Promote
    min(T1 t1, T2 t2)
    {
        return (t1 < t2) ? vigra::PromoteTraits<T1, T2>::toPromote(t1) : 
                           vigra::PromoteTraits<T1, T2>::toPromote(t2);
    }    
    \endcode
    
    In addition, the traits class provide static functions to cast the
    arguments to the promote type. For example, if <TT>T1</TT> were <TT>int</TT> and 
    <TT>T2</TT> were <TT>float</TT>, the <TT>Promote</TT> type would be <TT>float</TT>. 
    The following members are defined in 
    <b> <TT>PromoteTraits<ArithmeticType1, ArithmeticType2></TT></b>:
    
    <table>
    <tr>
    <td>
    <b> <TT>typedef ... Promote;</TT></b>
    </td><td>
            promote type 
    </td></tr>
    <tr><td>
    <b> <TT>static Promote toPromote(ArithmeticType1 v);</TT></b> 
    
    <b> <TT>static Promote toPromote(ArithmeticType2 v);</TT></b>
    </td><td>
        convert to <TT>Promote</TT> type 
    </td></tr>
    </table>
    
    PromoteTraits for the built-in types are defined in <b>\#include</b> 
    \<vigra/numerictraits.hxx\>
    
    Namespace: vigra
*/

/** \page SquareRootTraits template<> struct SquareRootTraits<ArithmeticType>

    Unary traits for the calculation of the square root of arithmetic objects.
    
    <b>\#include</b> 
    \<vigra/numerictraits.hxx\>

    This traits class is used to determine appropriate argument and result types
    for the function sqrt(). These traits are typically used like this:
    
    \code
    ArithmeticType t = ...;
    SquareRootTraits<ArithmeticType>::SquareRootResult r = 
          sqrt((SquareRootTraits<ArithmeticType>::SquareRootArgument)t);
    \endcode
    
    This approach avoids 'ambiguous overload errors' when taking the square root of 
    an integer type. It also takes care of determining the proper result of the
    sqrt() function of \ref vigra::FixedPoint and of the norm() function, when
    it is implemented via sqrt(squaredNorm(x)).
    The following members are defined in <b> <TT>SquareRootTraits<ArithmeticType></TT></b>:
    
    <table>
    <tr><td>
    <b> <TT>typedef ArithmeticType Type;</TT></b>
    </td><td>
            the type itself
    </td></tr>
    <tr><td>
    <b> <TT>typedef ... SquareRootArgument;</TT></b>
    </td><td>
            required argument type for srqt(), i.e. <tt>sqrt((SquareRootArgument)x)</tt>
    </td></tr>
    <tr><td>
    <b> <TT>typedef ... SquareRootResult;</TT></b>
    </td><td>
            result of <tt>sqrt((SquareRootArgument)x)</tt>
    </td></tr>
    </table>
    
    NormTraits for the built-in types are defined in <b>\#include</b> 
    \<vigra/numerictraits.hxx\>
    
    Namespace: vigra
*/

/** \page NormTraits template<> struct NormTraits<ArithmeticType>

    Unary traits for the calculation of the norm and squared norm of arithmetic objects.
    
    <b>\#include</b> 
    \<vigra/numerictraits.hxx\>

    This traits class is used to determine appropriate result types
    for the functions norm() and squaredNorm(). These functions are always 
    declared like this (where <tt>ArithmeticType</tt> is a type that supports a norm):
    
    \code
    NormTraits<ArithmeticType>::NormType        norm(ArithmeticType const & t);
    NormTraits<ArithmeticType>::SquaredNormType squaredNorm(ArithmeticType const & t);
    \endcode
    
    The following members are defined in <b> <TT>NormTraits<ArithmeticType></TT></b>:
    
    <table>
    <tr><td>
    <b> <TT>typedef ArithmeticType Type;</TT></b>
    </td><td>
            the type itself
    </td></tr>
    <tr><td>
    <b> <TT>typedef ... SquaredNormType;</TT></b>
    </td><td>
            result of <tt>squaredNorm(ArithmeticType)</tt>
    </td></tr>
    <tr><td>
    <b> <TT>typedef ... NormType;</TT></b>
    </td><td>
            result of <tt>norm(ArithmeticType)</tt><br>
            Usually equal to <tt>SquareRootTraits<SquaredNormType>::SquareRootResult</tt>
    </td></tr>
    </table>
    
    NormTraits for the built-in types are defined in <b>\#include</b> 
    \<vigra/numerictraits.hxx\>
    
    Namespace: vigra
*/

namespace vigra {

struct Error_NumericTraits_not_specialized_for_this_case { };
struct Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char { };

template<class A>
struct NumericTraits
{
    typedef Error_NumericTraits_not_specialized_for_this_case Type;
    typedef Error_NumericTraits_not_specialized_for_this_case Promote;
    typedef Error_NumericTraits_not_specialized_for_this_case UnsignedPromote;
    typedef Error_NumericTraits_not_specialized_for_this_case RealPromote;
    typedef Error_NumericTraits_not_specialized_for_this_case ComplexPromote;
    typedef Error_NumericTraits_not_specialized_for_this_case ValueType;

    typedef Error_NumericTraits_not_specialized_for_this_case isScalar;
    typedef Error_NumericTraits_not_specialized_for_this_case isIntegral;
    typedef Error_NumericTraits_not_specialized_for_this_case isSigned;
    typedef Error_NumericTraits_not_specialized_for_this_case isOrdered;
    typedef Error_NumericTraits_not_specialized_for_this_case isComplex;
};

template<>
struct NumericTraits<char>
{
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char Type;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char Promote;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char UnsignedPromote;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char RealPromote;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char ComplexPromote;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char ValueType;

    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char isScalar;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char isIntegral;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char isSigned;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char isOrdered;
    typedef Error_NumericTraits_char_is_not_a_numeric_type__use_signed_char_or_unsigned_char isComplex;
};

#ifndef NO_BOOL
template<>
struct NumericTraits<bool>
{
    typedef bool Type;
    typedef int Promote;
    typedef unsigned int UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraFalseType isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;
    
    static bool zero() { return false; }
    static bool one() { return true; }
    static bool nonZero() { return true; }
    static bool min() { return false; }
    static bool max() { return true; }
    
#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { minConst = false , maxConst = true };
#else
    static const bool minConst = false;
    static const bool maxConst = true;
#endif
    
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
    typedef unsigned int UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;
    
    static signed char zero() { return 0; }
    static signed char one() { return 1; }
    static signed char nonZero() { return 1; }
    static signed char min() { return SCHAR_MIN; }
    static signed char max() { return SCHAR_MAX; }
    
#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { minConst = SCHAR_MIN, maxConst = SCHAR_MIN };
#else
    static const signed char minConst = SCHAR_MIN;
    static const signed char maxConst = SCHAR_MIN;
#endif
    
    static Promote toPromote(signed char v) { return v; }
    static RealPromote toRealPromote(signed char v) { return v; }
    static signed char fromPromote(Promote v) { 
        return ((v < SCHAR_MIN) ? SCHAR_MIN : (v > SCHAR_MAX) ? SCHAR_MAX : v); 
    }
    static signed char fromRealPromote(RealPromote v) {
        return ((v < 0.0) 
                   ? ((v < (RealPromote)SCHAR_MIN) 
                       ? SCHAR_MIN 
                       : static_cast<signed char>(v - 0.5)) 
                   : (v > (RealPromote)SCHAR_MAX) 
                       ? SCHAR_MAX 
                       : static_cast<signed char>(v + 0.5)); 
    }
};

template<>
struct NumericTraits<unsigned char>
{
    typedef unsigned char Type;
    typedef int Promote;
    typedef unsigned int UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraFalseType isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;
    
    static unsigned char zero() { return 0; }
    static unsigned char one() { return 1; }
    static unsigned char nonZero() { return 1; }
    static unsigned char min() { return 0; }
    static unsigned char max() { return UCHAR_MAX; }
    
#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { minConst = 0, maxConst = UCHAR_MAX };
#else
    static const unsigned char minConst = 0;
    static const unsigned char maxConst = UCHAR_MAX;
#endif
    
    static Promote toPromote(unsigned char v) { return v; }
    static RealPromote toRealPromote(unsigned char v) { return v; }
    static unsigned char fromPromote(Promote const & v) { 
        return Type((v < 0) 
             ? 0 
             : (v > (Promote)UCHAR_MAX) 
                    ? UCHAR_MAX
                    : v); 
    }
    static unsigned char fromRealPromote(RealPromote const & v) {
            return Type((v < 0.0) 
                     ? 0 
                     : ((v > (RealPromote)UCHAR_MAX) 
                         ? UCHAR_MAX 
                         : v + 0.5));
    }
};

template<>
struct NumericTraits<short int>
{
    typedef short int Type;
    typedef int Promote;
    typedef unsigned int UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;
    
    static short int zero() { return 0; }
    static short int one() { return 1; }
    static short int nonZero() { return 1; }
    static short int min() { return SHRT_MIN; }
    static short int max() { return SHRT_MAX; }
    
#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { minConst = SHRT_MIN, maxConst = SHRT_MAX };
#else
    static const short int minConst = SHRT_MIN;
    static const short int maxConst = SHRT_MAX;
#endif
    
    static Promote toPromote(short int v) { return v; }
    static RealPromote toRealPromote(short int v) { return v; }
    static short int fromPromote(Promote v) { 
        return ((v < SHRT_MIN) ? SHRT_MIN : 
                (v > SHRT_MAX) ? SHRT_MAX : v); 
    }
    static short int fromRealPromote(RealPromote v) {
        return ((v < 0.0) 
                 ? ((v < (RealPromote)SHRT_MIN) 
                     ? SHRT_MIN 
                     : static_cast<short int>(v - 0.5)) 
                 : ((v > (RealPromote)SHRT_MAX) 
                     ? SHRT_MAX 
                     : static_cast<short int>(v + 0.5))); 
    }
};

template<>
struct NumericTraits<short unsigned int>
{
    typedef short unsigned int Type;
    typedef int Promote;
    typedef unsigned int UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraFalseType isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;

    static short unsigned int zero() { return 0; }
    static short unsigned int one() { return 1; }
    static short unsigned int nonZero() { return 1; }
    static short unsigned int min() { return 0; }
    static short unsigned int max() { return USHRT_MAX; }
    
#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { minConst = 0, maxConst = USHRT_MAX };
#else
    static const short unsigned int minConst = 0;
    static const short unsigned int maxConst = USHRT_MAX;
#endif

    static Promote toPromote(short unsigned int v) { return v; }
    static RealPromote toRealPromote(short unsigned int v) { return v; }
    static short unsigned int fromPromote(Promote v) { 
        return Type((v < 0) 
              ? 0 
              : (v > USHRT_MAX) 
                   ? USHRT_MAX 
                   : v); 
    }
    static short unsigned int fromRealPromote(RealPromote v) {
            return Type((v < 0.0) 
                     ? 0 
                     : ((v > (RealPromote)USHRT_MAX) 
                         ? USHRT_MAX 
                         : v + 0.5));
    }
};

template<>
struct NumericTraits<int>
{
    typedef int Type;
    typedef int Promote;
    typedef unsigned int UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;

    static int zero() { return 0; }
    static int one() { return 1; }
    static int nonZero() { return 1; }
    static int min() { return INT_MIN; }
    static int max() { return INT_MAX; }
    
#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { minConst = INT_MIN, maxConst = INT_MAX };
#else
    static const int minConst = INT_MIN;
    static const int maxConst = INT_MAX;
#endif

    static Promote toPromote(int v) { return v; }
    static RealPromote toRealPromote(int v) { return v; }
    static int fromPromote(Promote v) { return v; }
    static int fromRealPromote(RealPromote v) {
        return ((v < 0.0) 
                 ? ((v < (RealPromote)INT_MIN) 
                     ? INT_MIN 
                     : static_cast<int>(v - 0.5)) 
                 : ((v > (RealPromote)INT_MAX) 
                     ? INT_MAX 
                     : static_cast<int>(v + 0.5))); 
    }
};

template<>
struct NumericTraits<unsigned int>
{
    typedef unsigned int Type;
    typedef unsigned int Promote;
    typedef unsigned int UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraFalseType isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;
    
    static unsigned int zero() { return 0; }
    static unsigned int one() { return 1; }
    static unsigned int nonZero() { return 1; }
    static unsigned int min() { return 0; }
    static unsigned int max() { return UINT_MAX; }
    
#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { minConst = 0, maxConst = UINT_MAX };
#else
    static const unsigned int minConst = 0;
    static const unsigned int maxConst = UINT_MAX;
#endif

    static Promote toPromote(unsigned int v) { return v; }
    static RealPromote toRealPromote(unsigned int v) { return v; }
    static unsigned int fromPromote(Promote v) { return v; }
    static unsigned int fromRealPromote(RealPromote v) {
            return ((v < 0.0) 
                     ? 0 
                     : ((v > (RealPromote)UINT_MAX) 
                         ? UINT_MAX 
                         : static_cast<unsigned int>(v + 0.5)));
    }
};

template<>
struct NumericTraits<long>
{
    typedef long Type;
    typedef long Promote;
    typedef unsigned long UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;
    
    static long zero() { return 0; }
    static long one() { return 1; }
    static long nonZero() { return 1; }
    static long min() { return LONG_MIN; }
    static long max() { return LONG_MAX; }
    
#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { minConst = LONG_MIN, maxConst = LONG_MAX };
#else
    static const long minConst = LONG_MIN;
    static const long maxConst = LONG_MAX;
#endif

    static Promote toPromote(long v) { return v; }
    static RealPromote toRealPromote(long v) { return v; }
    static long fromPromote(Promote v) { return v; }
    static long fromRealPromote(RealPromote v) {
        return ((v < 0.0) 
                 ? ((v < (RealPromote)LONG_MIN) 
                     ? LONG_MIN 
                     : static_cast<long>(v - 0.5)) 
                 : ((v > (RealPromote)LONG_MAX) 
                     ? LONG_MAX 
                     : static_cast<long>(v + 0.5))); 
    }
};

template<>
struct NumericTraits<unsigned long>
{
    typedef unsigned long Type;
    typedef unsigned long Promote;
    typedef unsigned long UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraFalseType isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;
    
    static unsigned long zero() { return 0; }
    static unsigned long one() { return 1; }
    static unsigned long nonZero() { return 1; }
    static unsigned long min() { return 0; }
    static unsigned long max() { return ULONG_MAX; }
    
#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { minConst = 0, maxConst = ULONG_MAX };
#else
    static const unsigned long minConst = 0;
    static const unsigned long maxConst = ULONG_MAX;
#endif

    static Promote toPromote(unsigned long v) { return v; }
    static RealPromote toRealPromote(unsigned long v) { return v; }
    static unsigned long fromPromote(Promote v) { return v; }
    static unsigned long fromRealPromote(RealPromote v) {
            return ((v < 0.0) 
                     ? 0 
                     : ((v > (RealPromote)ULONG_MAX) 
                         ? ULONG_MAX 
                         : static_cast<unsigned long>(v + 0.5)));
    }
};

#ifdef LLONG_MAX
template<>
struct NumericTraits<long long>
{
    typedef long long Type;
    typedef long long Promote;
    typedef unsigned long long UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;
    
    static long long zero() { return 0; }
    static long long one() { return 1; }
    static long long nonZero() { return 1; }
    static long long min() { return LLONG_MIN; }
    static long long max() { return LLONG_MAX; }
    
#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { minConst = LLONG_MIN, maxConst = LLONG_MAX };
#else
    static const long long minConst = LLONG_MIN;
    static const long long maxConst = LLONG_MAX;
#endif

    static Promote toPromote(long long v) { return v; }
    static RealPromote toRealPromote(long long v) { return (RealPromote)v; }
    static long long fromPromote(Promote v) { return v; }
    static long long fromRealPromote(RealPromote v) {
        return ((v < 0.0) 
                 ? ((v < (RealPromote)LLONG_MIN) 
                     ? LLONG_MIN 
                     : static_cast<long long>(v - 0.5)) 
                 : ((v > (RealPromote)LLONG_MAX) 
                     ? LLONG_MAX 
                     : static_cast<long long>(v + 0.5))); 
    }
};

template<>
struct NumericTraits<unsigned long long>
{
    typedef unsigned long long Type;
    typedef unsigned long long Promote;
    typedef unsigned long long UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    typedef VigraTrueType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraFalseType isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;
    
    static unsigned long long zero() { return 0; }
    static unsigned long long one() { return 1; }
    static unsigned long long nonZero() { return 1; }
    static unsigned long long min() { return 0; }
    static unsigned long long max() { return ULLONG_MAX; }
    
#ifdef NO_INLINE_STATIC_CONST_DEFINITION
    enum { minConst = 0, maxConst = ULLONG_MAX };
#else
    static const unsigned long long minConst = 0;
    static const unsigned long long maxConst = ULLONG_MAX;
#endif

    static Promote toPromote(unsigned long long v) { return v; }
    static RealPromote toRealPromote(unsigned long long v) { return (RealPromote)v; }
    static unsigned long long fromPromote(Promote v) { return v; }
    static unsigned long long fromRealPromote(RealPromote v) {
            return ((v < 0.0) 
                     ? 0 
                     : ((v > (RealPromote)ULLONG_MAX) 
                         ? ULONG_MAX 
                         : static_cast<unsigned long long>(v + 0.5)));
    }
};
#endif // LLONG_MAX

template<>
struct NumericTraits<float>
{
    typedef float Type;
    typedef float Promote;
    typedef float UnsignedPromote;
    typedef float RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;
    
    typedef VigraFalseType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;
    
    static float zero() { return 0.0; }
    static float one() { return 1.0; }
    static float nonZero() { return 1.0; }
    static float epsilon() { return FLT_EPSILON; }
    static float smallestPositive() { return FLT_MIN; }
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
    typedef double UnsignedPromote;
    typedef double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    typedef VigraFalseType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;
    
    static double zero() { return 0.0; }
    static double one() { return 1.0; }
    static double nonZero() { return 1.0; }
    static double epsilon() { return DBL_EPSILON; }
    static double smallestPositive() { return DBL_MIN; }
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
    typedef long double UnsignedPromote;
    typedef long double RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef Type ValueType;

    typedef VigraFalseType isIntegral;
    typedef VigraTrueType isScalar;
    typedef VigraTrueType isSigned;
    typedef VigraTrueType isOrdered;
    typedef VigraFalseType isComplex;
    
    static long double zero() { return 0.0; }
    static long double one() { return 1.0; }
    static long double nonZero() { return 1.0; }
    static long double epsilon() { return LDBL_EPSILON; }
    static long double smallestPositive() { return LDBL_MIN; }
    static long double min() { return -LDBL_MAX; }
    static long double max() { return LDBL_MAX; }

    static Promote toPromote(long double v) { return v; }
    static RealPromote toRealPromote(long double v) { return v; }
    static long double fromPromote(Promote v) { return v; }
    static long double fromRealPromote(RealPromote v) { return v; }
};

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template<class T>
struct NumericTraits<std::complex<T> >
{
    typedef std::complex<T> Type;
    typedef std::complex<typename NumericTraits<T>::Promote> Promote;
    typedef std::complex<typename NumericTraits<T>::UnsignedPromote> UnsignedPromote;
    typedef std::complex<typename NumericTraits<T>::RealPromote> RealPromote;
    typedef std::complex<RealPromote> ComplexPromote;
    typedef T ValueType;

    typedef VigraFalseType isIntegral;
    typedef VigraFalseType isScalar;
    typedef typename NumericTraits<T>::isSigned isSigned;
    typedef VigraFalseType isOrdered;
    typedef VigraTrueType isComplex;
    
    static Type zero() { return Type(0.0); }
    static Type one() { return Type(1.0); }
    static Type nonZero() { return one(); }
    static Type epsilon() { return Type(NumericTraits<T>::epsilon()); }
    static Type smallestPositive() { return Type(NumericTraits<T>::smallestPositive()); }

    static Promote toPromote(Type const & v) { return v; }
    static Type fromPromote(Promote const & v) { return v; }
    static Type fromRealPromote(RealPromote v) { return Type(v); }
};

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

/********************************************************/
/*                                                      */
/*                    SquareRootTraits                  */
/*                                                      */
/********************************************************/

template<class T>
struct SquareRootTraits
{
    typedef T                                                    Type;
    typedef typename NumericTraits<T>::RealPromote               SquareRootResult;
    typedef typename NumericTraits<T>::RealPromote               SquareRootArgument;
};


/********************************************************/
/*                                                      */
/*                       NormTraits                     */
/*                                                      */
/********************************************************/

struct Error_NormTraits_not_specialized_for_this_case { };

template<class T>
struct NormTraits
{
    typedef T                                                Type;
    typedef Error_NormTraits_not_specialized_for_this_case   SquaredNormType;
    typedef Error_NormTraits_not_specialized_for_this_case   NormType;
};

#define VIGRA_DEFINE_NORM_TRAITS(T) \
    template <> struct NormTraits<T> { \
        typedef T Type; \
        typedef NumericTraits<T>::Promote SquaredNormType; \
        typedef T NormType; \
    };

VIGRA_DEFINE_NORM_TRAITS(bool)
VIGRA_DEFINE_NORM_TRAITS(signed char)
VIGRA_DEFINE_NORM_TRAITS(unsigned char)
VIGRA_DEFINE_NORM_TRAITS(short)
VIGRA_DEFINE_NORM_TRAITS(unsigned short)
VIGRA_DEFINE_NORM_TRAITS(int)
VIGRA_DEFINE_NORM_TRAITS(unsigned int)
VIGRA_DEFINE_NORM_TRAITS(long)
VIGRA_DEFINE_NORM_TRAITS(unsigned long)
VIGRA_DEFINE_NORM_TRAITS(float)
VIGRA_DEFINE_NORM_TRAITS(double)
VIGRA_DEFINE_NORM_TRAITS(long double)

#ifdef LLONG_MAX
VIGRA_DEFINE_NORM_TRAITS(long long)
VIGRA_DEFINE_NORM_TRAITS(unsigned long long)
#endif // LLONG_MAX

#undef VIGRA_DEFINE_NORM_TRAITS

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template<class T>
struct NormTraits<std::complex<T> >
{
    typedef std::complex<T>                                              Type;
    typedef typename NormTraits<T>::SquaredNormType                      SquaredNormType;
    typedef typename SquareRootTraits<SquaredNormType>::SquareRootResult NormType;
};

#endif // NO_PARTIAL_TEMPLATE_SPECIALIZATION

/********************************************************/
/*                                                      */
/*                      PromoteTraits                   */
/*                                                      */
/********************************************************/

namespace detail {

template <class T, class U>
struct PromoteType
{
    static T & t();
    static U & u();
    // let C++ figure out the promote type by adding a T and an U
    typedef typename SizeToType<sizeof(*typeToSize(t() + u()))>::result Promote;
    static Promote toPromote(T t) { return Promote(t); }
    static Promote toPromote(U u) { return Promote(u); }
};


template <class T>
struct PromoteType<T, T>
{
    static T & t();
    // let C++ figure out the promote type by adding two Ts
    typedef typename SizeToType<sizeof(*typeToSize(t() + t()))>::result Promote;
    static Promote toPromote(T t) { return Promote(t); }
};

} // namespace detail

struct Error_PromoteTraits_not_specialized_for_this_case { };

template<class A, class B>
struct PromoteTraits
{
    typedef Error_PromoteTraits_not_specialized_for_this_case Promote;
};

#include "promote_traits.hxx"

#ifndef NO_PARTIAL_TEMPLATE_SPECIALIZATION

template <class T>
struct PromoteTraits<std::complex<T>, std::complex<T> >
{
    typedef std::complex<typename PromoteTraits<T, T>::Promote> Promote;
    static Promote toPromote(std::complex<T> const & v) { return v; }
};

template <class T1, class T2>
struct PromoteTraits<std::complex<T1>, std::complex<T2> >
{
    typedef std::complex<typename PromoteTraits<T1, T2>::Promote> Promote;
    static Promote toPromote(std::complex<T1> const & v) { return v; }
    static Promote toPromote(std::complex<T2> const & v) { return v; }
};

template <class T1, class T2>
struct PromoteTraits<std::complex<T1>, T2 >
{
    typedef std::complex<typename PromoteTraits<T1, T2>::Promote> Promote;
    static Promote toPromote(std::complex<T1> const & v) { return v; }
    static Promote toPromote(T2 const & v) { return Promote(v); }
};

template <class T1, class T2>
struct PromoteTraits<T1, std::complex<T2> >
{
    typedef std::complex<typename PromoteTraits<T1, T2>::Promote> Promote;
    static Promote toPromote(T1 const & v) { return Promote(v); }
    static Promote toPromote(std::complex<T2> const & v) { return v; }
};

#endif

namespace detail {

template <class T>
struct RequiresExplicitCast {
    template <class U>
    static U const & cast(U const & v)
        { return v; }
};

#if !defined(_MSC_VER) || _MSC_VER >= 1300
#  define VIGRA_SPECIALIZED_CAST(type) \
    template <> \
    struct RequiresExplicitCast<type> { \
        static type cast(float v) \
            { return NumericTraits<type>::fromRealPromote(v); } \
        static type cast(double v) \
            { return NumericTraits<type>::fromRealPromote(v); } \
        static type cast(type v) \
            { return v; } \
        template <class U> \
        static type cast(U v) \
            { return static_cast<type>(v); } \
 \
    };
#else
#  define VIGRA_SPECIALIZED_CAST(type) \
    template <> \
    struct RequiresExplicitCast<type> { \
        static type cast(float v) \
            { return NumericTraits<type>::fromRealPromote(v); } \
        static type cast(double v) \
            { return NumericTraits<type>::fromRealPromote(v); } \
        static type cast(signed char v) \
            { return v; } \
        static type cast(unsigned char v) \
            { return v; } \
        static type cast(short v) \
            { return v; } \
        static type cast(unsigned short v) \
            { return v; } \
        static type cast(int v) \
            { return v; } \
        static type cast(unsigned int v) \
            { return v; } \
        static type cast(long v) \
            { return v; } \
        static type cast(unsigned long v) \
            { return v; } \
    };
#endif


VIGRA_SPECIALIZED_CAST(signed char)
VIGRA_SPECIALIZED_CAST(unsigned char)
VIGRA_SPECIALIZED_CAST(short)
VIGRA_SPECIALIZED_CAST(unsigned short)
VIGRA_SPECIALIZED_CAST(int)
VIGRA_SPECIALIZED_CAST(unsigned int)
VIGRA_SPECIALIZED_CAST(long)
VIGRA_SPECIALIZED_CAST(unsigned long)

template <>
struct RequiresExplicitCast<bool> {
    template <class U>
    static bool cast(U v)
    { return v == NumericTraits<U>::zero()
                ? false
                : true; }
};

template <>
struct RequiresExplicitCast<float> {
    static float cast(int v)
        { return (float)v; }

    static float cast(unsigned int v)
        { return (float)v; }

    static float cast(long v)
        { return (float)v; }

    static float cast(unsigned long v)
        { return (float)v; }

    static float cast(long long v)
        { return (float)v; }

    static float cast(unsigned long long v)
        { return (float)v; }

    static float cast(double v)
        { return (float)v; }

    static float cast(long double v)
        { return (float)v; }

    template <class U>
    static U cast(U v)
        { return v; }
};

template <>
struct RequiresExplicitCast<double> {
    static double cast(Int64 v)
        { return (double)v; }

    static double cast(UInt64 v)
        { return (double)v; }

    template <class U>
    static U cast(U v)
        { return v; }
};

#undef VIGRA_SPECIALIZED_CAST

} // namespace detail



} // namespace vigra

#endif // VIGRA_NUMERICTRAITS_HXX

