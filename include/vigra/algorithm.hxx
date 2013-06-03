/************************************************************************/
/*                                                                      */
/*               Copyright 2010-2011 by Ullrich Koethe                  */
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

#ifndef VIGRA_ALGORITHM_HXX
#define VIGRA_ALGORITHM_HXX

#include "sized_int.hxx"
#include "numerictraits.hxx"
#include "inspector_passes.hxx"
#include <algorithm>
#include <functional>
#include <iterator>

namespace vigra {

/** \addtogroup MathFunctions
*/
//@{
    /** \brief Find the minimum element in a sequence.
    
        The function returns the iterator referring to the minimum element.
        This is identical to the function <tt>std::min_element()</tt>.
        
        <b>Required Interface:</b>
        
        \code
        Iterator is a standard forward iterator.
        
        bool f = *first < NumericTraits<typename std::iterator_traits<Iterator>::value_type>::max();
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class Iterator>
Iterator argMin(Iterator first, Iterator last)
{
    if(first == last)
        return last;
    Iterator best = first;
    for(++first; first != last; ++first)
        if(*first < *best)
            best = first;
    return best;
}

    /** \brief Find the maximum element in a sequence.
    
        The function returns the iterator referring to the maximum element.
        This is identical to the function <tt>std::max_element()</tt>.
        
        <b>Required Interface:</b>
        
        \code
        Iterator is a standard forward iterator.
        
        bool f = NumericTraits<typename std::iterator_traits<Iterator>::value_type>::min() < *first;
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class Iterator>
Iterator argMax(Iterator first, Iterator last)
{
    if(first == last)
        return last;
    Iterator best = first;
    for(++first; first != last; ++first)
        if(*best < *first)
            best = first;
    return best;
}

    /** \brief Find the minimum element in a sequence conforming to a condition.
    
        The function returns the iterator referring to the minimum element,
        where only elements conforming to the condition (i.e. where 
        <tt>condition(*iterator)</tt> evaluates to <tt>true</tt>) are considered.
        If no element conforms to the condition, or the sequence is empty,
        the end iterator \a last is returned.
        
        <b>Required Interface:</b>
        
        \code
        Iterator is a standard forward iterator.
        
        bool c = condition(*first);
        
        bool f = *first < NumericTraits<typename std::iterator_traits<Iterator>::value_type>::max();
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class Iterator, class UnaryFunctor>
Iterator argMinIf(Iterator first, Iterator last, UnaryFunctor condition)
{
    for(; first != last; ++first)
        if(condition(*first))
            break;
    if(first == last)
        return last;
    Iterator best = first;
    for(++first; first != last; ++first)
        if(condition(*first) && *first < *best)
            best = first;
    return best;
}

    /** \brief Find the maximum element in a sequence conforming to a condition.
    
        The function returns the iterator referring to the maximum element,
        where only elements conforming to the condition (i.e. where 
        <tt>condition(*iterator)</tt> evaluates to <tt>true</tt>) are considered.
        If no element conforms to the condition, or the sequence is empty,
        the end iterator \a last is returned.
        
        <b>Required Interface:</b>
        
        \code
        Iterator is a standard forward iterator.
        
        bool c = condition(*first);
        
        bool f = NumericTraits<typename std::iterator_traits<Iterator>::value_type>::min() < *first;
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class Iterator, class UnaryFunctor>
Iterator argMaxIf(Iterator first, Iterator last, UnaryFunctor condition)
{
    for(; first != last; ++first)
        if(condition(*first))
            break;
    if(first == last)
        return last;
    Iterator best = first;
    for(++first; first != last; ++first)
        if(condition(*first) && *best < *first)
            best = first;
    return best;
}

    /** \brief Fill an array with a sequence of numbers.
    
        The sequence starts at \a start and is incremented with \a step. Default start
        and stepsize are 0 and 1 respectively.
        
        <b> Declaration:</b>

        \code
        namespace vigra {
            template <class Iterator, class Value>
            void linearSequence(Iterator first, Iterator last, 
                          Value const & start = 0, Value const & step = 1);
        }
        \endcode
        
        <b>Required Interface:</b>
        
        \code
        Iterator is a standard forward iterator.
        
        *first = start;
        start += step;
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class Iterator, class Value>
void linearSequence(Iterator first, Iterator last, Value start, Value step)
{
    for(; first != last; ++first, start += step)
        *first = start;
}

template <class Iterator, class Value>
void linearSequence(Iterator first, Iterator last, Value start)
{
    for(; first != last; ++first, ++start)
        *first = start;
}

template <class Iterator>
void linearSequence(Iterator first, Iterator last)
{
    typedef typename std::iterator_traits<Iterator>::value_type Value;
    
    linearSequence(first, last, NumericTraits<Value>::zero());
}

/** \brief Call an analyzing functor at every element of a sequence.

    This function can be used to collect statistics of the sequence
    <tt>[first, last)</tt> defined by these two input interators.
    The results must be stored in the functor, which serves as a return
    value.

    <b> Declarations:</b>

    \code
    namespace vigra {
        template <class InputIterator, class Functor>
        void
        inspectSequence(InputIterator first, InputIterator last, Functor & f);
    }
    \endcode

    <b> Usage:</b>

    <b>\#include</b> \<vigra/algorithm.hxx\><br>
    Namespace: vigra

    \code
    std::vector array(100);

    // init functor
    vigra::FindMinMax<int> minmax;

    vigra::inspectSequence(array.begin(), array.end(), minmax);

    cout << "Min: " << minmax.min << " Max: " << minmax.max;

    \endcode

*/
doxygen_overloaded_function(template <...> void inspectSequence)

namespace detail {

template <class InputIterator>
struct inspectSequence_binder
{
    InputIterator first;
    InputIterator last;
    inspectSequence_binder(InputIterator first_, InputIterator last_)
        : first(first_), last(last_) {}
    template <class Functor>
    void operator()(Functor & f)
    {
        for (InputIterator i = first; i != last; ++i)
            f(*i);
    }
};

} // namespace detail

template <class InputIterator, class Functor>
inline void
inspectSequence(InputIterator first, InputIterator last, Functor & f)
{
    detail::inspectSequence_binder<InputIterator> g(first, last);
    detail::extra_passes_select(g, f);
}
   
namespace detail {

template <class Iterator, class Compare>
struct IndexCompare
{
    Iterator i_;
    Compare c_;
    
    IndexCompare(Iterator i, Compare c)
    : i_(i),
      c_(c)
    {}

    template <class Index>
    bool operator()(Index const & l, Index const & r) const
    {
        return c_(i_[l], i_[r]);
    }
};

} // namespace detail

    /** \brief Return the index permutation that would sort the input array.
    
        To actually sort an array according to the ordering thus determined, use 
        \ref applyPermutation().
        
        <b> Declarations:</b>

        \code
        namespace vigra {
            // compare using std::less
            template <class Iterator, class IndexIterator>
            void indexSort(Iterator first, Iterator last, IndexIterator index_first);

            // compare using functor Compare
            template <class Iterator, class IndexIterator, class Compare>
            void indexSort(Iterator first, Iterator last, IndexIterator index_first, Compare compare);
        }
        \endcode
        
        <b>Required Interface:</b>
        
        \code
        Iterator and IndexIterators are random access iterators.
        
        bool res = compare(first[*index_first], first[*index_first]);
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class Iterator, class IndexIterator, class Compare>
void indexSort(Iterator first, Iterator last, IndexIterator index_first, Compare c)
{
    int size = last - first;
    linearSequence(index_first, index_first+size);
    std::sort(index_first, index_first+size, 
              detail::IndexCompare<Iterator, Compare>(first, c));
}

template <class Iterator, class IndexIterator>
void indexSort(Iterator first, Iterator last, IndexIterator index_first)
{
    typedef typename std::iterator_traits<Iterator>::value_type Value;
    indexSort(first, last, index_first, std::less<Value>());
}

    /** \brief Sort an array according to the given index permutation.
    
        The iterators \a in and \a out may not refer to the same array, as
        this would overwrite the input prematurely.
        
        <b> Declaration:</b>

        \code
        namespace vigra {
            template <class IndexIterator, class InIterator, class OutIterator>
            void applyPermutation(IndexIterator index_first, IndexIterator index_last, 
                                  InIterator in, OutIterator out);
        }
        \endcode
        
        <b>Required Interface:</b>
        
        \code
        OutIterator and IndexIterators are forward iterators.
        InIterator is a random access iterator.
        
        *out = in[*index_first];
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class IndexIterator, class InIterator, class OutIterator>
void applyPermutation(IndexIterator index_first, IndexIterator index_last, 
                      InIterator in, OutIterator out)
{
    for(; index_first != index_last; ++index_first, ++out)
        *out = in[*index_first];
}


    /** \brief Compute the inverse of a given permutation.
    
        This is just another name for \ref indexSort(), referring to
        another semantics.
        
        <b> Declaration:</b>

        \code
        namespace vigra {
            template <class InIterator, class OutIterator>
            void inversePermutation(InIterator first, InIterator last, 
                                    OutIterator out);
        }
        \endcode
        
        <b>Required Interface:</b>
        
        \code
        InIterator and OutIterator are random access iterators.
        
        *out = in[*index_first];
        \endcode

        <b>\#include</b> \<vigra/algorithm.hxx\><br>
        Namespace: vigra
    */
template <class InIterator, class OutIterator>
void inversePermutation(InIterator first, InIterator last, 
                        OutIterator out)
{
    indexSort(first, last, out);
}

namespace detail {

static bool isLittleEndian()
{
    static const UIntBiggest testint = 0x01;
    return ((UInt8 *)&testint)[0] == 0x01;
}

template <class INT>
struct ChecksumImpl
{
    static UInt32 table0[256];
    static UInt32 table1[256];
    static UInt32 table2[256];
    static UInt32 table3[256];
    
    template <class InIterator>
    static UInt32 exec(InIterator i, unsigned int size, UInt32 crc = 0xFFFFFFFF);
};

template <class INT>
UInt32 ChecksumImpl<INT>::table0[256] = {
    0x0U, 0x77073096U, 0xee0e612cU, 0x990951baU, 0x76dc419U, 0x706af48fU, 
    0xe963a535U, 0x9e6495a3U, 0xedb8832U, 0x79dcb8a4U, 0xe0d5e91eU, 0x97d2d988U, 
    0x9b64c2bU, 0x7eb17cbdU, 0xe7b82d07U, 0x90bf1d91U, 0x1db71064U, 0x6ab020f2U, 
    0xf3b97148U, 0x84be41deU, 0x1adad47dU, 0x6ddde4ebU, 0xf4d4b551U, 0x83d385c7U, 
    0x136c9856U, 0x646ba8c0U, 0xfd62f97aU, 0x8a65c9ecU, 0x14015c4fU, 0x63066cd9U, 
    0xfa0f3d63U, 0x8d080df5U, 0x3b6e20c8U, 0x4c69105eU, 0xd56041e4U, 0xa2677172U, 
    0x3c03e4d1U, 0x4b04d447U, 0xd20d85fdU, 0xa50ab56bU, 0x35b5a8faU, 0x42b2986cU, 
    0xdbbbc9d6U, 0xacbcf940U, 0x32d86ce3U, 0x45df5c75U, 0xdcd60dcfU, 0xabd13d59U, 
    0x26d930acU, 0x51de003aU, 0xc8d75180U, 0xbfd06116U, 0x21b4f4b5U, 0x56b3c423U, 
    0xcfba9599U, 0xb8bda50fU, 0x2802b89eU, 0x5f058808U, 0xc60cd9b2U, 0xb10be924U, 
    0x2f6f7c87U, 0x58684c11U, 0xc1611dabU, 0xb6662d3dU, 0x76dc4190U, 0x1db7106U, 
    0x98d220bcU, 0xefd5102aU, 0x71b18589U, 0x6b6b51fU, 0x9fbfe4a5U, 0xe8b8d433U, 
    0x7807c9a2U, 0xf00f934U, 0x9609a88eU, 0xe10e9818U, 0x7f6a0dbbU, 0x86d3d2dU, 
    0x91646c97U, 0xe6635c01U, 0x6b6b51f4U, 0x1c6c6162U, 0x856530d8U, 0xf262004eU, 
    0x6c0695edU, 0x1b01a57bU, 0x8208f4c1U, 0xf50fc457U, 0x65b0d9c6U, 0x12b7e950U, 
    0x8bbeb8eaU, 0xfcb9887cU, 0x62dd1ddfU, 0x15da2d49U, 0x8cd37cf3U, 0xfbd44c65U, 
    0x4db26158U, 0x3ab551ceU, 0xa3bc0074U, 0xd4bb30e2U, 0x4adfa541U, 0x3dd895d7U, 
    0xa4d1c46dU, 0xd3d6f4fbU, 0x4369e96aU, 0x346ed9fcU, 0xad678846U, 0xda60b8d0U, 
    0x44042d73U, 0x33031de5U, 0xaa0a4c5fU, 0xdd0d7cc9U, 0x5005713cU, 0x270241aaU, 
    0xbe0b1010U, 0xc90c2086U, 0x5768b525U, 0x206f85b3U, 0xb966d409U, 0xce61e49fU, 
    0x5edef90eU, 0x29d9c998U, 0xb0d09822U, 0xc7d7a8b4U, 0x59b33d17U, 0x2eb40d81U, 
    0xb7bd5c3bU, 0xc0ba6cadU, 0xedb88320U, 0x9abfb3b6U, 0x3b6e20cU, 0x74b1d29aU, 
    0xead54739U, 0x9dd277afU, 0x4db2615U, 0x73dc1683U, 0xe3630b12U, 0x94643b84U, 
    0xd6d6a3eU, 0x7a6a5aa8U, 0xe40ecf0bU, 0x9309ff9dU, 0xa00ae27U, 0x7d079eb1U, 
    0xf00f9344U, 0x8708a3d2U, 0x1e01f268U, 0x6906c2feU, 0xf762575dU, 0x806567cbU, 
    0x196c3671U, 0x6e6b06e7U, 0xfed41b76U, 0x89d32be0U, 0x10da7a5aU, 0x67dd4accU, 
    0xf9b9df6fU, 0x8ebeeff9U, 0x17b7be43U, 0x60b08ed5U, 0xd6d6a3e8U, 0xa1d1937eU, 
    0x38d8c2c4U, 0x4fdff252U, 0xd1bb67f1U, 0xa6bc5767U, 0x3fb506ddU, 0x48b2364bU, 
    0xd80d2bdaU, 0xaf0a1b4cU, 0x36034af6U, 0x41047a60U, 0xdf60efc3U, 0xa867df55U, 
    0x316e8eefU, 0x4669be79U, 0xcb61b38cU, 0xbc66831aU, 0x256fd2a0U, 0x5268e236U, 
    0xcc0c7795U, 0xbb0b4703U, 0x220216b9U, 0x5505262fU, 0xc5ba3bbeU, 0xb2bd0b28U, 
    0x2bb45a92U, 0x5cb36a04U, 0xc2d7ffa7U, 0xb5d0cf31U, 0x2cd99e8bU, 0x5bdeae1dU, 
    0x9b64c2b0U, 0xec63f226U, 0x756aa39cU, 0x26d930aU, 0x9c0906a9U, 0xeb0e363fU, 
    0x72076785U, 0x5005713U, 0x95bf4a82U, 0xe2b87a14U, 0x7bb12baeU, 0xcb61b38U, 
    0x92d28e9bU, 0xe5d5be0dU, 0x7cdcefb7U, 0xbdbdf21U, 0x86d3d2d4U, 0xf1d4e242U, 
    0x68ddb3f8U, 0x1fda836eU, 0x81be16cdU, 0xf6b9265bU, 0x6fb077e1U, 0x18b74777U, 
    0x88085ae6U, 0xff0f6a70U, 0x66063bcaU, 0x11010b5cU, 0x8f659effU, 0xf862ae69U, 
    0x616bffd3U, 0x166ccf45U, 0xa00ae278U, 0xd70dd2eeU, 0x4e048354U, 0x3903b3c2U, 
    0xa7672661U, 0xd06016f7U, 0x4969474dU, 0x3e6e77dbU, 0xaed16a4aU, 0xd9d65adcU, 
    0x40df0b66U, 0x37d83bf0U, 0xa9bcae53U, 0xdebb9ec5U, 0x47b2cf7fU, 0x30b5ffe9U, 
    0xbdbdf21cU, 0xcabac28aU, 0x53b39330U, 0x24b4a3a6U, 0xbad03605U, 0xcdd70693U, 
    0x54de5729U, 0x23d967bfU, 0xb3667a2eU, 0xc4614ab8U, 0x5d681b02U, 0x2a6f2b94U, 
    0xb40bbe37U, 0xc30c8ea1U, 0x5a05df1bU, 0x2d02ef8dU }; 

template <class INT>
UInt32 ChecksumImpl<INT>::table1[256] = {
    0x00000000U, 0x191b3141U, 0x32366282U, 0x2b2d53c3U, 0x646cc504U,
    0x7d77f445U, 0x565aa786U, 0x4f4196c7U, 0xc8d98a08U, 0xd1c2bb49U,
    0xfaefe88aU, 0xe3f4d9cbU, 0xacb54f0cU, 0xb5ae7e4dU, 0x9e832d8eU,
    0x87981ccfU, 0x4ac21251U, 0x53d92310U, 0x78f470d3U, 0x61ef4192U,
    0x2eaed755U, 0x37b5e614U, 0x1c98b5d7U, 0x05838496U, 0x821b9859U,
    0x9b00a918U, 0xb02dfadbU, 0xa936cb9aU, 0xe6775d5dU, 0xff6c6c1cU,
    0xd4413fdfU, 0xcd5a0e9eU, 0x958424a2U, 0x8c9f15e3U, 0xa7b24620U,
    0xbea97761U, 0xf1e8e1a6U, 0xe8f3d0e7U, 0xc3de8324U, 0xdac5b265U,
    0x5d5daeaaU, 0x44469febU, 0x6f6bcc28U, 0x7670fd69U, 0x39316baeU,
    0x202a5aefU, 0x0b07092cU, 0x121c386dU, 0xdf4636f3U, 0xc65d07b2U,
    0xed705471U, 0xf46b6530U, 0xbb2af3f7U, 0xa231c2b6U, 0x891c9175U,
    0x9007a034U, 0x179fbcfbU, 0x0e848dbaU, 0x25a9de79U, 0x3cb2ef38U,
    0x73f379ffU, 0x6ae848beU, 0x41c51b7dU, 0x58de2a3cU, 0xf0794f05U,
    0xe9627e44U, 0xc24f2d87U, 0xdb541cc6U, 0x94158a01U, 0x8d0ebb40U,
    0xa623e883U, 0xbf38d9c2U, 0x38a0c50dU, 0x21bbf44cU, 0x0a96a78fU,
    0x138d96ceU, 0x5ccc0009U, 0x45d73148U, 0x6efa628bU, 0x77e153caU,
    0xbabb5d54U, 0xa3a06c15U, 0x888d3fd6U, 0x91960e97U, 0xded79850U,
    0xc7cca911U, 0xece1fad2U, 0xf5facb93U, 0x7262d75cU, 0x6b79e61dU,
    0x4054b5deU, 0x594f849fU, 0x160e1258U, 0x0f152319U, 0x243870daU,
    0x3d23419bU, 0x65fd6ba7U, 0x7ce65ae6U, 0x57cb0925U, 0x4ed03864U,
    0x0191aea3U, 0x188a9fe2U, 0x33a7cc21U, 0x2abcfd60U, 0xad24e1afU,
    0xb43fd0eeU, 0x9f12832dU, 0x8609b26cU, 0xc94824abU, 0xd05315eaU,
    0xfb7e4629U, 0xe2657768U, 0x2f3f79f6U, 0x362448b7U, 0x1d091b74U,
    0x04122a35U, 0x4b53bcf2U, 0x52488db3U, 0x7965de70U, 0x607eef31U,
    0xe7e6f3feU, 0xfefdc2bfU, 0xd5d0917cU, 0xcccba03dU, 0x838a36faU,
    0x9a9107bbU, 0xb1bc5478U, 0xa8a76539U, 0x3b83984bU, 0x2298a90aU,
    0x09b5fac9U, 0x10aecb88U, 0x5fef5d4fU, 0x46f46c0eU, 0x6dd93fcdU,
    0x74c20e8cU, 0xf35a1243U, 0xea412302U, 0xc16c70c1U, 0xd8774180U,
    0x9736d747U, 0x8e2de606U, 0xa500b5c5U, 0xbc1b8484U, 0x71418a1aU,
    0x685abb5bU, 0x4377e898U, 0x5a6cd9d9U, 0x152d4f1eU, 0x0c367e5fU,
    0x271b2d9cU, 0x3e001cddU, 0xb9980012U, 0xa0833153U, 0x8bae6290U,
    0x92b553d1U, 0xddf4c516U, 0xc4eff457U, 0xefc2a794U, 0xf6d996d5U,
    0xae07bce9U, 0xb71c8da8U, 0x9c31de6bU, 0x852aef2aU, 0xca6b79edU,
    0xd37048acU, 0xf85d1b6fU, 0xe1462a2eU, 0x66de36e1U, 0x7fc507a0U,
    0x54e85463U, 0x4df36522U, 0x02b2f3e5U, 0x1ba9c2a4U, 0x30849167U,
    0x299fa026U, 0xe4c5aeb8U, 0xfdde9ff9U, 0xd6f3cc3aU, 0xcfe8fd7bU,
    0x80a96bbcU, 0x99b25afdU, 0xb29f093eU, 0xab84387fU, 0x2c1c24b0U,
    0x350715f1U, 0x1e2a4632U, 0x07317773U, 0x4870e1b4U, 0x516bd0f5U,
    0x7a468336U, 0x635db277U, 0xcbfad74eU, 0xd2e1e60fU, 0xf9ccb5ccU,
    0xe0d7848dU, 0xaf96124aU, 0xb68d230bU, 0x9da070c8U, 0x84bb4189U,
    0x03235d46U, 0x1a386c07U, 0x31153fc4U, 0x280e0e85U, 0x674f9842U,
    0x7e54a903U, 0x5579fac0U, 0x4c62cb81U, 0x8138c51fU, 0x9823f45eU,
    0xb30ea79dU, 0xaa1596dcU, 0xe554001bU, 0xfc4f315aU, 0xd7626299U,
    0xce7953d8U, 0x49e14f17U, 0x50fa7e56U, 0x7bd72d95U, 0x62cc1cd4U,
    0x2d8d8a13U, 0x3496bb52U, 0x1fbbe891U, 0x06a0d9d0U, 0x5e7ef3ecU,
    0x4765c2adU, 0x6c48916eU, 0x7553a02fU, 0x3a1236e8U, 0x230907a9U,
    0x0824546aU, 0x113f652bU, 0x96a779e4U, 0x8fbc48a5U, 0xa4911b66U,
    0xbd8a2a27U, 0xf2cbbce0U, 0xebd08da1U, 0xc0fdde62U, 0xd9e6ef23U,
    0x14bce1bdU, 0x0da7d0fcU, 0x268a833fU, 0x3f91b27eU, 0x70d024b9U,
    0x69cb15f8U, 0x42e6463bU, 0x5bfd777aU, 0xdc656bb5U, 0xc57e5af4U,
    0xee530937U, 0xf7483876U, 0xb809aeb1U, 0xa1129ff0U, 0x8a3fcc33U,
    0x9324fd72U };
    
template <class INT>
UInt32 ChecksumImpl<INT>::table2[256] = {
    0x00000000U, 0x01c26a37U, 0x0384d46eU, 0x0246be59U, 0x0709a8dcU,
    0x06cbc2ebU, 0x048d7cb2U, 0x054f1685U, 0x0e1351b8U, 0x0fd13b8fU,
    0x0d9785d6U, 0x0c55efe1U, 0x091af964U, 0x08d89353U, 0x0a9e2d0aU,
    0x0b5c473dU, 0x1c26a370U, 0x1de4c947U, 0x1fa2771eU, 0x1e601d29U,
    0x1b2f0bacU, 0x1aed619bU, 0x18abdfc2U, 0x1969b5f5U, 0x1235f2c8U,
    0x13f798ffU, 0x11b126a6U, 0x10734c91U, 0x153c5a14U, 0x14fe3023U,
    0x16b88e7aU, 0x177ae44dU, 0x384d46e0U, 0x398f2cd7U, 0x3bc9928eU,
    0x3a0bf8b9U, 0x3f44ee3cU, 0x3e86840bU, 0x3cc03a52U, 0x3d025065U,
    0x365e1758U, 0x379c7d6fU, 0x35dac336U, 0x3418a901U, 0x3157bf84U,
    0x3095d5b3U, 0x32d36beaU, 0x331101ddU, 0x246be590U, 0x25a98fa7U,
    0x27ef31feU, 0x262d5bc9U, 0x23624d4cU, 0x22a0277bU, 0x20e69922U,
    0x2124f315U, 0x2a78b428U, 0x2bbade1fU, 0x29fc6046U, 0x283e0a71U,
    0x2d711cf4U, 0x2cb376c3U, 0x2ef5c89aU, 0x2f37a2adU, 0x709a8dc0U,
    0x7158e7f7U, 0x731e59aeU, 0x72dc3399U, 0x7793251cU, 0x76514f2bU,
    0x7417f172U, 0x75d59b45U, 0x7e89dc78U, 0x7f4bb64fU, 0x7d0d0816U,
    0x7ccf6221U, 0x798074a4U, 0x78421e93U, 0x7a04a0caU, 0x7bc6cafdU,
    0x6cbc2eb0U, 0x6d7e4487U, 0x6f38fadeU, 0x6efa90e9U, 0x6bb5866cU,
    0x6a77ec5bU, 0x68315202U, 0x69f33835U, 0x62af7f08U, 0x636d153fU,
    0x612bab66U, 0x60e9c151U, 0x65a6d7d4U, 0x6464bde3U, 0x662203baU,
    0x67e0698dU, 0x48d7cb20U, 0x4915a117U, 0x4b531f4eU, 0x4a917579U,
    0x4fde63fcU, 0x4e1c09cbU, 0x4c5ab792U, 0x4d98dda5U, 0x46c49a98U,
    0x4706f0afU, 0x45404ef6U, 0x448224c1U, 0x41cd3244U, 0x400f5873U,
    0x4249e62aU, 0x438b8c1dU, 0x54f16850U, 0x55330267U, 0x5775bc3eU,
    0x56b7d609U, 0x53f8c08cU, 0x523aaabbU, 0x507c14e2U, 0x51be7ed5U,
    0x5ae239e8U, 0x5b2053dfU, 0x5966ed86U, 0x58a487b1U, 0x5deb9134U,
    0x5c29fb03U, 0x5e6f455aU, 0x5fad2f6dU, 0xe1351b80U, 0xe0f771b7U,
    0xe2b1cfeeU, 0xe373a5d9U, 0xe63cb35cU, 0xe7fed96bU, 0xe5b86732U,
    0xe47a0d05U, 0xef264a38U, 0xeee4200fU, 0xeca29e56U, 0xed60f461U,
    0xe82fe2e4U, 0xe9ed88d3U, 0xebab368aU, 0xea695cbdU, 0xfd13b8f0U,
    0xfcd1d2c7U, 0xfe976c9eU, 0xff5506a9U, 0xfa1a102cU, 0xfbd87a1bU,
    0xf99ec442U, 0xf85cae75U, 0xf300e948U, 0xf2c2837fU, 0xf0843d26U,
    0xf1465711U, 0xf4094194U, 0xf5cb2ba3U, 0xf78d95faU, 0xf64fffcdU,
    0xd9785d60U, 0xd8ba3757U, 0xdafc890eU, 0xdb3ee339U, 0xde71f5bcU,
    0xdfb39f8bU, 0xddf521d2U, 0xdc374be5U, 0xd76b0cd8U, 0xd6a966efU,
    0xd4efd8b6U, 0xd52db281U, 0xd062a404U, 0xd1a0ce33U, 0xd3e6706aU,
    0xd2241a5dU, 0xc55efe10U, 0xc49c9427U, 0xc6da2a7eU, 0xc7184049U,
    0xc25756ccU, 0xc3953cfbU, 0xc1d382a2U, 0xc011e895U, 0xcb4dafa8U,
    0xca8fc59fU, 0xc8c97bc6U, 0xc90b11f1U, 0xcc440774U, 0xcd866d43U,
    0xcfc0d31aU, 0xce02b92dU, 0x91af9640U, 0x906dfc77U, 0x922b422eU,
    0x93e92819U, 0x96a63e9cU, 0x976454abU, 0x9522eaf2U, 0x94e080c5U,
    0x9fbcc7f8U, 0x9e7eadcfU, 0x9c381396U, 0x9dfa79a1U, 0x98b56f24U,
    0x99770513U, 0x9b31bb4aU, 0x9af3d17dU, 0x8d893530U, 0x8c4b5f07U,
    0x8e0de15eU, 0x8fcf8b69U, 0x8a809decU, 0x8b42f7dbU, 0x89044982U,
    0x88c623b5U, 0x839a6488U, 0x82580ebfU, 0x801eb0e6U, 0x81dcdad1U,
    0x8493cc54U, 0x8551a663U, 0x8717183aU, 0x86d5720dU, 0xa9e2d0a0U,
    0xa820ba97U, 0xaa6604ceU, 0xaba46ef9U, 0xaeeb787cU, 0xaf29124bU,
    0xad6fac12U, 0xacadc625U, 0xa7f18118U, 0xa633eb2fU, 0xa4755576U,
    0xa5b73f41U, 0xa0f829c4U, 0xa13a43f3U, 0xa37cfdaaU, 0xa2be979dU,
    0xb5c473d0U, 0xb40619e7U, 0xb640a7beU, 0xb782cd89U, 0xb2cddb0cU,
    0xb30fb13bU, 0xb1490f62U, 0xb08b6555U, 0xbbd72268U, 0xba15485fU,
    0xb853f606U, 0xb9919c31U, 0xbcde8ab4U, 0xbd1ce083U, 0xbf5a5edaU,
    0xbe9834edU };
    
template <class INT>
UInt32 ChecksumImpl<INT>::table3[256] = {
    0x00000000U, 0xb8bc6765U, 0xaa09c88bU, 0x12b5afeeU, 0x8f629757U,
    0x37def032U, 0x256b5fdcU, 0x9dd738b9U, 0xc5b428efU, 0x7d084f8aU,
    0x6fbde064U, 0xd7018701U, 0x4ad6bfb8U, 0xf26ad8ddU, 0xe0df7733U,
    0x58631056U, 0x5019579fU, 0xe8a530faU, 0xfa109f14U, 0x42acf871U,
    0xdf7bc0c8U, 0x67c7a7adU, 0x75720843U, 0xcdce6f26U, 0x95ad7f70U,
    0x2d111815U, 0x3fa4b7fbU, 0x8718d09eU, 0x1acfe827U, 0xa2738f42U,
    0xb0c620acU, 0x087a47c9U, 0xa032af3eU, 0x188ec85bU, 0x0a3b67b5U,
    0xb28700d0U, 0x2f503869U, 0x97ec5f0cU, 0x8559f0e2U, 0x3de59787U,
    0x658687d1U, 0xdd3ae0b4U, 0xcf8f4f5aU, 0x7733283fU, 0xeae41086U,
    0x525877e3U, 0x40edd80dU, 0xf851bf68U, 0xf02bf8a1U, 0x48979fc4U,
    0x5a22302aU, 0xe29e574fU, 0x7f496ff6U, 0xc7f50893U, 0xd540a77dU,
    0x6dfcc018U, 0x359fd04eU, 0x8d23b72bU, 0x9f9618c5U, 0x272a7fa0U,
    0xbafd4719U, 0x0241207cU, 0x10f48f92U, 0xa848e8f7U, 0x9b14583dU,
    0x23a83f58U, 0x311d90b6U, 0x89a1f7d3U, 0x1476cf6aU, 0xaccaa80fU,
    0xbe7f07e1U, 0x06c36084U, 0x5ea070d2U, 0xe61c17b7U, 0xf4a9b859U,
    0x4c15df3cU, 0xd1c2e785U, 0x697e80e0U, 0x7bcb2f0eU, 0xc377486bU,
    0xcb0d0fa2U, 0x73b168c7U, 0x6104c729U, 0xd9b8a04cU, 0x446f98f5U,
    0xfcd3ff90U, 0xee66507eU, 0x56da371bU, 0x0eb9274dU, 0xb6054028U,
    0xa4b0efc6U, 0x1c0c88a3U, 0x81dbb01aU, 0x3967d77fU, 0x2bd27891U,
    0x936e1ff4U, 0x3b26f703U, 0x839a9066U, 0x912f3f88U, 0x299358edU,
    0xb4446054U, 0x0cf80731U, 0x1e4da8dfU, 0xa6f1cfbaU, 0xfe92dfecU,
    0x462eb889U, 0x549b1767U, 0xec277002U, 0x71f048bbU, 0xc94c2fdeU,
    0xdbf98030U, 0x6345e755U, 0x6b3fa09cU, 0xd383c7f9U, 0xc1366817U,
    0x798a0f72U, 0xe45d37cbU, 0x5ce150aeU, 0x4e54ff40U, 0xf6e89825U,
    0xae8b8873U, 0x1637ef16U, 0x048240f8U, 0xbc3e279dU, 0x21e91f24U,
    0x99557841U, 0x8be0d7afU, 0x335cb0caU, 0xed59b63bU, 0x55e5d15eU,
    0x47507eb0U, 0xffec19d5U, 0x623b216cU, 0xda874609U, 0xc832e9e7U,
    0x708e8e82U, 0x28ed9ed4U, 0x9051f9b1U, 0x82e4565fU, 0x3a58313aU,
    0xa78f0983U, 0x1f336ee6U, 0x0d86c108U, 0xb53aa66dU, 0xbd40e1a4U,
    0x05fc86c1U, 0x1749292fU, 0xaff54e4aU, 0x322276f3U, 0x8a9e1196U,
    0x982bbe78U, 0x2097d91dU, 0x78f4c94bU, 0xc048ae2eU, 0xd2fd01c0U,
    0x6a4166a5U, 0xf7965e1cU, 0x4f2a3979U, 0x5d9f9697U, 0xe523f1f2U,
    0x4d6b1905U, 0xf5d77e60U, 0xe762d18eU, 0x5fdeb6ebU, 0xc2098e52U,
    0x7ab5e937U, 0x680046d9U, 0xd0bc21bcU, 0x88df31eaU, 0x3063568fU,
    0x22d6f961U, 0x9a6a9e04U, 0x07bda6bdU, 0xbf01c1d8U, 0xadb46e36U,
    0x15080953U, 0x1d724e9aU, 0xa5ce29ffU, 0xb77b8611U, 0x0fc7e174U,
    0x9210d9cdU, 0x2aacbea8U, 0x38191146U, 0x80a57623U, 0xd8c66675U,
    0x607a0110U, 0x72cfaefeU, 0xca73c99bU, 0x57a4f122U, 0xef189647U,
    0xfdad39a9U, 0x45115eccU, 0x764dee06U, 0xcef18963U, 0xdc44268dU,
    0x64f841e8U, 0xf92f7951U, 0x41931e34U, 0x5326b1daU, 0xeb9ad6bfU,
    0xb3f9c6e9U, 0x0b45a18cU, 0x19f00e62U, 0xa14c6907U, 0x3c9b51beU,
    0x842736dbU, 0x96929935U, 0x2e2efe50U, 0x2654b999U, 0x9ee8defcU,
    0x8c5d7112U, 0x34e11677U, 0xa9362eceU, 0x118a49abU, 0x033fe645U,
    0xbb838120U, 0xe3e09176U, 0x5b5cf613U, 0x49e959fdU, 0xf1553e98U,
    0x6c820621U, 0xd43e6144U, 0xc68bceaaU, 0x7e37a9cfU, 0xd67f4138U,
    0x6ec3265dU, 0x7c7689b3U, 0xc4caeed6U, 0x591dd66fU, 0xe1a1b10aU,
    0xf3141ee4U, 0x4ba87981U, 0x13cb69d7U, 0xab770eb2U, 0xb9c2a15cU,
    0x017ec639U, 0x9ca9fe80U, 0x241599e5U, 0x36a0360bU, 0x8e1c516eU,
    0x866616a7U, 0x3eda71c2U, 0x2c6fde2cU, 0x94d3b949U, 0x090481f0U,
    0xb1b8e695U, 0xa30d497bU, 0x1bb12e1eU, 0x43d23e48U, 0xfb6e592dU,
    0xe9dbf6c3U, 0x516791a6U, 0xccb0a91fU, 0x740cce7aU, 0x66b96194U,
    0xde0506f1U };


template <class INT>
template <class InIterator>
UInt32 ChecksumImpl<INT>::exec(InIterator i, unsigned int size, UInt32 crc)
{
    InIterator end = i + size;
    
    if(isLittleEndian() && size > 3)
    {
        // take care of alignment
        for(; (std::size_t)i % 4 != 0; ++i)
        {
            crc = (crc >> 8) ^ table0[(crc ^ *i) & 0xFF];
        }
        for(; i < end-3; i+=4)
        {
            crc ^= *((UInt32 *)i);
            crc = table3[crc & 0xFF] ^
                  table2[(crc >> 8) & 0xFF] ^
                  table1[(crc >> 16) & 0xFF] ^
                  table0[crc >> 24];
        }
    }
    for(; i < end; ++i)
    {
        crc = (crc >> 8) ^ table0[(crc ^ *i) & 0xFF];
    }
    return ~crc;
}

} // namespace detail

    /** \brief Compute the CRC-32 checksum of a byte array.
    
        Implementation note: This function is slower on big-endian machines
        because the "4 bytes at a time" optimization is only implemented for 
        little-endian.
    */
inline UInt32 checksum(const char * data, unsigned int size)
{
    return detail::ChecksumImpl<UInt32>::exec(data, size);
}

    /** Concatenate a byte array to an existing CRC-32 checksum.
    */
inline UInt32 concatenateChecksum(UInt32 checksum, const char * data, unsigned int size)
{
    
    return detail::ChecksumImpl<UInt32>::exec(data, size, ~checksum);
}

template <class T>
void updateMin(T & x, const T & y)
{
    using std::min;
    x = min(x, y);
}

template <class T>
void updateMax(T & x, const T & y)
{
    using std::max;
    x = max(x, y);
}


//@}

} // namespace vigra

#endif /* VIGRA_ALGORITHM_HXX */
