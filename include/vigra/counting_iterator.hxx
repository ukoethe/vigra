/************************************************************************/
/*                                                                      */
/*                 Copyright 2015 by Thorsten Beier                     */
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

#ifndef VIGRA_COUNTING_ITERATOR_HXX
#define VIGRA_COUNTING_ITERATOR_HXX


#include <iterator>

namespace vigra{


    template<class T_INTEGER>
    class CountingIterator: 
    public std::iterator<std::random_access_iterator_tag,T_INTEGER, std::ptrdiff_t, T_INTEGER *, T_INTEGER>
    {
    public:
        CountingIterator(): count_(0) {}
        CountingIterator(const T_INTEGER count): count_(count) {}
        CountingIterator(const CountingIterator& other): count_(other.count_) {}
        const CountingIterator& operator=(const CountingIterator& other) {count_ = other.count_; return other;}

        CountingIterator& operator++()    {count_++; return *this;} // prefix++
        CountingIterator  operator++(int)const{CountingIterator tmp(*this); ++(*this); return tmp;} // postfix++
        CountingIterator& operator--()    {count_--; return *this;} // prefix--
        CountingIterator  operator--(int)const{CountingIterator tmp(*this); --(*this); return tmp;} // postfix--

        void     operator+=(const std::size_t& n)  {count_ += n;}
        void     operator+=(const CountingIterator& other) {count_ += other.count_;}
        CountingIterator operator+ (const std::size_t& n) const{CountingIterator tmp(*this); tmp += n; return tmp;}
        CountingIterator operator+ (const CountingIterator& other) {CountingIterator tmp(*this); tmp += other; return tmp;}

        void        operator-=(const std::size_t& n)  {count_ -= n;}
        void        operator-=(const CountingIterator& other) {count_ -= other.count_;}
        CountingIterator    operator- (const std::size_t& n) const{CountingIterator tmp(*this); tmp -= n; return tmp;}
        std::size_t operator- (const CountingIterator& other) {return count_ - other.count_;}

        bool operator< (const CountingIterator& other)const{return (count_-other.count_)< 0;}
        bool operator<=(const CountingIterator& other)const{return (count_-other.count_)<=0;}
        bool operator> (const CountingIterator& other)const{return (count_-other.count_)> 0;}
        bool operator>=(const CountingIterator& other)const{return (count_-other.count_)>=0;}
        bool operator==(const CountingIterator& other)const{return  count_ == other.count_; }
        bool operator!=(const CountingIterator& other)const{return  count_ != other.count_; }

        T_INTEGER   operator[](const int& n) const {
            return count_ + n;
        }
        T_INTEGER   operator[](const int& n) {
            return count_ + n;
        }
        T_INTEGER   operator*()const{
            return  count_;
        }
        T_INTEGER   operator*() {
            return  count_;
        }
        T_INTEGER * operator->(){
            return &count_;
        }
    private:
        T_INTEGER  count_;
    };

}

#endif
