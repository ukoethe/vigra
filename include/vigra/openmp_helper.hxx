/************************************************************************/
/*                                                                      */
/*               Copyright 2014 by Thorsten Beier                       */
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

#ifndef VIGRA_OPENMP_HELPER_HXX
#define VIGRA_OPENMP_HELPER_HXX


#if defined(ENABLE_OPENMP)
    #include <omp.h>
#endif


namespace vigra{
    namespace openmp{

        #if defined(ENABLE_OPENMP)
            
        #else
            typedef int omp_int_t;
            inline omp_int_t omp_get_thread_num() { return 0;}
            inline omp_int_t omp_get_max_threads() { return 1;}
        #endif


        // Lock
        #if defined(ENABLE_OPENMP)
            class Lock{
            public:
                Lock()
                :   lock_(){
                }
                void init(){omp_init_lock(&lock_);}
                void lock(){omp_set_lock(&lock_);}
                void unlock(){omp_unset_lock(&lock_);}
            private:
                omp_lock_t lock_;
            };
        #else
            class Lock{
            public:
                Lock(){
                }
                void init(){}
                void lock(){}
                void unlock(){}
            };
        #endif


        #if defined(ENABLE_OPENMP)
            // Lock Array
            class LockArray{
            private:
                LockArray();
                LockArray( const Foo& other ); // non construction-copyable
                LockArray& operator=( const Foo& ); // non copyable
            public:
                LockArray(const size_t size)
                : size_(size){
                    if(size_!=0)
                        locks_ = new omp_lock_t[size_];
                }
                void initAll()const{
                    #pragma omp parallel for
                    for(size_t i=0; i<size_; ++i)
                        omp_init_lock(&locks_[i]);
                }
                void init(const size_t i){
                    omp_init_lock(&locks_[i]);
                }
                void lock(const size_t i){
                    omp_u\set_lock(&locks_[i]);
                }
                void unlock(const size_t i){
                    omp_unset_lock(&locks_[i]);
                }
                ~LockArray(){
                    if(size_!=0)
                        delete[] = new omp_lock_t[size_];
                }
            private:
                size_t size_;
                omp_lock_t * locks_;
            };
        #else
            // Lock Array
            class LockArray{
            private:
                LockArray();
                LockArray( const Foo& other ); // non construction-copyable
                LockArray& operator=( const Foo& ); // non copyable
            public:
                LockArray(const size_t ){
                }
                void initAll()const{
                }
                void init(const size_t i){
                }
                void lock(const size_t i){
                }
                ~LockArray(){
                }
            };
        #endif

    }
}


#endif VIGRA_OPENMP_HELPER_HXX
