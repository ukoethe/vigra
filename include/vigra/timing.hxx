/************************************************************************/
/*                                                                      */
/*               Copyright 2008-2009 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
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


#ifndef VIGRA_TIMING_HXX
#define VIGRA_TIMING_HXX
#ifndef VIGRA_NO_TIMING

#include <iostream>
#include <sstream>
#ifdef MULTI_TICTOC
    #include <vector>
#endif
// usage:
// void time_it()
// {
//     USETICTOC;
//     TIC;
//      ...
//     std::cerr << TOC << " for time_it\n";
// }

#ifdef WIN32

    #include "windows.h"

    namespace {

    inline double queryTimerUnit()
    {
        LARGE_INTEGER frequency;
        QueryPerformanceFrequency(&frequency);
        return 1000.0 / frequency.QuadPart;
    }

    inline double tic_toc_diff_num(LARGE_INTEGER const & tic)
    {
        LARGE_INTEGER toc;
        QueryPerformanceCounter(&toc);
        static double unit = queryTimerUnit();
        return ((toc.QuadPart - tic.QuadPart) * unit);
    }

    inline std::string tic_toc_diff_string(LARGE_INTEGER const & tic)
    {
        double diff = tic_toc_diff_num(tic); 
        std::stringstream s;
        s << diff << " msec";
        return s.str();
    }

    inline void tic_toc_diff(LARGE_INTEGER const & tic)
    {
        std::cerr << tic_toc_diff_string(tic) <<std::endl;
    }

    } // unnamed namespace
    
#ifndef MULTI_TICTOC
    #define USETICTOC LARGE_INTEGER tic_timer;
    #define TIC QueryPerformanceCounter(&tic_timer);
    #define TOC  tic_toc_diff       (tic_timer);
    #define TOCN tic_toc_diff_num   (tic_timer);
    #define TOCS tic_toc_diff_string(tic_timer);
#else
    #define USETICTOC std::vector<LARGE_INTEGER> tic_timer;
    #define TIC tic_timer.push_back(LARGE_INTEGER());\
                QueryPerformanceCounter(&(tic_timer.back()));
    #define TOC  tic_toc_diff       (tic_timer.back());\
                 tic_timer.pop_back();
    #define TOCN tic_toc_diff_num   (tic_timer.back());\
                 tic_timer.pop_back();
    #define TOCS tic_toc_diff_string(tic_timer.back());\
                 tic_timer.pop_back();
#endif

#else

    #if defined(VIGRA_HIRES_TIMING) && !defined(__CYGWIN__)
        // requires linking against librt
    
        #include <time.h>

        namespace {

        inline double tic_toc_diff_num(timespec const & tic)
        {
            timespec toc;
            clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &toc);
            return ((toc.tv_sec*1000.0 + toc.tv_nsec/1000000.0) -
                  (tic.tv_sec*1000.0 + tic.tv_nsec/1000000.0));
        }

        inline std::string tic_toc_diff_string(timeval const & tic)
        {
            double diff = tic_toc_diff_num(tic); 
            std::stringstream s;
            s << diff << " msec";
            return s.str();
        }

        inline void tic_toc_diff(timeval const & tic)
        {
            std::cerr << tic_toc_diff_string(tic) << std::endl;
        }
        } // unnamed namespace

#ifndef MULTI_TICTOC
        #define USETICTOC timespec tic_timer;
        #define TIC clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tic_timer);
        #define TOC  tic_toc_diff       (tic_timer);
        #define TOCN tic_toc_diff_num   (tic_timer);
        #define TOCS tic_toc_diff_string(tic_timer);
#else

        #define USETICTOC std::vector<timespec> tic_timer;
        #define TIC tic_timer.push_back(timespec());\
                    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(tic_timer.back()));
        #define TOC  tic_toc_diff       (tic_timer.back());\
                     tic_timer.pop_back();
        #define TOCN tic_toc_diff_num   (tic_timer.back());\
                     tic_timer.pop_back();
        #define TOCS tic_toc_diff_string(tic_timer.back());\
                     tic_timer.pop_back();
#endif
    #else
    
        #include <sys/time.h>

        namespace {

        inline double tic_toc_diff_num(timeval const & tic)
        {
            timeval toc;
            gettimeofday(&toc, NULL);
            return  ((toc.tv_sec*1000.0 + toc.tv_usec/1000.0) -
                        (tic.tv_sec*1000.0 + tic.tv_usec/1000.0));
        }
        inline std::string tic_toc_diff_string(timeval const & tic)
        {
            double diff = tic_toc_diff_num(tic); 
            std::stringstream s;
            s << diff << " msec";
            return s.str();
        }
        inline void tic_toc_diff(timeval const & tic)
        {
            std::cerr << tic_toc_diff_string(tic)<< std::endl;
        }

        } // unnamed namespace

#ifndef MULTI_TICTOC
        #define USETICTOC timeval tic_timer;
        #define TIC  gettimeofday       (&tic_timer, NULL);
        #define TOC  tic_toc_diff       (tic_timer);
        #define TOCN tic_toc_diff_num   (tic_timer);
        #define TOCS tic_toc_diff_string(tic_timer);
#else

        #define USETICTOC std::vector<timeval> tic_timer;
        #define TIC tic_timer.push_back(timeval());\
                    gettimeofday(&(tic_timer.back()), NULL);
        #define TOC  tic_toc_diff       (tic_timer.back());\
                     tic_timer.pop_back();
        #define TOCN tic_toc_diff_num   (tic_timer.back());\
                     tic_timer.pop_back();
        #define TOCS tic_toc_diff_string(tic_timer.back());\
                     tic_timer.pop_back();
#endif

    #endif // VIGRA_HIRES_TIMING

#endif // WIN32




#else // NDEBUG

#define USETICTOC 
#define TIC
#define TOC
#define TOCN
#define TICS
#endif // NDEBUG

#endif // VIGRA_TIMING_HXX
