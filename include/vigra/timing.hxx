/************************************************************************/
/*                                                                      */
/*               Copyright 2008-2011 by Ullrich Koethe                  */
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
#include <vector>

// Usage:
//
// void time_it()
// {
//     USETICTOC
//
//     TIC
//      ... // code to be timed
//     TOC
//      ... // untimed code
//     TIC
//      ... // other code to be timed
//     TOC
// }
//
// Intead of TOC which outputs the time difference to std::cerr, 
// you may use TOCN (the time difference in msec as a double)
// or TOCS (the time difference as a std::string). 
//
// Alternatively, you can performe nested timing like so:
//
// void time_it()
// {
//     USE_NESTED_TICTOC
//
//     TICPUSH
//      ...      // code to be timed
//     TICPUSH
//      ...      // nested code to be timed
//     TOC       // print time for nested code
//      ...      // more code to be timed
//     TOC       // print total time
// }
//
// Timings below 1 msec are generally subject to round-off errors. Under
// LINUX, you can #define VIGRA_HIRES_TIMING to get better
// accuracy, but this requires linking against librt.

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

    inline double tic_toc_diff_num(std::vector<LARGE_INTEGER> & tic)
    {
        double res = tic_toc_diff_num(tic.back());
        tic.pop_back();
        return res;
    }

    inline std::string tic_toc_diff_string(std::vector<LARGE_INTEGER> & tic)
    {
        std::string res = tic_toc_diff_string(tic.back());
        tic.pop_back();
        return res;
    }

    inline void tic_toc_diff(std::vector<LARGE_INTEGER> & tic)
    {
        tic_toc_diff(tic.back());
        tic.pop_back();
    }

    } // unnamed namespace
    
    #define USETICTOC LARGE_INTEGER tic_timer;
    #define USE_NESTED_TICTOC std::vector<LARGE_INTEGER> tic_timer;
    #define TIC QueryPerformanceCounter(&tic_timer);
    #define TICPUSH tic_timer.push_back(LARGE_INTEGER());\
                    QueryPerformanceCounter(&(tic_timer.back()));
    #define TOC  tic_toc_diff       (tic_timer);
    #define TOCN tic_toc_diff_num   (tic_timer)
    #define TOCS tic_toc_diff_string(tic_timer)

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

        inline std::string tic_toc_diff_string(timespec const & tic)
        {
            double diff = tic_toc_diff_num(tic); 
            std::stringstream s;
            s << diff << " msec";
            return s.str();
        }

        inline void tic_toc_diff(timespec const & tic)
        {
            std::cerr << tic_toc_diff_string(tic) << std::endl;
        }
        
        inline double tic_toc_diff_num(std::vector<timespec> & tic)
        {
            double res = tic_toc_diff_num(tic.back());
            tic.pop_back();
            return res;
        }

        inline std::string tic_toc_diff_string(std::vector<timespec> & tic)
        {
            std::string res = tic_toc_diff_string(tic.back());
            tic.pop_back();
            return res;
        }

        inline void tic_toc_diff(std::vector<timespec> & tic)
        {
            tic_toc_diff(tic.back());
            tic.pop_back();
        }

        } // unnamed namespace

        #define USETICTOC timespec tic_timer;
        #define TIC clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &tic_timer);
        #define TOC  tic_toc_diff       (tic_timer);
        #define TOCN tic_toc_diff_num   (tic_timer)
        #define TOCS tic_toc_diff_string(tic_timer)
        #define USE_NESTED_TICTOC std::vector<timespec> tic_timer;
        #define TICPUSH tic_timer.push_back(timespec());\
                        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &(tic_timer.back()));

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

        inline double tic_toc_diff_num(std::vector<timeval> & tic)
        {
            double res = tic_toc_diff_num(tic.back());
            tic.pop_back();
            return res;
        }

        inline std::string tic_toc_diff_string(std::vector<timeval> & tic)
        {
            std::string res = tic_toc_diff_string(tic.back());
            tic.pop_back();
            return res;
        }

        inline void tic_toc_diff(std::vector<timeval> & tic)
        {
            tic_toc_diff(tic.back());
            tic.pop_back();
        }

        } // unnamed namespace

        #define USETICTOC timeval tic_timer;
        #define TIC  gettimeofday       (&tic_timer, NULL);
        #define TOC  tic_toc_diff       (tic_timer);
        #define TOCN tic_toc_diff_num   (tic_timer)
        #define TOCS tic_toc_diff_string(tic_timer)
        #define USE_NESTED_TICTOC std::vector<timeval> tic_timer;
        #define TICPUSH tic_timer.push_back(timeval());\
                        gettimeofday(&(tic_timer.back()), NULL);

    #endif // VIGRA_HIRES_TIMING

#endif // WIN32

// TICTOCLOOP runs the body inner_repetitions times, and minimizes the result over a number of outer_repetitions runs,
//  outputting the final minimal average to std::cerr
#define TICTOCLOOP_BEGIN(inner_repetitions,outer_repetitions) \
    { \
	USETICTOC \
	    double tictoc_best_, tictoc_inner_repetitions_=inner_repetitions; size_t tictoc_outer_repetitions_=outer_repetitions; \
	    for (size_t tictoccounter_=0; tictoccounter_<tictoc_outer_repetitions_; ++tictoccounter_) { \
		TIC \
		for (size_t tictocinnercounter_=0; tictocinnercounter_<inner_repetitions; ++tictocinnercounter_) { \

		
#define TICTOCLOOP_END \
                } \
		const double tictoc_cur_ = TOCN; \
                if ((tictoccounter_==0) || (tictoc_cur_ < tictoc_best_)) \
		    tictoc_best_ = tictoc_cur_; \
	    } \
	    std::cerr << tictoc_best_/tictoc_inner_repetitions_ \
			 << " msec (best-of-" << tictoc_outer_repetitions_ << ")" << std::endl; \
    }\



#else // NDEBUG

#define USETICTOC 
#define TIC
#define TOC
#define TOCN 0.0
#define TICS ""
#define USE_NESTED_TICTOC
#define TICPUSH
#define TICTOCLOOP_BEGIN {
#define TICTOCLOOP_END }
#endif // NDEBUG



#endif // VIGRA_TIMING_HXX
