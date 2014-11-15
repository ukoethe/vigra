/************************************************************************/
/*                                                                      */
/*     Copyright 2013-2014 by Martin Bidlingmaier and Ullrich Koethe    */
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

#ifndef VIGRA_PRINT_BACKTRACE_HXX
#define VIGRA_PRINT_BACKTRACE_HXX

/* Quick-and-dirty way to print a backtrace upon a signal in Linux.

   Especially useful if you can't use a debugger (e.g. on TravisCI).
   
   Usage:
   
   Make sure to compile in debug mode.
   Have "addr2line" installed (was already present on TravisCI and our local machines).
   
   #include <vigra/print_backtrace.hxx>
   
   int main(int argc, char** argv)
   {
       program_name = argv[0];
       signal(SIGSEGV, &vigra_print_backtrace);  // catch the desired signal
       
       run_buggy_code();
   }   
*/

#include <execinfo.h>
#include <stdio.h>
#include <stdlib.h>
	

static char * program_name;

static int vigra_addr2line(void const * const addr)
{
    char addr2line_cmd[512] = {0};
    sprintf(addr2line_cmd,"addr2line -C -f -p -i -e %.256s %p", program_name, addr);
    return system(addr2line_cmd);
}

static void vigra_print_backtrace(int sig)
{
    int i, trace_size = 0;
    char **messages = (char **)NULL;
    static const int BACKTRACE_SIZE = 100;
    void *stack_traces[BACKTRACE_SIZE];
    
    fprintf(stderr, "caught signal %d, printing backtrace\n\n", sig);

    trace_size = backtrace(stack_traces, BACKTRACE_SIZE);
    messages = backtrace_symbols(stack_traces, trace_size);

    for (i = 0; i < trace_size; ++i)
    {
        if (vigra_addr2line(stack_traces[i]) != 0)
        {
            fprintf(stderr, "  error determining line # for: %sn", messages[i]);
        }
    }
    if (messages) { free(messages); }
    exit(1);
}
    
#endif // VIGRA_PRINT_BACKTRACE_HXX
