/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
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
 

/*                                                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%  This file contains source code adapted from ImageMagick                    %
%                                                                             %
%  ImageMagick is Copyright 1998 E. I. du Pont de Nemours and Company         %
%                                                                             %
%  Permission is hereby granted, free of charge, to any person obtaining a    %
%  copy of this software and associated documentation files ("ImageMagick"),  %
%  to deal in ImageMagick without restriction, including without limitation   %
%  the rights to use, copy, modify, merge, publish, distribute, sublicense,   %
%  and/or sell copies of ImageMagick, and to permit persons to whom the       %
%  ImageMagick is furnished to do so, subject to the following conditions:    %
%                                                                             %
%  The above copyright notice and this permission notice shall be included in %
%  all copies or substantial portions of ImageMagick.                         %
%                                                                             %
%  The software is provided "as is", without warranty of any kind, express or %
%  implied, including but not limited to the warranties of merchantability,   %
%  fitness for a particular purpose and noninfringement.  In no event shall   %
%  E. I. du Pont de Nemours and Company be liable for any claim, damages or   %
%  other liability, whether in an action of contract, tort or otherwise,      %
%  arising from, out of or in connection with ImageMagick or the use or other %
%  dealings in ImageMagick.                                                   %
%                                                                             %
%  Except as contained in this notice, the name of the E. I. du Pont de       %
%  Nemours and Company shall not be used in advertising or otherwise to       %
%  promote the sale, use or other dealings in ImageMagick without prior       %
%  written authorization from the E. I. du Pont de Nemours and Company.       %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*
  Include declarations.
*/
#include <stdio.h>
#include <stdlib.h>
#if defined(_MSC_VER)
#    include <direct.h>
#else
#    include <unistd.h>
#endif
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "vigra/impex.h"
#include <errno.h>
#include "error.h"
#include "utility.h"

/*
  Define declarations.
*/
#if defined(sun) && !defined(SVR4)
#if !defined(strerror)
#define strerror(n) \
  (((n) >= 0 && (n) < sys_nerr) ? sys_errlist[n] : "unknown error")

extern char
  *sys_errlist[];

extern int
  sys_nerr;
#endif
#endif

/*
  Forward declaraations.
*/
static void
  DefaultErrorHandler(const char *,const char *),
  DefaultWarningHandler(const char *,const char *);

/*
  Global declarations.
*/
static ErrorHandler
  error_handler = DefaultErrorHandler,
  warning_handler = DefaultWarningHandler;

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   D e f a u l t E r r o r H a n d l e r                                     %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function DefaultErrorHandler displays an error message and then terminates
%  the program.
%
%  The format of the DefaultErrorHandler routine is:
%
%      DefaultErrorHandler(message,qualifier)
%
%  A description of each parameter follows:
%
%    o message: Specifies the message to display before terminating the
%      program.
%
%    o qualifier: Specifies any qualifier to the message.
%
%
*/
static void DefaultErrorHandler(const char *message,const char *qualifier)
{
  int
    magick_status;

  magick_status=vigraImpexSetErrorStatus(0);
  if (message == (char *) NULL)
    exit(magick_status % 100);
  (void) fprintf(stderr,"%s: %s",vigraImpexSetClientName((char *) NULL),message);
  if (qualifier != (char *) NULL)
    (void) fprintf(stderr," (%s)",qualifier);
  if (errno)
    (void) fprintf(stderr," [%s]",strerror(errno));
  (void) fprintf(stderr,".\n");
  exit(magick_status % 100);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   D e f a u l t W a r n i n g H a n d l e r                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function DefaultWarningHandler displays a warning message.
%
%  The format of the DefaultWarningHandler routine is:
%
+      DefaultWarningHandler(message,qualifier)
%
%  A description of each parameter follows:
%
%    o message: Specifies the message to display before terminating the
%      program.
%
%    o qualifier: Specifies any qualifier to the message.
%
%
*/
static void DefaultWarningHandler(const char *message,const char *qualifier)
{
  if (message == (char *) NULL)
    return;
  (void) fprintf(stderr,"%s: %s",vigraImpexSetClientName((char *) NULL),message);
  if (qualifier != (char *) NULL)
    (void) fprintf(stderr," (%s)",qualifier);
  (void) fprintf(stderr,".\n");
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   M a g i c k E r r o r                                                     %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexMagickError calls the error handler routines with an error message.
%
%  The format of the vigraImpexMagickError routine is:
%
%      vigraImpexMagickError(error,message,qualifier)
%
%  A description of each parameter follows:
%
%    o error: Specifies the numeric error category.
%
%    o message: Specifies the message to display before terminating the
%      program.
%
%    o qualifier: Specifies any qualifier to the message.
%
%
*/
Export void vigraImpexMagickError(const int error,const char *message,
  const char *qualifier)
{
  (void) vigraImpexSetErrorStatus(error);
  if (error_handler != (ErrorHandler) NULL)
    (*error_handler)(message,qualifier);
  errno=0;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   M a g i c k W a r n i n g                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexMagickWarning calls the warning handler routines with an warning
%  message.
%
%  The format of the vigraImpexMagickWarning routine is:
%
%      vigraImpexMagickWarning(warning,message,qualifier)
%
%  A description of each parameter follows:
%
%    o warning: Specifies the numeric warning category.
%
%    o message: Specifies the message to display before terminating the
%      program.
%
%    o qualifier: Specifies any qualifier to the message.
%
%
*/
Export void vigraImpexMagickWarning(const int warning,const char *message,
  const char *qualifier)
{
  (void) vigraImpexSetErrorStatus(warning);
  if (warning_handler != (ErrorHandler) NULL)
    (*warning_handler)(message,qualifier);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   S e t E r r o r H a n d l e r                                             %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexSetErrorHandler sets the error handler to the specified routine
%  and returns the previous error handler.
%
%  The format of the vigraImpexSetErrorHandler routine is:
%
%      previous_handler=vigraImpexSetErrorHandler(handler)
%
%  A description of each parameter follows:
%
%    o handler: Specifies a pointer to a routine to handle errors.
%
%
*/
Export ErrorHandler vigraImpexSetErrorHandler(ErrorHandler handler)
{
  ErrorHandler
    previous_handler;

  previous_handler=error_handler;
  error_handler=handler;
  return(previous_handler);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   S e t E r r o r S t a t u s                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexSetErrorStatus sets the error handler status and returns the
%  previous status.
%
%  The format of the vigraImpexSetErrorStatus routine is:
%
%      previous_status=vigraImpexSetErrorStatus(status)
%
%  A description of each parameter follows:
%
%    o previous_status: Function vigraImpexSetErrorStatus returns the current error
%      status.
%
%    o status: Specifies the new status value.
%
%
*/
Export int vigraImpexSetErrorStatus(int status)
{
  int
    previous_status;

  static int
    magick_status = 0;

  previous_status=magick_status;
  magick_status=status;
  return(previous_status);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   S e t W a r n i n g H a n d l e r                                         %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexSetWarningHandler sets the warning handler to the specified routine
%  and returns the previous warning handler.
%
%  The format of the vigraImpexSetWarningHandler routine is:
%
%      previous_handler=vigraImpexSetWarningHandler(handler)
%
%  A description of each parameter follows:
%
%    o handler: Specifies a pointer to a routine to handle warnings.
%
%
*/
Export ErrorHandler vigraImpexSetWarningHandler(ErrorHandler handler)
{
  ErrorHandler
    previous_handler;

  previous_handler=warning_handler;
  warning_handler=handler;
  return(previous_handler);
}
