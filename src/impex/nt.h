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
  Windows NT specific include declarations.
*/
#ifndef _NT_H
#define _NT_H

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#include <windows.h>
#include <errno.h>

/*
  Define declarations.
*/
#if !defined(_MSC_VER)
#define S_IREAD  00400
#define S_IWRITE  00200
#endif

/*
  Typedef declarations.
*/
#if !defined(XS_VERSION)
typedef struct _DIR
{
  HANDLE
    hSearch;

  WIN32_FIND_DATA
    Win32FindData;
} DIR;

struct dirent
{
  char
     d_name[2048];
 
  int
    d_namlen;
};
#endif

/*
  NT utilities routines.
*/
extern __declspec(dllexport) char
  *vigraImpexSetClientName(const char *);


#if !defined(XS_VERSION)
extern __declspec(dllexport) DIR
  *vigraImpexopendir(char *);
 
extern __declspec(dllexport) int
  vigraImpexNTTemporaryFilename(char *);

extern __declspec(dllexport) long
  vigraImpextelldir(DIR *);

extern __declspec(dllexport) struct dirent
  *vigraImpexreaddir(DIR *);
 
extern __declspec(dllexport) void
  vigraImpexclosedir(DIR *),
  vigraImpexNTErrorHandler(const char *,const char *),
  vigraImpexNTWarningHandler(const char *,const char *),
  vigraImpexseekdir(DIR *,long);
#endif

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif
