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
#include "vigra/impex.h"

#if !defined(_MSC_VER)
#if HAVE_SYS_NDIR_H || HAVE_SYS_DIR_H || HAVE_NDIR_H
# define dirent direct
# define NAMLEN(dirent) (dirent)->d_namlen
# if HAVE_SYS_NDIR_H
#  include <sys/ndir.h>
# endif
# if HAVE_SYS_DIR_H
#  include <sys/dir.h>
# endif
# if HAVE_NDIR_H
#  include <ndir.h>
# endif
#else
# include <dirent.h>
# define NAMLEN(dirent) Extent((dirent)->d_name)
#endif
#include <pwd.h>
#else
#if defined(_MSC_VER)
#include "nt.h"
#endif
#endif

/*
  ImageMagick include declarations.
*/
#if !defined(_MSC_VER)
#define Export
#else
#define Export  __declspec(dllexport)
#pragma warning( disable : 4018 )
#pragma warning( disable : 4244 )
#endif

#if !defined(S_ISDIR)
#define S_ISDIR(mode) (((mode) & S_IFMT) == S_IFDIR)
#endif


/*
  Define declarations.
*/
#define AbsoluteValue(x)  ((x) < 0 ? -(x) : (x))
#define DownShift(x) (((int) ((x)+(1L << 13))) >> 14)
#define Extent(string)  ((int) strlen(string))
#define False  0
#define Max(x,y)  (((x) > (y)) ? (x) : (y))
#define Min(x,y)  (((x) < (y)) ? (x) : (y))
#if !defined(M_PI)
#define M_PI  3.14159265358979323846
#endif
#define QuantumTick(i,image) \
  (((i+1) == image->packets) || ((i % image->rows) == 0))
#define Swap(x,y) ((x)^=(y), (y)^=(x), (x)^=(y))
#if !defined(STDIN_FILENO)
#define STDIN_FILENO  0
#endif
#define True  1
#define UpShift(x) ((int) (x) << 14)
#define UpShifted(x) ((int) ((x)*(1L << 14)+0.5))
#define DefaultInterlace  VigraImpexNoInterlace
#define DefaultImageQuality  "75"
#define DefaultPointSize  "12"
#define ReadBinaryType  "rb"
#define WriteBinaryType  "wb"

/*
  Image define declarations.
*/
#define ColorMatch(color,target,delta) \
  ((((int) ((color).red)-delta) <= (int) ((target).red)) && \
    ((int) ((target).red) <= ((int) ((color).red)+delta)) && \
   (((int) ((color).green)-delta) <= (int) ((target).green)) && \
    ((int) ((target).green) <= ((int) ((color).green)+delta)) && \
   (((int) ((color).blue)-delta) <= (int) ((target).blue)) && \
    ((int) ((target).blue) <= ((int) ((color).blue)+delta)))
#define DegreesToRadians(x) ((x)*M_PI/180.0)
#define Intensity(color)  \
  ((unsigned int) ((color).red*77+(color).green*150+(color).blue*29) >> 8)
#define IsFaxImage(color)  \
  (vigraImpexIsMonochromeImage(image) && ((image)->columns <= 2560))
#define IsGray(color)  \
  (((color).red == (color).green) && ((color).green == (color).blue))
#define MatteMatch(color,target,delta) \
  (ColorMatch(color,target,delta) && ((color).index == (target).index))
#define MaxColormapSize  65535L
#define MaxStacksize  (1 << 15)
#define MaxTextExtent  1664
#define PixelOffset(x,y) image->pixels+((y)*image->columns+(x))
#define Push(up,left,right,delta) \
  if ((p < (segment_stack+MaxStacksize)) && (((up)+(delta)) >= 0) && \
      (((up)+(delta)) < image->rows)) \
    { \
      p->y1=(up); \
      p->x1=(left); \
      p->x2=(right); \
      p->y2=(delta); \
      p++; \
    }
#define RadiansToDegrees(x) ((x)*180/M_PI)
#define ReadQuantum(quantum,p)  \
{  \
  if (image->depth == 8) \
    quantum=UpScale(*p++); \
  else \
    { \
      value=(*p++) << 8;  \
      value|=(*p++);  \
      quantum=value >> (image->depth-QuantumDepth); \
    } \
}
#define ReadQuantumFile(quantum)  \
{  \
  if (image->depth == 8) \
    quantum=UpScale(fgetc(image->file)); \
  else \
    quantum=vigraImpexMSBFirstReadShort(image->file) >> (image->depth-QuantumDepth); \
}
#define SharpenFactor  60.0
#define Transparent  0
#define UncompressImage  vigraImpexUncondenseImage
#define WriteQuantum(quantum,q)  \
{  \
  if (image->depth == 8) \
    *q++=DownScale(quantum); \
  else \
    { \
      value=(quantum); \
      if ((QuantumDepth-image->depth) > 0) \
        value*=257; \
      *q++=value >> 8; \
      *q++=value; \
    } \
}
#define WriteQuantumFile(quantum)  \
{  \
  if (image->depth == 8) \
    (void) fputc(DownScale(quantum),image->file); \
  else \
    if ((QuantumDepth-image->depth) > 0) \
      vigraImpexMSBFirstWriteShort((quantum)*257,image->file); \
    else \
      vigraImpexMSBFirstWriteShort(quantum,image->file); \
}

#if defined(QuantumLeap)
/*
  Color quantum is [0..65535].
*/
#define DownScale(quantum)  (((unsigned int) (quantum)) >> 8)
#define HexColorFormat "#%04x%04x%04x"
#define MaxRGB  65535L
#define MaxRunlength  65535L
#define Opaque  65535L
#define QuantumDepth  16
#define UpScale(quantum)  (((unsigned int) (quantum))*257)
#define XDownScale(color)  ((unsigned int) (color))
#define XUpScale(color)  ((unsigned int) (color))

typedef unsigned short Quantum;
#else
/*
  Color quantum is [0..255].
*/
#define DownScale(quantum)  ((unsigned int) (quantum))
#define HexColorFormat "#%02x%02x%02x"
#define MaxRGB  255
#define MaxRunlength  255
#define Opaque  255
#define QuantumDepth  8
#define UpScale(quantum)  ((unsigned int) (quantum))
#define XDownScale(color)  (((unsigned int) (color)) >> 8)
#define XUpScale(color)  (((unsigned int) (color))*257)

typedef unsigned char Quantum;
#endif
/*
  Utility define declarations.
*/
#if !defined(vms)
#define IsGlob(text) \
  ((strchr(text,'*') != (char *) NULL) || \
   (strchr(text,'?') != (char *) NULL) || \
   (strchr(text,'{') != (char *) NULL) || \
   (strchr(text,'}') != (char *) NULL) || \
   (strchr(text,'[') != (char *) NULL) || \
   (strchr(text,']') != (char *) NULL))
#else
#define IsGlob(text) \
  ((strchr(text,'*') != (char *) NULL) || \
   (strchr(text,'?') != (char *) NULL) || \
   (strchr(text,'{') != (char *) NULL) || \
   (strchr(text,'}') != (char *) NULL))
#endif
#if !defined(vms) && !defined(macintosh) && !defined(WIN32)
#define BasenameSeparator  "/"
#define DirectorySeparator  "/"
#define SystemCommand(command)  system(command)
#define TemporaryTemplate  "%s/magickXXXXXX"
#else
#if defined(vms)
#define BasenameSeparator  "]"
#define DirectorySeparator  ""
#define SystemCommand(command)  (!system(command))
#endif
#if defined(macintosh)
#define BasenameSeparator  ":"
#define DirectorySeparator  ":"
#define SystemCommand(command)  MACSystemCommand(command)
#endif
#if defined(WIN32)
#define BasenameSeparator  "/"
#define DirectorySeparator  "/"
#define SystemCommand(command)  NTSystemCommand(command)
#endif
#endif

/*
  Utilities routines.
*/
extern Export char
  *vigraImpexBaseFilename(const char *),
  **ListColors(const char *,int *),
  **ListFiles(char *,const char *,int *),
  *vigraImpexPostscriptGeometry(const char *),
  *vigraImpexSetClientName(const char *),
  **vigraImpexStringToArgv(char *,int *),
  **vigraImpexStringToList(char *);

extern Export int
  GlobExpression(char *,const char *),
  vigraImpexMultilineCensus(const char *),
  vigraImpexReadDataBlock(char *,FILE *);

extern Export unsigned int
  vigraImpexIsAccessible(const char *),
  vigraImpexIsDirectory(const char *),
  vigraImpexReadData(char *,const unsigned int,const unsigned int,FILE *);

extern Export unsigned long
  vigraImpexLSBFirstReadLong(FILE *),
  vigraImpexMSBFirstReadLong(FILE *);

extern Export unsigned short
  vigraImpexLSBFirstReadShort(FILE *),
  vigraImpexMSBFirstReadShort(FILE *);

extern Export void
  vigraImpexAppendImageFormat(const char *,char *),
  vigraImpexExpandFilename(char *),
  vigraImpexExpandFilenames(int *,char ***),
  LocaleFilename(char *),
  vigraImpexLSBFirstWriteLong(const unsigned long,FILE *),
  vigraImpexLSBFirstWriteShort(const unsigned int,FILE *),
  vigraImpexMSBFirstOrderLong(char *,const unsigned int),
  vigraImpexMSBFirstOrderShort(char *,const unsigned int),
  vigraImpexMSBFirstWriteLong(const unsigned long,FILE *),
  vigraImpexMSBFirstWriteShort(const unsigned int,FILE *),
  vigraImpexStrip(char *),
  TemporaryFilename(char *);

extern Export unsigned int vigraImpexIsPseudoClass(VigraImpexImage *image);
extern Export void vigraImpexGetQuantizeInfo(VigraImpexQuantizeInfo *quantize_info);
extern Export void vigraImpexQuantizeImage(VigraImpexQuantizeInfo *quantize_info,VigraImpexImage *image);
extern Export void vigraImpexNumberColors(VigraImpexImage *image,FILE *file);

