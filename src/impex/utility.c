/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2001 by Ullrich Koethe                  */
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
#include "error.h"
#include "utility.h"

#define TemporaryDirectory  "/tmp"
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  A p p e n d I m a g e F o r m a t                                          %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexAppendImageFormat appends the image format type to the filename.
%  If an extension to the file already exists, it is first removed.
%
%  The format of the vigraImpexAppendImageFormat routine is:
%
%       vigraImpexAppendImageFormat(format,filename)
%
%  A description of each parameter follows.
%
%   o  format:  Specifies a pointer to an array of characters.  This is the
%      format of the image.
%
%   o  filename:  Specifies a pointer to an array of characters.  The unique
%      file name is returned in this array.
%
%
*/
Export void vigraImpexAppendImageFormat(const char *format,char *filename)
{
  char
    staging[MaxTextExtent];

  register char
    *p;

  assert(format != (char *) NULL);
  assert(filename != (char *) NULL);
  if ((*format == '\0') || (*filename == '\0'))
    return;
  if (strcmp(filename,"-") == 0)
    {
      (void) sprintf(staging,"%s:%s",format,filename);
      (void) strcpy(filename,staging);
      return;
    }
  p=filename+Extent(filename)-1;
  while ((p > filename) && (*p != *BasenameSeparator))
  {
    if (*(p-1) == '.')
      {
        (void) strcpy(p,format);
        return;
      }
    p--;
  }
  (void) strcat(filename,".");
  (void) strcat(filename,format);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   B a s e F i l e n a m e                                                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexBaseFilename removes the path name component and any extensions.
%
%  The format of the vigraImpexBaseFilename function is:
%
%      vigraImpexBaseFilename(name)
%
%  A description of each parameter follows:
%
%    o name: Specifies a pointer to an character array that contains the
%      name.
%
%
*/
Export char *vigraImpexBaseFilename(const char *name)
{
  register char
    *p;

  static char
    basename[MaxTextExtent];

  /*
    Get basename of client.
  */
  assert(name != (char *) NULL);
  (void) strcpy(basename,name);
  p=basename+(Extent(basename)-1);
  while (p > basename)
  {
    if (*p == *BasenameSeparator)
      {
        (void) strcpy(basename,p+1);
        break;
      }
    p--;
  }
  /*
    Delete any extension.
  */
  p=basename+(Extent(basename)-1);
  while (p > basename)
  {
    if (*p == '.')
      {
        *p='\0';
        break;
      }
    p--;
  }
  return(basename);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   E x p a n d F i l e n a m e                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexExpandFilename expands '~' in a filename.
%
%  The format of the vigraImpexExpandFilename function is:
%
%      vigraImpexExpandFilename(filename)
%
%  A description of each parameter follows:
%
%    o filename: Specifies a pointer to an character array that contains the
%      filename.
%
%
*/
Export void vigraImpexExpandFilename(char *filename)
{
  char
    expanded_filename[MaxTextExtent];

  register char
    *p;

  if (filename == (char *) NULL)
    return;
  if (*filename != '~')
    return;
  (void) strcpy(expanded_filename,filename);
  if (*(filename+1) == '/')
    {
      /*
        Substitute ~ with $HOME.
      */
      p=(char *) getenv("HOME");
      if (p == (char *) NULL)
        p=".";
      (void) strcpy(expanded_filename,p);
      (void) strcat(expanded_filename,filename+1);
    }
  else
    {
#if !defined(vms) && !defined(macintosh) && !defined(WIN32)
      char
        username[MaxTextExtent];

      struct passwd
        *entry;

      /*
        Substitute ~ with home directory from password file.
      */
      (void) strcpy(username,filename+1);
      p=strchr(username,'/');
      if (p != (char *) NULL)
        *p='\0';
      entry=getpwnam(username);
      if (entry == (struct passwd *) NULL)
        return;
      (void) strcpy(expanded_filename,entry->pw_dir);
      if (p != (char *) NULL)
        {
          (void) strcat(expanded_filename,"/");
          (void) strcat(expanded_filename,p+1);
        }
#endif
    }
  (void) strcpy(filename,expanded_filename);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  I s A c c e s s i b l e                                                    %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexIsAccessible returns True if the file as defined by filename is
%  accessible.
%
%  The format of the vigraImpexIsAccessible routine is:
%
%       status=vigraImpexIsAccessible(filename)
%
%  A description of each parameter follows.
%
%   o  status:  Function vigraImpexIsAccessible returns True is the file as defined by
%      filename is accessible, otherwise False is returned.
%
%   o  filename:  Specifies a pointer to an array of characters.  The unique
%      file name is returned in this array.
%
%
*/
Export unsigned int vigraImpexIsAccessible(const char *filename)
{
  FILE
    *file;

  /*
    Return False if the file cannot be opened.
  */
  assert(filename != (char *) NULL);
  file=fopen(filename,ReadBinaryType);
  if (file == (FILE *) NULL)
    return(False);
  (void) fclose(file);
  return(True);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  I s D i r e c t o r y                                                      %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexIsDirectory returns True if the file as defined by filename is
%  a directory.  Once MetroWerks write a stat(2) function, we can remove the
%  chdir(2) function.
%
%  The format of the vigraImpexIsAccessible routine is:
%
%       status=vigraImpexIsDirectory(filename)
%
%  A description of each parameter follows.
%
%   o  status:  Function vigraImpexIsDirectory returns True is the file as defined by
%      filename is a directory, otherwise False is returned.
%
%   o  filename:  Specifies a pointer to an array of characters.  The unique
%      file name is returned in this array.
%
%
*/
Export unsigned int vigraImpexIsDirectory(const char *filename)
{
  int
    status;

#if !defined(WIN32)
  struct stat
    file_info;

  status=stat(filename,&file_info);
  if (status != 0)
    return(False);
  return(S_ISDIR(file_info.st_mode));
#else
  char
    current_directory[MaxTextExtent];

  (void) getcwd(current_directory,MaxTextExtent-1);
  status=chdir(filename);
  if (status == 0)
    (void) chdir(current_directory);
  return(status == 0);
#endif
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  L S B F i r s t R e a d L o n g                                            %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexLSBFirstReadLong reads a long value as a 32 bit quantity in
%  least-significant byte first order.
%
%  The format of the vigraImpexLSBFirstReadLong routine is:
%
%       value=vigraImpexLSBFirstReadLong(file)
%
%  A description of each parameter follows.
%
%    o value:  Function vigraImpexLSBFirstReadLong returns an unsigned long read from
%      the file.
%
%   o  file:  Specifies the file to read the data from.
%
%
*/
Export unsigned long vigraImpexLSBFirstReadLong(FILE *file)
{
  unsigned char
    buffer[4];

  unsigned int
    status;

  unsigned long
    value;

  assert(file != (FILE *) NULL);
  status=vigraImpexReadData((char *) buffer,1,4,file);
  if (status == False)
    return((unsigned long) ~0);
  value=(unsigned int) (buffer[3] << 24);
  value|=(unsigned int) (buffer[2] << 16);
  value|=(unsigned int) (buffer[1] << 8);
  value|=(unsigned int) (buffer[0]);
  return(value);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  L S B F i r s t R e a d S h o r t                                          %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexLSBFirstReadShort reads a short value as a 16 bit quantity in
%  least-significant byte first order.
%
%  The format of the vigraImpexLSBFirstReadShort routine is:
%
%       value=vigraImpexLSBFirstReadShort(file)
%
%  A description of each parameter follows.
%
%    o value:  Function vigraImpexLSBFirstReadShort returns an unsigned short read from
%      the file.
%
%   o  file:  Specifies the file to read the data from.
%
%
*/
Export unsigned short vigraImpexLSBFirstReadShort(FILE *file)
{
  unsigned char
    buffer[2];

  unsigned int
    status;

  unsigned short
    value;

  assert(file != (FILE *) NULL);
  status=vigraImpexReadData((char *) buffer,1,2,file);
  if (status == False)
    return((unsigned short) ~0);
  value=(unsigned short) (buffer[1] << 8);
  value|=(unsigned short) (buffer[0]);
  return(value);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  L S B F i r s t W r i t e L o n g                                          %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexLSBFirstWriteLong writes a long value as a 32 bit quantity in
%  least-significant byte first order.
%
%  The format of the vigraImpexLSBFirstWriteLong routine is:
%
%       vigraImpexLSBFirstWriteLong(value,file)
%
%  A description of each parameter follows.
%
%   o  value:  Specifies the value to write.
%
%   o  file:  Specifies the file to write the data to.
%
%
*/
Export void vigraImpexLSBFirstWriteLong(const unsigned long value,FILE *file)
{
  unsigned char
    buffer[4];

  assert(file != (FILE *) NULL);
  buffer[0]=(unsigned char) (value);
  buffer[1]=(unsigned char) ((value) >> 8);
  buffer[2]=(unsigned char) ((value) >> 16);
  buffer[3]=(unsigned char) ((value) >> 24);
  (void) fwrite((char *) buffer,1,4,file);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  L S B F i r s t W r i t e S h o r t                                        %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexLSBFirstWriteShort writes a long value as a 16 bit quantity in
%  least-significant byte first order.
%
%  The format of the vigraImpexLSBFirstWriteShort routine is:
%
%       vigraImpexLSBFirstWriteShort(value,file)
%
%  A description of each parameter follows.
%
%   o  value:  Specifies the value to write.
%
%   o  file:  Specifies the file to write the data to.
%
%
*/
Export void vigraImpexLSBFirstWriteShort(const unsigned int value,FILE *file)
{
  unsigned char
    buffer[2];

  assert(file != (FILE *) NULL);
  buffer[0]=(unsigned char) (value);
  buffer[1]=(unsigned char) ((value) >> 8);
  (void) fwrite((char *) buffer,1,2,file);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  M S B F i r s t O r d e r L o n g                                          %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexMSBFirstOrderLong converts a least-significant byte first buffer
%  of integers to most-significant byte first.
%
%  The format of the vigraImpexMSBFirstOrderLong routine is:
%
%       vigraImpexMSBFirstOrderLong(p,length);
%
%  A description of each parameter follows.
%
%   o  p:  Specifies a pointer to a buffer of integers.
%
%   o  length:  Specifies the length of the buffer.
%
%
*/
Export void vigraImpexMSBFirstOrderLong(register char *p,const unsigned int length)
{
  register char
    c,
    *q,
    *sp;

  assert(p != (char *) NULL);
  q=p+length;
  while (p < q)
  {
    sp=p+3;
    c=(*sp);
    *sp=(*p);
    *p++=c;
    sp=p+1;
    c=(*sp);
    *sp=(*p);
    *p++=c;
    p+=2;
  }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  M S B F i r s t O r d e r S h o r t                                        %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexMSBFirstOrderShort converts a least-significant byte first buffer
%  of integers to most-significant byte first.
%
%  The format of the vigraImpexMSBFirstOrderShort routine is:
%
%       vigraImpexMSBFirstOrderLongShort(p,length);
%
%  A description of each parameter follows.
%
%   o  p:  Specifies a pointer to a buffer of integers.
%
%   o  length:  Specifies the length of the buffer.
%
%
*/
Export void vigraImpexMSBFirstOrderShort(register char *p,const unsigned int length)
{
  register char
    c,
    *q;

  assert(p != (char *) NULL);
  q=p+length;
  while (p < q)
  {
    c=(*p);
    *p=(*(p+1));
    p++;
    *p++=c;
  }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  M S B F i r s t R e a d S h o r t                                          %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexMSBFirstReadShort reads a short value as a 16 bit quantity in
%  most-significant byte first order.
%
%  The format of the vigraImpexMSBFirstReadShort routine is:
%
%       value=vigraImpexMSBFirstReadShort(file)
%
%  A description of each parameter follows.
%
%    o value:  Function vigraImpexMSBFirstReadShort returns an unsigned short read from
%      the file.
%
%   o  file:  Specifies the file to read the data from.
%
%
*/
Export unsigned short vigraImpexMSBFirstReadShort(FILE *file)
{
  unsigned char
    buffer[2];

  unsigned int
    status;

  unsigned short
    value;

  assert(file != (FILE *) NULL);
  status=vigraImpexReadData((char *) buffer,1,2,file);
  if (status == False)
    return((unsigned short) ~0);
  value=(unsigned int) (buffer[0] << 8);
  value|=(unsigned int) (buffer[1]);
  return(value);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  M S B F i r s t R e a d L o n g                                            %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexMSBFirstReadLong reads a long value as a 32 bit quantity in
%  most-significant byte first order.
%
%  The format of the vigraImpexMSBFirstReadLong routine is:
%
%       value=vigraImpexMSBFirstReadLong(file)
%
%  A description of each parameter follows.
%
%    o value:  Function vigraImpexMSBFirstReadLong returns an unsigned long read from
%      the file.
%
%   o  file:  Specifies the file to read the data from.
%
%
*/
Export unsigned long vigraImpexMSBFirstReadLong(FILE *file)
{
  unsigned char
    buffer[4];

  unsigned int
    status;

  unsigned long
    value;

  assert(file != (FILE *) NULL);
  status=vigraImpexReadData((char *) buffer,1,4,file);
  if (status == False)
    return((unsigned long) ~0);
  value=(unsigned int) (buffer[0] << 24);
  value|=(unsigned int) (buffer[1] << 16);
  value|=(unsigned int) (buffer[2] << 8);
  value|=(unsigned int) (buffer[3]);
  return(value);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  M S B F i r s t W r i t e L o n g                                          %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexMSBFirstWriteLong writes a long value as a 32 bit quantity in
%  most-significant byte first order.
%
%  The format of the vigraImpexMSBFirstWriteLong routine is:
%
%       vigraImpexMSBFirstWriteLong(value,file)
%
%  A description of each parameter follows.
%
%   o  value:  Specifies the value to write.
%
%   o  file:  Specifies the file to write the data to.
%
%
*/
Export void vigraImpexMSBFirstWriteLong(const unsigned long value,FILE *file)
{
  unsigned char
    buffer[4];

  assert(file != (FILE *) NULL);
  buffer[0]=(unsigned char) ((value) >> 24);
  buffer[1]=(unsigned char) ((value) >> 16);
  buffer[2]=(unsigned char) ((value) >> 8);
  buffer[3]=(unsigned char) (value);
  (void) fwrite((char *) buffer,1,4,file);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  M S B F i r s t W r i t e S h o r t                                        %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexMSBFirstWriteShort writes a long value as a 16 bit quantity in
%  most-significant byte first order.
%
%  The format of the vigraImpexMSBFirstWriteShort routine is:
%
%       vigraImpexMSBFirstWriteShort(value,file)
%
%  A description of each parameter follows.
%
%   o  value:  Specifies the value to write.
%
%   o  file:  Specifies the file to write the data to.
%
%
*/
Export void vigraImpexMSBFirstWriteShort(const unsigned int value,FILE *file)
{
  unsigned char
    buffer[2];

  assert(file != (FILE *) NULL);
  buffer[0]=(unsigned char) ((value) >> 8);
  buffer[1]=(unsigned char) (value);
  (void) fwrite((char *) buffer,1,2,file);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  M u l t i l i n e C e n s u s                                              %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexMultilineCensus returns the number of lines within a label.  A line
%  is represented by a \n character.
%
%  The format of the MultilineCenus routine is:
%
%       MultilineCenus(label)
%
%  A description of each parameter follows.
%
%   o  label:  This character string is the label.
%
%
*/
Export int vigraImpexMultilineCensus(const char *label)
{
  int
    number_lines;

  /*
    Determine the number of lines within this label.
  */
  if (label == (char *) NULL)
    return(0);
  for (number_lines=1; *label != '\0'; label++)
    if (*label == '\n')
      number_lines++;
  return(number_lines);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  P o s t s c r i p t G e o m e t r y                                        %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexPostscriptGeometry replaces any page mneumonic with the equivalent
%  size in picas.
%
%  The format of the vigraImpexPostscriptGeometry routine is:
%
%       geometry=vigraImpexPostscriptGeometry(page)
%
%  A description of each parameter follows.
%
%   o  page:  Specifies a pointer to an array of characters.  The string is
%      either a Postscript page name (e.g. A4) or a postscript page geometry
%      (e.g. 612x792+36+36).
%
%
*/
Export char *vigraImpexPostscriptGeometry(const char *page)
{
  static char
    *PageSizes[][2]=
    {
      { "LETTER", "612x792" },
      { "TABLOID", "792x1224" },
      { "LEDGER", "1224x792" },
      { "LEGAL", " 612x1008" },
      { "STATEMENT", "396x612" },
      { "EXECUTIVE", "540x720" },
      { "A3", "842x1191" },
      { "A4", "595x842" },
      { "A5", "421x595" },
      { "B4", "729x1032" },
      { "B5", "516x729" },
      { "FOLIO", "612x936" },
      { "QUARTO", "610x780" },
      { "10x14", "720x1008" },
      { (char *) NULL, (char *) NULL }
    };

  char
    c,
    *geometry;

  register char
    *p;

  register int
    i;

  /*
    Allocate page geometry memory.
  */
  geometry=(char *) malloc((Extent(page)+MaxTextExtent)*sizeof(char));
  if (geometry == (char *) NULL)
    {
      vigraImpexMagickWarning(ResourceLimitWarning,"Unable to translate page geometry",
        "Memory allocation failed");
      return((char *) NULL);
    }
  *geometry='\0';
  if (page == (char *) NULL)
    return(geometry);
  /*
    Comparison is case insensitive.
  */
  (void) strcpy(geometry,page);
  if (!isdigit((int) (*geometry)))
    for (p=geometry; *p != '\0'; p++)
    {
      c=(*p);
      if (islower((int) c))
        *p=toupper(c);
    }
  /*
    Comparison is case insensitive.
  */
  for (i=0; *PageSizes[i] != (char *) NULL; i++)
    if (strncmp(PageSizes[i][0],geometry,Extent(PageSizes[i][0])) == 0)
      {
        /*
          Replace mneumonic with the equivalent size in dots-per-inch.
        */
        (void) strcpy(geometry,PageSizes[i][1]);
        (void) strcat(geometry,page+Extent(PageSizes[i][0]));
        break;
      }
  return(geometry);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  R e a d D a t a                                                            %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexReadData reads data from the image file and returns it.  If it
%  cannot read the requested number of items, False is returned indicating
%  an error.
%
%  The format of the vigraImpexReadData routine is:
%
%      status=vigraImpexReadData(data,size,number_items,file)
%
%  A description of each parameter follows:
%
%    o status:  Function vigraImpexReadData returns True if all the data requested
%      is obtained without error, otherwise False.
%
%    o data:  Specifies an area to place the information reuested from
%      the file.
%
%    o size:  Specifies an integer representing the length of an
%      individual item to be read from the file.
%
%    o number_items:  Specifies an integer representing the number of items
%      to read from the file.
%
%    o file:  Specifies a file to read the data.
%
%
*/
Export unsigned int vigraImpexReadData(char *data,const unsigned int size,
  const unsigned int number_items,FILE *file)
{
  long
    bytes,
    count;

  assert(data != (char *) NULL);
  assert(file != (FILE *) NULL);
  count=0;
  for (bytes=size*number_items; bytes > 0; bytes-=count)
  {
    count=(long) fread(data,1,bytes,file);
    if (count <= 0)
      return(False);
    data+=count;
  }
  return(True);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  R e a d D a t a B l o c k                                                  %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexReadDataBlock reads data from the image file and returns it.  The
%  amount of data is determined by first reading a count byte.  If
%  vigraImpexReadDataBlock cannot read the requested number of items, `-1' is returned
%  indicating an error.
%
%  The format of the vigraImpexReadData routine is:
%
%      status=vigraImpexReadData(data,file)
%
%  A description of each parameter follows:
%
%    o status:  Function vigraImpexReadData returns the number of characters read
%      unless there is an error, otherwise `-1'.
%
%    o data:  Specifies an area to place the information reuested from
%      the file.
%
%    o file:  Specifies a file to read the data.
%
%
*/
Export int vigraImpexReadDataBlock(char *data,FILE *file)
{
  unsigned char
    count;

  unsigned int
    status;

  assert(data != (char *) NULL);
  assert(file != (FILE *) NULL);
  status=vigraImpexReadData((char *) &count,1,1,file);
  if (status == False)
    return(-1);
  if (count == 0)
    return(0);
  status=vigraImpexReadData(data,1,(unsigned int) count,file);
  if (status == False)
    return(-1);
  return(count);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   S e t C l i e n t N a m e                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexSetClientName sets the client name if the name is specified.
%  Otherwise the current client name is returned.
%
%  The format of the vigraImpexSetClientName routine is:
%
%      client_name=vigraImpexSetClientName(name)
%
%  A description of each parameter follows:
%
%    o client_name: Function vigraImpexSetClientName returns the current client name.
%
%    o status: Specifies the new client name.
%
%
*/
Export char *vigraImpexSetClientName(const char *name)
{
  static char
    client_name[MaxTextExtent] = "Magick";

  if (name != (char *) NULL)
    (void) strcpy(client_name,vigraImpexBaseFilename(name));
  return(client_name);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  S t r i n g T o A r g v                                                    %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexStringToArgv converts a text string into command line arguments.
%
%  The format of the vigraImpexStringToArgv routine is:
%
%      argv=vigraImpexStringToArgv(text,argc)
%
%  A description of each parameter follows:
%
%    o argv:  Function vigraImpexStringToArgv returns the string list unless an error
%      occurs, otherwise NULL.
%
%    o text:  Specifies the string to segment into a list.
%
%    o argc:  This integer pointer returns the number of arguments in the
%      list.
%
%
*/
Export char **vigraImpexStringToArgv(char *text,int *argc)
{
  char
    **argv;

  register char
    *p,
    *q;

  register int
    i;

  *argc=0;
  if (text == (char *) NULL)
    return((char **) NULL);
  /*
    Determine the number of arguments.
  */
  for (p=text; *p != '\0'; )
  {
    while (isspace((int) (*p)))
      p++;
    (*argc)++;
    if (*p == '"')
      for (p++; (*p != '"') && (*p != '\0'); p++);
    if (*p == '\'')
      for (p++; (*p != '\'') && (*p != '\0'); p++);
    while (!isspace((int) (*p)) && (*p != '\0'))
      p++;
  }
  (*argc)++;
  argv=(char **) malloc((*argc+1)*sizeof(char *));
  if (argv == (char **) NULL)
    {
      vigraImpexMagickWarning(ResourceLimitWarning,"Unable to convert text",
        "Memory allocation failed");
      return((char **) NULL);
    }
  /*
    Convert string to an ASCII list.
  */
  argv[0]="magick";
  p=text;
  for (i=1; i < *argc; i++)
  {
    while (isspace((int) (*p)))
      p++;
    q=p;
    if (*q == '"')
      {
        p++;
        for (q++; (*q != '"') && (*q != '\0'); q++);
      }
    else
      if (*p == '\'')
        {
          p++;
          for (q++; (*q != '\'') && (*q != '\0'); q++);
        }
      else
        while (!isspace((int) (*q)) && (*q != '\0'))
          q++;
    argv[i]=(char *) malloc((q-p+1)*sizeof(char));
    if (argv[i] == (char *) NULL)
      {
        vigraImpexMagickWarning(ResourceLimitWarning,"Unable to convert text",
          "Memory allocation failed");
        return((char **) NULL);
      }
    (void) strncpy(argv[i],p,q-p);
    argv[i][q-p]='\0';
    p=q;
    while (!isspace((int) (*p)) && (*p != '\0'))
      p++;
  }
  argv[i]=(char *) NULL;
  return(argv);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  S t r i n g T o L i s t                                                    %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexStringToList converts a text string into a list by segmenting the
%  text string at each carriage return discovered.  The list is converted to
%  HEX characters if any control characters are discovered within the text
%  string.
%
%  The format of the vigraImpexStringToList routine is:
%
%      list=vigraImpexStringToList(text)
%
%  A description of each parameter follows:
%
%    o list:  Function vigraImpexStringToList returns the string list unless an error
%      occurs, otherwise NULL.
%
%    o text:  Specifies the string to segment into a list.
%
%
*/
Export char **vigraImpexStringToList(char *text)
{
  char
    **textlist;

  register char
    *p,
    *q;

  register int
    i;

  unsigned int
    lines;

  if (text == (char *) NULL)
    return((char **) NULL);
  for (p=text; *p != '\0'; p++)
    if (((unsigned char) *p < 32) && !isspace((int) (*p)))
      break;
  if (*p == '\0')
    {
      /*
        Convert string to an ASCII list.
      */
      lines=1;
      for (p=text; *p != '\0'; p++)
        if (*p == '\n')
          lines++;
      textlist=(char **) malloc((lines+1)*sizeof(char *));
      if (textlist == (char **) NULL)
        {
          vigraImpexMagickWarning(ResourceLimitWarning,"Unable to convert text",
            "Memory allocation failed");
          return((char **) NULL);
        }
      p=text;
      for (i=0; i < lines; i++)
      {
        for (q=p; *q != '\0'; q++)
          if ((*q == '\r') || (*q == '\n'))
            break;
        textlist[i]=(char *) malloc((q-p+1)*sizeof(char));
        if (textlist[i] == (char *) NULL)
          {
            vigraImpexMagickWarning(ResourceLimitWarning,"Unable to convert text",
              "Memory allocation failed");
            return((char **) NULL);
          }
        (void) strncpy(textlist[i],p,q-p);
        textlist[i][q-p]='\0';
        if (*q == '\r')
          q++;
        p=q+1;
      }
    }
  else
    {
      char
        hex_string[MaxTextExtent];

      register int
        j;

      /*
        Convert string to a HEX list.
      */
      lines=(Extent(text)/0x14)+1;
      textlist=(char **) malloc((lines+1)*sizeof(char *));
      if (textlist == (char **) NULL)
        {
          vigraImpexMagickWarning(ResourceLimitWarning,"Unable to convert text",
            "Memory allocation failed");
          return((char **) NULL);
        }
      p=text;
      for (i=0; i < lines; i++)
      {
        textlist[i]=(char *) malloc(900*sizeof(char));
        if (textlist[i] == (char *) NULL)
          {
            vigraImpexMagickWarning(ResourceLimitWarning,"Unable to convert text",
              "Memory allocation failed");
            return((char **) NULL);
          }
        (void) sprintf(textlist[i],"0x%08x: ",(unsigned int) (i*0x14));
        q=textlist[i]+Extent(textlist[i]);
        for (j=1; j <= Min(Extent(p),0x14); j++)
        {
          (void) sprintf(hex_string,"%02x",(unsigned int) (*(p+j)));
          (void) strcpy(q,hex_string);
          q+=2;
          if ((j % 0x04) == 0)
            *q++=' ';
        }
        for (; j <= 0x14; j++)
        {
          *q++=' ';
          *q++=' ';
          if ((j % 0x04) == 0)
            *q++=' ';
        }
        *q++=' ';
        for (j=1; j <= Min(Extent(p),0x14); j++)
        {
          if (isprint((int) (*p)))
            *q++=(*p);
          else
            *q++='-';
          p++;
        }
        *q='\0';
      }
    }
  textlist[i]=(char *) NULL;
  return(textlist);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   S t r i p                                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexStrip strips the whitespace from the beginning and end of a string
%  of characters.
%
%  The format of the vigraImpexStrip routine is:
%
%     vigraImpexStrip(data)
%
%  A description of each parameter follows:
%
%    o data: Specifies an array of characters.
%
%
*/
Export void vigraImpexStrip(char *data)
{
  long
    count;

  register char
    *p,
    *q;

  register int
    i;

  assert(data != (char *) NULL);
  if (*data == '\0')
    return;
  p=data;
  while (isspace((int) (*p)))
    p++;
  q=data+Extent(data)-1;
  while (isspace((int) (*q)) && (q > p))
    q--;
  count=q-p+1;
  q=data;
  for (i=0; i < count; i++)
    *q++=(*p++);
  *q='\0';
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  T e m p o r a r y F i l e n a m e                                          %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function TemporaryFilename replaces the contents of the string pointed to
%  by filename by a unique file name.
%
%  The format of the TemporaryFilename routine is:
%
%       TemporaryFilename(filename)
%
%  A description of each parameter follows.
%
%   o  filename:  Specifies a pointer to an array of characters.  The unique
%      file name is returned in this array.
%
%
*/
Export void TemporaryFilename(char *filename)
{
#if !defined(_MSC_VER)
  char
    *directory;
#endif

  assert(filename != (char *) NULL);
  *filename='\0';
#if !defined(_MSC_VER)
  directory=(char *) getenv("TMPDIR");
  if (directory == (char *) NULL)
    directory=(char *) getenv("TEMP");
  if (directory == (char *) NULL)
    directory=TemporaryDirectory;
  (void) sprintf(filename,TemporaryTemplate,directory);
  (void) mkstemp(filename);
#else
  (void) vigraImpexNTTemporaryFilename(filename);
#endif
}
