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

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   A l l o c a t e I m a g e                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexAllocateImage allocates an VigraImpexImage structure and initializes each
%  field to a default value.
%
%  The format of the vigraImpexAllocateImage routine is:
%
%      allocated_image=vigraImpexAllocateImage(image_info)
%
%  A description of each parameter follows:
%
%    o allocated_image: Function vigraImpexAllocateImage returns a pointer to an image
%      structure initialized to default values.  A null image is returned if
%      there is a memory shortage.
%
%    o image_info: Specifies a pointer to a VigraImpexImageInfo structure.
%
%
*/
Export VigraImpexImage *vigraImpexAllocateImage(const VigraImpexImageInfo *image_info)
{
  VigraImpexImage
    *allocated_image;


#if 0
  XColor
    color;
#endif /* #if 0 */


  /*
    Allocate image structure.
  */
  allocated_image=(VigraImpexImage *) malloc(sizeof(VigraImpexImage));
  if (allocated_image == (VigraImpexImage *) NULL)
    {
      vigraImpexMagickWarning(ResourceLimitWarning,"Unable to allocate image",
        "Memory allocation failed");
      return((VigraImpexImage *) NULL);
    }
  /*
    Initialize VigraImpexImage structure.
  */
  allocated_image->file=(FILE *) NULL;
  allocated_image->status=False;
  allocated_image->temporary=False;
  *allocated_image->filename='\0';
  allocated_image->filesize=0;
  allocated_image->pipe=False;
  (void) strcpy(allocated_image->magick,"MIFF");
  allocated_image->comments=(char *) NULL;
  allocated_image->label=(char *) NULL;
  allocated_image->text=(char *) NULL;
  allocated_image->id=VigraImpexUndefinedId;
  allocated_image->c_class=VigraImpexDirectClass;
  allocated_image->matte=False;
  allocated_image->compression=VigraImpexRunlengthEncodedCompression;
  allocated_image->columns=0;
  allocated_image->rows=0;
  allocated_image->depth=QuantumDepth;
  allocated_image->interlace=DefaultInterlace;
  allocated_image->scene=0;
  allocated_image->number_scenes=1;
  allocated_image->units=VigraImpexUndefinedResolution;
  allocated_image->x_resolution=0.0;
  allocated_image->y_resolution=0.0;
  allocated_image->montage=(char *) NULL;
  allocated_image->directory=(char *) NULL;
  allocated_image->colormap=(VigraImpexColorPacket *) NULL;
  allocated_image->colors=0;
  allocated_image->gamma=0.0;
  allocated_image->normalized_maximum_error=0.0;
  allocated_image->normalized_mean_error=0.0;
  allocated_image->mean_error_per_pixel=0;
  allocated_image->total_colors=0;
  allocated_image->signature=(char *) NULL;
  allocated_image->pixels=(VigraImpexRunlengthPacket *) NULL;
  allocated_image->packet=(VigraImpexRunlengthPacket *) NULL;
  allocated_image->packets=0;
  allocated_image->packet_size=0;
  allocated_image->packed_pixels=(unsigned char *) NULL;
  *allocated_image->magick_filename='\0';
  allocated_image->magick_columns=0;
  allocated_image->magick_rows=0;
  allocated_image->magick_time=time((time_t *) NULL);
  allocated_image->geometry=(char *) NULL;
  allocated_image->page=(char *) NULL;
  allocated_image->dispose=0;
  allocated_image->delay=0;
  allocated_image->iterations=1;

#if 0
  (void) XQueryColorDatabase(BackgroundColor,&color);
  allocated_image->background_color.red=XDownScale(color.red);
  allocated_image->background_color.green=XDownScale(color.green);
  allocated_image->background_color.blue=XDownScale(color.blue);
  allocated_image->background_color.index=0;
  (void) XQueryColorDatabase(BorderColor,&color);
  allocated_image->border_color.red=XDownScale(color.red);
  allocated_image->border_color.green=XDownScale(color.green);
  allocated_image->border_color.blue=XDownScale(color.blue);
  allocated_image->border_color.index=0;
  (void) XQueryColorDatabase(MatteColor,&color);
  allocated_image->matte_color.red=XDownScale(color.red);
  allocated_image->matte_color.green=XDownScale(color.green);
  allocated_image->matte_color.blue=XDownScale(color.blue);
  allocated_image->matte_color.index=0;
#endif /* #if 0 */

  if (image_info != (VigraImpexImageInfo *) NULL)
    {
      (void) strcpy(allocated_image->filename,image_info->filename);
      (void) strcpy(allocated_image->magick,image_info->magick);
    }
  allocated_image->orphan=False;
  allocated_image->previous=(VigraImpexImage *) NULL;
  allocated_image->list=(VigraImpexImage *) NULL;
  allocated_image->next=(VigraImpexImage *) NULL;
  return(allocated_image);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   A l l o c a t e N e x t I m a g e                                         %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexAllocateNextImage allocates an VigraImpexImage structure and initializes each
%  field to a default value.
%
%  The format of the vigraImpexAllocateNextImage routine is:
%
%      vigraImpexAllocateImage(image_info,image)
%
%  A description of each parameter follows:
%
%    o image_info: Specifies a pointer to a VigraImpexImageInfo structure.
%
%    o image: The address of a structure of type VigraImpexImage.
%
%
*/
Export void vigraImpexAllocateNextImage(const VigraImpexImageInfo *image_info,VigraImpexImage *image)
{
  /*
    Allocate image structure.
  */
  assert(image != (VigraImpexImage *) NULL);
  if (image->packets == (image->columns*image->rows))
    vigraImpexCondenseImage(image);
  image->next=vigraImpexAllocateImage(image_info);
  if (image->next == (VigraImpexImage *) NULL)
    return;
  (void) strcpy(image->next->filename,image->filename);
  if (image_info != (VigraImpexImageInfo *) NULL)
    (void) strcpy(image->next->filename,image_info->filename);
  image->next->file=image->file;
  image->next->filesize=image->filesize;
  image->next->scene=image->scene+1;
  image->next->previous=image;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   C l o s e I m a g e                                                       %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexCloseImage closes a file associated with the image.  If the
%  filename prefix is '|', the file is a pipe and is closed with PipeClose.
%
%  The format of the vigraImpexCloseImage routine is:
%
%      vigraImpexCloseImage(image)
%
%  A description of each parameter follows:
%
%    o image: The address of a structure of type VigraImpexImage.
%
%
*/
Export void vigraImpexCloseImage(VigraImpexImage *image)
{
  /*
    Close image file.
  */
  assert(image != (VigraImpexImage *) NULL);
  if (image->file == (FILE *) NULL)
    return;
  image->status=ferror(image->file);
#if !defined(_MSC_VER)
  if (image->pipe)
    (void) pclose(image->file);
  else
#endif
    if ((image->file != stdin) && (image->file != stdout))
      (void) fclose(image->file);
  image->file=(FILE *) NULL;
  if (!image->orphan)
    do
    {
      image->file=(FILE *) NULL;
      image=image->next;
    }
    while (image != (VigraImpexImage *) NULL);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   C o n d e n s e I m a g e                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexCondenseImage compresses an image to the minimum number of
%  runlength-encoded packets.
%
%  The format of the vigraImpexCondenseImage routine is:
%
%      vigraImpexCondenseImage(image)
%
%  A description of each parameter follows:
%
%    o image: The address of a structure of type VigraImpexImage.
%
%
*/
Export void vigraImpexCondenseImage(VigraImpexImage *image)
{
  register int
    i;

  register VigraImpexRunlengthPacket
    *p,
    *q;

  /*
    Compress image.
  */
  assert(image != (VigraImpexImage *) NULL);
  if (image == (VigraImpexImage *) NULL)
    return;
  p=image->pixels;
  image->runlength=p->length+1;
  image->packets=0;
  q=image->pixels;
  q->length=MaxRunlength;
  if (image->matte)
    for (i=0; i < (image->columns*image->rows); i++)
    {
      if (image->runlength != 0)
        image->runlength--;
      else
        {
          p++;
          image->runlength=p->length;
        }
      if ((p->red == q->red) && (p->green == q->green) &&
          (p->blue == q->blue) && (p->index == q->index) &&
          ((int) q->length < MaxRunlength))
        q->length++;
      else
        {
          if (image->packets != 0)
            q++;
          image->packets++;
          *q=(*p);
          q->length=0;
        }
    }
  else
    for (i=0; i < (image->columns*image->rows); i++)
    {
      if (image->runlength != 0)
        image->runlength--;
      else
        {
          p++;
          image->runlength=p->length;
        }
      if ((p->red == q->red) && (p->green == q->green) &&
          (p->blue == q->blue) && ((int) q->length < MaxRunlength))
        q->length++;
      else
        {
          if (image->packets != 0)
            q++;
          image->packets++;
          *q=(*p);
          q->length=0;
        }
    }
  image->pixels=(VigraImpexRunlengthPacket *)
    realloc((char *) image->pixels,image->packets*sizeof(VigraImpexRunlengthPacket));
  /*
    Runlength-encode only if it takes up less space than no compression.
  */
  if (image->c_class == VigraImpexDirectClass)
    {
      if (image->packets >= ((image->columns*image->rows*3) >> 2))
        image->compression=VigraImpexNoCompression;
      return;
    }
  if (image->packets >= ((image->columns*image->rows) >> 1))
    image->compression=VigraImpexNoCompression;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   C l o n e I m a g e                                                       %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexCloneImage returns a copy of all fields of the input image.  The
%  the pixel memory is allocated but the pixel data is not copied.
%
%  The format of the vigraImpexCloneImage routine is:
%
%      clone_image=vigraImpexCloneImage(image,columns,rows,clone_pixels)
%
%  A description of each parameter follows:
%
%    o clone_image: Function vigraImpexCloneImage returns a pointer to the image after
%      copying.  A null image is returned if there is a memory shortage.
%
%    o image: The address of a structure of type VigraImpexImage.
%
%    o columns: An integer that specifies the number of columns in the copied
%      image.
%
%    o rows: An integer that specifies the number of rows in the copied
%      image.
%
%    o clone_pixels: Specifies whether the pixel data is copied.  Must be
%      either True or False;
%
%
*/
Export VigraImpexImage *vigraImpexCloneImage(VigraImpexImage *image,const unsigned int columns,
  const unsigned int rows,const unsigned int clone_pixels)
{
  VigraImpexImage
    *clone_image;

  register int
    i;

  /*
    Allocate image structure.
  */
  assert(image != (VigraImpexImage *) NULL);
  clone_image=(VigraImpexImage *) malloc(sizeof(VigraImpexImage));
  if (clone_image == (VigraImpexImage *) NULL)
    return((VigraImpexImage *) NULL);
  /*
    Allocate the image pixels.
  */
  *clone_image=(*image);
  clone_image->columns=columns;
  clone_image->rows=rows;
  if (clone_pixels)
    clone_image->pixels=(VigraImpexRunlengthPacket *)
      malloc((unsigned int) image->packets*sizeof(VigraImpexRunlengthPacket));
  else
    {
      clone_image->packets=clone_image->columns*clone_image->rows;
      clone_image->pixels=(VigraImpexRunlengthPacket *)
        malloc((unsigned int) clone_image->packets*sizeof(VigraImpexRunlengthPacket));
    }
  if (clone_image->pixels == (VigraImpexRunlengthPacket *) NULL)
    return((VigraImpexImage *) NULL);
  if (clone_pixels)
    {
      register VigraImpexRunlengthPacket
        *p,
        *q;

      /*
        Copy image pixels.
      */
      p=image->pixels;
      q=clone_image->pixels;
      for (i=0; i < image->packets; i++)
      {
        *q=(*p);
        p++;
        q++;
      }
    }
  clone_image->packed_pixels=(unsigned char *) NULL;
  if (image->colormap != (VigraImpexColorPacket *) NULL)
    {
      /*
        Allocate and copy the image colormap.
      */
      clone_image->colormap=(VigraImpexColorPacket *)
        malloc(image->colors*sizeof(VigraImpexColorPacket));
      if (clone_image->colormap == (VigraImpexColorPacket *) NULL)
        return((VigraImpexImage *) NULL);
      for (i=0; i < image->colors; i++)
        clone_image->colormap[i]=image->colormap[i];
    }
  if (image->comments != (char *) NULL)
    {
      /*
        Allocate and copy the image comments.
      */
      clone_image->comments=(char *)
        malloc((unsigned int) Extent(image->comments)+1);
      if (clone_image->comments == (char *) NULL)
        return((VigraImpexImage *) NULL);
      (void) strcpy(clone_image->comments,image->comments);
    }
  if (image->label != (char *) NULL)
    {
      /*
        Allocate and copy the image label.
      */
      clone_image->label=(char *) malloc((unsigned int) Extent(image->label)+1);
      if (clone_image->label == (char *) NULL)
        return((VigraImpexImage *) NULL);
      (void) strcpy(clone_image->label,image->label);
    }
  if (image->signature != (char *) NULL)
    {
      /*
        Allocate and copy the image signature.
      */
      clone_image->signature=(char *)
        malloc((unsigned int) Extent(image->signature)+1);
      if (clone_image->signature == (char *) NULL)
        return((VigraImpexImage *) NULL);
      (void) strcpy(clone_image->signature,image->signature);
    }
  clone_image->page=(char *) NULL;
  if ((image->columns == columns) && (image->rows == rows))
    if (image->page != (char *) NULL)
      {
        /*
          Allocate and copy the image page.
        */
        clone_image->page=(char *)
          malloc((unsigned int) Extent(image->page)+1);
        if (clone_image->page == (char *) NULL)
          return((VigraImpexImage *) NULL);
        (void) strcpy(clone_image->page,image->page);
      }
  clone_image->montage=(char *) NULL;
  if ((image->columns == columns) && (image->rows == rows))
    if (image->montage != (char *) NULL)
      {
        /*
          Allocate and copy the image montage.
        */
        clone_image->montage=(char *)
          malloc((unsigned int) Extent(image->montage)+1);
        if (clone_image->montage == (char *) NULL)
          return((VigraImpexImage *) NULL);
        (void) strcpy(clone_image->montage,image->montage);
      }
  clone_image->directory=(char *) NULL;
  if ((image->columns == columns) && (image->rows == rows))
    if (image->directory != (char *) NULL)
      {
        /*
          Allocate and copy the image directory.
        */
        clone_image->directory=(char *)
          malloc((unsigned int) Extent(image->directory)+1);
        if (clone_image->directory == (char *) NULL)
          return((VigraImpexImage *) NULL);
        (void) strcpy(clone_image->directory,image->directory);
      }
  if (image->orphan)
    {
      clone_image->file=(FILE *) NULL;
      clone_image->previous=(VigraImpexImage *) NULL;
      clone_image->next=(VigraImpexImage *) NULL;
    }
  else
    {
      /*
        Link image into image list.
      */
      if (clone_image->previous != (VigraImpexImage *) NULL)
        clone_image->previous->next=clone_image;
      if (clone_image->next != (VigraImpexImage *) NULL)
        clone_image->next->previous=clone_image;
    }
  clone_image->orphan=False;
  return(clone_image);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   D e s t r o y I m a g e                                                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexDestroyImage deallocates memory associated with an image.
%
%  The format of the vigraImpexDestroyImage routine is:
%
%      vigraImpexDestroyImage(image)
%
%  A description of each parameter follows:
%
%    o image: The address of a structure of type VigraImpexImage.
%
%
*/
Export void vigraImpexDestroyImage(VigraImpexImage *image)
{
  /*
    Close image.
  */
  assert(image != (VigraImpexImage *) NULL);
  if (image->file != (FILE *) NULL)
    {
      vigraImpexCloseImage(image);
      if (image->temporary)
        (void) remove(image->filename);
    }
  /*
    Deallocate the image comments.
  */
  if (image->comments != (char *) NULL)
    free((char *) image->comments);
  /*
    Deallocate the image label.
  */
  if (image->label != (char *) NULL)
    free((char *) image->label);
  /*
    Deallocate the image montage directory.
  */
  if (image->montage != (char *) NULL)
    free((char *) image->montage);
  if (image->directory != (char *) NULL)
    free((char *) image->directory);
  /*
    Deallocate the image colormap.
  */
  if (image->colormap != (VigraImpexColorPacket *) NULL)
    free((char *) image->colormap);
  /*
    Deallocate the image signature.
  */
  if (image->signature != (char *) NULL)
    free((char *) image->signature);
  /*
    Deallocate the image pixels.
  */
  if (image->pixels != (VigraImpexRunlengthPacket *) NULL)
    free((char *) image->pixels);
  if (image->packed_pixels != (unsigned char *) NULL)
    free((char *) image->packed_pixels);
  /*
    Deallocate the image page geometry.
  */
  if (image->page != (char *) NULL)
    free((char *) image->page);
  if (!image->orphan)
    {
      /*
        Unlink from linked list.
      */
      if (image->previous != (VigraImpexImage *) NULL)
      {
        if (image->next != (VigraImpexImage *) NULL)
          image->previous->next=image->next;
        else
          image->previous->next=(VigraImpexImage *) NULL;
      }
      if (image->next != (VigraImpexImage *) NULL)
      {
        if (image->previous != (VigraImpexImage *) NULL)
          image->next->previous=image->previous;
        else
          image->next->previous=(VigraImpexImage *) NULL;
      }
    }
  /*
    Deallocate the image structure.
  */
  free((char *) image);
  image=(VigraImpexImage *) NULL;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   D e s t r o y I m a g e I n f o                                           %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexDestroyImageInfo deallocates memory associated with an VigraImpexImageInfo
%  structure.
%
%  The format of the vigraImpexDestroyImageInfo routine is:
%
%      vigraImpexDestroyImageInfo(image_info)
%
%  A description of each parameter follows:
%
%    o image_info: Specifies a pointer to a VigraImpexImageInfo structure.
%
%
*/
Export void vigraImpexDestroyImageInfo(VigraImpexImageInfo *image_info)
{
  assert(image_info != (VigraImpexImageInfo *) NULL);
  free((char *) image_info->filename);
  image_info->filename=(char *) NULL;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   D e s t r o y I m a g e s                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexDestroyImage deallocates memory associated with a linked list
%  of images.
%
%  The format of the vigraImpexDestroyImage routine is:
%
%      vigraImpexDestroyImage(image)
%
%  A description of each parameter follows:
%
%    o image: The address of a structure of type VigraImpexImage.
%
%
*/
Export void vigraImpexDestroyImages(VigraImpexImage *image)
{
  VigraImpexImage
    *next_image;

  /*
    Proceed to the top of the image list.
  */
  assert(image != (VigraImpexImage *) NULL);
  while (image->previous != (VigraImpexImage *) NULL)
    image=image->previous;
  do
  {
    /*
      Destroy this image.
    */
    next_image=image->next;
    if (next_image != (VigraImpexImage *)NULL)
      next_image->previous=(VigraImpexImage *)NULL;
    vigraImpexDestroyImage(image);
    image=next_image;
  } while (image != (VigraImpexImage *) NULL);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G e t I m a g e I n f o                                                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexGetImageInfo initializes the VigraImpexImageInfo structure.
%
%  The format of the vigraImpexGetImageInfo routine is:
%
%      vigraImpexGetImageInfo(image_info)
%
%  A description of each parameter follows:
%
%    o image_info: Specifies a pointer to a VigraImpexImageInfo structure.
%
%
*/
Export void vigraImpexGetImageInfo(VigraImpexImageInfo *image_info)
{
  assert(image_info != (VigraImpexImageInfo *) NULL);
  *image_info->magick='\0';
  image_info->filename=(char *) malloc(MaxTextExtent);
  if (image_info->filename == (char *) NULL)
    vigraImpexMagickError(ResourceLimitError,"Unable to get image info",
      "Memory allocation failed");
  *image_info->filename='\0';
  image_info->affirm=False;
  image_info->subimage=0;
  image_info->subrange=0;
  image_info->server_name=(char *) NULL;
  image_info->font=(char *) NULL;
  image_info->pen=(char *) NULL;
  image_info->box=(char *) NULL;
  image_info->size=(char *) NULL;
  image_info->tile=(char *) NULL;
  image_info->density=(char *) NULL;
  image_info->page=(char *) NULL;
  image_info->dispose=(char *) NULL;
  image_info->delay=(char *) NULL;
  image_info->iterations=(char *) NULL;
  image_info->texture=(char *) NULL;
  image_info->view=(char *) NULL;
  image_info->adjoin=True;
  image_info->colorspace=VigraImpexRGBColorspace;
  image_info->compression=VigraImpexUndefinedCompression;
  image_info->dither=True;
  image_info->interlace=DefaultInterlace;
  image_info->monochrome=False;
  image_info->pointsize=atoi(DefaultPointSize);
  image_info->quality=atoi(DefaultImageQuality);
  image_info->verbose=False;
  image_info->preview_type=VigraImpexGammaPreview;
  image_info->undercolor=(char *) NULL;
  image_info->ping=False;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%     I s G r a y I m a g e                                                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexIsGrayImage returns True if the image is grayscale otherwise
%  False is returned.  If the image is VigraImpexDirectClass and grayscale, it is demoted
%  to VigraImpexPseudoClass.
%
%  The format of the vigraImpexIsGrayImage routine is:
%
%      status=vigraImpexIsGrayImage(image)
%
%  A description of each parameter follows:
%
%    o status: Function vigraImpexIsGrayImage returns True if the image is grayscale
%      otherwise False is returned.
%
%    o image: The address of a structure of type VigraImpexImage;  returned from
%      ReadImage.
%
%
*/
Export unsigned int vigraImpexIsGrayImage(VigraImpexImage *image)
{
  register int
    i;

  unsigned int
    gray_scale;

  /*
    Determine if image is grayscale.
  */
  assert(image != (VigraImpexImage *) NULL);
  gray_scale=True;
  switch (image->c_class)
  {
    case VigraImpexDirectClass:
    default:
    {
      register VigraImpexRunlengthPacket
        *p;

      if (image->matte)
        return(False);
      p=image->pixels;
      for (i=0; i < image->packets; i++)
      {
        if (!IsGray(*p))
          {
            gray_scale=False;
            break;
          }
        p++;
      }
      if (gray_scale)
        {
          VigraImpexQuantizeInfo
            quantize_info;

          vigraImpexGetQuantizeInfo(&quantize_info);
          quantize_info.colorspace=VigraImpexGRAYColorspace;
          vigraImpexQuantizeImage(&quantize_info,image);
          vigraImpexSyncImage(image);
        }
      break;
    }
    case VigraImpexPseudoClass:
    {
      for (i=0; i < image->colors; i++)
        if (!IsGray(image->colormap[i]))
          {
            gray_scale=False;
            break;
          }
      break;
    }
  }
  return(gray_scale);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%     I s M o n o c h r o m e I m a g e                                       %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexIsMonochromeImage returns True if the image is monochrome otherwise
%  False is returned.  If the image is VigraImpexDirectClass and monochrome, it is
%  demoted to VigraImpexPseudoClass.
%
%  The format of the vigraImpexIsMonochromeImage routine is:
%
%      status=vigraImpexIsMonochromeImage(image)
%
%  A description of each parameter follows:
%
%    o status: Function vigraImpexIsMonochromeImage returns True if the image is
%      monochrome otherwise False is returned.
%
%    o image: The address of a structure of type VigraImpexImage;  returned from
%      ReadImage.
%
%
*/
Export unsigned int vigraImpexIsMonochromeImage(VigraImpexImage *image)
{
  /*
    Determine if image is monochrome.
  */
  assert(image != (VigraImpexImage *) NULL);
  if (image->pixels == (VigraImpexRunlengthPacket *) NULL)
    return(False);
  if (!vigraImpexIsGrayImage(image))
    return(False);
  if (image->colors > 2)
    return(False);
  if ((Intensity(image->colormap[0]) != 0) &&
      (Intensity(image->colormap[0]) != MaxRGB))
    return(False);
  if (image->colors == 2)
    if ((Intensity(image->colormap[1]) != 0) &&
        (Intensity(image->colormap[1]) != MaxRGB))
      return(False);
  return(True);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%     N o r m a l i z e I m a g e                                             %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexNormalizeImage normalizes the pixel values to span the full
%  range of color values.  This is a contrast enhancement technique.
%
%  The format of the vigraImpexNormalizeImage routine is:
%
%      vigraImpexNormalizeImage(image)
%
%  A description of each parameter follows:
%
%    o image: The address of a structure of type VigraImpexImage;  returned from
%      ReadImage.
%
%
*/
Export void vigraImpexNormalizeImage(VigraImpexImage *image)
{
#define vigraImpexNormalizeImageText  "  Normalizing image...  "

  int
    histogram[MaxRGB+1],
    threshold_intensity;

  Quantum
    gray_value,
    normalize_map[MaxRGB+1];

  register int
    i,
    intensity;

  register VigraImpexRunlengthPacket
    *p;

  unsigned int
    high,
    low;

  /*
    Form histogram.
  */
  assert(image != (VigraImpexImage *) NULL);
  for (i=0; i <= MaxRGB; i++)
    histogram[i]=0;
  p=image->pixels;
  for (i=0; i < image->packets; i++)
  {
    gray_value=Intensity(*p);
    histogram[gray_value]+=p->length+1;
    p++;
  }
  /*
    Find the histogram boundaries by locating the 1 percent levels.
  */
  threshold_intensity=(image->columns*image->rows)/100;
  intensity=0;
  for (low=0; low < MaxRGB; low++)
  {
    intensity+=histogram[low];
    if (intensity > threshold_intensity)
      break;
  }
  intensity=0;
  for (high=MaxRGB; high != 0; high--)
  {
    intensity+=histogram[high];
    if (intensity > threshold_intensity)
      break;
  }
  if (low == high)
    {
      /*
        Unreasonable contrast;  use zero threshold to determine boundaries.
      */
      threshold_intensity=0;
      intensity=0;
      for (low=0; low < MaxRGB; low++)
      {
        intensity+=histogram[low];
        if (intensity > threshold_intensity)
          break;
      }
      intensity=0;
      for (high=MaxRGB; high != 0; high--)
      {
        intensity+=histogram[high];
        if (intensity > threshold_intensity)
          break;
      }
      if (low == high)
        return;  /* zero span bound */
    }
  /*
    Stretch the histogram to create the normalized image mapping.
  */
  for (i=0; i <= MaxRGB; i++)
    if (i < (int) low)
      normalize_map[i]=0;
    else
      if (i > (int) high)
        normalize_map[i]=MaxRGB;
      else
        normalize_map[i]=(MaxRGB-1)*(i-low)/(high-low);
  /*
    Normalize the image.
  */
  switch (image->c_class)
  {
    case VigraImpexDirectClass:
    default:
    {
      /*
        Normalize VigraImpexDirectClass image.
      */
      p=image->pixels;
      for (i=0; i < image->packets; i++)
      {
        p->red=normalize_map[p->red];
        p->green=normalize_map[p->green];
        p->blue=normalize_map[p->blue];
        p++;
      }
      break;
    }
    case VigraImpexPseudoClass:
    {
      /*
        Normalize VigraImpexPseudoClass image.
      */
      for (i=0; i < image->colors; i++)
      {
        image->colormap[i].red=normalize_map[image->colormap[i].red];
        image->colormap[i].green=normalize_map[image->colormap[i].green];
        image->colormap[i].blue=normalize_map[image->colormap[i].blue];
      }
      vigraImpexSyncImage(image);
      break;
    }
  }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   O p e n I m a g e                                                         %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexOpenImage open a file associated with the image.  A file name of
%  '-' sets the file to stdin for type 'r' and stdout for type 'w'.  If the
%  filename suffix is '.gz' or '.Z', the image is decompressed for type 'r'
%  and compressed for type 'w'.  If the filename prefix is '|', it is piped
%  to or from a system command.
%
%  The format of the vigraImpexOpenImage routine is:
%
%      vigraImpexOpenImage(image_info,image,type)
%
%  A description of each parameter follows:
%
%    o image_info: Specifies a pointer to a VigraImpexImageInfo structure.
%
%    o image: The address of a structure of type VigraImpexImage.
%
%    o type: 'r' for reading; 'w' for writing.
%
*/
Export void vigraImpexOpenImage(const VigraImpexImageInfo *image_info,VigraImpexImage *image,const char *type)
{
  char
    filename[MaxTextExtent];

  assert(image_info != (VigraImpexImageInfo *) NULL);
  assert(image != (VigraImpexImage *) NULL);
  assert(type != (char *) NULL);
  (void) strcpy(filename,image->filename);

#if 0
  if (*filename != '|')
    if ((Extent(filename) > 4) &&
        (strcmp(filename+Extent(filename)-4,".pgp") == 0))
      {
        /*
          Decrypt image file with PGP encryption utilities.
        */
        if (*type == 'r')
          (void) sprintf(filename,PgpvCommand,image->filename);
      }
    else
      if ((Extent(filename) > 4) &&
          (strcmp(filename+Extent(filename)-4,".bz2") == 0))
        {
          /*
            Uncompress/compress image file with BZIP compress utilities.
          */
          if (*type == 'r')
            (void) sprintf(filename,BunzipCommand,image->filename);
          else
            (void) sprintf(filename,BzipCommand,image->filename);
        }
      else
        if ((Extent(filename) > 3) &&
            (strcmp(filename+Extent(filename)-3,".gz") == 0))
          {
            /*
              Uncompress/compress image file with GNU compress utilities.
            */
            if (*type == 'r')
              (void) sprintf(filename,GunzipCommand,image->filename);
            else
              (void) sprintf(filename,GzipCommand,image->filename);
          }
        else
          if ((Extent(filename) > 2) &&
              (strcmp(filename+Extent(filename)-2,".Z") == 0))
            {
              /*
                Uncompress/compress image file with UNIX compress utilities.
              */
              if (*type == 'r')
                (void) sprintf(filename,UncompressCommand,image->filename);
              else
                (void) sprintf(filename,CompressCommand,image->filename);
            }
#endif /* #if 0 */

  /*
    Open image file.
  */
  image->pipe=False;
  if (strcmp(filename,"-") == 0)
    image->file=(*type == 'r') ? stdin : stdout;
  else
#if !defined(_MSC_VER)
    if (*filename == '|')
      {
        char
          mode[MaxTextExtent];

        /*
          Pipe image to or from a system command.
        */
        if (*type == 'w')
          (void) signal(SIGPIPE,SIG_IGN);
        (void) strncpy(mode,type,1);
        image->file=(FILE *) popen(filename+1,mode);
        image->pipe=True;
      }
    else
#endif
      {
        if (*type == 'w')
          {
            /*
              Form filename for multi-part images.
            */
            (void) sprintf(filename,image->filename,image->scene);
            if (!image_info->adjoin)
              if ((image->previous != (VigraImpexImage *) NULL) ||
                  (image->next != (VigraImpexImage *) NULL))
                {
                  if ((strcmp(filename,image->filename) == 0) ||
                      (strchr(filename,'%') != (char *) NULL))
                    (void) sprintf(filename,"%s.%u",filename,image->scene);
                  if (image->next != (VigraImpexImage *) NULL)
                    (void) strcpy(image->next->magick,image->magick);
                }
            (void) strcpy(image->filename,filename);
          }
#if defined(macintosh)
        if (*type == 'w')
          SetApplicationType(filename,image_info->magick,'8BIM');
#endif
        image->file=(FILE *) fopen(filename,type);
        if (image->file != (FILE *) NULL)
          {
            (void) fseek(image->file,0L,SEEK_END);
            image->filesize=ftell(image->file);
            (void) fseek(image->file,0L,SEEK_SET);
          }
      }
  image->status=False;
  if (*type == 'r')
    {
      image->next=(VigraImpexImage *) NULL;
      image->previous=(VigraImpexImage *) NULL;
    }
  return;
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%     R G B T r a n s f o r m I m a g e                                       %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexRGBTransformImage converts the reference image from RGB to
%  an alternate colorspace.  The transformation matrices are not the standard
%  ones: the weights are rescaled to normalized the range of the transformed
%  values to be [0..MaxRGB].
%
%  The format of the vigraImpexRGBTransformImage routine is:
%
%      vigraImpexRGBTransformImage(image,colorspace)
%
%  A description of each parameter follows:
%
%    o image: The address of a structure of type VigraImpexImage;  returned from
%      ReadImage.
%
%    o colorspace: An unsigned integer value that indicates which colorspace
%      to transform the image.
%
%
*/
Export void vigraImpexRGBTransformImage(VigraImpexImage *image,const unsigned int colorspace)
{
#define vigraImpexRGBTransformImageText  "  Transforming image colors...  "
#define X 0
#define Y (MaxRGB+1)
#define Z (MaxRGB+1)*2

  long
    tx,
    ty,
    tz,
    *x,
    *y,
    *z;

  Quantum
    *range_table;

  register int
    blue,
    green,
    i,
    red;

  register Quantum
    *range_limit;

  register VigraImpexRunlengthPacket
    *p;

  assert(image != (VigraImpexImage *) NULL);
  if ((colorspace == VigraImpexRGBColorspace) || (colorspace == VigraImpexTransparentColorspace) ||
      (colorspace == VigraImpexCMYKColorspace))
    return;
  if (colorspace == VigraImpexGRAYColorspace)
    {
      /*
        Return if the image is already gray_scale.
      */
      p=image->pixels;
      for (i=0; i < image->packets; i++)
      {
        if ((p->red != p->green) || (p->green != p->blue))
          break;
        p++;
      }
      if (i == image->packets)
        return;
    }
  /*
    Allocate the tables.
  */
  x=(long *) malloc(3*(MaxRGB+1)*sizeof(long));
  y=(long *) malloc(3*(MaxRGB+1)*sizeof(long));
  z=(long *) malloc(3*(MaxRGB+1)*sizeof(long));
  range_table=(Quantum *) malloc(4*(MaxRGB+1)*sizeof(Quantum));
  if ((x == (long *) NULL) || (y == (long *) NULL) ||
      (z == (long *) NULL) || (range_table == (Quantum *) NULL))
    {
      vigraImpexMagickWarning(ResourceLimitWarning,"Unable to transform color space",
        "Memory allocation failed");
      return;
    }
  /*
    Pre-compute conversion tables.
  */
  for (i=0; i <= MaxRGB; i++)
  {
    range_table[i]=0;
    range_table[i+(MaxRGB+1)]=(Quantum) i;
    range_table[i+(MaxRGB+1)*2]=MaxRGB;
  }
  for (i=0; i <= MaxRGB; i++)
    range_table[i+(MaxRGB+1)*3]=MaxRGB;
  range_limit=range_table+(MaxRGB+1);
  tx=0;
  ty=0;
  tz=0;
  switch (colorspace)
  {
    case VigraImpexGRAYColorspace:
    {
      /*
        Initialize GRAY tables:

          G = 0.29900*R+0.58700*G+0.11400*B
      */
      for (i=0; i <= MaxRGB; i++)
      {
        x[i+X]=UpShifted(0.29900)*i;
        y[i+X]=UpShifted(0.58700)*i;
        z[i+X]=UpShifted(0.11400)*i;
        x[i+Y]=UpShifted(0.29900)*i;
        y[i+Y]=UpShifted(0.58700)*i;
        z[i+Y]=UpShifted(0.11400)*i;
        x[i+Z]=UpShifted(0.29900)*i;
        y[i+Z]=UpShifted(0.58700)*i;
        z[i+Z]=UpShifted(0.11400)*i;
      }
      break;
    }
    case VigraImpexOHTAColorspace:
    {
      /*
        Initialize OHTA tables:

          I1 = 0.33333*R+0.33334*G+0.33333*B
          I2 = 0.50000*R+0.00000*G-0.50000*B
          I3 =-0.25000*R+0.50000*G-0.25000*B

        I and Q, normally -0.5 through 0.5, are normalized to the range 0
        through MaxRGB.
      */
      ty=UpShifted((MaxRGB+1) >> 1);
      tz=UpShifted((MaxRGB+1) >> 1);
      for (i=0; i <= MaxRGB; i++)
      {
        x[i+X]=UpShifted(0.33333)*i;
        y[i+X]=UpShifted(0.33334)*i;
        z[i+X]=UpShifted(0.33333)*i;
        x[i+Y]=UpShifted(0.50000)*i;
        y[i+Y]=0;
        z[i+Y]=(-UpShifted(0.50000))*i;
        x[i+Z]=(-UpShifted(0.25000))*i;
        y[i+Z]=UpShifted(0.50000)*i;
        z[i+Z]=(-UpShifted(0.25000))*i;
      }
      break;
    }
    case VigraImpexXYZColorspace:
    {
      /*
        Initialize CIE XYZ tables:

          X = 0.412453*X+0.357580*Y+0.180423*Z
          Y = 0.212671*X+0.715160*Y+0.072169*Z
          Z = 0.019334*X+0.119193*Y+0.950227*Z
      */
      for (i=0; i <= MaxRGB; i++)
      {
        x[i+X]=UpShifted(0.412453)*i;
        y[i+X]=UpShifted(0.357580)*i;
        z[i+X]=UpShifted(0.180423)*i;
        x[i+Y]=UpShifted(0.212671)*i;
        y[i+Y]=UpShifted(0.715160)*i;
        z[i+Y]=UpShifted(0.072169)*i;
        x[i+Z]=UpShifted(0.019334)*i;
        y[i+Z]=UpShifted(0.119193)*i;
        z[i+Z]=UpShifted(0.950227)*i;
      }
      break;
    }
    case VigraImpexYCbCrColorspace:
    {
      /*
        Initialize YCbCr tables:

          Y =  0.299000*R+0.587000*G+0.114000*B
          Cb= -0.172586*R-0.338828*G+0.511414*B
          Cr=  0.511414*R-0.428246*G-0.083168*B

        Cb and Cr, normally -0.5 through 0.5, are normalized to the range 0
        through MaxRGB.
      */
      ty=UpShifted((MaxRGB+1) >> 1);
      tz=UpShifted((MaxRGB+1) >> 1);
      for (i=0; i <= MaxRGB; i++)
      {
        x[i+X]=UpShifted(0.299000)*i;
        y[i+X]=UpShifted(0.587000)*i;
        z[i+X]=UpShifted(0.114000)*i;
        x[i+Y]=(-UpShifted(0.172586))*i;
        y[i+Y]=(-UpShifted(0.338828))*i;
        z[i+Y]=UpShifted(0.511414)*i;
        x[i+Z]=UpShifted(0.511414)*i;
        y[i+Z]=(-UpShifted(0.428246))*i;
        z[i+Z]=(-UpShifted(0.083168))*i;
      }
      break;
    }
    case VigraImpexYCCColorspace:
    {
      /*
        Initialize YCC tables:

          Y =  0.29900*R+0.58700*G+0.11400*B
          C1= -0.29900*R-0.58700*G+0.88600*B
          C2=  0.70100*R-0.58700*G-0.11400*B

        YCC is scaled by 1.3584.  C1 zero is 156 and C2 is at 137.
      */
      ty=UpShifted((unsigned int) UpScale(156));
      tz=UpShifted((unsigned int) UpScale(137));
      for (i=0; i <= (int) (0.018*MaxRGB); i++)
      {
        x[i+X]=(long) (UpShifted(0.29900/1.3584)*0.018*MaxRGB*i);
        y[i+X]=(long) (UpShifted(0.58700/1.3584)*0.018*MaxRGB*i);
        z[i+X]=(long) (UpShifted(0.11400/1.3584)*0.018*MaxRGB*i);
        x[i+Y]=(long) ((-UpShifted(0.29900/2.2179))*0.018*MaxRGB*i);
        y[i+Y]=(long) ((-UpShifted(0.58700/2.2179))*0.018*MaxRGB*i);
        z[i+Y]=(long) (UpShifted(0.88600/2.2179)*0.018*MaxRGB*i);
        x[i+Z]=(long) (UpShifted(0.70100/1.8215)*0.018*MaxRGB*i);
        y[i+Z]=(long) ((-UpShifted(0.58700/1.8215))*0.018*MaxRGB*i);
        z[i+Z]=(long) ((-UpShifted(0.11400/1.8215))*0.018*MaxRGB*i);
      }
      for ( ; i <= MaxRGB; i++)
      {
        x[i+X]=(long) (UpShifted(0.29900/1.3584)*(1.099*i-0.099));
        y[i+X]=(long) (UpShifted(0.58700/1.3584)*(1.099*i-0.099));
        z[i+X]=(long) (UpShifted(0.11400/1.3584)*(1.099*i-0.099));
        x[i+Y]=(long) ((-UpShifted(0.29900/2.2179))*(1.099*i-0.099));
        y[i+Y]=(long) ((-UpShifted(0.58700/2.2179))*(1.099*i-0.099));
        z[i+Y]=(long) (UpShifted(0.88600/2.2179)*(1.099*i-0.099));
        x[i+Z]=(long) (UpShifted(0.70100/1.8215)*(1.099*i-0.099));
        y[i+Z]=(long) ((-UpShifted(0.58700/1.8215))*(1.099*i-0.099));
        z[i+Z]=(long) ((-UpShifted(0.11400/1.8215))*(1.099*i-0.099));
      }
      break;
    }
    case VigraImpexYIQColorspace:
    {
      /*
        Initialize YIQ tables:

          Y = 0.29900*R+0.58700*G+0.11400*B
          I = 0.50000*R-0.23000*G-0.27000*B
          Q = 0.20200*R-0.50000*G+0.29800*B

        I and Q, normally -0.5 through 0.5, are normalized to the range 0
        through MaxRGB.
      */
      ty=UpShifted((MaxRGB+1) >> 1);
      tz=UpShifted((MaxRGB+1) >> 1);
      for (i=0; i <= MaxRGB; i++)
      {
        x[i+X]=UpShifted(0.29900)*i;
        y[i+X]=UpShifted(0.58700)*i;
        z[i+X]=UpShifted(0.11400)*i;
        x[i+Y]=UpShifted(0.50000)*i;
        y[i+Y]=(-UpShifted(0.23000))*i;
        z[i+Y]=(-UpShifted(0.27000))*i;
        x[i+Z]=UpShifted(0.20200)*i;
        y[i+Z]=(-UpShifted(0.50000))*i;
        z[i+Z]=UpShifted(0.29800)*i;
      }
      break;
    }
    case VigraImpexYPbPrColorspace:
    {
      /*
        Initialize YPbPr tables:

          Y =  0.299000*R+0.587000*G+0.114000*B
          Pb= -0.168736*R-0.331264*G+0.500000*B
          Pr=  0.500000*R-0.418688*G-0.081312*B

        Pb and Pr, normally -0.5 through 0.5, are normalized to the range 0
        through MaxRGB.
      */
      ty=UpShifted((MaxRGB+1) >> 1);
      tz=UpShifted((MaxRGB+1) >> 1);
      for (i=0; i <= MaxRGB; i++)
      {
        x[i+X]=UpShifted(0.299000)*i;
        y[i+X]=UpShifted(0.587000)*i;
        z[i+X]=UpShifted(0.114000)*i;
        x[i+Y]=(-UpShifted(0.168736))*i;
        y[i+Y]=(-UpShifted(0.331264))*i;
        z[i+Y]=UpShifted(0.500000)*i;
        x[i+Z]=UpShifted(0.500000)*i;
        y[i+Z]=(-UpShifted(0.418688))*i;
        z[i+Z]=(-UpShifted(0.081312))*i;
      }
      break;
    }
    case VigraImpexYUVColorspace:
    default:
    {
      /*
        Initialize YUV tables:

          Y =  0.29900*R+0.58700*G+0.11400*B
          U = -0.14740*R-0.28950*G+0.43690*B
          V =  0.61500*R-0.51500*G-0.10000*B

        U and V, normally -0.5 through 0.5, are normalized to the range 0
        through MaxRGB.  Note that U = 0.493*(B-Y), V = 0.877*(R-Y).
      */
      ty=UpShifted((MaxRGB+1) >> 1);
      tz=UpShifted((MaxRGB+1) >> 1);
      for (i=0; i <= MaxRGB; i++)
      {
        x[i+X]=UpShifted(0.29900)*i;
        y[i+X]=UpShifted(0.58700)*i;
        z[i+X]=UpShifted(0.11400)*i;
        x[i+Y]=(-UpShifted(0.14740))*i;
        y[i+Y]=(-UpShifted(0.28950))*i;
        z[i+Y]=UpShifted(0.43690)*i;
        x[i+Z]=UpShifted(0.61500)*i;
        y[i+Z]=(-UpShifted(0.51500))*i;
        z[i+Z]=(-UpShifted(0.10000))*i;
      }
      break;
    }
  }
  /*
    Convert from RGB.
  */
  switch (image->c_class)
  {
    case VigraImpexDirectClass:
    default:
    {
      /*
        Convert VigraImpexDirectClass image.
      */
      p=image->pixels;
      for (i=0; i < image->packets; i++)
      {
        red=p->red;
        green=p->green;
        blue=p->blue;
        p->red=range_limit[DownShift(x[red+X]+y[green+X]+z[blue+X]+tx)];
        p->green=range_limit[DownShift(x[red+Y]+y[green+Y]+z[blue+Y]+ty)];
        p->blue=range_limit[DownShift(x[red+Z]+y[green+Z]+z[blue+Z]+tz)];
        p++;
      }
      break;
    }
    case VigraImpexPseudoClass:
    {
      /*
        Convert VigraImpexPseudoClass image.
      */
      for (i=0; i < image->colors; i++)
      {
        red=image->colormap[i].red;
        green=image->colormap[i].green;
        blue=image->colormap[i].blue;
        image->colormap[i].red=
          range_limit[DownShift(x[red+X]+y[green+X]+z[blue+X]+tx)];
        image->colormap[i].green=
          range_limit[DownShift(x[red+Y]+y[green+Y]+z[blue+Y]+ty)];
        image->colormap[i].blue=
          range_limit[DownShift(x[red+Z]+y[green+Z]+z[blue+Z]+tz)];
      }
      vigraImpexSyncImage(image);
      break;
    }
  }
  /*
    Free allocated memory.
  */
  free((char *) range_table);
  free((char *) z);
  free((char *) y);
  free((char *) x);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   S y n c I m a g e                                                         %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexSyncImage initializes the red, green, and blue intensities of each
%  pixel as defined by the colormap index.
%
%  The format of the vigraImpexSyncImage routine is:
%
%      vigraImpexSyncImage(image)
%
%  A description of each parameter follows:
%
%    o image: The address of a structure of type VigraImpexImage.
%
%
*/
Export void vigraImpexSyncImage(VigraImpexImage *image)
{
  register int
    i;

  register VigraImpexRunlengthPacket
    *p;

  register unsigned short
    index;

  assert(image != (VigraImpexImage *) NULL);
  if (image->c_class == VigraImpexDirectClass)
    return;
  for (i=0; i < image->colors; i++)
  {
    image->colormap[i].index=0;
    image->colormap[i].flags=0;
  }
  p=image->pixels;
  for (i=0; i < image->packets; i++)
  {
    index=p->index;
    p->red=image->colormap[index].red;
    p->green=image->colormap[index].green;
    p->blue=image->colormap[index].blue;
    p++;
  }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%     T r a n s f o r m R G B I m a g e                                       %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexTransformRGBImage converts the reference image from an alternate
%  colorspace.  The transformation matrices are not the standard ones:  the
%  weights are rescaled to normalized the range of the transformed values to
%  be [0..MaxRGB].
%
%  The format of the vigraImpexTransformRGBImage routine is:
%
%      vigraImpexTransformRGBImage(image,colorspace)
%
%  A description of each parameter follows:
%
%    o image: The address of a structure of type VigraImpexImage;  returned from
%      ReadImage.
%
%    o colorspace: An unsigned integer value that indicates the colorspace
%      the image is currently in.  On return the image is in the RGB
%      color space.
%
%
*/
Export void vigraImpexTransformRGBImage(VigraImpexImage *image,const unsigned int colorspace)
{
#define B (MaxRGB+1)*2
#define G (MaxRGB+1)
#define R 0
#define vigraImpexTransformRGBImageText  "  Transforming image colors...  "

  static Quantum
    PCDMap[348] =  /* Photo CD information beyond 100% white, Gamma 2.2 */
    {
        0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  11,  12,  13,  14,
       15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,
       29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
       43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  55,
       56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  66,  67,  68,
       69,  70,  71,  72,  73,  74,  75,  76,  76,  77,  78,  79,  80,  81,
       82,  83,  84,  84,  85,  86,  87,  88,  89,  90,  91,  92,  92,  93,
       94,  95,  96,  97,  98,  99,  99, 100, 101, 102, 103, 104, 105, 106,
      106, 107, 108, 109, 110, 111, 112, 113, 114, 114, 115, 116, 117, 118,
      119, 120, 121, 122, 122, 123, 124, 125, 126, 127, 128, 129, 129, 130,
      131, 132, 133, 134, 135, 136, 136, 137, 138, 139, 140, 141, 142, 142,
      143, 144, 145, 146, 147, 148, 148, 149, 150, 151, 152, 153, 153, 154,
      155, 156, 157, 158, 158, 159, 160, 161, 162, 163, 164, 165, 165, 166,
      167, 168, 169, 170, 171, 172, 173, 173, 174, 175, 176, 177, 178, 178,
      179, 180, 181, 182, 182, 183, 184, 185, 186, 186, 187, 188, 189, 190,
      191, 192, 193, 194, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203,
      204, 205, 205, 206, 207, 208, 209, 210, 210, 211, 212, 213, 214, 215,
      216, 216, 217, 218, 219, 220, 221, 221, 222, 223, 224, 225, 225, 226,
      227, 228, 228, 229, 230, 230, 231, 232, 233, 233, 234, 235, 235, 236,
      237, 237, 238, 239, 239, 240, 241, 241, 242, 242, 243, 243, 244, 244,
      245, 245, 245, 246, 246, 247, 247, 247, 248, 248, 248, 249, 249, 249,
      250, 250, 250, 250, 251, 251, 251, 251, 252, 252, 252, 252, 252, 252,
      253, 253, 253, 253, 253, 253, 253, 253, 253, 253, 254, 254, 254, 254,
      254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254, 254,
      254, 254, 254, 254, 254, 254, 255, 255, 255, 255, 255, 255
    };

  long
    *blue,
    *green,
    *red;

  Quantum
    *range_table;

  register int
    i,
    x,
    y,
    z;

  register Quantum
    *range_limit;

  register VigraImpexRunlengthPacket
    *p;

  assert(image != (VigraImpexImage *) NULL);
  if ((colorspace == VigraImpexRGBColorspace) || (colorspace == VigraImpexGRAYColorspace) ||
      (colorspace == VigraImpexTransparentColorspace) || (colorspace == VigraImpexCMYKColorspace))
    return;
  /*
    Allocate the tables.
  */
  red=(long *) malloc(3*(MaxRGB+1)*sizeof(long));
  green=(long *) malloc(3*(MaxRGB+1)*sizeof(long));
  blue=(long *) malloc(3*(MaxRGB+1)*sizeof(long));
  range_table=(Quantum *) malloc(4*(MaxRGB+1)*sizeof(Quantum));
  if ((red == (long *) NULL) || (green == (long *) NULL) ||
      (blue == (long *) NULL) || (range_table == (Quantum *) NULL))
    {
      vigraImpexMagickWarning(ResourceLimitWarning,"Unable to transform color space",
        "Memory allocation failed");
      return;
    }
  /*
    Initialize tables.
  */
  for (i=0; i <= MaxRGB; i++)
  {
    range_table[i]=0;
    range_table[i+(MaxRGB+1)]=(Quantum) i;
    range_table[i+(MaxRGB+1)*2]=MaxRGB;
  }
  for (i=0; i <= MaxRGB; i++)
    range_table[i+(MaxRGB+1)*3]=MaxRGB;
  range_limit=range_table+(MaxRGB+1);
  switch (colorspace)
  {
    case VigraImpexOHTAColorspace:
    {
      /*
        Initialize OHTA tables:

          R = I1+1.00000*I2-0.66668*I3
          G = I1+0.00000*I2+1.33333*I3
          B = I1-1.00000*I2-0.66668*I3

        I and Q, normally -0.5 through 0.5, must be normalized to the range 0
        through MaxRGB.
      */
      for (i=0; i <= MaxRGB; i++)
      {
        red[i+R]=UpShifted(1.00000)*i;
        green[i+R]=UpShifted(1.0000*0.5)*((i << 1)-MaxRGB);
        blue[i+R]=(-UpShifted(0.66668*0.5))*((i << 1)-MaxRGB);
        red[i+G]=UpShifted(1.00000)*i;
        green[i+G]=0;
        blue[i+G]=UpShifted(1.33333*0.5)*((i << 1)-MaxRGB);
        red[i+B]=UpShifted(1.00000)*i;
        green[i+B]=(-UpShifted(1.00000*0.5))*((i << 1)-MaxRGB);
        blue[i+B]=(-UpShifted(0.66668*0.5))*((i << 1)-MaxRGB);
      }
      break;
    }
    case VigraImpexXYZColorspace:
    {
      /*
        Initialize CIE XYZ tables:

          R =  3.240479*R-1.537150*G-0.498535*B
          G = -0.969256*R+1.875992*G+0.041556*B
          B =  0.055648*R-0.204043*G+1.057311*B
      */
      for (i=0; i <= MaxRGB; i++)
      {
        red[i+R]=UpShifted(3.240479)*i;
        green[i+R]=(-UpShifted(1.537150))*i;
        blue[i+R]=(-UpShifted(0.498535))*i;
        red[i+G]=(-UpShifted(0.969256))*i;
        green[i+G]=UpShifted(1.875992)*i;
        blue[i+G]=UpShifted(0.041556)*i;
        red[i+B]=UpShifted(0.055648)*i;
        green[i+B]=(-UpShifted(0.204043))*i;
        blue[i+B]=UpShifted(1.057311)*i;
      }
      break;
    }
    case VigraImpexYCbCrColorspace:
    {
      /*
        Initialize YCbCr tables:

          R = Y            +1.370707*Cr
          G = Y-0.336453*Cb-0.698195*Cr
          B = Y+1.732445*Cb

        Cb and Cr, normally -0.5 through 0.5, must be normalized to the range 0
        through MaxRGB.
      */
      for (i=0; i <= MaxRGB; i++)
      {
        red[i+R]=UpShifted(1.000000)*i;
        green[i+R]=0;
        blue[i+R]=UpShifted(1.370707*0.5)*((i << 1)-MaxRGB);
        red[i+G]=UpShifted(1.000000)*i;
        green[i+G]=(-UpShifted(0.336453*0.5))*((i << 1)-MaxRGB);
        blue[i+G]=(-UpShifted(0.698195*0.5))*((i << 1)-MaxRGB);
        red[i+B]=UpShifted(1.000000)*i;
        green[i+B]=UpShifted(1.732445*0.5)*((i << 1)-MaxRGB);
        blue[i+B]=0;
      }
      break;
    }
    case VigraImpexYCCColorspace:
    {
      /*
        Initialize YCC tables:

          R = Y            +1.340762*C2
          G = Y-0.317038*C1-0.682243*C2
          B = Y+1.632639*C1

        YCC is scaled by 1.3584.  C1 zero is 156 and C2 is at 137.
      */
      for (i=0; i <= MaxRGB; i++)
      {
        red[i+R]=UpShifted(1.3584)*i;
        green[i+R]=0;
        blue[i+R]=UpShifted(1.8215)*(i-UpScale(137));
        red[i+G]=UpShifted(1.3584)*i;
        green[i+G]=(-(int) (UpShifted(0.194*2.2179)*(i-UpScale(156))));
        blue[i+G]=(-(int) (UpShifted(0.509*1.8215)*(i-UpScale(137))));
        red[i+B]=UpShifted(1.3584)*i;
        green[i+B]=UpShifted(2.2179)*(i-UpScale(156));
        blue[i+B]=0;
        range_table[i+(MaxRGB+1)]=(Quantum) UpScale(PCDMap[DownScale(i)]);
      }
      for ( ; i < UpScale(348); i++)
        range_table[i+(MaxRGB+1)]=(Quantum) UpScale(PCDMap[DownScale(i)]);
      break;
    }
    case VigraImpexYIQColorspace:
    {
      /*
        Initialize YIQ tables:

          R = 0.97087*Y+1.17782*I+0.59800*Q
          G = 0.97087*Y-0.28626*I-0.72851*Q
          B = 0.97087*Y-1.27870*I+1.72801*Q

        I and Q, normally -0.5 through 0.5, must be normalized to the range 0
        through MaxRGB.
      */
      for (i=0; i <= MaxRGB; i++)
      {
        red[i+R]=UpShifted(0.97087)*i;
        green[i+R]=UpShifted(1.17782*0.5)*((i << 1)-MaxRGB);
        blue[i+R]=UpShifted(0.59800*0.5)*((i << 1)-MaxRGB);
        red[i+G]=UpShifted(0.97087)*i;
        green[i+G]=(-UpShifted(0.28626*0.5))*((i << 1)-MaxRGB);
        blue[i+G]=(-UpShifted(0.72851*0.5))*((i << 1)-MaxRGB);
        red[i+B]=UpShifted(0.97087)*i;
        green[i+B]=(-UpShifted(1.27870*0.5))*((i << 1)-MaxRGB);
        blue[i+B]=UpShifted(1.72801*0.5)*((i << 1)-MaxRGB);
      }
      break;
    }
    case VigraImpexYPbPrColorspace:
    {
      /*
        Initialize YPbPr tables:

          R = Y            +1.402000*C2
          G = Y-0.344136*C1+0.714136*C2
          B = Y+1.772000*C1

        Pb and Pr, normally -0.5 through 0.5, must be normalized to the range 0
        through MaxRGB.
      */
      for (i=0; i <= MaxRGB; i++)
      {
        red[i+R]=UpShifted(1.000000)*i;
        green[i+R]=0;
        blue[i+R]=UpShifted(1.402000*0.5)*((i << 1)-MaxRGB);
        red[i+G]=UpShifted(1.000000)*i;
        green[i+G]=(-UpShifted(0.344136*0.5))*((i << 1)-MaxRGB);
        blue[i+G]=UpShifted(0.714136*0.5)*((i << 1)-MaxRGB);
        red[i+B]=UpShifted(1.000000)*i;
        green[i+B]=UpShifted(1.772000*0.5)*((i << 1)-MaxRGB);
        blue[i+B]=0;
      }
      break;
    }
    case VigraImpexYUVColorspace:
    default:
    {
      /*
        Initialize YUV tables:

          R = Y          +1.13980*V
          G = Y-0.39380*U-0.58050*V
          B = Y+2.02790*U

        U and V, normally -0.5 through 0.5, must be normalized to the range 0
        through MaxRGB.
      */
      for (i=0; i <= MaxRGB; i++)
      {
        red[i+R]=UpShifted(1.00000)*i;
        green[i+R]=0;
        blue[i+R]=UpShifted(1.13980*0.5)*((i << 1)-MaxRGB);
        red[i+G]=UpShifted(1.00000)*i;
        green[i+G]=(-UpShifted(0.39380*0.5))*((i << 1)-MaxRGB);
        blue[i+G]=(-UpShifted(0.58050*0.5))*((i << 1)-MaxRGB);
        red[i+B]=UpShifted(1.00000)*i;
        green[i+B]=UpShifted(2.02790*0.5)*((i << 1)-MaxRGB);
        blue[i+B]=0;
      }
      break;
    }
  }
  /*
    Convert to RGB.
  */
  switch (image->c_class)
  {
    case VigraImpexDirectClass:
    default:
    {
      /*
        Convert VigraImpexDirectClass image.
      */
      p=image->pixels;
      for (i=0; i < image->packets; i++)
      {
        x=p->red;
        y=p->green;
        z=p->blue;
        p->red=range_limit[DownShift(red[x+R]+green[y+R]+blue[z+R])];
        p->green=range_limit[DownShift(red[x+G]+green[y+G]+blue[z+G])];
        p->blue=range_limit[DownShift(red[x+B]+green[y+B]+blue[z+B])];
        p++;
      }
      break;
    }
    case VigraImpexPseudoClass:
    {
      /*
        Convert VigraImpexPseudoClass image.
      */
      for (i=0; i < image->colors; i++)
      {
        x=image->colormap[i].red;
        y=image->colormap[i].green;
        z=image->colormap[i].blue;
        image->colormap[i].red=
          range_limit[DownShift(red[x+R]+green[y+R]+blue[z+R])];
        image->colormap[i].green=
          range_limit[DownShift(red[x+G]+green[y+G]+blue[z+G])];
        image->colormap[i].blue=
          range_limit[DownShift(red[x+B]+green[y+B]+blue[z+B])];
      }
      p=image->pixels;
      for (i=0; i < image->packets; i++)
      {
        x=p->red;
        y=p->green;
        z=p->blue;
        p->red=range_limit[DownShift(red[x+R]+green[y+R]+blue[z+R])];
        p->green=range_limit[DownShift(red[x+G]+green[y+G]+blue[z+G])];
        p->blue=range_limit[DownShift(red[x+B]+green[y+B]+blue[z+B])];
        p++;
      }
      break;
    }
  }
  /*
    Free allocated memory.
  */
  free((char *) range_table);
  free((char *) blue);
  free((char *) green);
  free((char *) red);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   U n c o m p r e s s I m a g e                                             %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexUncondenseImage uncompresses runlength-encoded pixels packets to
%  a rectangular array of pixels.
%
%  The format of the vigraImpexUncondenseImage routine is:
%
%      status=vigraImpexUncondenseImage(image)
%
%  A description of each parameter follows:
%
%    o status: Function vigraImpexUncondenseImage returns True if the image is
%      uncompressed otherwise False.
%
%    o image: The address of a structure of type VigraImpexImage.
%
%
*/
Export unsigned int vigraImpexUncondenseImage(VigraImpexImage *image)
{
  int
    length;

  register int
    i,
    j;

  register VigraImpexRunlengthPacket
    *p,
    *q;

  VigraImpexRunlengthPacket
    *uncompressed_pixels;

  assert(image != (VigraImpexImage *) NULL);
  if (image->packets == (image->columns*image->rows))
    return(True);
  /*
    Uncompress runlength-encoded packets.
  */
  uncompressed_pixels=(VigraImpexRunlengthPacket *) realloc((char *) image->pixels,
    image->columns*image->rows*sizeof(VigraImpexRunlengthPacket));
  if (uncompressed_pixels == (VigraImpexRunlengthPacket *) NULL)
    {
      vigraImpexMagickWarning(ResourceLimitWarning,"Unable to uncompress image",
        "Memory allocation failed");
      return(False);
    }
  p=uncompressed_pixels+(image->packets-1);
  q=uncompressed_pixels+(image->columns*image->rows-1);
  for (i=0; i < image->packets; i++)
  {
    length=p->length;
    for (j=0; j <= length; j++)
    {
      *q=(*p);
      q->length=0;
      q--;
    }
    p--;
  }
  image->packets=image->columns*image->rows;
  image->pixels=uncompressed_pixels;
  return(True);
}
