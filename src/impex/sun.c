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
#include <stdio.h>
#include <stdlib.h>
#if defined(_MSC_VER)
#include <direct.h>
#else
#include <unistd.h>
#endif
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "vigra/impex.h"
#include "utility.h"
#include "sun.h"


/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   S U N D e c o d e I m a g e                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexSUNDecodeImage unpacks the packed image pixels into
%  runlength-encoded pixel packets.
%
%  The format of the vigraImpexSUNDecodeImage routine is:
%
%      status=vigraImpexSUNDecodeImage(compressed_pixels,pixels,number_columns,
%        number_rows)
%
%  A description of each parameter follows:
%
%    o status:  Function vigraImpexSUNDecodeImage returns True if all the pixels are
%      uncompressed without error, otherwise False.
%
%    o compressed_pixels:  The address of a byte (8 bits) array of compressed
%      pixel data.
%
%    o pixels:  The address of a byte (8 bits) array of pixel data created by
%      the uncompression process.  The number of bytes in this array
%      must be at least equal to the number columns times the number of rows
%      of the source pixels.
%
%    o number_columns:  An integer value that is the number of columns or
%      width in pixels of your source image.
%
%    o number_rows:  An integer value that is the number of rows or
%      heigth in pixels of your source image.
%
%
*/
Export unsigned int vigraImpexSUNDecodeImage(unsigned char *compressed_pixels,
  unsigned char *pixels,const unsigned int number_columns,
  const unsigned int number_rows)
{
  register int
    count;

  register unsigned char
    *p,
    *q;

  unsigned char
    byte;

  assert(compressed_pixels != (unsigned char *) NULL);
  assert(pixels != (unsigned char *) NULL);
  p=compressed_pixels;
  q=pixels;
  while ((q-pixels) <= (number_columns*number_rows))
  {
    byte=(*p++);
    if (byte != 128)
      *q++=byte;
    else
      {
        /*
          Runlength-encoded packet: <count><byte>
        */
        count=(*p++);
        if (count > 0)
          byte=(*p++);
        while (count >= 0)
        {
          *q++=byte;
          count--;
        }
     }
  }
  return(True);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   R e a d S U N I m a g e                                                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexReadSUNImage reads a SUN image file and returns it.  It allocates
%  the memory necessary for the new VigraImpexImage structure and returns a pointer to
%  the new image.
%
%  The format of the vigraImpexReadSUNImage routine is:
%
%      image=vigraImpexReadSUNImage(image_info)
%
%  A description of each parameter follows:
%
%    o image:  Function vigraImpexReadSUNImage returns a pointer to the image after
%      reading.  A null image is returned if there is a a memory shortage or
%      if the image cannot be read.
%
%    o image_info: Specifies a pointer to an VigraImpexImageInfo structure.
%
%
*/
VigraImpexImage *vigraImpexReadSUNImage(VigraImpexImageInfo *image_info)
{
  VigraImpexImage
    *image;

  register int
    bit,
    i,
    x,
    y;

  register VigraImpexRunlengthPacket
    *q;

  register unsigned char
    *p;

  SUNHeader
    sun_header;

  unsigned char
    *sun_data,
    *sun_pixels;

  unsigned int
    bytes_per_line,
    status;

  /*
    Allocate image structure.
  */
  image=vigraImpexAllocateImage(image_info);
  if (image == (VigraImpexImage *) NULL)
    return((VigraImpexImage *) NULL);
  /*
    Open image file.
  */
  vigraImpexOpenImage(image_info,image,ReadBinaryType);
  if (image->file == (FILE *) NULL)
  {
    fprintf(stderr, "vigraImpexReadSUNImage(): Unable to open file %s\n",image->filename);
    vigraImpexDestroyImage(image);
    return 0;
  }
  /*
    Read SUN raster header.
  */
  sun_header.magic=vigraImpexMSBFirstReadLong(image->file);
  do
  {
    /*
      Verify SUN identifier.
    */
    if (sun_header.magic != 0x59a66a95)
    {
        fprintf(stderr, "vigraImpexReadSUNImage(): %s isn't a SUN raster image\n",image->filename);
        vigraImpexDestroyImage(image);
        return 0;
    }
    sun_header.width=vigraImpexMSBFirstReadLong(image->file);
    sun_header.height=vigraImpexMSBFirstReadLong(image->file);
    sun_header.depth=vigraImpexMSBFirstReadLong(image->file);
    sun_header.length=vigraImpexMSBFirstReadLong(image->file);
    sun_header.type=vigraImpexMSBFirstReadLong(image->file);
    sun_header.maptype=vigraImpexMSBFirstReadLong(image->file);
    sun_header.maplength=vigraImpexMSBFirstReadLong(image->file);
    image->columns=(unsigned int) sun_header.width;
    image->rows=(unsigned int) sun_header.height;
    if (sun_header.depth < 24)
      {
        image->c_class=VigraImpexPseudoClass;
        image->colors=sun_header.maplength;
        if (sun_header.maptype == RMT_NONE)
          image->colors=1 << sun_header.depth;
        if (sun_header.maptype == RMT_EQUAL_RGB)
          image->colors=(unsigned int) sun_header.maplength/3;
      }
    if (image_info->ping)
      {
        vigraImpexCloseImage(image);
        return(image);
      }
    switch ((int)sun_header.maptype)
    {
      case RMT_NONE:
      {
        if (sun_header.depth < 24)
          {
            /*
              Create linear color ramp.
            */
            image->colormap=(VigraImpexColorPacket *)
              malloc(image->colors*sizeof(VigraImpexColorPacket));
            if (image->colormap == (VigraImpexColorPacket *) NULL)
            {
                fprintf(stderr, "vigraImpexReadSUNImage(): Memory allocation failed\n");
                vigraImpexDestroyImage(image);
                return 0;
            }
            for (i=0; i < image->colors; i++)
            {
              image->colormap[i].red=(MaxRGB*i)/(image->colors-1);
              image->colormap[i].green=(MaxRGB*i)/(image->colors-1);
              image->colormap[i].blue=(MaxRGB*i)/(image->colors-1);
            }
          }
        break;
      }
      case RMT_EQUAL_RGB:
      {
        unsigned char
          *sun_colormap;

        /*
          Read SUN raster colormap.
        */
        image->colormap=(VigraImpexColorPacket *)
          malloc(image->colors*sizeof(VigraImpexColorPacket));
        sun_colormap=(unsigned char *)
          malloc(image->colors*sizeof(unsigned char));
        if ((image->colormap == (VigraImpexColorPacket *) NULL) ||
            (sun_colormap == (unsigned char *) NULL))
        {
            fprintf(stderr, "vigraImpexReadSUNImage(): Memory allocation failed\n");
            vigraImpexDestroyImage(image);
            return 0;
        }
        (void) vigraImpexReadData((char *) sun_colormap,1,image->colors,image->file);
        for (i=0; i < image->colors; i++)
          image->colormap[i].red=UpScale(sun_colormap[i]);
        (void) vigraImpexReadData((char *) sun_colormap,1,image->colors,image->file);
        for (i=0; i < image->colors; i++)
          image->colormap[i].green=UpScale(sun_colormap[i]);
        (void) vigraImpexReadData((char *) sun_colormap,1,image->colors,image->file);
        for (i=0; i < image->colors; i++)
          image->colormap[i].blue=UpScale(sun_colormap[i]);
        free((char *) sun_colormap);
        break;
      }
      case RMT_RAW:
      {
        unsigned char
          *sun_colormap;

        /*
          Read SUN raster colormap.
        */
        sun_colormap=(unsigned char *)
          malloc(sun_header.maplength*sizeof(unsigned char));
        if (sun_colormap == (unsigned char *) NULL)
        {
            fprintf(stderr, "vigraImpexReadSUNImage(): Memory allocation failed\n");
            vigraImpexDestroyImage(image);
            return 0;
        }
        (void) vigraImpexReadData((char *) sun_colormap,1,(unsigned int)
          sun_header.maplength,image->file);
        free((char *) sun_colormap);
        break;
      }
      default:
            fprintf(stderr, "vigraImpexReadSUNImage(): Colormap type is not supported\n");
            vigraImpexDestroyImage(image);
            return 0;
    }
    sun_data=(unsigned char *) malloc(sun_header.length*sizeof(unsigned char));
    if (sun_data == (unsigned char *) NULL)
    {
        fprintf(stderr, "vigraImpexReadSUNImage(): Memory allocation failed\n");
        vigraImpexDestroyImage(image);
        return 0;
    }
    status=vigraImpexReadData((char *) sun_data,1,(unsigned int) sun_header.length,
      image->file);
    if ((status == False) && (sun_header.type != RT_ENCODED))
    {
        fprintf(stderr, "vigraImpexReadSUNImage(): Unable to read image data\n");
        vigraImpexDestroyImage(image);
        return 0;
    }
    sun_pixels=sun_data;
    if (sun_header.type == RT_ENCODED)
      {
        unsigned int
          height;

        /*
          Read run-length encoded raster pixels.
        */
        height=(unsigned int) sun_header.height;
        bytes_per_line=(2*sun_header.width*sun_header.depth+15)/16;
        sun_pixels=(unsigned char *)
          malloc(bytes_per_line*height*sizeof(unsigned char));
        if (sun_pixels == (unsigned char *) NULL)
        {
            fprintf(stderr, "vigraImpexReadSUNImage(): Memory allocation failed\n");
            vigraImpexDestroyImage(image);
            return 0;
        }
        (void) vigraImpexSUNDecodeImage(sun_data,sun_pixels,bytes_per_line,height);
        free((char *) sun_data);
      }
    /*
      Initialize image structure.
    */
    image->matte=(sun_header.depth == 32);
    image->columns=(unsigned int) sun_header.width;
    image->rows=(unsigned int) sun_header.height;
    image->packets=image->columns*image->rows;
    image->pixels=(VigraImpexRunlengthPacket *)
      malloc(image->packets*sizeof(VigraImpexRunlengthPacket));
    if (image->pixels == (VigraImpexRunlengthPacket *) NULL)
    {
        fprintf(stderr, "vigraImpexReadSUNImage(): Memory allocation failed\n");
        vigraImpexDestroyImage(image);
        return 0;
    }
    /*
      Convert SUN raster image to runlength-encoded packets.
    */
    p=sun_pixels;
    q=image->pixels;
    if (sun_header.depth == 1)
      for (y=0; y < image->rows; y++)
      {
        /*
          Convert bitmap scanline to runlength-encoded color packets.
        */
        for (x=0; x < (image->columns >> 3); x++)
        {
          for (bit=7; bit >= 0; bit--)
          {
            q->index=((*p) & (0x01 << bit) ? 0x00 : 0x01);
            q->length=0;
            q++;
          }
          p++;
        }
        if ((image->columns % 8) != 0)
          {
            for (bit=7; bit >= (8-(image->columns % 8)); bit--)
            {
              q->index=((*p) & (0x01 << bit) ? 0x00 : 0x01);
              q->length=0;
              q++;
            }
            p++;
          }
        if ((((image->columns/8)+(image->columns % 8 ? 1 : 0)) % 2) != 0)
          p++;
      }
    else
      if (image->c_class == VigraImpexPseudoClass)
        for (y=0; y < image->rows; y++)
        {
          /*
            Convert PseudoColor scanline to runlength-encoded color packets.
          */
          for (x=0; x < image->columns; x++)
          {
            q->index=(*p++);
            q->length=0;
            q++;
          }
          if ((image->columns % 2) != 0)
            p++;
        }
      else
        for (y=0; y < image->rows; y++)
        {
          /*
            Convert DirectColor scanline to runlength-encoded color packets.
          */
          for (x=0; x < image->columns; x++)
          {
            q->index=0;
            if (image->matte)
              q->index=UpScale(*p++);
            if (sun_header.type == RT_STANDARD)
              {
                q->blue=UpScale(*p++);
                q->green=UpScale(*p++);
                q->red=UpScale(*p++);
              }
            else
              {
                q->red=UpScale(*p++);
                q->green=UpScale(*p++);
                q->blue=UpScale(*p++);
              }
            if (image->colors != 0)
              {
                q->red=image->colormap[q->red].red;
                q->green=image->colormap[q->green].green;
                q->blue=image->colormap[q->blue].blue;
              }
            q->length=0;
            q++;
          }
          if (((image->columns % 2) != 0) && (image->matte == False))
            p++;
        }
    free((char *) sun_pixels);
    if (image->c_class == VigraImpexPseudoClass)
      vigraImpexSyncImage(image);
    /*
      Proceed to next image.
    */
    if (image_info->subrange != 0)
      if (image->scene >= (image_info->subimage+image_info->subrange-1))
        break;
    sun_header.magic=vigraImpexMSBFirstReadLong(image->file);
    if (sun_header.magic == 0x59a66a95)
      {
        /*
          Allocate next image structure.
        */
        vigraImpexAllocateNextImage(image_info,image);
        if (image->next == (VigraImpexImage *) NULL)
          {
            vigraImpexDestroyImage(image);
            return((VigraImpexImage *) NULL);
          }
        image=image->next;
      }
  } while (sun_header.magic == 0x59a66a95);
  vigraImpexCondenseImage(image);
  while (image->previous != (VigraImpexImage *) NULL)
    image=image->previous;
  vigraImpexCloseImage(image);
  return(image);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   W r i t e S U N I m a g e                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexWriteSUNImage writes an image in the SUN rasterfile format.
%
%  The format of the vigraImpexWriteSUNImage routine is:
%
%      status=vigraImpexWriteSUNImage(image_info,image)
%
%  A description of each parameter follows.
%
%    o status: Function vigraImpexWriteSUNImage return True if the image is written.
%      False is returned is there is a memory shortage or if the image file
%      fails to write.
%
%    o image_info: Specifies a pointer to an VigraImpexImageInfo structure.
%
%    o image:  A pointer to a VigraImpexImage structure.
%
%
*/
unsigned int vigraImpexWriteSUNImage(VigraImpexImageInfo *image_info,VigraImpexImage *image)
{
  register int
    i,
    j,
    x;

  register VigraImpexRunlengthPacket
    *p;

  SUNHeader
    sun_header;

  unsigned int
    scene;

  /*
    Open output image file.
  */
  vigraImpexOpenImage(image_info,image,WriteBinaryType);
  if (image->file == (FILE *) NULL)
  {
    fprintf(stderr,"vigraImpexWriteSUNImage(): Unable to open file %s\n",image->filename);
    return False;
  }
  scene=0;
  do
  {
    /*
      Initialize SUN raster file header.
    */
    sun_header.magic=0x59a66a95;
    sun_header.width=image->columns;
    sun_header.height=image->rows;
    sun_header.type=(image->c_class == VigraImpexDirectClass ? RT_FORMAT_RGB : RT_STANDARD);
    sun_header.maptype=RMT_NONE;
    sun_header.maplength=0;
    if (!vigraImpexIsPseudoClass(image) && !vigraImpexIsGrayImage(image))
      {
        /*
          Full color SUN raster.
        */
        sun_header.depth=(image->matte ? 32 : 24);
        sun_header.length=image->columns*image->rows*(image->matte ? 4 : 3);
        sun_header.length+=image->columns % 2 ? image->rows : 0;
      }
    else
      if (vigraImpexIsMonochromeImage(image))
        {
          /*
            Monochrome SUN raster.
          */
          sun_header.depth=1;
          sun_header.length=((image->columns+7) >> 3)*image->rows;
          sun_header.length+=((image->columns/8)+(image->columns % 8 ? 1 : 0)) %
            2 ? image->rows : 0;
        }
      else
        {
          /*
            Colormapped SUN raster.
          */
          sun_header.depth=8;
          sun_header.length=image->columns*image->rows;
          sun_header.length+=image->columns % 2 ? image->rows : 0;
          sun_header.maptype=RMT_EQUAL_RGB;
          sun_header.maplength=image->colors*3;
        }
    /*
      Write SUN header.
    */
    vigraImpexMSBFirstWriteLong(sun_header.magic,image->file);
    vigraImpexMSBFirstWriteLong(sun_header.width,image->file);
    vigraImpexMSBFirstWriteLong(sun_header.height,image->file);
    vigraImpexMSBFirstWriteLong(sun_header.depth,image->file);
    vigraImpexMSBFirstWriteLong(sun_header.length,image->file);
    vigraImpexMSBFirstWriteLong(sun_header.type,image->file);
    vigraImpexMSBFirstWriteLong(sun_header.maptype,image->file);
    vigraImpexMSBFirstWriteLong(sun_header.maplength,image->file);
    /*
      Convert MIFF to SUN raster pixels.
    */
    p=image->pixels;
    x=0;
    if (!vigraImpexIsPseudoClass(image) && !vigraImpexIsGrayImage(image))
      {
        /*
          Convert VigraImpexDirectClass packet to SUN RGB pixel.
        */
        for (i=0; i < image->packets; i++)
        {
          for (j=0; j <= ((int) p->length); j++)
          {
            if (image->matte)
              (void) fputc(DownScale(p->index),image->file);
            (void) fputc(DownScale(p->red),image->file);
            (void) fputc(DownScale(p->green),image->file);
            (void) fputc(DownScale(p->blue),image->file);
            x++;
            if (x == image->columns)
              {
                if ((image->columns % 2) != 0)
                  (void) fputc(0,image->file); /* pad scanline */
                x=0;
              }
          }
          p++;
        }
      }
    else
      if (vigraImpexIsMonochromeImage(image))
        {
          register unsigned char
            bit,
            byte,
            polarity;

          /*
            Convert VigraImpexPseudoClass image to a SUN monochrome image.
          */
          polarity=
            Intensity(image->colormap[0]) > Intensity(image->colormap[1]);
          bit=0;
          byte=0;
          for (i=0; i < image->packets; i++)
          {
            for (j=0; j <= ((int) p->length); j++)
            {
              byte<<=1;
              if (p->index == polarity)
                byte|=0x01;
              bit++;
              if (bit == 8)
                {
                  (void) fputc(byte,image->file);
                  bit=0;
                  byte=0;
                }
              x++;
              if (x == image->columns)
                {
                  /*
                    Advance to the next scanline.
                  */
                  if (bit != 0)
                    (void) fputc(byte << (8-bit),image->file);
                  if ((((image->columns/8)+
                      (image->columns % 8 ? 1 : 0)) % 2) != 0)
                    (void) fputc(0,image->file);  /* pad scanline */
                  bit=0;
                  byte=0;
                  x=0;
               }
            }
            p++;
          }
        }
      else
        {
          /*
            Dump colormap to file.
          */
          for (i=0; i < image->colors; i++)
            (void) fputc(DownScale(image->colormap[i].red),image->file);
          for (i=0; i < image->colors; i++)
            (void) fputc(DownScale(image->colormap[i].green),image->file);
          for (i=0; i < image->colors; i++)
            (void) fputc(DownScale(image->colormap[i].blue),image->file);
          /*
            Convert VigraImpexPseudoClass packet to SUN colormapped pixel.
          */
          for (i=0; i < image->packets; i++)
          {
            for (j=0; j <= ((int) p->length); j++)
            {
              (void) fputc(p->index,image->file);
              x++;
              if (x == image->columns)
                {
                  if ((image->columns % 2) != 0)
                    (void) fputc(0,image->file);  /* pad scanline */
                  x=0;
                }
            }
            p++;
          }
        }
    if (image->next == (VigraImpexImage *) NULL)
      break;
    image->next->file=image->file;
    image=image->next;
  } while (image_info->adjoin);
  vigraImpexCloseImage(image);
  return(True);
}
