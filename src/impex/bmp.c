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

/*
  Include declarations.
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
#include "bmp.h"
#include "utility.h"

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   B M P E n c o d e I m a g e                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexBMPEncodeImage compresses pixels using a runlength encoded format.
%
%  The format of the vigraImpexBMPEncodeImage routine is:
%
%      status=vigraImpexBMPEncodeImage(pixels,compressed_pixels,number_columns,
%        number_rows)
%
%  A description of each parameter follows:
%
%    o status:  Function vigraImpexBMPEncodeImage returns the number of bytes in the
%      runlength encoded compress_pixels array.
%
%    o pixels:  The address of a byte (8 bits) array of pixel data created by
%      the compression process.
%
%    o compressed_pixels:  The address of a byte (8 bits) array of compressed
%      pixel data.
%
%    o number_columns:  An integer value that is the number of columns or
%      width in pixels of your source image.
%
%    o number_rows:  An integer value that is the number of rows or
%      heigth in pixels of your source image.
%
%
*/
Export unsigned int vigraImpexBMPEncodeImage(unsigned char *pixels,
  unsigned char *compressed_pixels,const unsigned int number_columns,
  const unsigned int number_rows)
{
  register int
    i,
    x,
    y;

  register unsigned char
    *p,
    *q;

  /*
    Runlength encode pixels.
  */
  assert(pixels != (unsigned char *) NULL);
  assert(compressed_pixels != (unsigned char *) NULL);
  p=pixels;
  q=compressed_pixels;
  i=0;
  for (y=0; y < number_rows; y++)
  {
    for (x=0; x < number_columns; x+=i)
    {
      /*
        Determine runlength.
      */
      for (i=1; ((x+i) < number_columns); i++)
        if ((*(p+i) != *p) || (i == 255))
          break;
      *q++=i;
      *q++=(*p);
      p+=i;
    }
    /*
      End of line.
    */
    *q++=0;
    *q++=0x00;
  }
  /*
    End of bitmap.
  */
  *q++=0;
  *q++=0x01;
  return((unsigned int) (q-compressed_pixels));
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   B M P D e c o d e I m a g e                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexBMPDecodeImage unpacks the packed image pixels into
%  runlength-encoded pixel packets.
%
%  The format of the vigraImpexBMPDecodeImage routine is:
%
%      status=vigraImpexBMPDecodeImage(file,pixels,compression,number_columns,number_rows)
%
%  A description of each parameter follows:
%
%    o status:  Function vigraImpexBMPDecodeImage returns True if all the pixels are
%      uncompressed without error, otherwise False.
%
%    o file: The address of a structure of type FILE.  BMP encoded pixels
%      are read from this file.
%
%    o pixels:  The address of a byte (8 bits) array of pixel data created by
%      the decoding process.
%
%    o compression:  A value of 1 means the compressed pixels are runlength
%      encoded for a 256-color bitmap.  A value of 2 means a 16-color bitmap.
%
%    o number_columns:  An integer value that is the number of columns or
%      width in pixels of your source image.
%
%    o number_rows:  An integer value that is the number of rows or
%      heigth in pixels of your source image.
%
%
*/
Export unsigned int vigraImpexBMPDecodeImage(FILE *file,unsigned char *pixels,
  const unsigned int compression,const unsigned int number_columns,
  const unsigned int number_rows)
{
  int
    byte,
    count;

  register int
    i,
    x,
    y;

  register unsigned char
    *q;

  assert(file != (FILE *) NULL);
  assert(pixels != (unsigned char *) NULL);
  for (i=0; i < (number_columns*number_rows); i++)
    pixels[i]=0;
  byte=0;
  x=0;
  q=pixels;
  for (y=0; y < number_rows; )
  {
    count=fgetc(file);
    if (count != 0)
      {
        /*
          Encoded mode.
        */
        byte=fgetc(file);
        for (i=0; i < count; i++)
        {
          if (compression == 1)
            *q++=byte;
          else
            *q++=(i & 0x01) ? (byte & 0x0f) : ((byte >> 4) & 0x0f);
          x++;
        }
      }
    else
      {
        /*
          Escape mode.
        */
        count=fgetc(file);
        if (count == 0x01)
          return(True);
        switch (count)
        {
          case 0x00:
          {
            /*
              End of line.
            */
            x=0;
            y++;
            q=pixels+y*number_columns;
            break;
          }
          case 0x02:
          {
            /*
              Delta mode.
            */
            x+=fgetc(file);
            y+=fgetc(file);
            q=pixels+y*number_columns+x;
            break;
          }
          default:
          {
            /*
              Absolute mode.
            */
            for (i=0; i < count; i++)
            {
              if (compression == 1)
                *q++=fgetc(file);
              else
                {
                  if ((i & 0x01) == 0)
                    byte=fgetc(file);
                  *q++=(i & 0x01) ? (byte & 0x0f) : ((byte >> 4) & 0x0f);
                }
              x++;
            }
            /*
              Read pad byte.
            */
            if (compression == 1)
              {
                if (count & 0x01)
                  (void) fgetc(file);
              }
            else
              if (((count & 0x03) == 1) || ((count & 0x03) == 2))
                (void) fgetc(file);
            break;
          }
        }
      }
  }
  return(False);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   R e a d B M P I m a g e                                                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexReadBMPImage reads a Microsoft Windows bitmap image file and
%  returns it.  It allocates the memory necessary for the new VigraImpexImage structure
%  and returns a pointer to the new image.
%
%  The format of the vigraImpexReadBMPImage routine is:
%
%      image=vigraImpexReadBMPImage(image_info)
%
%  A description of each parameter follows:
%
%    o image:  Function vigraImpexReadBMPImage returns a pointer to the image after
%      reading.  A null image is returned if there is a a memory shortage or
%      if the image cannot be read.
%
%    o image_info: Specifies a pointer to an VigraImpexImageInfo structure.
%
%
*/
VigraImpexImage *vigraImpexReadBMPImage(VigraImpexImageInfo *image_info)
{

  BMPHeader
    bmp_header;

  VigraImpexImage
    *image;

  long
    start_position;

  register int
    bit,
    i,
    x,
    y;

  register VigraImpexRunlengthPacket
    *q;

  register unsigned char
    *p;

  unsigned char
    *bmp_pixels,
    magick[12];

  unsigned int
    bytes_per_line,
    image_size,
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
      fprintf(stderr, "vigraImpexReadBMPImage(): Unable to open file %s\n", image->filename);
      vigraImpexDestroyImage(image);
      return 0;
  }
  /*
    Determine if this is a BMP file.
  */
  status=vigraImpexReadData((char *) magick,1,2,image->file);
  do
  {
    /*
      Verify BMP identifier.
    */
    start_position=ftell(image->file)-2;
    if ((status == False) || (strncmp((char *) magick,"BM",2) != 0))
    {
      fprintf(stderr, "vigraImpexReadBMPImage(): %s isn't a BMP image file\n", image->filename);
      vigraImpexDestroyImage(image);
      return 0;
    }
    bmp_header.file_size=vigraImpexLSBFirstReadLong(image->file);
    bmp_header.reserved[0]=vigraImpexLSBFirstReadShort(image->file);
    bmp_header.reserved[1]=vigraImpexLSBFirstReadShort(image->file);
    bmp_header.offset_bits=vigraImpexLSBFirstReadLong(image->file);
    bmp_header.size=vigraImpexLSBFirstReadLong(image->file);
    if (bmp_header.size == 12)
      {
        /*
          OS/2 BMP image file.
        */
        bmp_header.width=(unsigned long) vigraImpexLSBFirstReadShort(image->file);
        bmp_header.height=(unsigned long) vigraImpexLSBFirstReadShort(image->file);
        bmp_header.planes=vigraImpexLSBFirstReadShort(image->file);
        bmp_header.bits_per_pixel=vigraImpexLSBFirstReadShort(image->file);
        bmp_header.x_pixels=0;
        bmp_header.y_pixels=0;
        bmp_header.number_colors=0;
        bmp_header.compression=0;
        bmp_header.image_size=0;
      }
    else
      {
        /*
          Microsoft Windows BMP image file.
        */
        bmp_header.width=vigraImpexLSBFirstReadLong(image->file);
        bmp_header.height=vigraImpexLSBFirstReadLong(image->file);
        bmp_header.planes=vigraImpexLSBFirstReadShort(image->file);
        bmp_header.bits_per_pixel=vigraImpexLSBFirstReadShort(image->file);
        bmp_header.compression=vigraImpexLSBFirstReadLong(image->file);
        bmp_header.image_size=vigraImpexLSBFirstReadLong(image->file);
        bmp_header.x_pixels=vigraImpexLSBFirstReadLong(image->file);
        bmp_header.y_pixels=vigraImpexLSBFirstReadLong(image->file);
        bmp_header.number_colors=vigraImpexLSBFirstReadLong(image->file);
        bmp_header.colors_important=vigraImpexLSBFirstReadLong(image->file);
        for (i=0; i < ((int) bmp_header.size-40); i++)
          (void) fgetc(image->file);
      }
    image->columns=(unsigned int) bmp_header.width;
    image->rows=(unsigned int) bmp_header.height;
    if ((bmp_header.number_colors != 0) || (bmp_header.bits_per_pixel < 16))
      {
        image->c_class=VigraImpexPseudoClass;
        image->colors=(unsigned int) bmp_header.number_colors;
        if (image->colors == 0)
          image->colors=1 << bmp_header.bits_per_pixel;
      }
    if (image_info->ping)
      {
        vigraImpexCloseImage(image);
        return(image);
      }
    if (image->c_class == VigraImpexPseudoClass)
      {
        /*
          Allocate image colormap.
        */
        image->colormap=(VigraImpexColorPacket *)
          malloc(image->colors*sizeof(VigraImpexColorPacket));
        if (image->colormap == (VigraImpexColorPacket *) NULL)
        {
          fprintf(stderr, "vigraImpexReadBMPImage(): Memory allocation failed\n");
          vigraImpexDestroyImage(image);
          return 0;
        }
        for (i=0; i < image->colors; i++)
        {
          image->colormap[i].red=(Quantum)
            ((long) (MaxRGB*i)/(image->colors-1));
          image->colormap[i].green=(Quantum)
            ((long) (MaxRGB*i)/(image->colors-1));
          image->colormap[i].blue=(Quantum)
            ((long) (MaxRGB*i)/(image->colors-1));
        }
        if (bmp_header.bits_per_pixel <= 8)
          {
            unsigned char
              *bmp_colormap;

            unsigned int
              packet_size;

            /*
              Read BMP raster colormap.
            */
            bmp_colormap=(unsigned char *)
              malloc(4*image->colors*sizeof(unsigned char));
            if (bmp_colormap == (unsigned char *) NULL)
            {
              fprintf(stderr, "vigraImpexReadBMPImage(): Memory allocation failed\n");
              vigraImpexDestroyImage(image);
              return 0;
            }
            packet_size=4;
            if (bmp_header.size == 12)
              packet_size=3;
            (void) vigraImpexReadData((char *) bmp_colormap,packet_size,image->colors,
              image->file);
            p=bmp_colormap;
            for (i=0; i < image->colors; i++)
            {
              image->colormap[i].blue=UpScale(*p++);
              image->colormap[i].green=UpScale(*p++);
              image->colormap[i].red=UpScale(*p++);
              if (bmp_header.size != 12)
                p++;
            }
            free((char *) bmp_colormap);
          }
      }
    /*
      Read image data.
    */
    while (ftell(image->file) < (start_position+bmp_header.offset_bits))
      (void) fgetc(image->file);
    if (bmp_header.compression == 2)
      bmp_header.bits_per_pixel<<=1;
    bytes_per_line=((image->columns*bmp_header.bits_per_pixel+31)/32)*4;
    image_size=bytes_per_line*bmp_header.height;
    bmp_pixels=(unsigned char *) malloc(image_size*sizeof(unsigned char));
    if (bmp_pixels == (unsigned char *) NULL)
    {
      fprintf(stderr, "vigraImpexReadBMPImage(): Memory allocation failed\n");
      vigraImpexDestroyImage(image);
      return 0;
    }
    if (bmp_header.compression == 0)
      (void) vigraImpexReadData((char *) bmp_pixels,1,(unsigned int) image_size,
        image->file);
    else
      {
        /*
          Convert run-length encoded raster pixels.
        */
        (void) vigraImpexBMPDecodeImage(image->file,bmp_pixels,
          (unsigned int) bmp_header.compression,(unsigned int) bmp_header.width,
          (unsigned int) bmp_header.height);
      }
    /*
      Initialize image structure.
    */
    image->columns=(unsigned int) bmp_header.width;
    image->rows=(unsigned int) bmp_header.height;
    image->units=VigraImpexPixelsPerCentimeterResolution;
    image->x_resolution=bmp_header.x_pixels/100.0;
    image->y_resolution=bmp_header.y_pixels/100.0;
    image->packets=image->columns*image->rows;
    image->pixels=(VigraImpexRunlengthPacket *)
      malloc(image->packets*sizeof(VigraImpexRunlengthPacket));
    if (image->pixels == (VigraImpexRunlengthPacket *) NULL)
    {
      fprintf(stderr, "vigraImpexReadBMPImage(): Memory allocation failed\n");
      vigraImpexDestroyImage(image);
      return 0;
    }
    /*
      Convert BMP raster image to runlength-encoded packets.
    */
    switch (bmp_header.bits_per_pixel)
    {
      case 1:
      {
        /*
          Convert bitmap scanline to runlength-encoded color packets.
        */
        for (y=image->rows-1; y >= 0; y--)
        {
          p=bmp_pixels+(image->rows-y-1)*bytes_per_line;
          q=image->pixels+(y*image->columns);
          for (x=0; x < ((int) image->columns-7); x+=8)
          {
            for (bit=0; bit < 8; bit++)
            {
              q->index=((*p) & (0x80 >> bit) ? 0x01 : 0x00);
              q->length=0;
              q++;
            }
            p++;
          }
          if ((image->columns % 8) != 0)
            {
              for (bit=0; bit < (image->columns % 8); bit++)
              {
                q->index=((*p) & (0x80 >> bit) ? 0x01 : 0x00);
                q->length=0;
                q++;
              }
              p++;
            }
        }
        vigraImpexSyncImage(image);
        break;
      }
      case 4:
      {
        /*
          Convert PseudoColor scanline to runlength-encoded color packets.
        */
        for (y=image->rows-1; y >= 0; y--)
        {
          p=bmp_pixels+(image->rows-y-1)*bytes_per_line;
          q=image->pixels+(y*image->columns);
          for (x=0; x < ((int) image->columns-1); x+=2)
          {
            q->index=(*p >> 4) & 0xf;
            q->length=0;
            q++;
            q->index=(*p) & 0xf;
            q->length=0;
            p++;
            q++;
          }
          if ((image->columns % 2) != 0)
            {
              q->index=(*p >> 4) & 0xf;
              q->length=0;
              q++;
              p++;
            }
        }
        vigraImpexSyncImage(image);
        break;
      }
      case 8:
      {
        /*
          Convert PseudoColor scanline to runlength-encoded color packets.
        */
        if ((bmp_header.compression == 1) || (bmp_header.compression == 2))
          bytes_per_line=image->columns;
        for (y=image->rows-1; y >= 0; y--)
        {
          p=bmp_pixels+(image->rows-y-1)*bytes_per_line;
          q=image->pixels+(y*image->columns);
          for (x=0; x < image->columns; x++)
          {
            q->index=(*p++);
            q->length=0;
            q++;
          }
        }
        vigraImpexSyncImage(image);
        break;
      }
      case 16:
      {
        /*
          Convert PseudoColor scanline to runlength-encoded color packets.
        */
        if (bmp_header.compression == 1)
          bytes_per_line=image->columns << 1;
        for (y=image->rows-1; y >= 0; y--)
        {
          p=bmp_pixels+(image->rows-y-1)*bytes_per_line;
          q=image->pixels+(y*image->columns);
          for (x=0; x < image->columns; x++)
          {
            q->index=(*p++);
            q->index|=(*p++) << 8;
            q->red=XDownScale(q->index);
            q->green=XDownScale(q->index);
            q->blue=XDownScale(q->index);
            q->length=0;
            q++;
          }
        }
        if (image->c_class == VigraImpexPseudoClass)
          vigraImpexSyncImage(image);
        break;
      }
      case 24:
      {
        /*
          Convert DirectColor scanline to runlength-encoded color packets.
        */
        for (y=image->rows-1; y >= 0; y--)
        {
          p=bmp_pixels+(image->rows-y-1)*bytes_per_line;
          q=image->pixels+(y*image->columns);
          for (x=0; x < image->columns; x++)
          {
            q->index=0;
            if (image->matte)
              q->index=UpScale(*p++);
            q->blue=UpScale(*p++);
            q->green=UpScale(*p++);
            q->red=UpScale(*p++);
            q->length=0;
            q++;
          }
        }
        break;
      }
      default:
        {
          fprintf(stderr, "vigraImpexReadBMPImage(): %s is not a BMP image file\n", image->filename);
          vigraImpexDestroyImage(image);
          return 0;
        }
    }
    free((char *) bmp_pixels);
    /*
      Proceed to next image.
    */
    if (image_info->subrange != 0)
      if (image->scene >= (image_info->subimage+image_info->subrange-1))
        break;
    status=vigraImpexReadData((char *) magick,1,2,image->file);
    if ((status == True) && (strncmp((char *) magick,"BM",2) == 0))
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
  } while ((status == True) && (strncmp((char *) magick,"BM",2) == 0));
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
+   W r i t e B M P I m a g e                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexWriteBMPImage writes an image in Microsoft Windows bitmap encoded
%  image format.
%
%  The format of the vigraImpexWriteBMPImage routine is:
%
%      status=vigraImpexWriteBMPImage(image_info,image)
%
%  A description of each parameter follows.
%
%    o status: Function vigraImpexWriteBMPImage return True if the image is written.
%      False is returned is there is a memory shortage or if the image file
%      fails to write.
%
%    o image_info: Specifies a pointer to an VigraImpexImageInfo structure.
%
%    o image:  A pointer to a VigraImpexImage structure.
%
%
*/
unsigned int vigraImpexWriteBMPImage(VigraImpexImageInfo *image_info,VigraImpexImage *image)
{
  BMPHeader
    bmp_header;

  register int
    i,
    j,
    x,
    y;

  register VigraImpexRunlengthPacket
    *p;

  register unsigned char
    *q;

  unsigned char
    *bmp_data,
    *bmp_pixels;

  unsigned int
    bytes_per_line,
    scene;

  /*
    Open output image file.
  */
  vigraImpexOpenImage(image_info,image,WriteBinaryType);
  if (image->file == (FILE *) NULL)
  {
      fprintf(stderr, "vigraImpexWriteBMPImage(): Unable to open file %s\n", image->filename);
      return False;
  }
  scene=0;
  do
  {
    /*
      Initialize BMP raster file header.
    */
    bmp_header.file_size=14+40;
    bmp_header.offset_bits=14+40;
    if ((strcmp(image_info->magick,"BMP24") == 0) ||
        (!vigraImpexIsPseudoClass(image) && !vigraImpexIsGrayImage(image)))
      {
        /*
          Full color BMP raster.
        */
        image->c_class=VigraImpexDirectClass;
        bmp_header.number_colors=0;
        bmp_header.bits_per_pixel=24;
        bytes_per_line=((image->columns*bmp_header.bits_per_pixel+31)/32)*4;
      }
    else
      {
        /*
          Colormapped BMP raster.
        */
        bmp_header.bits_per_pixel=8;
        bytes_per_line=((image->columns*bmp_header.bits_per_pixel+31)/32)*4;
        if (image_info->compression != VigraImpexNoCompression)
          bytes_per_line=image->columns;
        if (vigraImpexIsMonochromeImage(image))
          {
            bmp_header.bits_per_pixel=1;
            bytes_per_line=((image->columns*bmp_header.bits_per_pixel+31)/32)*4;
          }
        bmp_header.file_size+=4*(1 << bmp_header.bits_per_pixel);
        bmp_header.offset_bits+=4*(1 << bmp_header.bits_per_pixel);
        bmp_header.number_colors=1 << bmp_header.bits_per_pixel;
      }
    bmp_header.reserved[0]=0;
    bmp_header.reserved[1]=0;
    bmp_header.size=40;
    bmp_header.width=image->columns;
    bmp_header.height=image->rows;
    bmp_header.planes=1;
    bmp_header.compression=0;
    bmp_header.image_size=bytes_per_line*image->rows;
    bmp_header.file_size+=bmp_header.image_size;
    bmp_header.x_pixels=75*39;
    bmp_header.y_pixels=75*39;
    if (image->units == VigraImpexPixelsPerInchResolution)
      {
        bmp_header.x_pixels=(unsigned long) (100.0*image->x_resolution/2.54);
        bmp_header.y_pixels=(unsigned long) (100.0*image->y_resolution/2.54);
      }
    if (image->units == VigraImpexPixelsPerCentimeterResolution)
      {
        bmp_header.x_pixels=(unsigned long) (100.0*image->x_resolution);
        bmp_header.y_pixels=(unsigned long) (100.0*image->y_resolution);
      }
    bmp_header.colors_important=bmp_header.number_colors;
    /*
      Convert MIFF to BMP raster pixels.
    */
    bmp_pixels=(unsigned char *)
      malloc(bmp_header.image_size*sizeof(unsigned char));
    if (bmp_pixels == (unsigned char *) NULL)
      {
          fprintf(stderr, "vigraImpexWriteBMPImage(): Memory allocation failed\n");
          vigraImpexCloseImage(image);
          return False;
      }
    x=0;
    y=image->rows-1;
    switch (bmp_header.bits_per_pixel)
    {
      case 1:
      {
        register unsigned char
          bit,
          byte,
          polarity;

        /*
          Convert VigraImpexPseudoClass image to a BMP monochrome image.
        */
        p=image->pixels;
        polarity=0;
        if (image->colors == 2)
          polarity=
            Intensity(image->colormap[1]) > Intensity(image->colormap[0]);
        bit=0;
        byte=0;
        q=bmp_pixels+y*bytes_per_line;
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
                *q++=byte;
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
                  *q++=byte << (8-bit);
                bit=0;
                byte=0;
                x=0;
                y--;
                q=bmp_pixels+y*bytes_per_line;
             }
          }
          p++;
        }
        break;
      }
      case 8:
      {
        /*
          Convert VigraImpexPseudoClass packet to BMP pixel.
        */
        p=image->pixels;
        q=bmp_pixels+y*bytes_per_line;
        for (i=0; i < image->packets; i++)
        {
          for (j=0; j <= ((int) p->length); j++)
          {
            *q++=p->index;
            x++;
            if (x == image->columns)
              {
                x=0;
                y--;
                q=bmp_pixels+y*bytes_per_line;
              }
          }
          p++;
        }
        break;
      }
      case 24:
      {
        /*
          Convert VigraImpexDirectClass packet to BMP RGB pixel.
        */
        p=image->pixels;
        q=bmp_pixels+y*bytes_per_line;
        for (i=0; i < image->packets; i++)
        {
          for (j=0; j <= ((int) p->length); j++)
          {
            *q++=DownScale(p->blue);
            *q++=DownScale(p->green);
            *q++=DownScale(p->red);
            x++;
            if (x == image->columns)
              {
                x=0;
                y--;
                q=bmp_pixels+y*bytes_per_line;
              }
          }
          p++;
        }
        break;
      }
    }
    if (bmp_header.bits_per_pixel == 8)
      if (image_info->compression != VigraImpexNoCompression)
        {
          unsigned int
            packets;

          /*
            Convert run-length encoded raster pixels.
          */
          packets=(unsigned int)
            ((bytes_per_line*(bmp_header.height+2)+1) << 1);
          bmp_data=(unsigned char *) malloc(packets*sizeof(unsigned char));
          if (bmp_pixels == (unsigned char *) NULL)
          {
              fprintf(stderr, "vigraImpexWriteBMPImage(): Memory allocation failed\n");
              free((char *) bmp_pixels);
              vigraImpexCloseImage(image);
              return False;
          }
          bmp_header.image_size=
            vigraImpexBMPEncodeImage(bmp_pixels,bmp_data,image->columns,image->rows);
          free((char *) bmp_pixels);
          bmp_pixels=bmp_data;
          bmp_header.compression=1;
        }
    /*
      Write BMP header.
    */
    (void) fwrite("BM",1,2,image->file);
    vigraImpexLSBFirstWriteLong(bmp_header.file_size,image->file);
    vigraImpexLSBFirstWriteShort(bmp_header.reserved[0],image->file);
    vigraImpexLSBFirstWriteShort(bmp_header.reserved[1],image->file);
    vigraImpexLSBFirstWriteLong(bmp_header.offset_bits,image->file);
    vigraImpexLSBFirstWriteLong(bmp_header.size,image->file);
    vigraImpexLSBFirstWriteLong(bmp_header.width,image->file);
    vigraImpexLSBFirstWriteLong(bmp_header.height,image->file);
    vigraImpexLSBFirstWriteShort(bmp_header.planes,image->file);
    vigraImpexLSBFirstWriteShort(bmp_header.bits_per_pixel,image->file);
    vigraImpexLSBFirstWriteLong(bmp_header.compression,image->file);
    vigraImpexLSBFirstWriteLong(bmp_header.image_size,image->file);
    vigraImpexLSBFirstWriteLong(bmp_header.x_pixels,image->file);
    vigraImpexLSBFirstWriteLong(bmp_header.y_pixels,image->file);
    vigraImpexLSBFirstWriteLong(bmp_header.number_colors,image->file);
    vigraImpexLSBFirstWriteLong(bmp_header.colors_important,image->file);
    if (image->c_class == VigraImpexPseudoClass)
      {
        unsigned char
          *bmp_colormap;

        /*
          Dump colormap to file.
        */
        bmp_colormap=(unsigned char *)
          malloc(4*(1 << bmp_header.bits_per_pixel)*sizeof(unsigned char));
        if (bmp_colormap == (unsigned char *) NULL)
        {
          fprintf(stderr, "vigraImpexWriteBMPImage(): Memory allocation failed\n");
          vigraImpexCloseImage(image);
          return False;
        }
        q=bmp_colormap;
        for (i=0; i < image->colors; i++)
        {
          *q++=DownScale(image->colormap[i].blue);
          *q++=DownScale(image->colormap[i].green);
          *q++=DownScale(image->colormap[i].red);
          q++;
        }
        for ( ; i < (int) (1 << bmp_header.bits_per_pixel); i++)
        {
          *q++=(Quantum) 0x0;
          *q++=(Quantum) 0x0;
          *q++=(Quantum) 0x0;
          q++;
        }
        (void) fwrite((char *) bmp_colormap,4,1 << bmp_header.bits_per_pixel,
          image->file);
        free((char *) bmp_colormap);
      }
    (void) fwrite((char *) bmp_pixels,1,(int) bmp_header.image_size,
      image->file);
    free((char *) bmp_pixels);
    if (image->next == (VigraImpexImage *) NULL)
      break;
    image->next->file=image->file;
    image=image->next;
  } while (image_info->adjoin);
  vigraImpexCloseImage(image);
  return(True);
}
