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
#include <setjmp.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "jpeg.h"
#include "utility.h"

#if defined(HasJPEG)

#include "jpeglib.h"

static VigraImpexImage
  *image;

static jmp_buf
  error_recovery;
  
#define PrematureExit(message, image) \
{ \
    fprintf(stderr, "JPEG library: %s %s\n", (char *) message, image->filename); \
    vigraImpexDestroyImages(image); \
    return 0; \
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   R e a d J P E G I m a g e                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexReadJPEGImage reads a JPEG image file and returns it.  It allocates
%  the memory necessary for the new Image structure and returns a pointer to
%  the new image.
%
%  The format of the vigraImpexReadJPEGImage routine is:
%
%      image=vigraImpexReadJPEGImage(image_info)
%
%  A description of each parameter follows:
%
%    o image:  Function vigraImpexReadJPEGImage returns a pointer to the image after
%      reading.  A null image is returned if there is a a memory shortage or
%      if the image cannot be read.
%
%    o filename:  Specifies the name of the jpeg image to read.
%
%
*/

static unsigned int GetCharacter(j_decompress_ptr jpeg_info)
{
  struct jpeg_source_mgr
    *data;

  data=jpeg_info->src;
  if (data->bytes_in_buffer == 0)
    (*data->fill_input_buffer) (jpeg_info);
  data->bytes_in_buffer--;
  return(GETJOCTET(*data->next_input_byte++));
}

static boolean CommentHandler(j_decompress_ptr jpeg_info)
{
  long int
    length;

  register char
    *p;

  /*
    Determine length of comment.
  */
  length=GetCharacter(jpeg_info) << 8;
  length+=GetCharacter(jpeg_info);
  length-=2;
  if (image->comments != (char *) NULL)
    image->comments=(char *) realloc((char *) image->comments,
      (unsigned int) (Extent(image->comments)+length+1)*sizeof(char));
  else
    {
      image->comments=(char *)
        malloc((unsigned int) (length+1)*sizeof(char));
      if (image->comments != (char *) NULL)
        *image->comments='\0';
    }
  if (image->comments == (char *) NULL)
    {
      fprintf(stderr, "vigraImpexReadJPEGImage(): Memory allocation failed\n");
      return(0);
    }
  /*
    Read comment.
  */
  p=image->comments+Extent(image->comments);
  while (--length >= 0)
    *p++=GetCharacter(jpeg_info);
  *p='\0';
  return(True);
}

static void EmitMessage(j_common_ptr jpeg_info,int level)
{
  char
    message[JMSG_LENGTH_MAX];

  struct jpeg_error_mgr
    *jpeg_error;

  jpeg_error=jpeg_info->err;
  (jpeg_error->format_message) (jpeg_info,message);
  if (level < 0)
    {
      if ((jpeg_error->num_warnings == 0) || (jpeg_error->trace_level >= 3))
        fprintf(stderr, "vigraImpexReadJPEGImage(): %s %s\n", (char *) message,image->filename);
      jpeg_error->num_warnings++;
    }
  else
    if (jpeg_error->trace_level >= level)
        fprintf(stderr, "vigraImpexReadJPEGImage(): %s %s\n", (char *) message,image->filename);
}

static void ErrorExit(j_common_ptr jpeg_info)
{
  EmitMessage(jpeg_info,0);
  longjmp(error_recovery,1);
}

Export VigraImpexImage *vigraImpexReadJPEGImage( VigraImpexImageInfo *image_info)
{
  int
    x,
    y;

  JSAMPLE
    *jpeg_pixels;

  JSAMPROW
    scanline[1];

  Quantum
    blue,
    green,
    red;

  register int
    i;

  register JSAMPLE
    *p;

  register VigraImpexRunlengthPacket
    *q;

  struct jpeg_decompress_struct
    jpeg_info;

  struct jpeg_error_mgr
    jpeg_error;

  unsigned int
    black,
    cyan,
    magenta,
    packets,
    yellow;

  unsigned short
    index;

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
    PrematureExit("Unable to open file",image);
  /*
    Initialize image structure.
  */
  jpeg_info.err=jpeg_std_error(&jpeg_error);
  jpeg_info.err->emit_message=EmitMessage;
  jpeg_info.err->error_exit=ErrorExit;
  image->pixels=(VigraImpexRunlengthPacket *) NULL;
  jpeg_pixels=(JSAMPLE *) NULL;
  if (setjmp(error_recovery))
    {
      /*
        JPEG image is corrupt.
      */
      if (jpeg_pixels != (JSAMPLE *) NULL)
        free((char *) jpeg_pixels);
      jpeg_destroy_decompress(&jpeg_info);
      vigraImpexDestroyImage(image);
      return((VigraImpexImage *) NULL);
    }
  jpeg_create_decompress(&jpeg_info);
  jpeg_set_marker_processor(&jpeg_info,JPEG_COM,CommentHandler);
  jpeg_stdio_src(&jpeg_info,image->file);
  (void) jpeg_read_header(&jpeg_info,True);
  if (jpeg_info.saw_JFIF_marker)
    {
      /*
        Set image resolution.
      */
      image->x_resolution=jpeg_info.X_density;
      image->y_resolution=jpeg_info.Y_density;
      if (jpeg_info.density_unit == 1)
        image->units=VigraImpexPixelsPerInchResolution;
      if (jpeg_info.density_unit == 2)
        image->units=VigraImpexPixelsPerCentimeterResolution;
    }
  if ((image_info->subrange != 0) || (image_info->size != (char *) NULL))
    {
      /*
        Let the JPEG library subsample for us.
      */
      jpeg_calc_output_dimensions(&jpeg_info);
      image->magick_columns=jpeg_info.output_width;
      image->magick_rows=jpeg_info.output_height;
      jpeg_info.scale_denom=image_info->subimage;
      if (image_info->size != (char *) NULL)
        {
          unsigned int
            height,
            width;

          unsigned long
            scale_factor;

          width=jpeg_info.output_width;
          height=jpeg_info.output_height;
          x=0;
          y=0;
          if (width == 0)
            width=1;
          scale_factor=UpShift(jpeg_info.output_width)/width;
          if (height == 0)
            height=1;
          if (scale_factor > (UpShift(jpeg_info.output_height)/height))
            scale_factor=UpShift(jpeg_info.output_height)/height;
          jpeg_info.scale_denom=DownShift(scale_factor);
        }
      jpeg_calc_output_dimensions(&jpeg_info);
    }
#if (JPEG_LIB_VERSION >= 61)
  jpeg_info.dct_method=JDCT_FLOAT;
  image->interlace=jpeg_info.progressive_mode ? VigraImpexPlaneInterlace : VigraImpexNoInterlace;
#endif
  jpeg_start_decompress(&jpeg_info);
  image->columns=jpeg_info.output_width;
  image->rows=jpeg_info.output_height;
  if (image_info->ping)
    {
      jpeg_destroy_decompress(&jpeg_info);
      vigraImpexCloseImage(image);
      return(image);
    }
  image->packets=0;
  packets=Max((image->columns*image->rows+2) >> 2,1);
  image->pixels=(VigraImpexRunlengthPacket *) malloc(packets*sizeof(VigraImpexRunlengthPacket));
  jpeg_pixels=(JSAMPLE *)
    malloc(jpeg_info.output_components*image->columns*sizeof(JSAMPLE));
  if ((image->pixels == (VigraImpexRunlengthPacket *) NULL) ||
      (jpeg_pixels == (JSAMPLE *) NULL))
    PrematureExit("Memory allocation failed",image);
  if (jpeg_info.out_color_space == JCS_GRAYSCALE)
    {
      /*
        Initialize grayscale colormap.
      */
      image->c_class=VigraImpexPseudoClass;
      image->colors=MaxRGB+1;
      image->colormap=(VigraImpexColorPacket *) malloc(image->colors*sizeof(VigraImpexColorPacket));
      if (image->colormap == (VigraImpexColorPacket *) NULL)
        PrematureExit("Memory allocation failed",image);
      for (i=0; i < image->colors; i++)
      {
        image->colormap[i].red=UpScale(i);
        image->colormap[i].green=UpScale(i);
        image->colormap[i].blue=UpScale(i);
      }
    }
  /*
    Convert JPEG pixels to runlength-encoded packets.
  */
  red=0;
  green=0;
  blue=0;
  index=0;
  scanline[0]=(JSAMPROW) jpeg_pixels;
  q=image->pixels;
  q->length=MaxRunlength;
  for (y=0; y < image->rows; y++)
  {
    (void) jpeg_read_scanlines(&jpeg_info,scanline,1);
    p=jpeg_pixels;
    for (x=0; x < image->columns; x++)
    {
      if (jpeg_info.data_precision > QuantumDepth)
        {
          if (jpeg_info.out_color_space == JCS_GRAYSCALE)
            index=GETJSAMPLE(*p++) >> 4;
          else
            {
              red=(Quantum) (GETJSAMPLE(*p++) >> 4);
              green=(Quantum) (GETJSAMPLE(*p++) >> 4);
              blue=(Quantum) (GETJSAMPLE(*p++) >> 4);
              if (jpeg_info.out_color_space == JCS_CMYK)
                index=(Quantum) (GETJSAMPLE(*p++) >> 4);
            }
         }
       else
         if (jpeg_info.out_color_space == JCS_GRAYSCALE)
           index=GETJSAMPLE(*p++);
         else
           {
             red=(Quantum) UpScale(GETJSAMPLE(*p++));
             green=(Quantum) UpScale(GETJSAMPLE(*p++));
             blue=(Quantum) UpScale(GETJSAMPLE(*p++));
             if (jpeg_info.out_color_space == JCS_CMYK)
               index=(Quantum) UpScale(GETJSAMPLE(*p++));
           }
      if (jpeg_info.out_color_space == JCS_CMYK)
        {
          cyan=red;
          magenta=green;
          yellow=blue;
          black=index;
          if ((cyan+black) > MaxRGB)
            red=0;
          else
            red=MaxRGB-(cyan+black);
          if ((magenta+black) > MaxRGB)
            green=0;
          else
            green=MaxRGB-(magenta+black);
          if ((yellow+black) > MaxRGB)
            blue=0;
          else
            blue=MaxRGB-(yellow+black);
          index=0;
        }
      if ((red == q->red) && (green == q->green) && (blue == q->blue) &&
          (index == q->index) && ((int) q->length < MaxRunlength))
        q->length++;
      else
        {
          if (image->packets != 0)
            q++;
          image->packets++;
          if (image->packets == packets)
            {
              packets<<=1;
              image->pixels=(VigraImpexRunlengthPacket *)
                realloc((char *) image->pixels,packets*sizeof(VigraImpexRunlengthPacket));
              if (image->pixels == (VigraImpexRunlengthPacket *) NULL)
                {
                  free((char *) jpeg_pixels);
                  jpeg_destroy_decompress(&jpeg_info);
                  PrematureExit("Memory allocation failed",image);
                }
              q=image->pixels+image->packets-1;
            }
          q->red=red;
          q->green=green;
          q->blue=blue;
          q->index=index;
          q->length=0;
        }
    }
  }
  (void) jpeg_finish_decompress(&jpeg_info);
  image->pixels=(VigraImpexRunlengthPacket *)
    realloc((char *) image->pixels,image->packets*sizeof(VigraImpexRunlengthPacket));
  if (image->c_class == VigraImpexPseudoClass)
    vigraImpexSyncImage(image);
  /*
    Free memory.
  */
  jpeg_destroy_decompress(&jpeg_info);
  free((char *) jpeg_pixels);
  vigraImpexCloseImage(image);
  return(image);
}
#undef PrematureExit
#define PrematureExit(message, image) \
{ \
    fprintf(stderr, "JPEG library: %s %s\n", (char *) message, image->filename); \
    if(image->file) vigraImpexCloseImage(image); \
    return 0; \
}
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+  W r i t e J P E G I m a g e                                                %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexWriteJPEGImage writes a JPEG image file and returns it.  It
%  allocates the memory necessary for the new Image structure and returns a
%  pointer to the new image.
%
%  The format of the vigraImpexWriteJPEGImage routine is:
%
%      status=vigraImpexWriteJPEGImage(image_info,image)
%
%  A description of each parameter follows:
%
%    o status:  Function vigraImpexWriteJPEGImage return True if the image is written.
%      False is returned is there is of a memory shortage or if the image
%      file cannot be opened for writing.
%
%    o image_info: Specifies a pointer to an VigraImpexImageInfo structure.
%
%    o jpeg_image:  A pointer to a Image structure.
%
%
*/

Export unsigned int vigraImpexWriteJPEGImage( VigraImpexImageInfo *image_info,VigraImpexImage *image)
{
  float
     black_generation,
     undercolor;

   int
     black,
     cyan,
     magenta,
     yellow;

  JSAMPLE
    *jpeg_pixels;

  JSAMPROW
    scanline[1];

  register int
    i,
    j,
    x;

  register JSAMPLE
    *q;

  register VigraImpexRunlengthPacket
    *p;

  struct jpeg_compress_struct
    jpeg_info;

  struct jpeg_error_mgr
    jpeg_error;

  unsigned int
    packets;

  /*
    Open image file.
  */
  vigraImpexOpenImage(image_info,image,WriteBinaryType);
  if (image->file == (FILE *) NULL)
    PrematureExit("Unable to open file",image);
  /*
    Initialize JPEG parameters.
  */
  jpeg_info.err=jpeg_std_error(&jpeg_error);
  jpeg_info.err->emit_message=EmitMessage;
  jpeg_create_compress(&jpeg_info);
  jpeg_stdio_dest(&jpeg_info,image->file);
  jpeg_info.image_width=image->columns;
  jpeg_info.image_height=image->rows;
  jpeg_info.input_components=3;
  jpeg_info.in_color_space=JCS_RGB;
  if (vigraImpexIsGrayImage(image))
    {
      jpeg_info.input_components=1;
      jpeg_info.in_color_space=JCS_GRAYSCALE;
    }
  if (image_info->colorspace == VigraImpexCMYKColorspace)
    {
      jpeg_info.input_components=4;
      jpeg_info.in_color_space=JCS_CMYK;
      undercolor=1.0;
      black_generation=1.0;
      if (image_info->undercolor != (char *) NULL)
        {
          (void) sscanf(image_info->undercolor,"%fx%f",&undercolor,
            &black_generation);
          if (black_generation == 1.0)
            black_generation=undercolor;
        }
    }
  jpeg_set_defaults(&jpeg_info);
  jpeg_info.density_unit=0;
  jpeg_info.X_density=(short) image->x_resolution;
  jpeg_info.Y_density=(short) image->y_resolution;
  if (image->units == VigraImpexPixelsPerInchResolution)
    jpeg_info.density_unit=1;
  if (image->units == VigraImpexPixelsPerCentimeterResolution)
    jpeg_info.density_unit=2;
  for (i=0; i < MAX_COMPONENTS; i++)
  {
    jpeg_info.comp_info[i].h_samp_factor=1;
    jpeg_info.comp_info[i].v_samp_factor=1;
  }
  jpeg_set_quality(&jpeg_info,image_info->quality,True);
  jpeg_info.optimize_coding=True;
#if (JPEG_LIB_VERSION >= 61)
  jpeg_info.dct_method=JDCT_FLOAT;
  if (image_info->interlace != VigraImpexNoInterlace)
    jpeg_simple_progression(&jpeg_info);
#endif
  jpeg_start_compress(&jpeg_info,True);
  if (image->comments != (char *) NULL)
    for (i=0; i < Extent(image->comments); i+=65533)
      jpeg_write_marker(&jpeg_info,JPEG_COM,(unsigned char *) image->comments+i,
        (unsigned int) Min(Extent(image->comments+i),65533));
  /*
    Convert MIFF to JPEG raster pixels.
  */
  packets=jpeg_info.input_components*image->columns;
  jpeg_pixels=(JSAMPLE *) malloc(packets*sizeof(JSAMPLE));
  if (jpeg_pixels == (JSAMPLE *) NULL)
    PrematureExit("Memory allocation failed",image);
  p=image->pixels;
  q=jpeg_pixels;
  x=0;
  scanline[0]=(JSAMPROW) jpeg_pixels;
  if ((jpeg_info.data_precision > 8) && (QuantumDepth > 8))
    {
      if (jpeg_info.in_color_space == JCS_GRAYSCALE)
        for (i=0; i < image->packets; i++)
        {
          for (j=0; j <= ((int) p->length); j++)
          {
            *q++=(JSAMPLE) (Intensity(*p) >> 4);
            x++;
            if (x == image->columns)
              {
                (void) jpeg_write_scanlines(&jpeg_info,scanline,1);
                q=jpeg_pixels;
                x=0;
              }
          }
          p++;
        }
      else
        if (jpeg_info.in_color_space == JCS_RGB)
          for (i=0; i < image->packets; i++)
          {
            for (j=0; j <= ((int) p->length); j++)
            {
              *q++=(JSAMPLE) (p->red >> 4);
              *q++=(JSAMPLE) (p->green >> 4);
              *q++=(JSAMPLE) (p->blue >> 4);
              x++;
              if (x == image->columns)
                {
                  (void) jpeg_write_scanlines(&jpeg_info,scanline,1);
                  q=jpeg_pixels;
                  x=0;
                }
            }
            p++;
          }
        else
          for (i=0; i < image->packets; i++)
          {
            cyan=MaxRGB-p->red;
            magenta=MaxRGB-p->green;
            yellow=MaxRGB-p->blue;
            black=cyan;
            if (magenta < black)
              black=magenta;
            if (yellow < black)
              black=yellow;
            for (j=0; j <= ((int) p->length); j++)
            {
              /*
                Convert VigraImpexDirectClass packets to contiguous RGB scanlines.
              */
              *q++=(JSAMPLE) (cyan-undercolor*black) >> 4;
              *q++=(JSAMPLE) (magenta-undercolor*black) >> 4;
              *q++=(JSAMPLE) (yellow-undercolor*black) >> 4;
              *q++=(JSAMPLE) (black_generation*black) >> 4;
              x++;
              if (x == image->columns)
                {
                  (void) jpeg_write_scanlines(&jpeg_info,scanline,1);
                  q=jpeg_pixels;
                  x=0;
                }
            }
            p++;
          }
    }
  else
    if (jpeg_info.in_color_space == JCS_GRAYSCALE)
      for (i=0; i < image->packets; i++)
      {
        for (j=0; j <= ((int) p->length); j++)
        {
          *q++=(JSAMPLE) DownScale(Intensity(*p));
          x++;
          if (x == image->columns)
            {
              (void) jpeg_write_scanlines(&jpeg_info,scanline,1);
              q=jpeg_pixels;
              x=0;
            }
        }
        p++;
      }
    else
      if (jpeg_info.in_color_space == JCS_RGB)
        for (i=0; i < image->packets; i++)
        {
          for (j=0; j <= ((int) p->length); j++)
          {
            *q++=(JSAMPLE) DownScale(p->red);
            *q++=(JSAMPLE) DownScale(p->green);
            *q++=(JSAMPLE) DownScale(p->blue);
            x++;
            if (x == image->columns)
              {
                (void) jpeg_write_scanlines(&jpeg_info,scanline,1);
                q=jpeg_pixels;
                x=0;
              }
          }
          p++;
        }
      else
        for (i=0; i < image->packets; i++)
        {
          cyan=MaxRGB-p->red;
          magenta=MaxRGB-p->green;
          yellow=MaxRGB-p->blue;
          black=cyan;
          if (magenta < black)
            black=magenta;
          if (yellow < black)
            black=yellow;
          for (j=0; j <= ((int) p->length); j++)
          {
            /*
              Convert VigraImpexDirectClass packets to contiguous RGB scanlines.
            */
            *q++=(JSAMPLE) DownScale(cyan-undercolor*black);
            *q++=(JSAMPLE) DownScale(magenta-undercolor*black);
            *q++=(JSAMPLE) DownScale(yellow-undercolor*black);
            *q++=(JSAMPLE) DownScale(black_generation*black);
            x++;
            if (x == image->columns)
              {
                (void) jpeg_write_scanlines(&jpeg_info,scanline,1);
                q=jpeg_pixels;
                x=0;
              }
          }
          p++;
        }
  jpeg_finish_compress(&jpeg_info);
  /*
    Free memory.
  */
  jpeg_destroy_compress(&jpeg_info);
  free((char *) jpeg_pixels);
  vigraImpexCloseImage(image);
  return(True);
}
#else
Export VigraImpexImage *vigraImpexReadJPEGImage( VigraImpexImageInfo *image_info)
{
  fprintf(stderr, "JPEG library is not available");
  return 0;
}
Export unsigned int vigraImpexWriteJPEGImage( VigraImpexImageInfo *image_info,VigraImpexImage *image)
{
  fprintf(stderr, "JPEG library is not available");
  return 0;
}
#endif
