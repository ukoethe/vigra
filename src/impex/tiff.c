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
#include <setjmp.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "utility.h"
#include "vtiff.h"

#if defined(HasTIFF)

#include "tiffio.h"

#define PrematureExit(message, image) \
{ \
    fprintf(stderr, "TIFF library: %s %s\n", (char *) message, image->filename); \
    vigraImpexDestroyImage(image); \
    return 0; \
}
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   R e a d T I F F I m a g e                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function ReadTIFFImage reads a Tagged image file and returns it.  It
%  allocates the memory necessary for the new VigraImpexImage structure and returns a
%  pointer to the new image.
%
%  The format of the ReadTIFFImage routine is:
%
%      image=ReadTIFFImage(image_info)
%
%  A description of each parameter follows:
%
%    o image:  Function ReadTIFFImage returns a pointer to the image after
%      reading.  A null image is returned if there is a a memory shortage or
%      if the image cannot be read.
%
%    o image_info: Specifies a pointer to an VigraImpexImageInfo structure.
%
%
*/

static void TIFFWarningMessage(const char *module,const char *format,
  va_list warning)
{
  char
    message[MaxTextExtent];

  register char
    *p;

  p=message;
  if (module != (char *) NULL)
    {
      (void) sprintf(p,"%s: ",module);
      p+=Extent(message);
    }
  (void) vsprintf(p,format,warning);
  (void) strcat(p,".");
  fprintf(stderr, "TIFF library: %s\n",message);
}

Export VigraImpexImage *vigraImpexReadTIFFImage( VigraImpexImageInfo *image_info)
{
  char
    *text;

  VigraImpexImage
    *image;

  int
    range;

  Quantum
    blue,
    green,
    red;

  register int
    i,
    x,
    y;

  register VigraImpexRunlengthPacket
    *q;

  register unsigned char
    *p;

  TIFF
    *tiff;

  uint16
    extra_samples,
    *sample_info;

  unsigned char
    *scanline;

  unsigned int
    height,
    method,
    packets,
    status,
    width;

  unsigned short
    bits_per_sample,
    index,
    interlace,
    max_sample_value,
    min_sample_value,
    pages,
    photometric,
    samples_per_pixel,
    units,
    value;

  /*
    Allocate image structure.
  */
  image=vigraImpexAllocateImage(image_info);
  if (image == (VigraImpexImage *) NULL)
    return((VigraImpexImage *) NULL);
  /*
    Open image.
  */
  vigraImpexOpenImage(image_info,image,ReadBinaryType);
  if (image->file == (FILE *) NULL)
    PrematureExit("Unable to open file",image);
  if ((image->file == stdin) || image->pipe)
    {
      FILE
        *file;

      int
        c;

      /*
        Copy standard input or pipe to temporary file.
      */
      TemporaryFilename(image_info->filename);
      file=fopen(image_info->filename,WriteBinaryType);
      if (file == (FILE *) NULL)
        PrematureExit("Unable to write file",image);
      c=fgetc(image->file);
      while (c != EOF)
      {
        (void) putc(c,file);
        c=fgetc(image->file);
      }
      (void) fclose(file);
      (void) strcpy(image->filename,image_info->filename);
      image->temporary=True;
    }
  vigraImpexCloseImage(image);
  TIFFSetErrorHandler(TIFFWarningMessage);
  TIFFSetWarningHandler(TIFFWarningMessage);
  tiff=TIFFOpen(image->filename,ReadBinaryType);
  if (tiff == (TIFF *) NULL)
    PrematureExit("Unable to open file",image);
  if (image_info->subrange != 0)
    while (image->scene < image_info->subimage)
    {
      /*
        Skip to next image.
      */
      image->scene++;
      status=TIFFReadDirectory(tiff);
      if (status == False)
        PrematureExit("Unable to read subimage",image);
    }
  do
  {
    if (image_info->verbose)
      TIFFPrintDirectory(tiff,stderr,False);
    TIFFGetField(tiff,TIFFTAG_IMAGEWIDTH,&width);
    TIFFGetField(tiff,TIFFTAG_IMAGELENGTH,&height);
    TIFFGetFieldDefaulted(tiff,TIFFTAG_PLANARCONFIG,&interlace);
    TIFFGetFieldDefaulted(tiff,TIFFTAG_BITSPERSAMPLE,&bits_per_sample);
    TIFFGetFieldDefaulted(tiff,TIFFTAG_MINSAMPLEVALUE,&min_sample_value);
    TIFFGetFieldDefaulted(tiff,TIFFTAG_MAXSAMPLEVALUE,&max_sample_value);
    TIFFGetFieldDefaulted(tiff,TIFFTAG_PHOTOMETRIC,&photometric);
    TIFFGetFieldDefaulted(tiff,TIFFTAG_SAMPLESPERPIXEL,&samples_per_pixel);
    TIFFGetFieldDefaulted(tiff,TIFFTAG_RESOLUTIONUNIT,&units);
    TIFFGetFieldDefaulted(tiff,TIFFTAG_XRESOLUTION,&image->x_resolution);
    TIFFGetFieldDefaulted(tiff,TIFFTAG_YRESOLUTION,&image->y_resolution);
    /*
      Allocate memory for the image and pixel buffer.
    */
    image->columns=width;
    image->rows=height;
    range=max_sample_value-min_sample_value;
    if ((samples_per_pixel == 1) && !TIFFIsTiled(tiff))
      {
        image->c_class=VigraImpexPseudoClass;
        image->colors=range+1;
        if (bits_per_sample > QuantumDepth)
          image->colors=MaxRGB+1;
      }
    if (image_info->ping)
      {
        TIFFClose(tiff);
        vigraImpexCloseImage(image);
        return(image);
      }
    if (units == RESUNIT_INCH)
      image->units=VigraImpexPixelsPerInchResolution;
    if (units == RESUNIT_CENTIMETER)
      image->units=VigraImpexPixelsPerCentimeterResolution;
    image->depth=bits_per_sample;
    if (bits_per_sample < 8)
      image->depth=8;
    image->packets=0;
    packets=Max((image->columns*image->rows+4) >> 3,1);
    if (bits_per_sample == 1)
      packets=Max((image->columns*image->rows+8) >> 4,1);
    image->pixels=(VigraImpexRunlengthPacket *) malloc(packets*sizeof(VigraImpexRunlengthPacket));
    if (image->pixels == (VigraImpexRunlengthPacket *) NULL)
      {
        TIFFClose(tiff);
        PrematureExit("Memory allocation failed",image);
      }
    TIFFGetFieldDefaulted(tiff,TIFFTAG_PAGENUMBER,&value,&pages);
    image->scene=value;
    text=(char *) NULL;
    TIFFGetField(tiff,TIFFTAG_PAGENAME,&text);
    if (text != (char *) NULL)
      {
        image->label=(char *)
          malloc((unsigned int) (Extent(text)+1)*sizeof(char));
        if (image->label == (char *) NULL)
          {
            TIFFClose(tiff);
            PrematureExit("Memory allocation failed",
              image);
          }
        (void) strcpy(image->label,text);
      }
    text=(char *) NULL;
    TIFFGetField(tiff,TIFFTAG_IMAGEDESCRIPTION,&text);
    if (text != (char *) NULL)
      {
        image->comments=(char *)
          malloc((unsigned int) (Extent(text)+1)*sizeof(char));
        if (image->comments == (char *) NULL)
          {
            TIFFClose(tiff);
            PrematureExit("Memory allocation failed",
              image);
          }
        (void) strcpy(image->comments,text);
      }
    if (range < 0)
      range=max_sample_value;
    q=image->pixels;
    q->length=MaxRunlength;
    method=0;
    if ((samples_per_pixel > 1) || TIFFIsTiled(tiff))
      {
        method=2;
        if ((samples_per_pixel >= 3) && (photometric == PHOTOMETRIC_RGB) &&
            (interlace == PLANARCONFIG_CONTIG))
          method=1;
      }
    switch (method)
    {
      case 0:
      {
        Quantum
          *quantum_scanline;

        register Quantum
          *r;

        /*
          Convert TIFF image to VigraImpexPseudoClass MIFF image.
        */
        image->colormap=(VigraImpexColorPacket *)
          malloc(image->colors*sizeof(VigraImpexColorPacket));
        quantum_scanline=(Quantum *) malloc(width*sizeof(Quantum));
        scanline=(unsigned char *) malloc(TIFFScanlineSize(tiff)+4);
        if ((image->colormap == (VigraImpexColorPacket *) NULL) ||
            (quantum_scanline == (Quantum *) NULL) ||
            (scanline == (unsigned char *) NULL))
          {
            TIFFClose(tiff);
            PrematureExit("Memory allocation failed",
              image);
          }
        /*
          Create colormap.
        */
        switch (photometric)
        {
          case PHOTOMETRIC_MINISBLACK:
          {
            for (i=0; i < image->colors; i++)
            {
              image->colormap[i].red=(MaxRGB*i)/(image->colors-1);
              image->colormap[i].green=(MaxRGB*i)/(image->colors-1);
              image->colormap[i].blue=(MaxRGB*i)/(image->colors-1);
            }
            break;
          }
          case PHOTOMETRIC_MINISWHITE:
          {
            unsigned int
              colors;

            colors=image->colors;
            for (i=0; i < image->colors; i++)
            {
              image->colormap[colors-i-1].red=(MaxRGB*i)/(image->colors-1);
              image->colormap[colors-i-1].green=(MaxRGB*i)/(image->colors-1);
              image->colormap[colors-i-1].blue=(MaxRGB*i)/(image->colors-1);
            }
            break;
          }
          case PHOTOMETRIC_PALETTE:
          {
            long
              range;

            unsigned short
              *blue_colormap,
              *green_colormap,
              *red_colormap;

            TIFFGetField(tiff,TIFFTAG_COLORMAP,&red_colormap,&green_colormap,
              &blue_colormap);
            range=256L;  /* might be old style 8-bit colormap */
            for (i=0; i < image->colors; i++)
              if ((red_colormap[i] >= 256) || (green_colormap[i] >= 256) ||
                  (blue_colormap[i] >= 256))
                {
                  range=65535L;
                  break;
                }
            for (i=0; i < image->colors; i++)
            {
              image->colormap[i].red=(Quantum)
                ((long) (MaxRGB*red_colormap[i])/range);
              image->colormap[i].green=(Quantum)
                ((long) (MaxRGB*green_colormap[i])/range);
              image->colormap[i].blue=(Quantum)
                ((long) (MaxRGB*blue_colormap[i])/range);
            }
            break;
          }
          default:
            break;
        }
        /*
          Convert image to VigraImpexPseudoClass runlength-encoded packets.
        */
        for (y=0; y < image->rows; y++)
        {
          TIFFReadScanline(tiff,(char *) scanline,y,0);
          p=scanline;
          r=quantum_scanline;
          switch (bits_per_sample)
          {
            case 1:
            {
              register int
                bit;

              for (x=0; x < ((int) width-7); x+=8)
              {
                for (bit=7; bit >= 0; bit--)
                  *r++=((*p) & (0x01 << bit) ? 0x01 : 0x00);
                p++;
              }
              if ((width % 8) != 0)
                {
                  for (bit=7; bit >= (8-(width % 8)); bit--)
                    *r++=((*p) & (0x01 << bit) ? 0x01 : 0x00);
                  p++;
                }
              break;
            }
            case 2:
            {
              for (x=0; x < ((int) width-3); x+=4)
              {
                *r++=(*p >> 6) & 0x3;
                *r++=(*p >> 4) & 0x3;
                *r++=(*p >> 2) & 0x3;
                *r++=(*p) & 0x3;
                p++;
              }
              if ((width % 4) != 0)
                {
                  for (i=3; i >= (4-(width % 4)); i--)
                    *r++=(*p >> (i*2)) & 0x03;
                  p++;
                }
              break;
            }
            case 4:
            {
              for (x=0; x < ((int) width-1); x+=2)
              {
                *r++=(*p >> 4) & 0xf;
                *r++=(*p) & 0xf;
                p++;
              }
              if ((width % 2) != 0)
                *r++=(*p++ >> 4) & 0xf;
              break;
            }
            case 8:
            {
              for (x=0; x < width; x++)
                *r++=(*p++);
              break;
            }
            case 16:
            {
              for (x=0; x < image->columns; x++)
              {
                ReadQuantum(*r,p);
                r++;
              }
              break;
            }
            default:
              break;
          }
          /*
            Transfer image scanline.
          */
          r=quantum_scanline;
          for (x=0; x < image->columns; x++)
          {
            index=(*r++);
            if ((index == q->index) && ((int) q->length < MaxRunlength))
              q->length++;
            else
              {
                if (image->packets != 0)
                  q++;
                image->packets++;
                if (image->packets == packets)
                  {
                    packets<<=1;
                    image->pixels=(VigraImpexRunlengthPacket *) realloc((char *)
                      image->pixels,packets*sizeof(VigraImpexRunlengthPacket));
                    if (image->pixels == (VigraImpexRunlengthPacket *) NULL)
                      {
                        free((char *) scanline);
                        free((char *) quantum_scanline);
                        TIFFClose(tiff);
                        PrematureExit("Memory allocation failed",image);
                      }
                    q=image->pixels+image->packets-1;
                  }
                q->index=index;
                q->length=0;
              }
          }
        }
        free((char *) scanline);
        free((char *) quantum_scanline);
        if (image->c_class == VigraImpexPseudoClass)
          vigraImpexSyncImage(image);
        break;
      }
      case 1:
      {
        /*
          Convert TIFF image to VigraImpexDirectClass MIFF image.
        */
        scanline=(unsigned char *) malloc((TIFFScanlineSize(tiff) << 1)+4);
        if (scanline == (unsigned char *) NULL)
          {
            TIFFClose(tiff);
            PrematureExit("Memory allocation failed",
              image);
          }
        TIFFGetFieldDefaulted(tiff,TIFFTAG_EXTRASAMPLES,&extra_samples,
          &sample_info);
        image->matte=extra_samples == 1;
        for (y=0; y < image->rows; y++)
        {
          TIFFReadScanline(tiff,(char *) scanline,y,0);
          if (bits_per_sample == 4)
            {
              register unsigned char
                *r;

              width=TIFFScanlineSize(tiff);
              p=scanline+width-1;
              r=scanline+(width << 1)-1;
              for (x=0; x < (int) width; x++)
              {
                *r--=((*p) & 0xf) << 4;
                *r--=((*p >> 4) & 0xf) << 4;
                p--;
              }
            }
          p=scanline;
          for (x=0; x < image->columns; x++)
          {
            ReadQuantum(red,p);
            ReadQuantum(green,p);
            ReadQuantum(blue,p);
            index=0;
            if (samples_per_pixel == 4)
              ReadQuantum(index,p);
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
                    image->pixels=(VigraImpexRunlengthPacket *) realloc((char *)
                      image->pixels,packets*sizeof(VigraImpexRunlengthPacket));
                    if (image->pixels == (VigraImpexRunlengthPacket *) NULL)
                      {
                        TIFFClose(tiff);
                        free((char *) scanline);
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
        free((char *) scanline);
        break;
      }
      case 2:
      default:
      {
        register uint32
          *p,
          *pixels;

        /*
          Convert TIFF image to VigraImpexDirectClass MIFF image.
        */
        TIFFGetFieldDefaulted(tiff,TIFFTAG_EXTRASAMPLES,&extra_samples,
          &sample_info);
        image->matte=
          ((extra_samples == 1) && (sample_info[0] == EXTRASAMPLE_ASSOCALPHA));
        pixels=(uint32 *)
          malloc((image->columns*image->rows+image->columns)*sizeof(uint32));
        if (pixels == (uint32 *) NULL)
          {
            TIFFClose(tiff);
            PrematureExit("Memory allocation failed",
              image);
          }
        status=TIFFReadRGBAImage(tiff,image->columns,image->rows,pixels,0);
        if (status == False)
          {
            free((char *) pixels);
            TIFFClose(tiff);
            PrematureExit("Unable to read image",image);
          }
        /*
          Convert image to VigraImpexDirectClass runlength-encoded packets.
        */
        for (y=image->rows-1; y >= 0; y--)
        {
          p=pixels+y*image->columns;
          for (x=0; x < image->columns; x++)
          {
            red=UpScale(TIFFGetR(*p));
            green=UpScale(TIFFGetG(*p));
            blue=UpScale(TIFFGetB(*p));
            index=image->matte ? UpScale(TIFFGetA(*p)) : 0;
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
                    image->pixels=(VigraImpexRunlengthPacket *) realloc((char *)
                      image->pixels,packets*sizeof(VigraImpexRunlengthPacket));
                    if (image->pixels == (VigraImpexRunlengthPacket *) NULL)
                      {
                        free((char *) pixels);
                        TIFFClose(tiff);
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
            p++;
          }
        }
        free((char *) pixels);
        (void) vigraImpexIsPseudoClass(image);
        break;
      }
    }
    image->pixels=(VigraImpexRunlengthPacket *)
      realloc((char *) image->pixels,image->packets*sizeof(VigraImpexRunlengthPacket));
    /*
      Proceed to next image.
    */
    if (image_info->subrange != 0)
      if (image->scene >= (image_info->subimage+image_info->subrange-1))
        break;
    status=TIFFReadDirectory(tiff);
    if (status == True)
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
  } while (status == True);
  TIFFClose(tiff);
  if (image->temporary)
    {
      (void) remove(image->filename);
      image->temporary=False;
    }
  while (image->previous != (VigraImpexImage *) NULL)
    image=image->previous;
  return(image);
}
#undef PrematureExit
#define PrematureExit(message, image) \
{ \
    fprintf(stderr, "TIFF library: %s %s\n", (char *) message, image->filename); \
    if(tiff) TIFFClose(tiff); \
    return 0; \
}
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   W r i t e T I F F I m a g e                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function WriteTIFFImage writes an image in the Tagged image file format.
%
%  The format of the WriteTIFFImage routine is:
%
%      status=WriteTIFFImage(image_info,image)
%
%  A description of each parameter follows:
%
%    o status:  Function WriteTIFFImage return True if the image is written.
%      False is returned is there is of a memory shortage or if the image
%      file cannot be opened for writing.
%
%    o image_info: Specifies a pointer to an VigraImpexImageInfo structure.
%
%    o image:  A pointer to a VigraImpexImage structure.
%
%
*/
Export unsigned int vigraImpexWriteTIFFImage( VigraImpexImageInfo *image_info,VigraImpexImage *image)
{
#if !defined(TIFFDefaultvigraImpexStripSize)
#define TIFFDefaultvigraImpexStripSize(tiff,request)  ((8*1024)/TIFFScanlineSize(tiff))
#endif

  VigraImpexCompressionType
    compression;

  VigraImpexImage
    encode_image;


  register VigraImpexRunlengthPacket
    *p;

  register int
    i,
    j,
    x,
    y;

  register unsigned char
    *q;

  TIFF
    *tiff = 0;

  uint16
    compress_tag,
    photometric;

  unsigned char
    *scanline;

  unsigned int
    scene,
    x_resolution,
    y_resolution;

  unsigned short
    units,
    value;

  /*
    Open TIFF file.
  */
  vigraImpexOpenImage(image_info,image,WriteBinaryType);
  if (image->file == (FILE *) NULL)
    PrematureExit("Unable to open file",image);
  if ((image->file != stdout) && !image->pipe)
    (void) remove(image->filename);
  else
    {
      /*
        Write standard output or pipe to temporary file.
      */
      encode_image=(*image);
      TemporaryFilename(image->filename);
      image->temporary=True;
    }
  vigraImpexCloseImage(image);
  tiff=TIFFOpen(image->filename,WriteBinaryType);
  if (tiff == (TIFF *) NULL)
    PrematureExit("Unable to open file",image);
  compression=image_info->compression;
#if defined(HasLZW)
  if (compression == VigraImpexUndefinedCompression)
    compression=VigraImpexLZWCompression;
#endif
  scene=0;
  do
  {
    /*
      Initialize TIFF fields.
    */
    TIFFSetField(tiff,TIFFTAG_IMAGELENGTH,(uint32) image->rows);
    TIFFSetField(tiff,TIFFTAG_IMAGEWIDTH,(uint32) image->columns);
    TIFFSetField(tiff,TIFFTAG_BITSPERSAMPLE,8);
    if (image->depth == 16)
      TIFFSetField(tiff,TIFFTAG_BITSPERSAMPLE,16);
    compress_tag=COMPRESSION_NONE;
    if (compression == VigraImpexJPEGCompression)
      compress_tag=COMPRESSION_JPEG;
    if (compression == VigraImpexLZWCompression)
      compress_tag=COMPRESSION_LZW;
    if (compression == VigraImpexRunlengthEncodedCompression)
      compress_tag=COMPRESSION_PACKBITS;
    if (compression == VigraImpexZipCompression)
      compress_tag=COMPRESSION_DEFLATE;
    if (image_info->colorspace == VigraImpexCMYKColorspace)
      {
        photometric=PHOTOMETRIC_SEPARATED;
        TIFFSetField(tiff,TIFFTAG_SAMPLESPERPIXEL,4);
        TIFFSetField(tiff,TIFFTAG_INKSET,INKSET_CMYK);
      }
    else
      if ((strcmp(image_info->magick,"TIFF24") == 0) ||
          (!vigraImpexIsPseudoClass(image) && !vigraImpexIsGrayImage(image)))
        {
          /*
            Full color TIFF raster.
          */
          photometric=PHOTOMETRIC_RGB;
          TIFFSetField(tiff,TIFFTAG_SAMPLESPERPIXEL,(image->matte ? 4 : 3));
          if (image->matte)
            {
              uint16
                extra_samples,
                sample_info[1];

              /*
                TIFF has a matte channel.
              */
              extra_samples=1;
              sample_info[0]=EXTRASAMPLE_ASSOCALPHA;
              TIFFSetField(tiff,TIFFTAG_EXTRASAMPLES,extra_samples,
                &sample_info);
            }
        }
      else
        {
          /*
            Colormapped TIFF raster.
          */
          TIFFSetField(tiff,TIFFTAG_SAMPLESPERPIXEL,1);
          photometric=PHOTOMETRIC_PALETTE;
          if (image->colors <= 2)
            {
              if (vigraImpexIsMonochromeImage(image))
                photometric=PHOTOMETRIC_MINISWHITE;
              if (compression != VigraImpexNoCompression)
                compress_tag=COMPRESSION_CCITTFAX4;
              TIFFSetField(tiff,TIFFTAG_BITSPERSAMPLE,1);
            }
          else
            if (vigraImpexIsGrayImage(image))
              photometric=PHOTOMETRIC_MINISBLACK;
        }
    TIFFSetField(tiff,TIFFTAG_PHOTOMETRIC,photometric);
    TIFFSetField(tiff,TIFFTAG_COMPRESSION,compress_tag);
    TIFFSetField(tiff,TIFFTAG_FILLORDER,FILLORDER_MSB2LSB);
    TIFFSetField(tiff,TIFFTAG_ORIENTATION,ORIENTATION_TOPLEFT);
    TIFFSetField(tiff,TIFFTAG_PLANARCONFIG,PLANARCONFIG_CONTIG);
    TIFFSetField(tiff,TIFFTAG_ROWSPERSTRIP,image->rows);
    if (photometric == PHOTOMETRIC_RGB)
      if ((image_info->interlace == VigraImpexPlaneInterlace) ||
          (image_info->interlace == VigraImpexPartitionInterlace))
        TIFFSetField(tiff,TIFFTAG_PLANARCONFIG,PLANARCONFIG_SEPARATE);
    x_resolution=72;
    y_resolution=72;
    units=RESUNIT_NONE;
    if (image->units == VigraImpexPixelsPerInchResolution)
      units=RESUNIT_INCH;
    if (image->units == VigraImpexPixelsPerCentimeterResolution)
      units=RESUNIT_CENTIMETER;
    if ((image->x_resolution == 0.0) || (image->y_resolution == 0.0))
      {
        units=RESUNIT_NONE;
        image->x_resolution=image->columns;
        image->y_resolution=image->rows;
      }
    TIFFSetField(tiff,TIFFTAG_RESOLUTIONUNIT,(uint16) units);
    TIFFSetField(tiff,TIFFTAG_XRESOLUTION,image->x_resolution);
    TIFFSetField(tiff,TIFFTAG_YRESOLUTION,image->y_resolution);
    TIFFSetField(tiff,TIFFTAG_DOCUMENTNAME,image->filename);
    TIFFSetField(tiff,TIFFTAG_SOFTWARE,"");
    if (image->number_scenes > 1)
      {
        TIFFSetField(tiff,TIFFTAG_SUBFILETYPE,FILETYPE_PAGE);
        TIFFSetField(tiff,TIFFTAG_PAGENUMBER,(unsigned short) image->scene,
          image->number_scenes);
      }
    if (image->label != (char *) NULL)
      TIFFSetField(tiff,TIFFTAG_PAGENAME,image->label);
    if (image->comments != (char *) NULL)
      TIFFSetField(tiff,TIFFTAG_IMAGEDESCRIPTION,image->comments);
    /*
      Write image scanlines.
    */
    scanline=(unsigned char *) malloc(TIFFScanlineSize(tiff));
    if (scanline == (unsigned char *) NULL)
      PrematureExit("Memory allocation failed",image);
    p=image->pixels;
    q=scanline;
    x=0;
    y=0;
    switch (photometric)
    {
      case PHOTOMETRIC_RGB:
      {
        /*
          RGB TIFF image.
        */
        switch (image_info->interlace)
        {
          case VigraImpexNoInterlace:
          default:
          {
            for (i=0; i < image->packets; i++)
            {
              for (j=0; j <= ((int) p->length); j++)
              {
                /*
                  Convert VigraImpexDirectClass packets to contiguous RGB scanlines.
                */
                WriteQuantum(p->red,q);
                WriteQuantum(p->green,q);
                WriteQuantum(p->blue,q);
                if (image->matte)
                  WriteQuantum(p->index,q);
                x++;
                if (x == image->columns)
                  {
                    if (TIFFWriteScanline(tiff,(char *) scanline,y,0) < 0)
                      break;
                    q=scanline;
                    x=0;
                    y++;
                  }
              }
              p++;
            }
            break;
          }
          case VigraImpexPlaneInterlace:
          case VigraImpexPartitionInterlace:
          {
            /*
              Plane interlacing:  RRRRRR...GGGGGG...BBBBBB...
            */
            p=image->pixels;
            for (i=0; i < image->packets; i++)
            {
              for (j=0; j <= ((int) p->length); j++)
              {
                WriteQuantum(p->red,q);
                x++;
                if (x == image->columns)
                  {
                    if (TIFFWriteScanline(tiff,(char *) scanline,y,0) < 0)
                      break;
                    q=scanline;
                    x=0;
                    y++;
                  }
              }
              p++;
            }
            p=image->pixels;
            y=0;
            for (i=0; i < image->packets; i++)
            {
              for (j=0; j <= ((int) p->length); j++)
              {
                WriteQuantum(p->green,q);
                x++;
                if (x == image->columns)
                  {
                    if (TIFFWriteScanline(tiff,(char *) scanline,y,1) < 0)
                      break;
                    q=scanline;
                    x=0;
                    y++;
                  }
              }
              p++;
            }
            p=image->pixels;
            y=0;
            for (i=0; i < image->packets; i++)
            {
              for (j=0; j <= ((int) p->length); j++)
              {
                WriteQuantum(p->blue,q);
                x++;
                if (x == image->columns)
                  {
                    if (TIFFWriteScanline(tiff,(char *) scanline,y,2) < 0)
                      break;
                    q=scanline;
                    x=0;
                    y++;
                  }
              }
              p++;
            }
            p=image->pixels;
            y=0;
            if (image->matte)
              for (i=0; i < image->packets; i++)
              {
                for (j=0; j <= ((int) p->length); j++)
                {
                  WriteQuantum(p->index,q);
                  x++;
                  if (x == image->columns)
                    {
                      if (TIFFWriteScanline(tiff,(char *) scanline,y,3) < 0)
                        break;
                      q=scanline;
                      x=0;
                      y++;
                    }
                }
                p++;
              }
            break;
          }
        }
        break;
      }
      case PHOTOMETRIC_SEPARATED:
      {
        float
           black_generation,
           undercolor;

         int
           black,
           cyan,
           magenta,
           yellow;

        /*
          CMYK TIFF image.
        */
        undercolor=1.0;
        black_generation=1.0;
        if (image_info->undercolor != (char *) NULL)
          {
            (void) sscanf(image_info->undercolor,"%fx%f",&undercolor,
              &black_generation);
            if (black_generation == 1.0)
              black_generation=undercolor;
          }
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
            WriteQuantum((unsigned int) (cyan-undercolor*black),q);
            WriteQuantum((unsigned int) (magenta-undercolor*black),q);
            WriteQuantum((unsigned int) (yellow-undercolor*black),q);
            WriteQuantum((unsigned int) (black_generation*black),q);
            x++;
            if (x == image->columns)
              {
                if (TIFFWriteScanline(tiff,(char *) scanline,y,0) < 0)
                  break;
                q=scanline;
                x=0;
                y++;
              }
          }
          p++;
        }
        break;
      }
      case PHOTOMETRIC_PALETTE:
      {
        unsigned short
          *blue,
          *green,
          *red;

        /*
          Colormapped TIFF image.
        */
        blue=(unsigned short *)
          malloc((1 << image->depth)*sizeof(unsigned short));
        green=(unsigned short *)
          malloc((1 << image->depth)*sizeof(unsigned short));
        red=(unsigned short *)
          malloc((1 << image->depth)*sizeof(unsigned short));
        if ((blue == (unsigned short *) NULL) ||
            (green == (unsigned short *) NULL) ||
            (red == (unsigned short *) NULL))
          PrematureExit("Memory allocation failed",
            image);
        /*
          Initialize TIFF colormap.
        */
        for (i=0; i < image->colors; i++)
        {
          red[i]=(unsigned int) (image->colormap[i].red*65535L)/MaxRGB;
          green[i]=(unsigned int) (image->colormap[i].green*65535L)/MaxRGB;
          blue[i]=(unsigned int) (image->colormap[i].blue*65535L)/MaxRGB;
        }
        for ( ; i < (1 << image->depth); i++)
        {
          red[i]=0;
          green[i]=0;
          blue[i]=0;
        }
        TIFFSetField(tiff,TIFFTAG_COLORMAP,red,green,blue);
        free((char *) red);
        free((char *) green);
        free((char *) blue);
      }
      default:
      {
        register unsigned char
          bit,
          byte,
          polarity;

        if (image->colors > 2)
          {
            /*
              Convert VigraImpexPseudoClass packets to contiguous grayscale scanlines.
            */
            for (i=0; i < image->packets; i++)
            {
              for (j=0; j <= ((int) p->length); j++)
              {
                if (photometric == PHOTOMETRIC_PALETTE)
                  WriteQuantum(p->index,q)
                else
                  WriteQuantum(Intensity(*p),q);
                x++;
                if (x == image->columns)
                  {
                    if (TIFFWriteScanline(tiff,(char *) scanline,y,0) < 0)
                      break;
                    q=scanline;
                    x=0;
                    y++;
                  }
              }
              p++;
            }
            break;
          }
        /*
          Convert VigraImpexPseudoClass packets to contiguous monochrome scanlines.
        */
        polarity=0;
        if (photometric == PHOTOMETRIC_PALETTE)
          polarity=1;
        else
          if (image->colors == 2)
            {
              polarity=
                Intensity(image->colormap[0]) > Intensity(image->colormap[1]);
              if (photometric == PHOTOMETRIC_MINISBLACK)
                polarity=!polarity;
            }
        bit=0;
        byte=0;
        x=0;
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
                if (TIFFWriteScanline(tiff,(char *) scanline,y,0) < 0)
                  break;
                q=scanline;
                bit=0;
                byte=0;
                x=0;
                y++;
             }
          }
          p++;
        }
        break;
      }
    }
    free((char *) scanline);
    if (image_info->verbose == True)
      TIFFPrintDirectory(tiff,stderr,False);
    TIFFWriteDirectory(tiff);
    if (image->next == (VigraImpexImage *) NULL)
      break;
    image->next->file=image->file;
    image=image->next;
  } while (image_info->adjoin);
  (void) TIFFClose(tiff);
  if (image->temporary)
    {
      FILE
        *file;

      int
        c;

      /*
        Copy temporary file to standard output or pipe.
      */
      file=fopen(image->filename,ReadBinaryType);
      if (file == (FILE *) NULL)
        PrematureExit("Unable to open file",image);
      c=fgetc(file);
      while (c != EOF)
      {
        (void) putc(c,encode_image.file);
        c=fgetc(file);
      }
      (void) fclose(file);
      (void) remove(image->filename);
      image->temporary=False;
      vigraImpexCloseImage(&encode_image);
    }
  return(True);
}
#else
Export VigraImpexImage *vigraImpexReadTIFFImage( VigraImpexImageInfo *image_info)
{
  fprintf(stderr, "TIFF library is not available");
  return 0;
}
Export unsigned int vigraImpexWriteTIFFImage( VigraImpexImageInfo *image_info,VigraImpexImage *image)
{
  fprintf(stderr, "TIFF library is not available");
  return 0;
}
#endif
