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
#include "error.h"

#define PrematureExit(arg, message, image) \
{ \
    fprintf(stderr, message " %s\n", image->filename); \
    vigraImpexDestroyImage(image); \
    return 0; \
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   R e a d V I F F I m a g e                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexReadVIFFImage reads a Khoros Visualization image file and returns
%  it.  It allocates the memory necessary for the new VigraImpexImage structure and
%  returns a pointer to the new image.
%
%  The format of the vigraImpexReadVIFFImage routine is:
%
%      image=vigraImpexReadVIFFImage(image_info)
%
%  A description of each parameter follows:
%
%    o image: Function vigraImpexReadVIFFImage returns a pointer to the image after
%      reading.  A null image is returned if there is a a memory shortage or if
%      the image cannot be read.
%
%    o filename: Specifies the name of the image to read.
%
%
*/
Export VigraImpexImage *vigraImpexReadVIFFImage(VigraImpexImageInfo *image_info)
{
#define VFF_CM_genericRGB  15
#define VFF_CM_ntscRGB  1
#define VFF_CM_NONE  0
#define VFF_DEP_DECORDER  0x4
#define VFF_DEP_NSORDER  0x8
#define VFF_DES_RAW  0
#define VFF_LOC_IMPLICIT  1
#define VFF_MAPTYP_NONE  0
#define VFF_MAPTYP_1_BYTE  1
#define VFF_MS_NONE  0
#define VFF_MS_ONEPERBAND  1
#define VFF_MS_SHARED  3
#define VFF_TYP_BIT  0
#define VFF_TYP_1_BYTE  1
#define VFF_TYP_2_BYTE  2
#define VFF_TYP_4_BYTE  4

  typedef struct _ViffHeader
  {
    unsigned char
      identifier,
      file_type,
      release,
      version,
      machine_dependency,
      reserve[3];

    char
      comment[512];

    unsigned long
      rows,
      columns,
      subrows;

    long
      x_offset,
      y_offset;

    float
      x_bits_per_pixel,
      y_bits_per_pixel;

    unsigned long
      location_type,
      location_dimension,
      number_of_images,
      number_data_bands,
      data_storage_type,
      data_encode_scheme,
      map_scheme,
      map_storage_type,
      map_rows,
      map_columns,
      map_subrows,
      map_enable,
      maps_per_cycle,
      color_space_model;
  } ViffHeader;

  VigraImpexImage
    *image;

  register int
    bit,
    i,
    x,
    y;

  register Quantum
    *p;

  register VigraImpexRunlengthPacket
    *q;

  unsigned char
    buffer[7],
    *viff_pixels;

  unsigned int
    bytes_per_pixel,
    status;

  unsigned long
    packets;

  ViffHeader
    viff_header;

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
    PrematureExit(FileOpenWarning,"Unable to open file",image);
  /*
    Read VIFF header (1024 bytes).
  */
  status=vigraImpexReadData((char *) &viff_header.identifier,1,1,image->file);
  do
  {
    /*
      Verify VIFF identifier.
    */
    if ((status == False) || ((unsigned char) viff_header.identifier != 0xab))
      PrematureExit(CorruptImageWarning,"Not a VIFF raster",image);
    /*
      Initialize VIFF image.
    */
    (void) vigraImpexReadData((char *) buffer,1,7,image->file);
    viff_header.file_type=buffer[0];
    viff_header.release=buffer[1];
    viff_header.version=buffer[2];
    viff_header.machine_dependency=buffer[3];
    (void) vigraImpexReadData((char *) viff_header.comment,1,512,image->file);
    viff_header.comment[511]='\0';
    if (Extent(viff_header.comment) > 4)
      {
        image->comments=(char *)
          malloc((unsigned int) (Extent(viff_header.comment)+1)*sizeof(char));
        if (image->comments == (char *) NULL)
          PrematureExit(ResourceLimitWarning,"Memory allocation failed",image);
        (void) strcpy(image->comments,viff_header.comment);
      }
    if ((viff_header.machine_dependency == VFF_DEP_DECORDER) ||
        (viff_header.machine_dependency == VFF_DEP_NSORDER))
      {
        viff_header.rows=vigraImpexLSBFirstReadLong(image->file);
        viff_header.columns=vigraImpexLSBFirstReadLong(image->file);
        viff_header.subrows=vigraImpexLSBFirstReadLong(image->file);
        viff_header.x_offset=(long) vigraImpexLSBFirstReadLong(image->file);
        viff_header.y_offset=(long) vigraImpexLSBFirstReadLong(image->file);
        viff_header.x_bits_per_pixel=(float) vigraImpexLSBFirstReadLong(image->file);
        viff_header.y_bits_per_pixel=(float) vigraImpexLSBFirstReadLong(image->file);
        viff_header.location_type=vigraImpexLSBFirstReadLong(image->file);
        viff_header.location_dimension=vigraImpexLSBFirstReadLong(image->file);
        viff_header.number_of_images=vigraImpexLSBFirstReadLong(image->file);
        viff_header.number_data_bands=vigraImpexLSBFirstReadLong(image->file);
        viff_header.data_storage_type=vigraImpexLSBFirstReadLong(image->file);
        viff_header.data_encode_scheme=vigraImpexLSBFirstReadLong(image->file);
        viff_header.map_scheme=vigraImpexLSBFirstReadLong(image->file);
        viff_header.map_storage_type=vigraImpexLSBFirstReadLong(image->file);
        viff_header.map_rows=vigraImpexLSBFirstReadLong(image->file);
        viff_header.map_columns=vigraImpexLSBFirstReadLong(image->file);
        viff_header.map_subrows=vigraImpexLSBFirstReadLong(image->file);
        viff_header.map_enable=vigraImpexLSBFirstReadLong(image->file);
        viff_header.maps_per_cycle=vigraImpexLSBFirstReadLong(image->file);
        viff_header.color_space_model=vigraImpexLSBFirstReadLong(image->file);
      }
    else
      {
        viff_header.rows=vigraImpexMSBFirstReadLong(image->file);
        viff_header.columns=vigraImpexMSBFirstReadLong(image->file);
        viff_header.subrows=vigraImpexMSBFirstReadLong(image->file);
        viff_header.x_offset=(long) vigraImpexMSBFirstReadLong(image->file);
        viff_header.y_offset=(long) vigraImpexMSBFirstReadLong(image->file);
        viff_header.x_bits_per_pixel=(float) vigraImpexMSBFirstReadLong(image->file);
        viff_header.y_bits_per_pixel=(float) vigraImpexMSBFirstReadLong(image->file);
        viff_header.location_type=vigraImpexMSBFirstReadLong(image->file);
        viff_header.location_dimension=vigraImpexMSBFirstReadLong(image->file);
        viff_header.number_of_images=vigraImpexMSBFirstReadLong(image->file);
        viff_header.number_data_bands=vigraImpexMSBFirstReadLong(image->file);
        viff_header.data_storage_type=vigraImpexMSBFirstReadLong(image->file);
        viff_header.data_encode_scheme=vigraImpexMSBFirstReadLong(image->file);
        viff_header.map_scheme=vigraImpexMSBFirstReadLong(image->file);
        viff_header.map_storage_type=vigraImpexMSBFirstReadLong(image->file);
        viff_header.map_rows=vigraImpexMSBFirstReadLong(image->file);
        viff_header.map_columns=vigraImpexMSBFirstReadLong(image->file);
        viff_header.map_subrows=vigraImpexMSBFirstReadLong(image->file);
        viff_header.map_enable=vigraImpexMSBFirstReadLong(image->file);
        viff_header.maps_per_cycle=vigraImpexMSBFirstReadLong(image->file);
        viff_header.color_space_model=vigraImpexMSBFirstReadLong(image->file);
      }
    for (i=0; i < 420; i++)
      (void) fgetc(image->file);
    image->columns=(unsigned int) viff_header.rows;
    image->rows=(unsigned int) viff_header.columns;
    /*
      Verify that we can read this VIFF image.
    */
    if ((viff_header.columns*viff_header.rows) == 0)
      PrematureExit(CorruptImageWarning,
        "VigraImpexImage column or row size is not supported",image);
    if ((viff_header.data_storage_type != VFF_TYP_BIT) &&
        (viff_header.data_storage_type != VFF_TYP_1_BYTE) &&
        (viff_header.data_storage_type != VFF_TYP_2_BYTE) &&
        (viff_header.data_storage_type != VFF_TYP_4_BYTE))
      PrematureExit(CorruptImageWarning,
        "Data storage type is not supported",image);
    if (viff_header.data_encode_scheme != VFF_DES_RAW)
      PrematureExit(CorruptImageWarning,
        "Data encoding scheme is not supported",image);
    if ((viff_header.map_storage_type != VFF_MAPTYP_NONE) &&
        (viff_header.map_storage_type != VFF_MAPTYP_1_BYTE))
      PrematureExit(CorruptImageWarning,
        "Map storage type is not supported",image);
    if ((viff_header.color_space_model != VFF_CM_NONE) &&
        (viff_header.color_space_model != VFF_CM_ntscRGB) &&
        (viff_header.color_space_model != VFF_CM_genericRGB))
      PrematureExit(CorruptImageWarning,
        "Colorspace model is not supported",image);
    if (viff_header.location_type != VFF_LOC_IMPLICIT)
      {
        vigraImpexMagickWarning(CorruptImageWarning,
          "Location type is not supported",image->filename);
        vigraImpexDestroyImage(image);
        return((VigraImpexImage *) NULL);
      }
    if (viff_header.number_of_images != 1)
      PrematureExit(CorruptImageWarning,
        "Number of images is not supported",image);
    switch ((int)viff_header.map_scheme)
    {
      case VFF_MS_NONE:
      {
        if (viff_header.number_data_bands < 3)
          {
            /*
              Create linear color ramp.
            */
            if (viff_header.data_storage_type == VFF_TYP_BIT)
              image->colors=2;
            else
              image->colors=1 << (viff_header.number_data_bands*QuantumDepth);
            image->colormap=(VigraImpexColorPacket *)
              malloc(image->colors*sizeof(VigraImpexColorPacket));
            if (image->colormap == (VigraImpexColorPacket *) NULL)
              PrematureExit(ResourceLimitWarning,"Memory allocation failed",
                image);
            for (i=0; i < image->colors; i++)
            {
              image->colormap[i].red=(MaxRGB*i)/(image->colors-1);
              image->colormap[i].green=(MaxRGB*i)/(image->colors-1);
              image->colormap[i].blue=(MaxRGB*i)/(image->colors-1);
            }
          }
        break;
      }
      case VFF_MS_ONEPERBAND:
      case VFF_MS_SHARED:
      {
        unsigned char
          *viff_colormap;

        /*
          Read VIFF raster colormap.
        */
        image->colors=(unsigned int) viff_header.map_columns;
        image->colormap=(VigraImpexColorPacket *)
          malloc(image->colors*sizeof(VigraImpexColorPacket));
        viff_colormap=(unsigned char *)
          malloc(image->colors*sizeof(unsigned char));
        if ((image->colormap == (VigraImpexColorPacket *) NULL) ||
            (viff_colormap == (unsigned char *) NULL))
          PrematureExit(ResourceLimitWarning,"Memory allocation failed",image);
        (void) vigraImpexReadData((char *) viff_colormap,1,image->colors,image->file);
        for (i=0; i < image->colors; i++)
        {
          image->colormap[i].red=UpScale(viff_colormap[i]);
          image->colormap[i].green=UpScale(viff_colormap[i]);
          image->colormap[i].blue=UpScale(viff_colormap[i]);
        }
        if (viff_header.map_rows > 1)
          {
            (void) vigraImpexReadData((char *) viff_colormap,1,image->colors,image->file);
            for (i=0; i < image->colors; i++)
              image->colormap[i].green=UpScale(viff_colormap[i]);
          }
        if (viff_header.map_rows > 2)
          {
            (void) vigraImpexReadData((char *) viff_colormap,1,image->colors,image->file);
            for (i=0; i < image->colors; i++)
              image->colormap[i].blue=UpScale(viff_colormap[i]);
          }
        free((char *) viff_colormap);
        break;
      }
      default:
        PrematureExit(CorruptImageWarning,"Colormap type is not supported",
          image);
    }
    /*
      Allocate VIFF pixels.
    */
    bytes_per_pixel=1;
    if (viff_header.data_storage_type == VFF_TYP_2_BYTE)
      bytes_per_pixel=2;
    if (viff_header.data_storage_type == VFF_TYP_4_BYTE)
      bytes_per_pixel=4;
    if (viff_header.data_storage_type == VFF_TYP_BIT)
      packets=((viff_header.columns+7) >> 3)*viff_header.rows;
    else
      packets=
        viff_header.columns*viff_header.rows*viff_header.number_data_bands;
    viff_pixels=(unsigned char *)
      malloc(bytes_per_pixel*packets*sizeof(Quantum));
    if (viff_pixels == (unsigned char *) NULL)
      PrematureExit(ResourceLimitWarning,"Memory allocation failed",image);
    (void) vigraImpexReadData((char *) viff_pixels,bytes_per_pixel,(unsigned int) packets,
      image->file);
    switch ((int)viff_header.data_storage_type)
    {
      int
        max_value,
        min_value,
        value;

      register Quantum
        *q;

      unsigned long
        scale_factor;

      case VFF_TYP_1_BYTE:
      {
        register unsigned char
          *p;

        if (QuantumDepth == 8)
          break;
        /*
          Scale integer pixels to [0..MaxRGB].
        */
        p=viff_pixels;
        q=(Quantum *) viff_pixels;
        p+=packets-1;
        q+=packets-1;
        for (i=0; i < packets; i++)
        {
          value=UpScale(*p);
          *q=(Quantum) value;
          p--;
          q--;
        }
        break;
      }
      case VFF_TYP_2_BYTE:
      {
        register short int
          *p;

        /*
          Ensure the header byte-order is most-significant byte first.
        */
        if ((viff_header.machine_dependency == VFF_DEP_DECORDER) ||
            (viff_header.machine_dependency == VFF_DEP_NSORDER))
          vigraImpexMSBFirstOrderShort((char *) &viff_header,
            (unsigned int) (bytes_per_pixel*packets));
        /*
          Determine scale factor.
        */
        p=(short int *) viff_pixels;
        max_value=(*p);
        min_value=(*p);
        for (i=0; i < packets; i++)
        {
          if (*p > max_value)
            max_value=(*p);
          else
            if (*p < min_value)
              min_value=(*p);
          p++;
        }
        if ((min_value == 0) && (max_value == 0))
          scale_factor=0;
        else
          if (min_value == max_value)
            {
              scale_factor=UpShift(MaxRGB)/min_value;
              min_value=0;
            }
          else
            scale_factor=UpShift(MaxRGB)/(max_value-min_value);
        /*
          Scale integer pixels to [0..MaxRGB].
        */
        p=(short int *) viff_pixels;
        q=(Quantum *) viff_pixels;
        for (i=0; i < packets; i++)
        {
          value=DownShift((*p-min_value)*scale_factor);
          if (value > MaxRGB)
            value=MaxRGB;
          else
            if (value < 0)
              value=0;
          *q=(Quantum) value;
          p++;
          q++;
        }
        break;
      }
      case VFF_TYP_4_BYTE:
      {
        register int
          *p;

        /*
          Ensure the header byte-order is most-significant byte first.
        */
        if ((viff_header.machine_dependency == VFF_DEP_DECORDER) ||
            (viff_header.machine_dependency == VFF_DEP_NSORDER))
          vigraImpexMSBFirstOrderLong((char *) &viff_header,
            (unsigned int) (bytes_per_pixel*packets));
        /*
          Determine scale factor.
        */
        p=(int *) viff_pixels;
        max_value=(*p);
        min_value=(*p);
        for (i=0; i < packets; i++)
        {
          if (*p > max_value)
            max_value=(*p);
          else
            if (*p < min_value)
              min_value=(*p);
          p++;
        }
        if ((min_value == 0) && (max_value == 0))
          scale_factor=0;
        else
          if (min_value == max_value)
            {
              scale_factor=UpShift(MaxRGB)/min_value;
              min_value=0;
            }
          else
            scale_factor=UpShift(MaxRGB)/(max_value-min_value);
        /*
          Scale integer pixels to [0..MaxRGB].
        */
        p=(int *) viff_pixels;
        q=(Quantum *) viff_pixels;
        for (i=0; i < packets; i++)
        {
          value=DownShift((*p-min_value)*scale_factor);
          if (value > MaxRGB)
            value=MaxRGB;
          else
            if (value < 0)
              value=0;
          *q=(unsigned char) value;
          p++;
          q++;
        }
        break;
      }
    }
    /*
      Initialize image structure.
    */
    image->matte=(viff_header.number_data_bands == 4);
    image->c_class=
      (viff_header.number_data_bands < 3 ? VigraImpexPseudoClass : VigraImpexDirectClass);
    image->columns=(unsigned int) viff_header.rows;
    image->rows=(unsigned int) viff_header.columns;
    if (image_info->ping)
      {
        vigraImpexCloseImage(image);
        return(image);
      }
    image->packets=image->columns*image->rows;
    image->pixels=(VigraImpexRunlengthPacket *)
      malloc(image->packets*sizeof(VigraImpexRunlengthPacket));
    if (image->pixels == (VigraImpexRunlengthPacket *) NULL)
    {
      PrematureExit(ResourceLimitWarning,"Memory allocation failed",image);
    }
    /*
      Convert VIFF raster image to runlength-encoded packets.
    */
    p=(Quantum *) viff_pixels;
    q=image->pixels;
    if (viff_header.data_storage_type == VFF_TYP_BIT)
      {
        unsigned int
          polarity;

        /*
          Convert bitmap scanline to runlength-encoded color packets.
        */
        polarity=(viff_header.machine_dependency == VFF_DEP_DECORDER) ||
          (viff_header.machine_dependency == VFF_DEP_NSORDER);
        for (y=0; y < image->rows; y++)
        {
          /*
            Convert bitmap scanline to runlength-encoded color packets.
          */
          for (x=0; x < (image->columns >> 3); x++)
          {
            for (bit=0; bit < 8; bit++)
            {
              q->index=
                ((*p) & (0x01 << bit) ? (int) polarity : (int) !polarity);
              q->length=0;
              q++;
            }
            p++;
          }
          if ((image->columns % 8) != 0)
            {
              for (bit=0; bit < (image->columns % 8); bit++)
              {
                q->index=
                  ((*p) & (0x01 << bit) ? (int) polarity : (int) !polarity);
                q->length=0;
                q++;
              }
              p++;
            }
        }
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
       }
      else
        {
          unsigned long
            offset;

          /*
            Convert DirectColor scanline to runlength-encoded color packets.
          */
          offset=image->columns*image->rows;
          for (y=0; y < image->rows; y++)
          {
            for (x=0; x < image->columns; x++)
            {
              q->red=(*p);
              q->green=(*(p+offset));
              q->blue=(*(p+offset*2));
              if (image->colors != 0)
                {
                  q->red=image->colormap[q->red].red;
                  q->green=image->colormap[q->green].green;
                  q->blue=image->colormap[q->blue].blue;
                }
              q->index=(unsigned short) (image->matte ? (*(p+offset*3)) : 0);
              q->length=0;
              p++;
              q++;
            }
          }
        }
    free((char *) viff_pixels);
    if (image->c_class == VigraImpexPseudoClass)
      vigraImpexSyncImage(image);
    /*
      Proceed to next image.
    */
    if (image_info->subrange != 0)
      if (image->scene >= (image_info->subimage+image_info->subrange-1))
        break;
    status=vigraImpexReadData((char *) &viff_header.identifier,1,1,image->file);
    if ((status == True) && (viff_header.identifier == 0xab))
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
  } while ((status == True) && (viff_header.identifier == 0xab));
  vigraImpexCondenseImage(image);
  while (image->previous != (VigraImpexImage *) NULL)
    image=image->previous;
  vigraImpexCloseImage(image);
  return(image);
}
