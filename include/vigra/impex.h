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
#ifndef VIGRA_VIGRAIMPEX_H
#define VIGRA_VIGRAIMPEX_H

/*
  System include declarations.
*/

#if defined(__hpux)
#define _HPUX_SOURCE  1
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


#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

#define MaxTextExtent  1664


/*
  Enumeration declarations.
*/
typedef enum
{
  VigraImpexUndefinedClass,
  VigraImpexDirectClass,
  VigraImpexPseudoClass
} VigraImpexClassType;

typedef enum
{
  VigraImpexUndefinedColorspace,
  VigraImpexRGBColorspace,
  VigraImpexGRAYColorspace,
  VigraImpexTransparentColorspace,
  VigraImpexOHTAColorspace,
  VigraImpexXYZColorspace,
  VigraImpexYCbCrColorspace,
  VigraImpexYCCColorspace,
  VigraImpexYIQColorspace,
  VigraImpexYPbPrColorspace,
  VigraImpexYUVColorspace,
  VigraImpexCMYKColorspace
} VigraImpexColorspaceType;

typedef enum
{
  VigraImpexUndefinedCompression,
  VigraImpexNoCompression,
  VigraImpexLZWCompression,
  VigraImpexJPEGCompression,
  VigraImpexRunlengthEncodedCompression,
  VigraImpexZipCompression
} VigraImpexCompressionType;

typedef enum
{
  VigraImpexUndefinedId,
  VigraImpexImageMagickId
} VigraImpexIdType;

typedef enum
{
  VigraImpexUndefinedInterlace,
  VigraImpexNoInterlace,
  VigraImpexLineInterlace,
  VigraImpexPlaneInterlace,
  VigraImpexPartitionInterlace
} VigraImpexInterlaceType;

typedef enum
{
  VigraImpexGammaPreview
} VigraImpexPreviewType;

typedef enum
{
  VigraImpexUndefinedResolution,
  VigraImpexPixelsPerInchResolution,
  VigraImpexPixelsPerCentimeterResolution
} VigraImpexResolutionType;

/*
  Typedef declarations.
*/
typedef struct _VigraImpexColorPacket
{
  unsigned short /* Quantum */
    red,
    green,
    blue;

  unsigned char
    flags;

  char
    key[3];

  unsigned short
    index;
} VigraImpexColorPacket;

typedef struct _VigraImpexContributionInfo
{
  int
    pixel;

  long
    weight;
} VigraImpexContributionInfo;

typedef struct _VigraImpexFrameInfo
{
  int
    x,
    y;

  unsigned int
    width,
    height;

  int
    inner_bevel,
    outer_bevel;
} VigraImpexFrameInfo;

typedef struct _VigraImpexImageInfo
{
  char
    *filename,
    magick[MaxTextExtent];

  unsigned int
    affirm,
    subimage,
    subrange;

  char
    *server_name,
    *font,
    *pen,
    *box,
    *size,
    *tile,
    *density,
    *page,
    *dispose,
    *delay,
    *iterations,
    *texture,
    *view;

  unsigned int
    adjoin;

  VigraImpexColorspaceType
    colorspace;

  VigraImpexCompressionType
    compression;

  unsigned int
    dither;

  VigraImpexInterlaceType
    interlace;

  unsigned int
    monochrome,
    pointsize,
    quality,
    verbose;

  VigraImpexPreviewType
    preview_type;

  char
    *undercolor;

  unsigned int
    ping;
} VigraImpexImageInfo;

typedef struct _VigraImpexPointInfo
{
  float
    x,
    y;
} VigraImpexPointInfo;

typedef struct _VigraImpexQuantizeInfo
{
  unsigned int
    number_colors,
    tree_depth,
    dither;

  VigraImpexColorspaceType
    colorspace;
} VigraImpexQuantizeInfo;

typedef struct _VigraImpexRectangleInfo
{
  unsigned int
    width,
    height;

  int
    x,
    y;
} VigraImpexRectangleInfo;

typedef struct _VigraImpexRunlengthPacket
{
  unsigned short /* Quantum */
    red,
    green,
    blue,
    length;

  unsigned short
    index;
} VigraImpexRunlengthPacket;

typedef struct _VigraImpexSegmentInfo
{
  int
    x1,
    y1,
    x2,
    y2;
} VigraImpexSegmentInfo;

struct _VigraImpexImage
{
  FILE
    *file;

  int
    status,
    temporary;

  char
    filename[MaxTextExtent];

  long int
    filesize;

  int
    pipe;

  char
    magick[MaxTextExtent],
    *comments,
    *label,
    *text;

  VigraImpexIdType
    id;

  VigraImpexClassType
    c_class;

  unsigned int
    matte;

  VigraImpexCompressionType
    compression;

  unsigned int
    columns,
    rows,
    depth;

  VigraImpexInterlaceType
    interlace;

  unsigned int
    scene,
    number_scenes;

  char
    *montage,
    *directory;

  VigraImpexColorPacket
    *colormap;

  unsigned int
    colors;

  double
    gamma;

  VigraImpexResolutionType
    units;

  float
    x_resolution,
    y_resolution;

  unsigned int
    mean_error_per_pixel;

  double
    normalized_mean_error,
    normalized_maximum_error;

  unsigned long
    total_colors;

  char
    *signature;

  VigraImpexRunlengthPacket
    *pixels,
    *packet;

  unsigned int
    packets,
    runlength,
    packet_size;

  unsigned char
    *packed_pixels;

  VigraImpexColorPacket
    background_color,
    border_color,
    matte_color;

  long int
    magick_time;

  char
    magick_filename[MaxTextExtent];

  unsigned int
    magick_columns,
    magick_rows;

  char
    *geometry,
    *page;

  unsigned int
    dispose,
    delay,
    iterations;

  unsigned int
    orphan;

  struct _VigraImpexImage
    *previous,
    *list,
    *next;
};

typedef struct _VigraImpexImage VigraImpexImage;


/*
  Image utilities routines.
*/
Export VigraImpexImage *vigraImpexAllocateImage(const VigraImpexImageInfo *image_info);
Export void vigraImpexAllocateNextImage(const VigraImpexImageInfo *image_info,VigraImpexImage *image);
Export void vigraImpexCloseImage(VigraImpexImage *image);
Export void vigraImpexCondenseImage(VigraImpexImage *image);
Export VigraImpexImage *vigraImpexCloneImage(VigraImpexImage *image,const unsigned int columns,
  const unsigned int rows,const unsigned int clone_pixels);
Export void vigraImpexDestroyImage(VigraImpexImage *image);
Export void vigraImpexDestroyImageInfo(VigraImpexImageInfo *image_info);
Export void vigraImpexDestroyImages(VigraImpexImage *image);
Export void vigraImpexGetImageInfo(VigraImpexImageInfo *image_info);
Export unsigned int vigraImpexIsGrayImage(VigraImpexImage *image);
Export unsigned int vigraImpexIsMonochromeImage(VigraImpexImage *image);
Export void vigraImpexNormalizeImage(VigraImpexImage *image);
Export void vigraImpexOpenImage(const VigraImpexImageInfo *image_info,VigraImpexImage *image,const char *type);
Export void vigraImpexRGBTransformImage(VigraImpexImage *image,const unsigned int colorspace);
unsigned int vigraImpexIsGeometry(char * m);
unsigned int vigraImpexIsSubimage(char * m, unsigned int i);
Export void vigraImpexSyncImage(VigraImpexImage *image);
Export void vigraImpexTransformRGBImage(VigraImpexImage *image,const unsigned int colorspace);
Export unsigned int vigraImpexUncondenseImage(VigraImpexImage *image);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif


#endif /* VIGRA_VIGRAIMPEX_H */
