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
#include "gif.h"
#include "utility.h"


/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G I F D e c o d e I m a g e                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexGIFDecodeImage uncompresses an image via GIF-coding.
%
%  The format of the vigraImpexGIFDecodeImage routine is:
%
%      status=vigraImpexGIFDecodeImage(image)
%
%  A description of each parameter follows:
%
%    o status:  Function vigraImpexGIFDecodeImage returns True if all the pixels are
%      uncompressed without error, otherwise False.
%
%    o image: The address of a structure of type VigraImpexImage.
%
%
*/
Export unsigned int vigraImpexGIFDecodeImage(VigraImpexImage *image)
{
#define MaxStackSize  4096
#define NullCode  (-1)

  int
    available,
    clear,
    code_mask,
    code_size,
    end_of_information,
    in_code,
    old_code,
    status;

  register int
    bits,
    code,
    count,
    i;

  register VigraImpexRunlengthPacket
    *p;

  register unsigned char
    *c;

  register unsigned int
    datum;

  short
    *prefix;

  unsigned char
    data_size,
    first,
    *packet,
    *pixel_stack,
    *suffix,
    *top_stack;

  /*
    Allocate decoder tables.
  */
  assert(image != (VigraImpexImage *) NULL);
  packet=(unsigned char *) malloc(256*sizeof(unsigned char));
  prefix=(short *) malloc(MaxStackSize*sizeof(short));
  suffix=(unsigned char *) malloc(MaxStackSize*sizeof(unsigned char));
  pixel_stack=(unsigned char *) malloc((MaxStackSize+1)*sizeof(unsigned char));
  if ((packet == (unsigned char *) NULL) ||
      (prefix == (short *) NULL) ||
      (suffix == (unsigned char *) NULL) ||
      (pixel_stack == (unsigned char *) NULL))
    return(False);
  /*
    Initialize GIF data stream decoder.
  */
  data_size=fgetc(image->file);
  clear=1 << data_size;
  end_of_information=clear+1;
  available=clear+2;
  old_code=NullCode;
  code_size=data_size+1;
  code_mask=(1 << code_size)-1;
  for (code=0; code < clear; code++)
  {
    prefix[code]=0;
    suffix[code]=code;
  }
  /*
    Decode GIF pixel stream.
  */
  datum=0;
  bits=0;
  c=0;
  count=0;
  first=0;
  top_stack=pixel_stack;
  p=image->pixels;
  for (i=0; i < image->packets; )
  {
    if (top_stack == pixel_stack)
      {
        if (bits < code_size)
          {
            /*
              Load bytes until there is enough bits for a code.
            */
            if (count == 0)
              {
                /*
                  Read a new data block.
                */
                count=vigraImpexReadDataBlock((char *) packet,image->file);
                if (count <= 0)
                  break;
                c=packet;
              }
            datum+=(*c) << bits;
            bits+=8;
            c++;
            count--;
            continue;
          }
        /*
          Get the next code.
        */
        code=datum & code_mask;
        datum>>=code_size;
        bits-=code_size;
        /*
          Interpret the code
        */
        if ((code > available) || (code == end_of_information))
          break;
        if (code == clear)
          {
            /*
              Reset decoder.
            */
            code_size=data_size+1;
            code_mask=(1 << code_size)-1;
            available=clear+2;
            old_code=NullCode;
            continue;
          }
        if (old_code == NullCode)
          {
            *top_stack++=suffix[code];
            old_code=code;
            first=code;
            continue;
          }
        in_code=code;
        if (code == available)
          {
            *top_stack++=first;
            code=old_code;
          }
        while (code > clear)
        {
          *top_stack++=suffix[code];
          code=prefix[code];
        }
        first=suffix[code];
        /*
          Add a new string to the string table,
        */
        if (available >= MaxStackSize)
          break;
        *top_stack++=first;
        prefix[available]=old_code;
        suffix[available]=first;
        available++;
        if (((available & code_mask) == 0) && (available < MaxStackSize))
          {
            code_size++;
            code_mask+=available;
          }
        old_code=in_code;
      }
    /*
      Pop a pixel off the pixel stack.
    */
    top_stack--;
    p->index=(unsigned short) *top_stack;
    p->length=0;
    p++;
    i++;
  }
  /*
    Initialize any remaining color packets to a known color.
  */
  status=i == image->packets;
  for ( ; i < image->packets; i++)
  {
    p->index=0;
    p->length=0;
    p++;
  }
  vigraImpexSyncImage(image);
  /*
    Free decoder memory.
  */
  free((char *) pixel_stack);
  free((char *) suffix);
  free((char *) prefix);
  free((char *) packet);
  return(status);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%   G I F E n c o d e I m a g e                                               %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexGIFEncodeImage compresses an image via GIF-coding.
%
%  The format of the vigraImpexGIFEncodeImage routine is:
%
%      status=vigraImpexGIFEncodeImage(image,data_size)
%
%  A description of each parameter follows:
%
%    o status:  Function vigraImpexGIFEncodeImage returns True if all the pixels are
%      compressed without error, otherwise False.
%
%    o image: The address of a structure of type VigraImpexImage.
%
%
*/
Export unsigned int vigraImpexGIFEncodeImage(VigraImpexImage *image,const unsigned int data_size)
{
#define MaxCode(number_bits)  ((1 << (number_bits))-1)
#define MaxHashTable  5003
#define MaxGIFBits  12
#if defined(HasLZW)
#define MaxGIFTable  (1 << MaxGIFBits)
#else
#define MaxGIFTable  max_code
#endif
#define GIFOutputCode(code) \
{ \
  /*  \
    Emit a code. \
  */ \
  if (bits > 0) \
    datum|=((long) code << bits); \
  else \
    datum=(long) code; \
  bits+=number_bits; \
  while (bits >= 8) \
  { \
    /*  \
      Add a character to current packet. \
    */ \
    packet[byte_count++]=(unsigned char) (datum & 0xff); \
    if (byte_count >= 254) \
      { \
        (void) fputc(byte_count,image->file); \
        (void) fwrite((char *) packet,1,byte_count,image->file); \
        byte_count=0; \
      } \
    datum>>=8; \
    bits-=8; \
  } \
  if (free_code > max_code)  \
    { \
      number_bits++; \
      if (number_bits == MaxGIFBits) \
        max_code=MaxGIFTable; \
      else \
        max_code=MaxCode(number_bits); \
    } \
}

  int
    bits,
    byte_count,
    i,
    number_bits;

  long
    datum;

  register int
    j,
    k;

  register VigraImpexRunlengthPacket
    *p;

  short
    clear_code,
    end_of_information_code,
    free_code,
    *hash_code,
    *hash_prefix,
    index,
    max_code,
    waiting_code;

  unsigned char
    *packet,
    *hash_suffix;

  /*
    Allocate encoder tables.
  */
  assert(image != (VigraImpexImage *) NULL);
  packet=(unsigned char *) malloc(256*sizeof(unsigned char));
  hash_code=(short *) malloc(MaxHashTable*sizeof(short));
  hash_prefix=(short *) malloc(MaxHashTable*sizeof(short));
  hash_suffix=(unsigned char *) malloc(MaxHashTable*sizeof(unsigned char));
  if ((packet == (unsigned char *) NULL) || (hash_code == (short *) NULL) ||
      (hash_prefix == (short *) NULL) ||
      (hash_suffix == (unsigned char *) NULL))
    return(False);
  /*
    Initialize GIF encoder.
  */
  number_bits=data_size;
  max_code=MaxCode(number_bits);
  clear_code=((short) 1 << (data_size-1));
  end_of_information_code=clear_code+1;
  free_code=clear_code+2;
  byte_count=0;
  datum=0;
  bits=0;
  for (i=0; i < MaxHashTable; i++)
    hash_code[i]=0;
  GIFOutputCode(clear_code);
  /*
    Encode pixels.
  */
  p=image->pixels;
  waiting_code=p->index;
  for (i=0; i < image->packets; i++)
  {
    for (j=(i == 0) ? 1 : 0; j <= ((int) p->length); j++)
    {
      /*
        Probe hash table.
      */
      index=p->index & 0xff;
      k=(int) ((int) index << (MaxGIFBits-8))+waiting_code;
      if (k >= MaxHashTable)
        k-=MaxHashTable;
#if defined(HasLZW)
      if (hash_code[k] > 0)
        {
          if ((hash_prefix[k] == waiting_code) && (hash_suffix[k] == index))
            {
              waiting_code=hash_code[k];
              continue;
            }
          if (k == 0)
            displacement=1;
          else
            displacement=MaxHashTable-k;
          next_pixel=False;
          for ( ; ; )
          {
            k-=displacement;
            if (k < 0)
              k+=MaxHashTable;
            if (hash_code[k] == 0)
              break;
            if ((hash_prefix[k] == waiting_code) && (hash_suffix[k] == index))
              {
                waiting_code=hash_code[k];
                next_pixel=True;
                break;
              }
          }
          if (next_pixel == True)
            continue;
        }
#endif
      GIFOutputCode(waiting_code);
      if (free_code < MaxGIFTable)
        {
          hash_code[k]=free_code++;
          hash_prefix[k]=waiting_code;
          hash_suffix[k]=index;
        }
      else
        {
          /*
            Fill the hash table with empty entries.
          */
          for (k=0; k < MaxHashTable; k++)
            hash_code[k]=0;
          /*
            Reset compressor and issue a clear code.
          */
          free_code=clear_code+2;
          GIFOutputCode(clear_code);
          number_bits=data_size;
          max_code=MaxCode(number_bits);
        }
      waiting_code=index;
    }
    p++;
  }
  /*
    Flush out the buffered code.
  */
  GIFOutputCode(waiting_code);
  GIFOutputCode(end_of_information_code);
  if (bits > 0)
    {
      /*
        Add a character to current packet.
      */
      packet[byte_count++]=(unsigned char) (datum & 0xff);
      if (byte_count >= 254)
        {
          (void) fputc(byte_count,image->file);
          (void) fwrite((char *) packet,1,byte_count,image->file);
          byte_count=0;
        }
    }
  /*
    Flush accumulated data.
  */
  if (byte_count > 0)
    {
      (void) fputc(byte_count,image->file);
      (void) fwrite((char *) packet,1,byte_count,image->file);
    }
  /*
    Free encoder memory.
  */
  free((char *) hash_suffix);
  free((char *) hash_prefix);
  free((char *) hash_code);
  free((char *) packet);
  if (i < image->packets)
    return(False);
  return(True);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   W r i t e G I F I m a g e                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexWriteGIFImage writes an image to a file in the Compuserve Graphics
%  image format.
%
%  The format of the vigraImpexWriteGIFImage routine is:
%
%      status=vigraImpexWriteGIFImage(image_info,image)
%
%  A description of each parameter follows.
%
%    o status: Function vigraImpexWriteGIFImage return True if the image is written.
%      False is returned is there is a memory shortage or if the image file
%      fails to write.
%
%    o image_info: Specifies a pointer to an VigraImpexImageInfo structure.
%
%    o image:  A pointer to a VigraImpexImage structure.
%
%
*/
Export unsigned int vigraImpexWriteGIFImage(VigraImpexImageInfo *image_info,VigraImpexImage *image)
{
  VigraImpexColorPacket
    transparent_pixel;

  VigraImpexImage
    *next_image;

  VigraImpexRectangleInfo
    page_info;

  register int
    i,
    x;

  register VigraImpexRunlengthPacket
    *p;

  register unsigned char
    *q;

  unsigned char
    bits_per_pixel,
    c,
    *colormap;

  unsigned int
    colors,
    global_colormap,
    interlace,
    scene,
    status;

  /*
    Open output image file.
  */
  vigraImpexOpenImage(image_info,image,WriteBinaryType);
  if (image->file == (FILE *) NULL)
  {
      fprintf(stderr, "vigraImpexWriteGIFImage(): Unable to open file %s\n", image->filename);
      return False;
  }
  /*
    Allocate colormap.
  */
  colormap=(unsigned char *) malloc(3*256*sizeof(unsigned char));
  if (colormap == (unsigned char *) NULL)
  {
      fprintf(stderr, "vigraImpexWriteGIFImage(): Memory allocation failed\n");
      vigraImpexCloseImage(image);
      return False;
  }
  /*
    Write GIF header.
  */
  if ((image->comments == (char *) NULL) && !image_info->adjoin &&
      !image->matte)
    (void) fwrite("GIF87a",1,6,image->file);
  else
    if (strcmp(image_info->magick,"GIF87") == 0)
      (void) fwrite("GIF87a",1,6,image->file);
    else
      (void) fwrite("GIF89a",1,6,image->file);
  /*
    Determine image bounding box and global colormap status.
  */
  page_info.width=image->columns;
  page_info.height=image->rows;
  page_info.x=0;
  page_info.y=0;
  global_colormap=image->c_class == VigraImpexPseudoClass;
  next_image=image;
  for ( ; next_image != (VigraImpexImage *) NULL; next_image=next_image->next)
  {
    if ((next_image->columns+page_info.x) > page_info.width)
      page_info.width=next_image->columns+page_info.x;
    if ((next_image->rows+page_info.y) > page_info.height)
      page_info.height=next_image->rows+page_info.y;
    if (!global_colormap)
      continue;
    if ((next_image->c_class == VigraImpexDirectClass) ||
        (next_image->colors != image->colors))
      {
        global_colormap=False;
        continue;
      }
    for (i=0; i < image->colors; i++)
      if (!ColorMatch(next_image->colormap[i],image->colormap[i],0))
        {
          global_colormap=False;
          break;
        }
  }

  vigraImpexLSBFirstWriteShort(page_info.width,image->file);
  vigraImpexLSBFirstWriteShort(page_info.height,image->file);
  /*
    Write images to file.
  */
  interlace=image_info->interlace;
  if (image_info->adjoin && (image->next != (VigraImpexImage *) NULL))
    interlace=VigraImpexNoInterlace;
  scene=0;
  do
  {
    transparent_pixel.flags=False;
    if (vigraImpexIsPseudoClass(image))
      colors=image->colors;
    else
      {
        VigraImpexQuantizeInfo
          quantize_info;

        unsigned char
          *matte_image;

        matte_image=(unsigned char *) NULL;
        if (image->matte)
          {
            /*
              Track all the transparent pixels.
            */
            if (!vigraImpexUncondenseImage(image))
              return(False);
            matte_image=(unsigned char *)
              malloc(image->packets*sizeof(unsigned char));
            if (matte_image == (unsigned char *) NULL)
              {
                  fprintf(stderr, "vigraImpexWriteGIFImage(): Memory allocation failed\n");
                  vigraImpexCloseImage(image);
                  return False;
              }
            p=image->pixels;
            for (i=0; i < image->packets; i++)
            {
              matte_image[i]=p->index == Transparent;
              if (p->index == Transparent)
                {
                  transparent_pixel.red=p->red;
                  transparent_pixel.green=p->green;
                  transparent_pixel.blue=p->blue;
                  transparent_pixel.flags=True;
                }
              p++;
            }
          }
        colors=transparent_pixel.flags ? 255 : 256;
        vigraImpexGetQuantizeInfo(&quantize_info);
        quantize_info.number_colors=colors;
        quantize_info.dither=image_info->dither;
        vigraImpexQuantizeImage(&quantize_info,image);
        vigraImpexSyncImage(image);
        colors=image->colors;
        if (transparent_pixel.flags)
          {
            /*
              Set the transparent pixel index.
            */
            image->c_class=VigraImpexDirectClass;
            if (!vigraImpexUncondenseImage(image))
              return(False);
            p=image->pixels;
            for (i=0; i < image->packets; i++)
            {
              if (matte_image[i])
                p->index=image->colors;
              p++;
            }
            colors++;
          }
        if (matte_image != (unsigned char *) NULL)
          free((char *) matte_image);
      }
    for (bits_per_pixel=1; bits_per_pixel < 8; bits_per_pixel++)
      if ((1 << bits_per_pixel) >= colors)
        break;
    q=colormap;
    for (i=0; i < image->colors; i++)
    {
      *q++=DownScale(image->colormap[i].red);
      *q++=DownScale(image->colormap[i].green);
      *q++=DownScale(image->colormap[i].blue);
    }
    if (transparent_pixel.flags)
      {
        *q++=DownScale(transparent_pixel.red);
        *q++=DownScale(transparent_pixel.green);
        *q++=DownScale(transparent_pixel.blue);
        i++;
      }
    for ( ; i < (int) (1 << bits_per_pixel); i++)
    {
      *q++=(Quantum) 0x0;
      *q++=(Quantum) 0x0;
      *q++=(Quantum) 0x0;
    }
    if (!image_info->adjoin || (image->previous == (VigraImpexImage *) NULL))
      {
        register int
          j;

        /*
          Write global colormap.
        */
        c=0x80;
        c|=(8-1) << 4;  /* color resolution */
        c|=(bits_per_pixel-1);   /* size of global colormap */
        (void) fputc((char) c,image->file);
        for (j=0; j < (image->colors-1); j++)
          if (ColorMatch(image->background_color,image->colormap[j],0))
            break;
        (void) fputc(j,image->file);  /* background color */
        (void) fputc(0x0,image->file);  /* reserved */
        (void) fwrite(colormap,1,3*(1 << bits_per_pixel),image->file);
      }
    if (strcmp(image_info->magick,"GIF87") != 0)
      {
        /*
          Write Graphics Control extension.
        */
        (void) fputc(0x21,image->file);
        (void) fputc(0xf9,image->file);
        (void) fputc(0x04,image->file);
        c=image->dispose << 2;
        if (transparent_pixel.flags)
          c|=0x01;
        (void) fputc(c,image->file);
        vigraImpexLSBFirstWriteShort(image->delay,image->file);
        (void) fputc((char) image->colors,image->file);
        (void) fputc(0x00,image->file);
        if (image->comments != (char *) NULL)
          {
            register char
              *p;

            register unsigned int
              count;

            /*
              Write Comment extension.
            */
            (void) fputc(0x21,image->file);
            (void) fputc(0xfe,image->file);
            p=image->comments;
            while (Extent(p) > 0)
            {
              count=Min(Extent(p),255);
              (void) fputc(count,image->file);
              for (i=0; i < count; i++)
                (void) fputc(*p++,image->file);
            }
            (void) fputc(0x0,image->file);
          }
        if ((image->previous == (VigraImpexImage *) NULL) &&
            (image->next != (VigraImpexImage *) NULL) && (image->iterations != 1))
          {
            /*
              Write Netscape Loop extension.
            */
            (void) fputc(0x21,image->file);
            (void) fputc(0xff,image->file);
            (void) fputc(0x0b,image->file);
            (void) fwrite("NETSCAPE2.0",1,11,image->file);
            (void) fputc(0x03,image->file);
            (void) fputc(0x01,image->file);
            vigraImpexLSBFirstWriteShort(image->iterations,image->file);
            (void) fputc(0x00,image->file);
          }
      }
    (void) fputc(',',image->file);  /* image separator */
    /*
      Write the image header.
    */

    vigraImpexLSBFirstWriteShort(page_info.x,image->file);
    vigraImpexLSBFirstWriteShort(page_info.y,image->file);
    vigraImpexLSBFirstWriteShort(image->columns,image->file);
    vigraImpexLSBFirstWriteShort(image->rows,image->file);
    c=0x00;
    if (interlace != VigraImpexNoInterlace)
      c|=0x40;  /* pixel data is interlaced */
    if (global_colormap)
      (void) fputc((char) c,image->file);
    else
      {
        c|=0x80;
        c|=(bits_per_pixel-1);   /* size of local colormap */
        (void) fputc((char) c,image->file);
        (void) fwrite(colormap,1,3*(1 << bits_per_pixel),image->file);
      }
    /*
      Write the image data.
    */
    c=Max(bits_per_pixel,2);
    (void) fputc((char) c,image->file);
    if (interlace == VigraImpexNoInterlace)
      status=vigraImpexGIFEncodeImage(image,Max(bits_per_pixel,2)+1);
    else
      {
        VigraImpexImage
          *interlaced_image;

        int
          pass,
          y;

        register VigraImpexRunlengthPacket
          *q;

        static int
          interlace_rate[4] = { 8, 8, 4, 2 },
          interlace_start[4] = { 0, 4, 2, 1 };

        /*
          Interlace image.
        */
        if (!vigraImpexUncondenseImage(image))
          return(False);
        image->orphan=True;
        interlaced_image=vigraImpexCloneImage(image,image->columns,image->rows,False);
        image->orphan=False;
        if (interlaced_image == (VigraImpexImage *) NULL)
          {
              fprintf(stderr, "vigraImpexWriteGIFImage(): Memory allocation failed\n");
              vigraImpexCloseImage(image);
              return False;
          }
        p=image->pixels;
        q=interlaced_image->pixels;
        for (pass=0; pass < 4; pass++)
        {
          y=interlace_start[pass];
          while (y < image->rows)
          {
            p=image->pixels+(y*image->columns);
            for (x=0; x < image->columns; x++)
            {
              *q=(*p);
              p++;
              q++;
            }
            y+=interlace_rate[pass];
          }
        }
        interlaced_image->file=image->file;
        status=vigraImpexGIFEncodeImage(interlaced_image,Max(bits_per_pixel,2)+1);
        interlaced_image->file=(FILE *) NULL;
        vigraImpexDestroyImage(interlaced_image);
      }
    if (status == False)
    {
      fprintf(stderr, "vigraImpexWriteGIFImage(): Memory allocation failed\n");
      vigraImpexCloseImage(image);
      return False;
    }
    (void) fputc(0x0,image->file);
    if (image->next == (VigraImpexImage *) NULL)
      break;
    image->next->file=image->file;
    image=image->next;
  } while (image_info->adjoin);
  (void) fputc(';',image->file); /* terminator */
  free((char *) colormap);
  vigraImpexCloseImage(image);
  return(True);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
+   R e a d G I F I m a g e                                                   %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Function vigraImpexReadGIFImage reads a Compuserve Graphics image file and returns it.
%  It allocates the memory necessary for the new VigraImpexImage structure and returns a
%  pointer to the new image.
%
%  The format of the vigraImpexReadGIFImage routine is:
%
%      image=vigraImpexReadGIFImage(image_info)
%
%  A description of each parameter follows:
%
%    o image:  Function vigraImpexReadGIFImage returns a pointer to the image after
%      reading.  A null image is returned if there is a a memory shortage or
%      an error occurs.
%
%    o image_info: Specifies a pointer to an VigraImpexImageInfo structure.
%
%
*/
Export VigraImpexImage *vigraImpexReadGIFImage(VigraImpexImageInfo *image_info)
{
#define BitSet(byte,bit)  (((byte) & (bit)) == (bit))
#define LSBFirstOrder(x,y)  (((y) << 8) | (x))

  char
    geometry[MaxTextExtent];

  VigraImpexImage
    *image;

  int
    status,
    x,
    y;

  VigraImpexRectangleInfo
    page_info;

  register int
    i;

  register VigraImpexRunlengthPacket
    *q;

  register unsigned char
    *p;

  short int
    transparency_index;

  unsigned char
    background,
    c,
    flag,
    *global_colormap,
    header[MaxTextExtent],
    magick[12];

  unsigned int
    delay,
    dispose,
    global_colors,
    image_count,
    iterations;

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
      fprintf(stderr, "vigraImpexReadGIFImage(): Unable to open file %s\n",image->filename);
      vigraImpexDestroyImage(image);
      return 0;
  }
  /*
    Determine if this is a GIF file.
  */
  status=vigraImpexReadData((char *) magick,1,6,image->file);
  if ((status == False) || ((strncmp((char *) magick,"GIF87",5) != 0) &&
      (strncmp((char *) magick,"GIF89",5) != 0)))
  {
      fprintf(stderr, "vigraImpexReadGIFImage(): %s is not a GIF image file\n",image->filename);
      vigraImpexDestroyImage(image);
      return 0;
  }
  global_colors=0;
  global_colormap=(unsigned char *) NULL;
  page_info.width=vigraImpexLSBFirstReadShort(image->file);
  page_info.height=vigraImpexLSBFirstReadShort(image->file);
  (void) vigraImpexReadData((char *) &flag,1,1,image->file);
  (void) vigraImpexReadData((char *) &background,1,1,image->file);
  (void) vigraImpexReadData((char *) &c,1,1,image->file);  /* reserved */
  if (BitSet(flag,0x80))
    {
      /*
        Read global colormap.
      */
      global_colors=1 << ((flag & 0x07)+1);
      global_colormap=(unsigned char *)
        malloc(3*global_colors*sizeof(unsigned char));
      if (global_colormap == (unsigned char *) NULL)
      {
          fprintf(stderr, "vigraImpexReadGIFImage(): Unable to read image colormap file\n");
          vigraImpexDestroyImage(image);
          return 0;
      }
      (void) vigraImpexReadData((char *) global_colormap,3,global_colors,image->file);
    }
  delay=0;
  dispose=0;
  iterations=1;
  transparency_index=(-1);
  image_count=0;
  for ( ; ; )
  {
    status=vigraImpexReadData((char *) &c,1,1,image->file);
    if (status == False)
      break;
    if (c == ';')
      break;  /* terminator */
    if (c == '!')
      {
        /*
          GIF Extension block.
        */
        status=vigraImpexReadData((char *) &c,1,1,image->file);
        if (status == False)
          {
              fprintf(stderr, "vigraImpexReadGIFImage(): Unable to read extention block\n");
              vigraImpexDestroyImage(image);
              return 0;
          }
        switch (c)
        {
          case 0xf9:
          {
            /*
              Read Graphics Control extension.
            */
            while (vigraImpexReadDataBlock((char *) header,image->file) > 0);
            dispose=header[0] >> 2;
            delay=(header[2] << 8) | header[1];
            if ((header[0] & 0x01) == 1)
              transparency_index=header[3];
            break;
          }
          case 0xfe:
          {
            int
              length;

            /*
              Read Comment extension.
            */
            for ( ; ; )
            {
              length=vigraImpexReadDataBlock((char *) header,image->file);
              if (length <= 0)
                break;
              if (image->comments != (char *) NULL)
                image->comments=(char *) realloc((char *) image->comments,
                  (Extent(image->comments)+length+1)*sizeof(char));
              else
                {
                  image->comments=(char *) malloc((length+1)*sizeof(char));
                  if (image->comments != (char *) NULL)
                    *image->comments='\0';
                }
              if (image->comments == (char *) NULL)
              {
                  fprintf(stderr, "vigraImpexReadGIFImage(): Memory allocation failed\n");
                  vigraImpexDestroyImage(image);
                  return 0;
              }
              header[length]='\0';
              (void) strcat(image->comments,(char *) header);
            }
            break;
          }
          case 0xff:
          {
            /*
              Read Netscape Loop extension.
            */
            while (vigraImpexReadDataBlock((char *) header,image->file) > 0);
            iterations=(header[2] << 8) | header[1];
            break;
          }
          default:
          {
            while (vigraImpexReadDataBlock((char *) header,image->file) > 0);
            break;
          }
        }
      }
    if (c != ',')
      continue;
    if (image_count != 0)
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
    image_count++;
    /*
      Read image attributes.
    */
    image->c_class=VigraImpexPseudoClass;
    page_info.x=vigraImpexLSBFirstReadShort(image->file);
    page_info.y=vigraImpexLSBFirstReadShort(image->file);
    image->columns=vigraImpexLSBFirstReadShort(image->file);
    image->rows=vigraImpexLSBFirstReadShort(image->file);
    (void) vigraImpexReadData((char *) &flag,1,1,image->file);
    image->interlace=BitSet(flag,0x40) ? VigraImpexPlaneInterlace : VigraImpexNoInterlace;
    image->colors=!BitSet(flag,0x80) ? global_colors : 1 << ((flag & 0x07)+1);
    (void) sprintf(geometry,"%ux%u%+d%+d",page_info.width,page_info.height,
      page_info.x,page_info.y);
    image->page=vigraImpexPostscriptGeometry(geometry);
    image->delay=delay;
    image->dispose=dispose;
    image->iterations=iterations;
    delay=0;
    dispose=0;
    iterations=1;
    if (image_info->ping)
      {
        if (transparency_index >= 0)
          image->c_class=VigraImpexDirectClass;
        vigraImpexCloseImage(image);
        return(image);
      }
    image->packets=image->columns*image->rows;
    if (image->packets == 0)
    {
      fprintf(stderr, "vigraImpexReadGIFImage(): Size of image %s is 0\n",image->filename);
      vigraImpexDestroyImage(image);
      return 0;
    }
    /*
      Allocate image.
    */
    if (image->pixels != (VigraImpexRunlengthPacket *) NULL)
      free((char *) image->pixels);
    image->pixels=(VigraImpexRunlengthPacket *)
      malloc(image->packets*sizeof(VigraImpexRunlengthPacket));
    if (image->pixels == (VigraImpexRunlengthPacket *) NULL)
    {
      fprintf(stderr, "vigraImpexReadGIFImage(): Memory allocation failed");
      vigraImpexDestroyImage(image);
      return 0;
    }
    /*
      Inititialize colormap.
    */
    image->colormap=(VigraImpexColorPacket *) malloc(image->colors*sizeof(VigraImpexColorPacket));
    if (image->colormap == (VigraImpexColorPacket *) NULL)
    {
      fprintf(stderr, "vigraImpexReadGIFImage(): Memory allocation failed");
      vigraImpexDestroyImage(image);
      return 0;
    }
    if (!BitSet(flag,0x80))
      {
        /*
          Use global colormap.
        */
        p=global_colormap;
        for (i=0; i < image->colors; i++)
        {
          image->colormap[i].red=UpScale(*p++);
          image->colormap[i].green=UpScale(*p++);
          image->colormap[i].blue=UpScale(*p++);
        }
        image->background_color=
          image->colormap[Min(background,image->colors-1)];
      }
    else
      {
        unsigned char
          *colormap;

        /*
          Read local colormap.
        */
        colormap=(unsigned char *)
          malloc(3*image->colors*sizeof(unsigned char));
        if (colormap == (unsigned char *) NULL)
        {
          fprintf(stderr, "vigraImpexReadGIFImage(): Memory allocation failed");
          vigraImpexDestroyImage(image);
          return 0;
        }
        (void) vigraImpexReadData((char *) colormap,3,image->colors,image->file);
        p=colormap;
        for (i=0; i < image->colors; i++)
        {
          image->colormap[i].red=UpScale(*p++);
          image->colormap[i].green=UpScale(*p++);
          image->colormap[i].blue=UpScale(*p++);
        }
        free((char *) colormap);
      }
    /*
      Decode image.
    */
    status=vigraImpexGIFDecodeImage(image);
    if (image->interlace != VigraImpexNoInterlace)
      {
        VigraImpexImage
          *interlaced_image;

        int
          pass;

        register VigraImpexRunlengthPacket
          *p;

        static int
          interlace_rate[4] = { 8, 8, 4, 2 },
          interlace_start[4] = { 0, 4, 2, 1 };

        /*
          Interlace image.
        */
        image_info->interlace=VigraImpexLineInterlace;
        image->orphan=True;
        interlaced_image=vigraImpexCloneImage(image,image->columns,image->rows,True);
        image->orphan=False;
        if (interlaced_image == (VigraImpexImage *) NULL)
        {
          fprintf(stderr, "vigraImpexReadGIFImage(): Memory allocation failed");
          vigraImpexDestroyImage(image);
          return 0;
        }
        p=interlaced_image->pixels;
        q=image->pixels;
        for (pass=0; pass < 4; pass++)
        {
          y=interlace_start[pass];
          while (y < image->rows)
          {
            q=image->pixels+(y*image->columns);
            for (x=0; x < image->columns; x++)
            {
              *q=(*p);
              p++;
              q++;
            }
            y+=interlace_rate[pass];
          }
        }
        vigraImpexDestroyImage(interlaced_image);
      }
    if (transparency_index >= 0)
      {
        /*
          Create matte channel.
        */
        q=image->pixels;
        for (i=0; i < image->packets; i++)
        {
          if (q->index != (unsigned short) transparency_index)
            q->index=Opaque;
          else
            q->index=Transparent;
          q++;
        }
        transparency_index=(-1);
        image->c_class=VigraImpexDirectClass;
        image->matte=True;
      }
    if (status == False)
      {
        break;
      }
    if (image_info->subrange != 0)
      {
          if (image->scene < image_info->subimage)
            {
              VigraImpexImage
                subimage;

              /*
                Destroy image.
              */
              subimage=(*image);
              image->file=(FILE *) NULL;
              vigraImpexDestroyImage(image);
              image=vigraImpexAllocateImage(image_info);
              if (image == (VigraImpexImage *) NULL)
                return((VigraImpexImage *) NULL);
              image->file=subimage.file;
              image->scene=subimage.scene+1;
              image_count=0;
            }
          else
            if (image->scene >= (image_info->subimage+image_info->subrange-1))
              break;
      }
  }
  if (global_colormap != (unsigned char *) NULL)
    free((char *) global_colormap);
  if (image->pixels == (VigraImpexRunlengthPacket *) NULL)
    {
      fprintf(stderr, "vigraImpexReadGIFImage(): Corrupt GIF image or subimage (%s)",
                      image->filename);
      vigraImpexDestroyImage(image);
      return 0;
    }
  vigraImpexCondenseImage(image);
  while (image->previous != (VigraImpexImage *) NULL)
    image=image->previous;
  vigraImpexCloseImage(image);
  return(image);
}
