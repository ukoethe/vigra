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
#include <stdarg.h>
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
#include "pnm.h"
#include "utility.h" 

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                             %
%                                                                             %
%                                                                             %
%  F o r m a t S t r i n g                                                    %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Method FormatString prints formatted output of a variable argument list.
%
%  The format of the FormatString method is:
%
%      void FormatString(char *string,const char *format,...)
%
%  A description of each parameter follows.
%
%   o  string:  Method FormatString returns the formatted string in this
%      character buffer.
%
%   o  format:  A string describing the format to use to write the remaining
%      arguments.
%
%
*/
void FormatString(char *string,const char *format, ...)
{
  va_list
    operands;

  va_start(operands,format);
#if !defined(HAVE_VSNPRINTF)
  (void) vsprintf(string,format,operands);
#else
  (void) vsnprintf(string,MaxTextExtent,format,operands);
#endif
  va_end(operands);
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%
%
%
%
%
%
%   R e a d P N M I m a g e
%
%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%
%  Method ReadPNMImage reads a Portable Anymap image file and returns it.
%  It allocates the memory necessary for the new Image structure and returns
%  a pointer to the new image.
%
%  The format of the ReadPNMImage method is:
%
%      Image *ReadPNMImage(const ImageInfo *image_info,ExceptionInfo
*exception)
%
%  A description of each parameter follows:
%
%    o image:  Method ReadPNMImage returns a pointer to the image after
%      reading.  A null image is returned if there is a memory shortage or
%      if the image cannot be read.
%
%    o image_info: Specifies a pointer to an ImageInfo structure.
%
%    o exception: return any errors or warnings in this structure.
%
%
*/


static unsigned int PNMInteger(VigraImpexImage *image,const unsigned int base)
{
#define P7Comment  "END_OF_COMMENTS\n"


  int
    c;


  unsigned int
    value;

  /*
    Skip any leading whitespace.
  */
  do
  {
    c=fgetc(image->file);

    if (c == EOF )
      return(0);
    if (c == '#')
      {
        char
          *comment;


        register char
          *p,
          *q;


        unsigned int
          length,
          offset;


        /*
          Read comment.
        */
        length=MaxTextExtent;
        comment=(char *) malloc(length+strlen(P7Comment)+1);
        p=comment;
        offset=p-comment;
        if (comment != (char *) NULL)
          for ( ; (c != EOF) && (c != '\n'); p++)
          {
            if ((p-comment) >= length)
              {
                length<<=1;
                length+=MaxTextExtent;
              
               /*MIKHAIL ReallocateMemory((void **)&comment,length+strlen(P7Comment)+1); */
                                
                comment=(char *) malloc(length+strlen(P7Comment)+1);


                if (comment == (char *) NULL)
                  break;
                p=comment+Extent(comment);
              }


            c=fgetc(image->file);            

            *p=c;
            *(p+1)='\0';
          }
        if (comment == (char *) NULL)
          return(0);
        q=comment+offset;
        if (strcmp(q,P7Comment) == 0)
          *q='\0';
        
		/*MIKHAIL (void) SetImageAttribute(image,"Comment",comment); */
        
        continue;
      }
  } while (!isdigit(c));
  if (base == 2)
    return((char)c-'0');
  /*
    Evaluate number.
  */
  
  value=0;
  do
  {
    value*=10;
    value+=(char)c-'0';
            
    c=fgetc(image->file);
        
    if (c == EOF)
      return(0);
  }
  while (isdigit(c));
  return(value);
}


VigraImpexImage *vigraImpexReadPNMImage( VigraImpexImageInfo *image_info)
{
  char
    format;


  VigraImpexImage
    *image;

  int
    y;


  register int
    i,
    x;


  register VigraImpexRunlengthPacket
    *q;


  register unsigned char
    *p;

  Quantum
    blue,
    green,
    red,
    *scale;

  unsigned char
    *pixels,
	magick[12];
	

  unsigned int
    max_value,
    packets,
    status;

  /*
    Open image file.
  */  
  image=vigraImpexAllocateImage(image_info);  
  if (image == (VigraImpexImage *) NULL)
    return((VigraImpexImage *) NULL);
  
  vigraImpexOpenImage(image_info,image,ReadBinaryType);
    
     
  if (image->file==(FILE *)NULL)
  {
          fprintf(stderr,"vigraImpexReadPNMImage(): Unable to open file %s \n", image->filename);
          vigraImpexDestroyImage(image);
          return 0;
  }
  /*
        Determine if this is a PNM file.
  */  
  /*MIKHAIL  consider reading into image->magick */
  status=vigraImpexReadData((char *)magick,1,1,image->file);
  status=vigraImpexReadData((char *)&format,1,1,image->file);                 
  magick[1]=format;  
  strcpy(image->magick,(char *)magick);
  do
  {
/*MIKHAIL       Verify PNM identifier. */
    /*
    */
        if ((status == False) || (strncmp((char *)magick,"P",1) != 0))
        {
          fprintf(stderr,
			  "vigraImpexReadPNMImage():%s isn't a PNM image file\n",
			  image->filename);
          vigraImpexDestroyImage(image);
          return 0;
        }
              
    /*
      Initialize image structure.
          format=ReadByte(image);
    */
		/*MIKHAIL (char *)		 */
    if (format=='7')
      (void) PNMInteger(image,10);


    image->columns=PNMInteger(image,10);
    image->rows=PNMInteger(image,10);
        
    if ((format == '1') || (format == '4'))
      max_value=1;  /* bitmap */
    else
      max_value=PNMInteger(image,10);
	
	/*MIKHAIL image->depth=max_value < 256 ? 8 : 16;		//QuantumDepth; */
	/*MIKHAIL #define DownScale(quantum)  (((unsigned int) (quantum)) >> 8) */
	/*MIKHAIL #define HexColorFormat "#%04x%04x%04x" */
	/*MIKHAIL #define MaxRGB  65535L */
	/*MIKHAIL #define MaxRunlength  65535L */
	/*MIKHAIL #define Opaque  65535L */
	/*MIKHAIL #define QuantumDepth  16 */
	/*MIKHAIL #define UpScale(quantum)  (((unsigned int) (quantum))*257) */
	/*MIKHAIL #define XDownScale(color)  ((unsigned int) (color)) */
	/*MIKHAIL #define XUpScale(color)  ((unsigned int) (color)) */

	/*if (maxvalue<256)
		image->depth=8;		
	else if (maxvalue<4096)
		image->depth=12;
	else if (maxvalue<16384)
		image->depth=14;
	else if (maxvalue<65536)
		image->depth=16;
	*/	
	image->depth=(max_value<256) ? 8 : QuantumDepth;

    
    if ((format != '3') && (format != '6'))
      {
        image->c_class=VigraImpexPseudoClass;        
        image->colors=max_value+1;
      }
    if (image_info->ping)
      {                 
            while(image->previous!=(VigraImpexImage *)NULL)
                    image=image->previous;
            vigraImpexCloseImage(image);
			return(image);
      }
    if ((image->columns*image->rows) == 0)
        {
		  fprintf(stderr,"vigraImpexReadPNMImage():Unable to read %s: image dimensions are zero \n", image->filename);
          vigraImpexDestroyImage(image);
          return 0;
        }
                
    scale=(Quantum *) NULL;
    if (image->c_class == VigraImpexPseudoClass)
      {
        /*
          Create colormap.
        */
/*MIKHAIL         if (!AllocateImageColormap(image,image->colors)) */
/*MIKHAIL           ThrowReaderException(ResourceLimitWarning,"Memory allocation failed", */
/*MIKHAIL             image); */
                
        image->colormap=(VigraImpexColorPacket *)
	          malloc(image->colors*sizeof(VigraImpexColorPacket));
        if (image->colormap == (VigraImpexColorPacket *) NULL)
        {
			fprintf(stderr, "vigraImpexReadPNMImage(): Memory allocation failed\n");
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

        if (format == '7')
          {
            /*
              Initialize 332 colormap.
            */
            i=0;
            for (red=0; red < 8; red++)
              for (green=0; green < 8; green++)
                for (blue=0; blue < 4; blue++)
                {
                  image->colormap[i].red=((unsigned long) (MaxRGB*red)/0x07);
                  image->colormap[i].green=
                    ((unsigned long) (green*MaxRGB)/0x07);
                  image->colormap[i].blue=((unsigned long) (MaxRGB*blue)/0x03);
                  i++;
                }
          }
      }
    else
      if (max_value != MaxRGB)
        {
          /*
            Compute pixel scaling table.
          */
          scale=(Quantum *)malloc((max_value+1)*sizeof(Quantum));
          if (scale == (Quantum *) NULL)
                  {
          fprintf(stderr, "vigraImpexReadPNMImage(): Memory allocation failed\n");
          vigraImpexDestroyImage(image);
          return 0;
                  }
          for (i=0; i <= (int) max_value; i++)
            scale[i]=(Quantum) ((i*MaxRGB+(max_value >> 1))/max_value);
        }
    /*
      Convert PNM pixels to runlength-encoded MIFF packets.
    */

    image->units=VigraImpexPixelsPerCentimeterResolution;
    image->packets=image->columns*image->rows;
    image->pixels=(VigraImpexRunlengthPacket *)
      malloc(image->packets*sizeof(VigraImpexRunlengthPacket));     
          
        switch (format)
    {
      case '1':
      {
        /*
          Convert PBM image to pixel packets.
        */
		image->colors=2;
		q=image->pixels;/*MIKHAIL +(y*image->columns);              */
        for (y=0; y < (int) image->rows; y++)
        {                           
		  
          for (x=0; x < (int) image->columns; x++)
          {
			/*MIKHAIL byte|=0x01; */
            /*MIKHAIL q->index=!PNMInteger(image,2); */
			q->index=(!PNMInteger(image,2)) ? 0x01:0x00;
			q->length=0x00;
			q++;			
            /*MIKHAIL q++; */
          }
/*MIKHAIL           if (!SyncImagePixels(image)) */
/*MIKHAIL             break; */
/*MIKHAIL           if (image->previous == (Image *) NULL) */
/*MIKHAIL             if (QuantumTick(y,image->rows)) */
/*MIKHAIL               ProgressMonitor(LoadImageText,y,image->rows); */
        }
		break;
      }
      case '2':
      {
        /*
          Convert PGM image to pixel packets.
        */
		q=image->pixels;/*MIKHAIL +(y*image->columns);              */
        for (y=0; y < (int) image->rows; y++)
        {          
          for (x=0; x < (int) image->columns; x++)
          {
            q->index=PNMInteger(image,10);
            q->length=0;
            q++;
          }
/*MIKHAIL           if (!SyncImagePixels(image)) */
/*MIKHAIL             break; */
/*MIKHAIL           if (image->previous == (VigraImpexImage *) NULL) */
/*MIKHAIL             if (QuantumTick(y,image->rows)) */
/*MIKHAIL               ProgressMonitor(LoadImageText,y,image->rows); */
        }
        break;
      }
      case '3':
      {
        /*
          Convert PNM image to pixel packets.
        */
        for (y=0; y < (int) image->rows; y++)
        {
          q=image->pixels+(y*image->columns);             
          for (x=0; x < (int) image->columns; x++)
          {
            red=PNMInteger(image,10);
            green=PNMInteger(image,10);
            blue=PNMInteger(image,10);
            if (scale != (Quantum *) NULL)
              {
                red=scale[red];
                green=scale[green];
                blue=scale[blue];
              }			
            q->red=red;
            q->green=green;
            q->blue=blue;
			q->length=0;
            q++;
          }
/*MIKHAIL           if (!SyncImagePixels(image)) */
/*MIKHAIL             break; */
/*MIKHAIL         if (image->previous == (Image *) NULL) */
/*MIKHAIL           if (QuantumTick(y,image->rows)) */
/*MIKHAIL             ProgressMonitor(LoadImageText,y,image->rows); */
        }
        break;
      }
      case '4':
      {
        unsigned char
          bit,
          byte;


        /*
          Convert PBM raw image to pixel packets.
        */
        bit=0;
        byte=0;

        q=image->pixels;/*MIKHAIL +(y*image->columns); */
        for (y=0; y < (int) image->rows; y++)
        {
          for (x=0; x < (int) image->columns; x++)
          {				  
            if (bit == 0)
			byte=fgetc(image->file);
			if (byte==0)
				q->index=1;
			else
				q->index=(byte & 0x80) ? 0 : 1;
            q->length=0;
            q++;
            bit++;
            if (bit == 8)
              bit=0;
            byte<<=1;
          }
/*MIKHAIL           if (!SyncImagePixels(image)) */
/*MIKHAIL             break; */
/*MIKHAIL           if (image->previous == (Image *) NULL) */
/*MIKHAIL             if (QuantumTick(y,image->rows)) */
/*MIKHAIL               ProgressMonitor(LoadImageText,y,image->rows); */
        }
        break;
      }
      case '5':
      case '7':
      {
        /*
          Convert PGM raw image to pixel packets.
        */
        packets=image->depth <= 8 ? 1 : 2;
        pixels=(unsigned char *) malloc(packets*image->columns);
        if (pixels == (unsigned char *) NULL)
		{
		  fprintf(stderr, "vigraImpexReadPNMImage(): Memory allocation failed\n");
          vigraImpexDestroyImage(image);
          return 0;
		}	
        for (y=0; y < (int) image->rows; y++)
        {			
		  vigraImpexReadData((char *)pixels,1,packets*image->columns,image->file);  
		  /*MIKHAIL status=ReadBlob(image,packets*image->columns,pixels); */
          if (status == False)
			{
			    fprintf(stderr, "vigraImpexReadPNMImage(): Unable to read image.\n");
				vigraImpexDestroyImage(image);
				return 0;
			}	

			  
		  p=pixels;
          q=image->pixels+(y*image->columns);             
          for (x=0; x < (int) image->columns; x++)
          {
            if (image->depth <= 8)
              q->index=(*p++);
            else
              {
                q->index=(*p++) << 8;
                q->index|=(*p++);
              }


            if (q->index >= image->colors)
			{
			    fprintf(stderr, "vigraImpexReadPNMImage(): Invalid colormap. \n");
				vigraImpexDestroyImage(image);
				return 0;
			}	
				
            q->length=0;
            q++;
          }
/*MIKHAIL           if (!SyncImagePixels(image)) */
/*MIKHAIL             break; */
/*MIKHAIL           if (image->previous == (VigraImpexImage *) NULL) */
/*MIKHAIL             if (QuantumTick(y,image->rows)) */
/*MIKHAIL               ProgressMonitor(LoadImageText,y,image->rows); */
        }        
        break;
      }
      case '6':
      {
        /*
          Convert PNM raster image to pixel packets.
        */
        packets=image->depth <= 8 ? 3 : 6;
        pixels=(unsigned char *)malloc(packets*image->columns);
        if (pixels == (unsigned char *) NULL)
			{
			    fprintf(stderr, "vigraImpexReadPNMImage(): Unable to allocate memory.\n");
				vigraImpexDestroyImage(image);
				return 0;
			}	
        for (y=0; y < (int) image->rows; y++)
        {
          
		  vigraImpexReadData((char *)pixels,1,packets*image->columns,image->file);  
		  /*MIKHAIL  is it really packets? */
		  /*MIKHAIL status=ReadBlob(image,packets*image->columns,pixels); */
          if (status == False)
			{
			    fprintf(stderr, "vigraImpexReadPNMImage(): Unable to read image.\n");
				vigraImpexDestroyImage(image);
				return 0;
			}	
		  p=pixels;
          q=image->pixels+(y*image->columns);		            
          for (x=0; x < (int) image->columns; x++)
          {
            if (image->depth <= 8)
              {
				/*MIKHAIL  is this RGB  ??? */
				
				blue=(*p++);				
				red=(*p++);
				green=(*p++);				
              }
            else
              {
                blue=(*p++) << 8;
                blue|=(*p++);
                red=(*p++) << 8;
                red|=(*p++);
                green=(*p++) << 8;
                green|=(*p++);
              }
            if (scale != (Quantum *) NULL)
              {
                red=scale[red];
                green=scale[green];
                blue=scale[blue];
              }
            q->red=red;
			q->green=green;
            q->blue=blue;
			q->length=0;
            q++;
          }
/*MIKHAIL           if (!SyncImagePixels(image)) */
/*MIKHAIL             break; */
/*MIKHAIL           if (image->previous == (Image *) NULL) */
/*MIKHAIL             if (QuantumTick(y,image->rows)) */
/*MIKHAIL               ProgressMonitor(LoadImageText,y,image->rows); */
        }        
/*MIKHAIL         handler=SetMonitorHandler((MonitorHandler) NULL); */
/*MIKHAIL         (void) SetMonitorHandler(handler); */
        break;
      }
      default:
			{
			    fprintf(stderr, "vigraImpexReadPNMImage(): Not a PNM image file.\n");
				vigraImpexDestroyImage(image);
				return 0;
			}	
		  
    }
    if (scale != (Quantum *) NULL)

		
/*MIKHAIL PUT THIS IN LATER */
/*MIKHAIL     if (EOFBlob(image)) */
/*MIKHAIL 			{ */
/*MIKHAIL 			    fprintf(stderr, "vigraImpexReadPNMImage():not enough pixels.\n"); */
/*MIKHAIL 				vigraImpexDestroyImage(image); */
/*MIKHAIL 				return 0; */
/*MIKHAIL 			}			 */
    /*
      Proceed to next image.
    */
    if (image_info->subrange != 0)
      if (image->scene >= (image_info->subimage+image_info->subrange-1))
        break;
    if ((format == '1') || (format == '2') || (format == '3'))
      do
      {
        /*
          Skip to end of line.
        */
		status=vigraImpexReadData((char *)&format,1,1,image->file);  
		/*MIKHAIL status=ReadBlob(image,1,&format); */
        if (status == False)
          break;
      } while (format != '\n' || format != EOF);
	
	  status=vigraImpexReadData((char *) &format,1,1,image->file);  		  
	  /*MIKHAIL status=ReadBlob(image,1,(char *) &format); */
    if ((status == True) && (format == 'P'))
      {
        /*
          Allocate next image structure.
        */
		image->next=vigraImpexAllocateImage(image_info);          
        if (image->next == (VigraImpexImage *) NULL)
          {
			vigraImpexDestroyImage(image);
            return((VigraImpexImage *) NULL);
          }
        image=image->next;
/*MIKHAIL         ProgressMonitor(LoadImagesText,TellBlob(image),image->filesize); */
      }
  } while ((status == True) && (format == 'P'));
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
%   W r i t e P N M I m a g e                                                 %
%                                                                             %
%                                                                             %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Procedure WritePNMImage writes an image to a file in the PNM rasterfile
%  format.
%
%  The format of the WritePNMImage method is:
%
%      unsigned int WritePNMImage(const ImageInfo *image_info,Image *image)
%
%  A description of each parameter follows.
%
%    o status: Method WritePNMImage return True if the image is written.
%      False is returned is there is a memory shortage or if the image file
%      fails to write.
%
%    o image_info: Specifies a pointer to an ImageInfo structure.
%
%    o image:  A pointer to a Image structure.
%
%
*/
/*MIKHAIL static unsigned int WritePNMImage(const ImageInfo *image_info,Image *image) */

unsigned int vigraImpexWritePNMImage(VigraImpexImageInfo *image_info,VigraImpexImage *image)
{
  char
    buffer[MaxTextExtent],
    *magick;

  int
    j,
    y;

  register int
    i,
    x;

  register VigraImpexRunlengthPacket
    *p;

  unsigned char
    format;

  unsigned int
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
      Promote/Demote image based on image type.
    */	 
    (void) vigraImpexIsPseudoClass(image);
    
	/*MIKHAIL if (LocaleCompare(image_info->magick,"PPM") == 0) */
	if (strcmp(image_info->magick,"PPM") == 0 )
		image->c_class=VigraImpexDirectClass;
	      
    magick=(char *) image_info->magick;

	if ((strcmp(image_info->magick,"PGM") == 0  && !vigraImpexIsGrayImage(image)) ||
		(strcmp(image_info->magick,"PBM") == 0  && !vigraImpexIsMonochromeImage(image)))

      {
        VigraImpexQuantizeInfo quantize_info;

        vigraImpexGetQuantizeInfo(&quantize_info);
        quantize_info.number_colors=MaxRGB+1;
        
		if (strcmp(image_info->magick,"PBM") == 0)
			quantize_info.number_colors=2;
        quantize_info.dither=image_info->dither;
        quantize_info.colorspace=VigraImpexGRAYColorspace;
        (void) vigraImpexQuantizeImage(&quantize_info,image);
      }	


    /*
      Write PNM file header.
    */
	image_info->compression = VigraImpexNoCompression;
    if (!vigraImpexIsPseudoClass(image) && !vigraImpexIsGrayImage(image))
      {
        /*
          Full color PNM image.
        */
        format='6';
        if ((image_info->compression == VigraImpexNoCompression) || (image->depth > 8))
          format='3';
      }
    else
      {
        /*
          Colormapped PNM image.
        */
        format='6';

        if ((image_info->compression == VigraImpexNoCompression) || (image->depth > 8))
          format='3';
        
		/*MIKHAIL if ((LocaleCompare(magick,"PPM") != 0) && IsGrayImage(image)) */
		if ((strcmp(image_info->magick,"PPM") == 0) && vigraImpexIsGrayImage(image))
        
          {
            /*
              Grayscale PNM image.
            */
            format='5';
            if ((image_info->compression == VigraImpexNoCompression) ||
                (image->depth > 8))
              format='2';
            if (strcmp(image_info->magick,"PGM") == 0 ) 			        			
              if (image->colors == 2)
                {
                  format='4';
                  if (image_info->compression == VigraImpexNoCompression)
                    format='1';
                }
          }
		if (strcmp(image_info->magick,"PBM") == 0  && vigraImpexIsGrayImage(image))
                {
			if ((image_info->compression == VigraImpexNoCompression) || (image->depth > 8))
				format='1';
			else
				format='4';
		}
                if (strcmp(image_info->magick,"PGM") == 0  && vigraImpexIsGrayImage(image))
		{
                	if ((image_info->compression == VigraImpexNoCompression) || (image->depth > 8))
				format='2';
			else
				format='5';
                }

      }
    /*MIKHAIL if (LocaleCompare(magick,"P7") == 0) */
	if (strcmp(image_info->magick,"P7") == 0)
      {
        format='7';
        (void) strcpy(buffer,"P7 332\n");
      }
    else
      FormatString(buffer,"P%c\n",format);
    /*MIKHAIL (void) WriteBlob(image,strlen(buffer),buffer); */
	(void) fwrite(buffer,1,strlen(buffer),image->file);
    
	/*MIKHAIL GetImageAttribute(image,"Label"); */
	/*MIKHAIL attribute=image->label; */
			
    if (image->label != (char *) NULL)
      {
        register char
          *p;

        /*
          Write comments to file.
        */
        /*MIKHAIL (void) WriteByte(image,'#'); */
		(void) fputc('#',image->file);
        /*MIKHAIL for (p=attribute->value; *p != '\0'; p++) */
		for (p=image->label; *p != '\0'; p++)
        {
          /*MIKHAIL (void) WriteByte(image,*p); */
		(void) fputc(*p,image->file);
          if ((*p == '\n') && (*(p+1) != '\0'))
            /*MIKHAIL (void) WriteByte(image,'#'); */
			(void) fputc('#',image->file);
        }
        /*MIKHAIL (void) WriteByte(image,'\n'); */
		(void) fputc('\n',image->file); \
      }
    if (format != '7')
      {
        FormatString(buffer,"%u %u\n",image->columns,image->rows);
        /*MIKHAIL (void) WriteBlob(image,strlen(buffer),buffer); */
		(void) fwrite(buffer,1,strlen(buffer),image->file);
      }
    /*
      Convert runlength encoded to PNM raster pixels.
    */
/*MIKHAIL 	(void) vigraImpexUncondenseImage(image); */

	switch (format)
    {
      case '1':
      {
        register unsigned char
          polarity;

        /*
          Convert image to a PBM image.
        */
        polarity=Intensity(image->colormap[0]) > (MaxRGB >> 1);		
        if (image->colors == 2)          
			polarity=
				Intensity(image->colormap[0]) > Intensity(image->colormap[1]);
        i=0;
		/*MIKHAIL p=GetImagePixels(image,0,y,image->columns,1); */
		p=image->pixels;
        if (p == (VigraImpexRunlengthPacket *) NULL)
            break;
         /*MIKHAIL indexes=GetIndexes(image); */
	    /*MIKHAIL q=image->pixels+(y*image->columns);		 		             */
        for (y=0; y < (int) image->packets; y++)
        {
		  for (x=0; x <= (int) p->length; x++)
          {
            
            /*MIKHAIL oid) WriteBlob(image,strlen(buffer),buffer); */
			/*MIKHAIL FormatString(buffer,"%d ",(int) (p->index == polarity)); */
			
			  (void) fputc(((p->index==polarity) ? '1' : '0'),image->file);
			  (void) fputc(' ',image->file);
			/*MIKHAIL (void) fwrite(p->index,1,strlen(buffer),image->file) */

	        /*MIKHAIL FormatString(buffer,"%lu\n",DownScale(MaxRGB)); */
		    /*MIKHAIL (void) WriteBlob(image,strlen(buffer),buffer); */
			/*MIKHAIL (void) fwrite(buffer,1,strlen(buffer),image->file); */

            i++;
            if (i == 36)
              {
				/*MIKHAIL (void) WriteByte(image,'\n'); */
			    (void) fputc(13,image->file);
                i=0;
              }
          }
/*MIKHAIL           if (image->previous == (VigraImpexImage *) NULL) */
/*MIKHAIL             if (QuantumTick(y,image->rows)) */
/*MIKHAIL               ProgressMonitor(SaveImageText,y,image->rows); */
		  p++;
        }
        break;
      }
      case '2':
      {
        /*
          Convert image to a PGM image.
        */
		FormatString(buffer,"%ld\n",MaxRGB);

		
		
		/*MIKHAIL (void) WriteBlob(image,2,"BM"); */

        /*(void) fwrite("BM",1,2,image->file);
		  (void) WriteBlob(image,strlen(buffer),buffer);
		*/
		(void) fwrite(buffer,1,strlen(buffer),image->file);


        i=0;		
		p=image->pixels;/*MIKHAIL ;+(y*image->columns); */
        if (p == (VigraImpexRunlengthPacket *) NULL)
            break;
        for (y=0; y < (int) image->packets; y++)
        {
        /*MIKHAIL 	p=GetImagePixels(image,0,y,image->columns,1); */
		  p->index=Intensity(*p);
		  FormatString(buffer,"%d ",p->index);            
          for (x=0; x <= (int) p->length; x++)
          {          
			/*MIKHAIL (void) WriteBlob(image,strlen(buffer),buffer);			   */
			(void) fwrite(buffer,1,strlen(buffer),image->file);

            i++;
            if (i == 12)
              {
                /*MIKHAIL (void) WriteByte(image,'\n'); */
				(void) fputc('\n',image->file); \
                i=0;
              }
		  
		  }
		  p++; 
/*MIKHAIL           if (image->previous == (Image *) NULL) */
/*MIKHAIL             if (QuantumTick(y,image->rows)) */
/*MIKHAIL               ProgressMonitor(SaveImageText,y,image->rows); */
        }
        break;
      }
      case '3':
      {
        /*
          Convert image to a PNM image.
        */
        FormatString(buffer,"%ld\n",MaxRGB);
        (void) fwrite(buffer,1,strlen(buffer),image->file);
		/*MIKHAIL (void) WriteBlob(image,strlen(buffer),buffer); */
        i=0;
		p=image->pixels;
		for (y=0; y < (int) image->packets; y++)
		{
          /*MIKHAIL p=GetImagePixels(image,0,y,image->columns,1);			 */
          /*MIKHAIL indexes=GetIndexes(image); */
		  if (p == (VigraImpexRunlengthPacket *) NULL)
			    break;
          FormatString(buffer,"%d %d %d ",p->red,p->green,p->blue);
		  for (x=0; x <= (int) p->length; x++)
          {            
            /*MIKHAIL (void) WriteBlob(image,strlen(buffer),buffer); */
			(void) fwrite(buffer,1,strlen(buffer),image->file);
            i++;
            if (i == 4)
              {
                /*MIKHAIL (void) WriteByte(image,'\n'); */
				(void) fputc('\n',image->file);
                i=0;
              }

          }
		  p++;

		}
        break;
      }
      case '4':
      {
        register unsigned char
          byte, polarity, bit;

        /*
          Convert image to a PBM image.
        */
        polarity=Intensity(image->colormap[0]) > (MaxRGB >> 1);
        if (image->colors == 2)
			polarity= Intensity(image->colormap[0]) > Intensity(image->colormap[1]);
        p=image->pixels;
		if (p == (VigraImpexRunlengthPacket *) NULL)
            break;
		bit=0;
		byte=0;

		for (y=0; y < (int) image->packets; y++)
		{
          /*MIKHAIL p=GetImagePixels(image,0,y,image->columns,1);			 */
          /*MIKHAIL indexes=GetIndexes(image); */
          for (x=0; x <= (int) p->length; x++)
          {
/*MIKHAIL 			p++; */
			  
			byte<<=1;			  
            if (p->index == polarity)			
               byte|=0x01;
            bit++;
            if (bit == 8)
            {
				/*MIKHAIL (void) WriteByte(image,byte); */
                (void) fputc(byte,image->file);											
                bit=0;
                byte=0;
            }
            
          }
		  p++;
          /*MIKHAIL if (bit != 0) */
            /*MIKHAIL (void) WriteByte(image,byte << (8-bit)); */
			/*MIKHAIL (void) fputc(byte << (8-bit) ,image->file); */
/*MIKHAIL           if (image->previous == (Image *) NULL) */
/*MIKHAIL             if (QuantumTick(y,image->rows)) */
/*MIKHAIL               ProgressMonitor(SaveImageText,y,image->rows); */
        }
        break;
      }
      case '5':
      {
        /*
          Convert image to a PGM image.
        */
        FormatString(buffer,"%lu\n",MaxRGB);
        /*MIKHAIL (void) WriteBlob(image,strlen(buffer),buffer); */
		/*MIKHAIL p=GetImagePixels(image,0,y,image->columns,1);	 */
		/*MIKHAIL (void) WriteByte(image,index); */
		(void) fwrite(buffer,1,strlen(buffer),image->file);
        p=image->pixels;        
		if (p == (VigraImpexRunlengthPacket *) NULL)
            break;				
		for (y=0; y < (int) image->packets; y++)
        {
		  p->index=Intensity(*p);
          for (x=0; x <= (int) p->length; x++)
          {          
			(void) fputc(p->index,image->file);
          }
		   p++;
/*MIKHAIL           if (image->previous == (Image *) NULL) */
/*MIKHAIL             if (QuantumTick(y,image->rows)) */
/*MIKHAIL               ProgressMonitor(SaveImageText,y,image->rows); */
        }
        break;
      }
      case '6':
      {
        register unsigned char
          *q;

        unsigned char
          *pixels;

        /*
          Allocate memory for pixels.
        */
        pixels=(unsigned char *)
          malloc(image->rows*image->columns*3*sizeof(VigraImpexRunlengthPacket));
        if (pixels == (unsigned char *) NULL)
        {
          fprintf(stderr, "vigraImpexWritePNMImage(): Memory allocation failed.");
          vigraImpexDestroyImage(image);
          return 0;
        }
			
		/*
          Convert image to a PNM image.
        */
        FormatString(buffer,"%lu\n",MaxRGB);
        /*MIKHAIL (void) WriteBlob(image,strlen(buffer),buffer); */
		(void) fwrite(buffer,1,strlen(buffer),image->file);
		/*MIKHAIL p=GetImagePixels(image,0,y,image->columns,1); */
		p=image->pixels;
		q=pixels;		
        for (y=0; y < (int) image->packets; y++)
	    { 
	        if (p == (VigraImpexRunlengthPacket *) NULL)
		       break;

			for (x=0; x <= (int) p->length; x++)
			{
				*q++=p->blue;
				*q++=p->red;			
				*q++=p->green;			
									
			}
            p++;                 
		  /*MIKHAIL (void) WriteBlob(image,q-pixels,(char *) pixels);		   */
		   
            
          /*MIKHAIL if (image->previous == (VigraImpexImage *) NULL) */
/*MIKHAIL             if (QuantumTick(y,image->rows)) */
              /*MIKHAIL ProgressMonitor(SaveImageText,y,image->rows); */
        }
        /*MIKHAIL FreeMemory((void **) &pixels); */
		(void) fwrite((char *)pixels,1,image->rows*image->columns*3,image->file);
		/*MIKHAIL 		free(pixels);								 */
        break;
      }
      case '7':
      {
        static const short int
          dither_red[2][16]=
          {
            {-16,  4, -1, 11,-14,  6, -3,  9,-15,  5, -2, 10,-13,  7, -4,  8},
            { 15, -5,  0,-12, 13, -7,  2,-10, 14, -6,  1,-11, 12, -8,  3, -9}
          },
          dither_green[2][16]=
          {
            { 11,-15,  7, -3,  8,-14,  4, -2, 10,-16,  6, -4,  9,-13,  5, -1},
            {-12, 14, -8,  2, -9, 13, -5,  1,-11, 15, -7,  3,-10, 12, -6,  0}
          },
          dither_blue[2][16]=
          {
            { -3,  9,-13,  7, -1, 11,-15,  5, -4,  8,-14,  6, -2, 10,-16,  4},
            {  2,-10, 12, -8,  0,-12, 14, -6,  3, -9, 13, -7,  1,-11, 15, -5}
          };

        int
          value;

        Quantum
          pixel;

        unsigned short
          *blue_map[2][16],
          *green_map[2][16],
          *red_map[2][16];

        /*
          Allocate and initialize dither maps.
        */
        for (i=0; i < 2; i++)
          for (j=0; j < 16; j++)
          {
            red_map[i][j]=(unsigned short *)
              malloc(256*sizeof(unsigned short));
            green_map[i][j]=(unsigned short *)
              malloc(256*sizeof(unsigned short));
            blue_map[i][j]=(unsigned short *)
              malloc(256*sizeof(unsigned short));
            if ((red_map[i][j] == (unsigned short *) NULL) ||
                (green_map[i][j] == (unsigned short *) NULL) ||
                (blue_map[i][j] == (unsigned short *) NULL))		
				{
					fprintf(stderr,"Memory Allocation failed");
					return 0;
				}

          }
        /*
          Initialize dither tables.
        */
        for (i=0; i < 2; i++)
          for (j=0; j < 16; j++)
            for (x=0; x < 256; x++)
            {
              value=x-16;
              if (x < 48)
                value=x/2+8;
              value+=dither_red[i][j];
              red_map[i][j][x]=(unsigned short)
                ((value < 0) ? 0 : (value > 255) ? 255 : value);
              value=x-16;
              if (x < 48)
                value=x/2+8;
              value+=dither_green[i][j];
              green_map[i][j][x]=(unsigned short)
                ((value < 0) ? 0 : (value > 255) ? 255 : value);
              value=x-32;
              if (x < 112)
                value=x/2+24;
              value+=2*dither_blue[i][j];
              blue_map[i][j][x]=(unsigned short)
                ((value < 0) ? 0 : (value > 255) ? 255 : value);
            }
        /*
          Convert image to a P7 image.
        */
        (void) strcpy(buffer,"#END_OF_COMMENTS\n");
        /*MIKHAIL (void) WriteBlob(image,strlen(buffer),buffer); */
		(void) fwrite(buffer,1,strlen(buffer),image->file);
        FormatString(buffer,"%u %u 255\n",image->columns,image->rows);
        /*MIKHAIL (void) WriteBlob(image,strlen(buffer),buffer); */
		(void) fwrite(buffer,1,strlen(buffer),image->file);
        i=0;
        j=0;
        for (y=0; y < (int) image->rows; y++)
        {
          /*MIKHAIL p=GetImagePixels(image,0,y,image->columns,1); */
			p=image->pixels+(y*image->columns);
          if (p == (VigraImpexRunlengthPacket *) NULL)
            break;
          for (x=0; x < (int) image->columns; x++)
          {
            if (!image_info->dither)
              pixel=((p->red & 0xe0) |
                ((p->green & 0xe0) >> 3) |
                ((p->blue & 0xc0) >> 6));
            else
              pixel=((red_map[i][j][p->red] & 0xe0) |
                ((green_map[i][j][p->green] & 0xe0) >> 3) |
                ((blue_map[i][j][p->blue] & 0xc0) >> 6));
            /*MIKHAIL (void) WriteByte(image,pixel); */
			(void) fputc(pixel,image->file);
            p++;
            j++;
            if (j == 16)
              j=0;
          }
          i++;
          if (i == 2)
            i=0;
/*MIKHAIL           if (QuantumTick(y,image->rows)) */
/*MIKHAIL             ProgressMonitor(SaveImageText,y,image->rows); */
        }
        /*
          Free allocated memory.
        */
        for (i=0; i < 2; i++)
          for (j=0; j < 16; j++)
          {
			free(green_map[i][j]);
			free(blue_map[i][j]);
			free(red_map[i][j]);
/*MIKHAIL             FreeMemory((void **) &green_map[i][j]); */
  /*MIKHAIL /          FreeMemory((void **) &blue_map[i][j]); */
     /*MIKHAIL        FreeMemory((void **) &red_map[i][j]); */
          }
        break;
      }
    }
    if (image->next == (VigraImpexImage *) NULL)
      break;
    image=image->next;
    /*MIKHAIL ProgressMonitor(SaveImagesText,scene++,GetNumberScenes(image)); */
  } while (image_info->adjoin);
  if (image_info->adjoin)
    while (image->previous != (VigraImpexImage *) NULL)
      image=image->previous;
  vigraImpexCloseImage(image);
  return(True);
}
