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
#if !defined(HasTIFF)

#include "vigra/error.hxx"
#include "vigra/tiff.h"

Export VigraImpexImage *vigraImpexReadTIFFImage( VigraImpexImageInfo *image_info)
{
  fail( "TIFF library is not available");
  return 0;
}
Export unsigned int vigraImpexWriteTIFFImage( VigraImpexImageInfo *image_info,VigraImpexImage *image)
{
  fail( "TIFF library is not available");
  return 0;
}
Export	void TIFFClose(TiffImage*)
{
  fail( "TIFF library is not available");
}
Export	TiffImage* TIFFOpen(const char*, const char*)
{
  fail( "TIFF library is not available");
  return 0;
}
Export	int TIFFGetField(TiffImage*, ttag_t, ...)
{
  fail( "TIFF library is not available");
  return 0;
}
Export	int TIFFSetField(TiffImage*, ttag_t, ...)
{
  fail( "TIFF library is not available");
  return 0;
}
Export	tsize_t TIFFScanlineSize(TiffImage*)
{
  fail( "TIFF library is not available");
  return 0;
}
Export	int TIFFReadScanline(TiffImage*, tdata_t, uint32, tsample_t = 0)
{
  fail( "TIFF library is not available");
  return 0;
}
Export	int TIFFWriteScanline(TiffImage*, tdata_t, uint32, tsample_t = 0)
{
  fail( "TIFF library is not available");
  return 0;
}
Export	int TIFFReadRGBAImage(TiffImage*, uint32, uint32, uint32*, int = 0)
{
  fail( "TIFF library is not available");
  return 0;
}
#endif
