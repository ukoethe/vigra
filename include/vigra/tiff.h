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
 
#ifndef _VIGRA_TIFF_H
#define	_VIGRA_TIFF_H

/*
 * Tag Image File Format (TIFF)
 *
 *    Highjacked from:
 *    lib TIFF Rev 6.0
 *    Developer's Desk
 *    Aldus Corporation
 *    411 First Ave. South
 *    Suite 200
 *    Seattle, WA  98104
 *    206-622-5500
 */
#define	TIFF_VERSION	42

#ifndef _TIFF_DATA_TYPEDEFS_
#define _TIFF_DATA_TYPEDEFS_

#ifdef __STDC__
typedef	signed char int8;	/* NB: non-ANSI compilers may not grok */
#else
typedef	char int8;
#endif
typedef	unsigned char uint8;
typedef	short int16;
typedef	unsigned short uint16;	/* sizeof (uint16) must == 2 */
#if defined(__alpha) || (defined(_MIPS_SZLONG) && _MIPS_SZLONG == 64)
typedef	int int32;
typedef	unsigned int uint32;	/* sizeof (uint32) must == 4 */
#else
typedef	long int32;
typedef	unsigned long uint32;	/* sizeof (uint32) must == 4 */
#endif
#endif /* _TIFF_DATA_TYPEDEFS_ */

typedef	uint32 ttag_t;		/* directory tag */
typedef	int32 tsize_t;		/* i/o size in bytes */
typedef	void* tdata_t;		/* image data ref */
typedef	uint16 tsample_t;	/* sample number */

struct TiffImage;

extern "C" void TIFFClose(TiffImage*);
extern "C" TiffImage* TIFFOpen(const char*, const char*);
extern "C" int TIFFGetField(TiffImage*, ttag_t, ...);
extern "C" int TIFFSetField(TiffImage*, ttag_t, ...);
extern "C" tsize_t TIFFScanlineSize(TiffImage*);
extern "C" int TIFFReadScanline(TiffImage*, tdata_t, uint32, tsample_t = 0);
extern "C" int TIFFWriteScanline(TiffImage*, tdata_t, uint32, tsample_t = 0);
extern "C" int TIFFReadRGBAImage(TiffImage*, uint32, uint32, uint32*, int = 0);

/*
 * Macros for extracting components from the
 * packed ABGR form returned by TIFFReadRGBAImage.
 */
#ifndef TIFFGetR

#define	TIFFGetR(abgr)	((abgr) & 0xff)
#define	TIFFGetG(abgr)	(((abgr) >> 8) & 0xff)
#define	TIFFGetB(abgr)	(((abgr) >> 16) & 0xff)
#define	TIFFGetA(abgr)	(((abgr) >> 24) & 0xff)

#endif /* TIFFGetR */

/*
 * TIFF Tag Definitions.
 */

#ifndef TIFFTAG_IMAGEWIDTH

#define	TIFFTAG_IMAGEWIDTH		256	/* image width in pixels */
#define	TIFFTAG_IMAGELENGTH		257	/* image height in pixels */
#define	TIFFTAG_BITSPERSAMPLE		258	/* bits per channel (sample) */
#define	TIFFTAG_COMPRESSION		259	/* data compression technique */
#define	    COMPRESSION_NONE		1	/* dump mode */
#define	    COMPRESSION_LZW		5       /* Lempel-Ziv  & Welch */
#define	    COMPRESSION_JPEG		7	/* %JPEG DCT compression */
#define	    COMPRESSION_PACKBITS	32773	/* Macintosh RLE */
#define	TIFFTAG_JPEGQUALITY		65537	/* Compression quality level */
#define	TIFFTAG_PHOTOMETRIC		262	/* photometric interpretation */
#define	    PHOTOMETRIC_MINISWHITE	0	/* min value is white */
#define	    PHOTOMETRIC_MINISBLACK	1	/* min value is black */
#define	    PHOTOMETRIC_RGB		2	/* RGB color model */
#define	    PHOTOMETRIC_PALETTE		3	/* color map indexed */
#define	TIFFTAG_FILLORDER		266	/* data order within a byte */
#define	    FILLORDER_MSB2LSB		1	/* most significant -> least */
#define	    FILLORDER_LSB2MSB		2	/* least significant -> most */
#define	TIFFTAG_SAMPLESPERPIXEL		277	/* samples per pixel */
#define	TIFFTAG_PLANARCONFIG		284	/* storage organization */
#define	    PLANARCONFIG_CONTIG		1	/* single image plane */
#define	    PLANARCONFIG_SEPARATE	2	/* separate planes of data */
#define	TIFFTAG_SAMPLEFORMAT		339	/* !data sample format */
#define	    SAMPLEFORMAT_UINT		1	/* !unsigned integer data */
#define	    SAMPLEFORMAT_INT		2	/* !signed integer data */
#define	    SAMPLEFORMAT_IEEEFP		3	/* !IEEE floating point data */
#define	    SAMPLEFORMAT_VOID		4	/* !untyped data */
#endif /* TIFFTAG_IMAGEWIDTH */

#endif /* _VIGRA_TIFF_H */



