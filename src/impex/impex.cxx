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
 

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "vigra/impex.hxx"
#include "gif.h"
#include "bmp.h"
#include "sun.h"
#include "jpeg.h"
#include "pnm.h"

namespace vigra {

static ImageFileTypeInfo imageFileTypeInfo[] = 
{
    {"GIF", ".gif", "GIF8", 4, &vigraImpexWriteGIFImage, &vigraImpexReadGIFImage}, 
#ifdef HasJPEG
    {"JPEG", ".jpg", "\377\330\377", 3, &vigraImpexWriteJPEGImage, &vigraImpexReadJPEGImage}, 
    {"JPEG", ".jpeg", "\377\330\377", 3, &vigraImpexWriteJPEGImage, &vigraImpexReadJPEGImage}, 
#endif
    {"BMP", ".bmp", "BM", 2, &vigraImpexWriteBMPImage, &vigraImpexReadBMPImage}, 
    {"SUN", ".ras", "\131\246\152\225", 4, &vigraImpexWriteSUNImage, &vigraImpexReadSUNImage}, 
    {"PBM", ".pbm", "P1", 2, &vigraImpexWritePNMImage, &vigraImpexReadPNMImage}, 
    {"PGM", ".pgm", "P2", 2, &vigraImpexWritePNMImage, &vigraImpexReadPNMImage}, 
    {"PPM", ".ppm", "P3", 2, &vigraImpexWritePNMImage, &vigraImpexReadPNMImage}, 
    {"PBM", ".pbm", "P4", 2, &vigraImpexWritePNMImage, &vigraImpexReadPNMImage}, 	
    {"PGM", ".pgm", "P5", 2, &vigraImpexWritePNMImage, &vigraImpexReadPNMImage}, 	
    {"PPM", ".ppm", "P6", 2, &vigraImpexWritePNMImage, &vigraImpexReadPNMImage},
    {"P7",  ".p7", "P7", 2,   &vigraImpexWritePNMImage, &vigraImpexReadPNMImage}, 	
    {"PNM", ".pnm", "P8", 2, &vigraImpexWritePNMImage, &vigraImpexReadPNMImage}, 	    
#ifdef HasTIFF
    {"TIFF", ".tif", "\115\115\000\052", 4, 0, 0}, 
    {"TIFF", ".tiff", "\111\111\052\000", 4, 0, 0}, 
#endif
    {"VIFF", ".xv", "\253\1", 2, 0, 0}, 
    {0, 0, 0, 0, 0, 0}
};

std::string impexListFormats()
{
    std::string res;
    
    ImageFileTypeInfo * f = imageFileTypeInfo;
    
    for(; f->typeTag != 0; ++f)
    {
        if(f[1].typeTag && strcmp(f->typeTag, f[1].typeTag) == 0) continue;
        res += f->typeTag;
        res += " ";
    }
    
    return res;
}

ImageImportInfo::ImageImportInfo(char const * filename)
: filename_(filename),
  filetype_(0),
  colorspace_(UNDEF),
  viff_(0),
  tiff_(0),
  impex_(0)
{
    loadImage(filename);
}

ImageImportInfo::~ImageImportInfo()
{
    deleteInternalImages();
}

void ImageImportInfo::deleteInternalImages()
{
    if(viff_) 
    {
        freeViffImage(viff_);
        viff_ = 0;
    }
    if(tiff_) 
    {
        if(strcmp(filename_.c_str(), "-") != 0) TIFFClose(tiff_);
        tiff_ = 0;
    }
    if(impex_) 
    {
        vigraImpexDestroyImage(impex_);
        impex_ = 0;
    }
    colorspace_ = UNDEF;
}

void ImageImportInfo::loadImage(char const * filename)
{
    deleteInternalImages();
    
    findFileTypeFromMagicString(filename);
    
    if(strcmp(filetype_->typeTag, "VIFF") == 0)
    {
        viff_ = readViffImage((char *)filename);
        postcondition(viff_ != 0, 
               "ImageImportInfo::loadImage(): Unable to read image");
               
        if(viff_->num_data_bands == 3 || viff_->map_scheme == VFF_MS_ONEPERBAND)
        {
            colorspace_ = RGB;
        }
        else
        {
            colorspace_ = GRAY;
        }
        width_ = viff_->row_size;
        height_ = viff_->col_size;
        switch(viff_->data_storage_type)
        {
            case VFF_TYP_1_BYTE: pixelType_ = UINT8;  break;
            case VFF_TYP_2_BYTE: pixelType_ = INT16;  break;
            case VFF_TYP_4_BYTE: pixelType_ = INT32;  break;
            case VFF_TYP_FLOAT:  pixelType_ = FLOAT;  break;
            case VFF_TYP_DOUBLE: pixelType_ = DOUBLE; break;
            default:
                fail("ImageImportInfo::loadImage(): unsupported pixel type.");
        }
    }
    else if(strcmp(filetype_->typeTag, "TIFF") == 0)
    {
        tiff_ = TIFFOpen((char *)filename, "r");
        postcondition(tiff_ != 0, 
               "ImageImportInfo::loadImage(): Unable to open image");
               
        uint16 sampleFormat = 1, bitsPerSample, photometric;
        uint32 w,h;
        TIFFGetField(tiff_, TIFFTAG_SAMPLEFORMAT, &sampleFormat);
        TIFFGetField(tiff_, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
        TIFFGetField(tiff_, TIFFTAG_IMAGEWIDTH, &w);
        TIFFGetField(tiff_, TIFFTAG_IMAGELENGTH, &h);
        TIFFGetField(tiff_, TIFFTAG_PHOTOMETRIC, &photometric);
        
        switch(photometric)
        {
            case PHOTOMETRIC_MINISWHITE:
            case PHOTOMETRIC_MINISBLACK:
                colorspace_ = GRAY;
                break;
            case PHOTOMETRIC_RGB:
            case PHOTOMETRIC_PALETTE:
                colorspace_ = RGB;
                break;
            default:
                fail("ImageImportInfo::loadImage(): unsupported color model.");
        }
        width_ = w;
        height_ = h;
        switch(sampleFormat)
        {
            case SAMPLEFORMAT_UINT:
            case SAMPLEFORMAT_INT:
                switch(bitsPerSample)
                {
                    case 8: pixelType_ = UINT8;  break;
                    case 16: pixelType_ = INT16;  break;
                    case 32: pixelType_ = INT32;  break;
                    default:
                        fail("ImageImportInfo::loadImage(): unsupported pixel type.");
                }
                break;
            case SAMPLEFORMAT_IEEEFP:
                switch(bitsPerSample)
                {
                    case 8*sizeof(float):  pixelType_ = FLOAT;  break;
                    case 8*sizeof(double): pixelType_ = DOUBLE; break;
                    default:
                        fail("ImageImportInfo::loadImage(): unsupported pixel type.");
                }
                break;
            default:
                fail("ImageImportInfo::loadImage(): unsupported pixel type.");
        }
    }
    else
    {
        VigraImpexImageInfo image_info;

        vigraImpexGetImageInfo(&image_info);

        strcpy(image_info.filename, filename);

        impex_ = (*(filetype_->importFunction))(&image_info);
        
        vigraImpexDestroyImageInfo(&image_info);
        
        postcondition(impex_ != 0, 
               "ImageImportInfo::loadImage(): Unable to read image");
               
        if(vigraImpexIsGrayImage(impex_))
        {
            colorspace_ = GRAY;
        }
        else
        {
            colorspace_ = RGB;
        }
        width_ = impex_->columns;
        height_ = impex_->rows;
        pixelType_ = UINT8;
    }
}

void ImageImportInfo::findFileTypeFromMagicString(char const * filename)
{
    FILE * file;
    
    if(strcmp(filename, "-") == 0) 
    {
        file = stdin;
    }
    else
    {
        file = fopen(filename, "r");
    }
    postcondition(file != 0, 
      "ImageImportInfo::findFileTypeFromMagicString(): Unable to open file");
    
    const int bufsize = 10;
    char magic_string[bufsize];
    for(int i=0; i<bufsize; ++i) magic_string[i] = getc(file);
    
    if(file == stdin) 
    {
        for(int i=bufsize-1; i>=0; --i) ungetc(magic_string[i], stdin);
    }
    else
    {
        fclose(file);
    }
    
    ImageFileTypeInfo * f = imageFileTypeInfo;
    
    for(; f->typeTag != 0; ++f)
    {
        if(strncmp(magic_string, 
                   f->fileMagicString, f->lengthOfMagicString) == 0) break;
    }
    postcondition(f->typeTag != 0, 
      "ImageImportInfo::findFileTypeFromMagicString(): Unknown file type");
      
    filetype_ = f;
}

ImageExportInfo::ImageExportInfo(char const * filename)
: filename_(filename), 
  filetype_(0),
  compression_(VigraImpexUndefinedCompression)
{
    guessFiletypeFromExtension(filename);
}

ImageExportInfo & ImageExportInfo::setFileType(char const * filetype)
{
    ImageFileTypeInfo * f = imageFileTypeInfo;
    
    for(; f->typeTag != 0; ++f)
    {
        if(strcmp(f->typeTag, filetype) == 0) break;
    }
    precondition(f->typeTag != 0, 
           "ImageExportInfo::setFileType(): Unknown file type");
           
    filetype_ = f;
    
    return *this;
}

ImageExportInfo & ImageExportInfo::setCompression(char const * compression)
{
    float quality = -1;
    sscanf(compression, "%f", &quality);
    
    if(0.0 < quality && quality <= 100.0)
    {
        compression_ = VigraImpexJPEGCompression;
        quality_ = (int)quality;
    }
    else if(strcmp(compression, "LZW") == 0)
    {
        compression_ = VigraImpexLZWCompression;
    }
    else if(strcmp(compression, "RunLength") == 0)
    {
        compression_ = VigraImpexRunlengthEncodedCompression;
    }
    else if(strcmp(compression, "None") == 0)
    {
        compression_ = VigraImpexNoCompression;
    }
    else
    {
        precondition(0, 
          "ImageExportInfo::setCompression(): Unknown compression type");
    }
    
    return *this;
}

void ImageExportInfo::guessFiletypeFromExtension(char const * filename)
{
    if(filetype_ != 0) return;
    
    if(filename == 0 || *filename == 0) return;

    char const * dot = filename + strlen(filename) - 1;
    
    for(; dot > filename; --dot)
    {
        if(*dot == '.') break;
    }
    if(*dot != '.') return;
    
    char extension[10];
    int i;
    
    for(i=0; dot[i] != 0 && i < 9; ++i)   extension[i] = tolower(dot[i]);
    extension[i] = 0;
    
    ImageFileTypeInfo * f = imageFileTypeInfo;
    
    for(; f->typeTag != 0; ++f)
    {
        if(strcmp(f->fileExtension, extension) == 0) break;
    }
    
    if(f->typeTag != 0) filetype_ = f;
}

void ImageExportInfo::initImageInfo(VigraImpexImageInfo & image_info) const
{
    vigraImpexGetImageInfo(&image_info);

    precondition(filetype_ != 0, "ImageExportInfo::initImageInfo(): No filetype specified");
    
    strcpy(image_info.filename, filename_.c_str());
    strcpy(image_info.magick, filetype_->typeTag);
    image_info.compression = compression_;
    if(compression_ == VigraImpexJPEGCompression) 
    {
        image_info.quality = quality_;
    }
}


ImageExportFunctionPointer ImageExportInfo::exportFunction() const
{
    if(filetype_ == 0) return 0;
    return filetype_->exportFunction;
}

} // namespace vigra
