#include "vigra/tiff_file.hxx"

namespace vigra {

TIFFFile::TIFFFile(const char *filename, const char *mode)
    : h_(NULL),
      fileType_(INVALID),
      hasDepth_(false)
{
    h_ = TIFFOpen(filename, mode);

    TIFFGetField(h_, TIFFTAG_IMAGEWIDTH,  &imageSize_[0]);
    TIFFGetField(h_, TIFFTAG_IMAGELENGTH, &imageSize_[1]);
    hasDepth_ = TIFFGetField(h_, TIFFTAG_IMAGEDEPTH,  &imageSize_[2]);

    bool
        hasTileWidth  = TIFFGetField(h_, TIFFTAG_TILEWIDTH,  &tileSize_[0]),
        hasTileLength = TIFFGetField(h_, TIFFTAG_TILELENGTH, &tileSize_[1]),
        hasTileDepth  = TIFFGetField(h_, TIFFTAG_TILEDEPTH,  &tileSize_[2]);

    if(hasTileWidth && hasTileLength && (hasTileDepth || !hasDepth_))
    {
        fileType_ = TILED;
    }
    else if(!hasTileWidth && !hasTileLength && !hasTileDepth)
    {
        fileType_ = STRIPED;
    }

    if(!hasDepth_)
    {
        imageSize_[2] = 1;
        if(fileType() == TILED)
            tileSize_[2] = 1;
    }
}

} // namespace vigra
