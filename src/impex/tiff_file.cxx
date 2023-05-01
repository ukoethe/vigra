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

tsize_t TIFFFile::readTile(tdata_t buf, uint32 x, uint32 y, uint32 z, tsample_t sample)
{
    return TIFFReadTile(h_, buf, x, y, z, sample);    
}

std::string TIFFFile::pixelType() const
{
    // get bits per sample, default to 8
    uint16 bits_per_sample;
    if(!TIFFGetField( h_, TIFFTAG_BITSPERSAMPLE, &bits_per_sample))
        bits_per_sample = 8;

    // get pixeltype
    if(bits_per_sample == 1)
        return "BILEVEL";

    // try the sampleformat tag
    uint16 sampleformat;
    if(TIFFGetField(h_, TIFFTAG_SAMPLEFORMAT, &sampleformat))
    {
        switch (sampleformat)
        {
        case SAMPLEFORMAT_UINT:
            // added by dangelo, support for UINT 16 & 32 bit
            switch(bits_per_sample)
            {
            case 8:
                return "UINT8";
            case 16:
                return "UINT16";
            case 32:
                return "UINT32";
            }
            break;

        case SAMPLEFORMAT_INT:
            switch (bits_per_sample)
            {
            case 8:
                return "INT8";
            case 16:
                return "INT16";
            case 32:
                return "INT32";
            }
            break;

        case SAMPLEFORMAT_IEEEFP:
            switch (bits_per_sample)
            {
            case 32:
                return "FLOAT";
            case 64:
                return "DOUBLE";
            }
            break;

        default:
            ;
        }
    }

    // fall back to the (obsolete) datatype tag
    uint16 datatype;
    if(TIFFGetField(h_, TIFFTAG_DATATYPE, &datatype))
    {
        // dangelo: correct parsing of INT/UINT (given in tiff.h)
        switch (datatype)
        {
        case TIFF_BYTE:
            return "UINT8";
        case TIFF_SBYTE:
            return "INT8";
        case TIFF_SHORT:
            return "UINT16";
        case TIFF_SSHORT:
            return "INT16";
        case TIFF_LONG:
            return "UINT32";
        case TIFF_SLONG:
            return "INT32";
        case TIFF_FLOAT:
            return "FLOAT";
        case TIFF_DOUBLE:
            return "DOUBLE";
        default:
            ;
        }
    }

    // ugly: no useable pixeltype found..
    // e.g. imagemagick writes files without it; try to guess a suitable one here:
    switch(bits_per_sample)
    {
    case 8:
        return "UINT8";
    case 16:
        return "UINT16";
    case 32:
        return "UINT32"; // prefer int over float
    case 64:
        return  "DOUBLE";
    default:
        ;
    }

    vigra_fail( "TIFFDecoderImpl::init(): Sampleformat or Datatype tag undefined and guessing sampletype from Bits per Sample failed." );
    return "<unreached>";
}

} // namespace vigra
