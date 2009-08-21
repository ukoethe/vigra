/************************************************************************/
/*                                                                      */
/*               Copyright 2002 by Gunnar Kedenburg                     */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */
/*                                                                      */
/************************************************************************/
/* Modifications by Pablo d'Angelo
 * updated to vigra 1.4 by Douglas Wilkins
 * as of 18 Febuary 2006:
 *  - Added UINT16 and UINT32 pixel types.
 *  - Added support for obtaining extra bands beyond RGB.
 *  - Added support for a position field that indicates the start of this
 *    image relative to some global origin.
 *  - Added support for x and y resolution fields.
 *  - Added support for ICC profiles
 */

#include <iostream>
#include <algorithm>
#include <functional>
#include <iterator>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <sys/types.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "vigra/imageinfo.hxx"
#include "codecmanager.hxx"
#include "vigra/multi_impex.hxx"
#include "vigra/hdf5impex.hxx"

#if defined(_WIN32)
#  include "vigra/windows.h"
#else
#  include <dirent.h>
#endif

namespace vigra
{

namespace detail
{

struct NumberCompare
{
    bool operator()(std::string const & l, std::string const & r) const
    {
        return atoi(l.c_str()) < atoi(r.c_str());
    }
};

bool splitString(
    const std::string &s, char separator, std::string &a, std::string &b,
    bool reverse = false)
{
    std::size_t splitPos = (reverse ? s.rfind(separator) : s.find(separator));
    if(splitPos >= s.size())
        return false;
    a = std::string(s.begin(), s.begin() + splitPos);
    b = std::string(s.begin() + splitPos + 1, s.end());
    return true;
}

std::string trimString(const std::string &s)
{
    unsigned int begin = 0;
    while(begin < s.size() && ((s[begin] == ' ') || (s[begin] == '\t')))
        ++begin;
    std::size_t end = s.size();
    while(end > 0 && ((s[end-1] == ' ') || (s[end-1] == '\t')))
        --end;
    return std::string(s.begin() + begin, s.begin() + end);
}

} // namespace detail


// find filenames matching the pattern "<path>/base[0-9]+ext"
#ifdef _WIN32
void splitPathFromFilename(const std::string &pathAndName,
                           std::string &path, std::string &name)
{
    // on Windows, both '/' and '\' are valid path separators
    // note: std::basic_string.rfind() may return 'unsigned int', so explicitely cast to 'int'
    int split = std::max(static_cast<int>(pathAndName.rfind('/')), static_cast<int>(pathAndName.rfind('\\')));
    if(split == static_cast<int>(std::string::npos))
    {
        path = ".";
        name = pathAndName;
    }
    else
    {
        for(int i=0; i<split; ++i)
        {
            if(pathAndName[i] == '/')
                path += '\\';
            else
                path += pathAndName[i];
        }
        name.append(pathAndName, split+1, pathAndName.size() - split - 1);
    }
}

VIGRA_EXPORT void findImageSequence(const std::string &name_base,
                       const std::string &name_ext,
                       std::vector<std::string> & numbers)
{
    // find out how many images we have
    BOOL            fFinished;
    HANDLE          hList;
    TCHAR           szDir[MAX_PATH+1];
    WIN32_FIND_DATA FileData;

    std::string path, base;
    splitPathFromFilename(name_base, path, base);

    std::vector<std::string> result;
    char numbuf[21], extbuf[1024];
    std::string pattern = base + "%20[0-9]%1023s";

    // Get the proper directory path
    sprintf(szDir, "%s\\%s*%s", path.c_str(), base.c_str(), name_ext.c_str());

    // Get the first file
    hList = FindFirstFile(szDir, &FileData);
    if (hList == INVALID_HANDLE_VALUE)
    {
        std::string message("importVolume(): No files matching '");
        message = message + szDir + "'.";
        vigra_fail(message.c_str());
    }
    else
    {
        // Traverse through the directory structure
        fFinished = FALSE;
        while (!fFinished)
        {
            if(sscanf(FileData.cFileName, pattern.c_str(), numbuf, extbuf) == 2)
            {
                if(strcmp(name_ext.c_str(), extbuf) == 0)
                {
                    std::string num(numbuf);
                    std::string name = name_base + num + name_ext;
                    // skip matching files names that are not images
                    if(isImage(name.c_str()))
                        result.push_back(num);
                }
            }
            if (!FindNextFile(hList, &FileData))
            {
                if (GetLastError() == ERROR_NO_MORE_FILES)
                {
                    fFinished = TRUE;
                }
            }
        }
    }

    FindClose(hList);

    std::sort(result.begin(), result.end(), detail::NumberCompare());
    numbers.swap(result);
}

#else // _WIN32

void splitPathFromFilename(const std::string &pathAndName,
                           std::string &path, std::string &name)
{
    int split = pathAndName.rfind('/');
    if(split == -1)
    {
        path = ".";
        name = pathAndName;
    }
    else
    {
        path.append(pathAndName, 0, split);
        name.append(pathAndName, split+1, pathAndName.size() - split - 1);
    }
}

void findImageSequence(const std::string &name_base,
                       const std::string &name_ext,
                       std::vector<std::string> & numbers)
{
    // find out how many images we have
    std::string path, base;
    splitPathFromFilename(name_base, path, base);

    DIR * dir = opendir(path.c_str());
    if(!dir)
    {
        std::string message("importVolume(): Unable to open directory '");
        message = message + path + "'.";
        vigra_fail(message.c_str());
    }

    std::vector<std::string> result;
    dirent * dp;
    errno = 0;
    char numbuf[21], extbuf[1024];
    std::string pattern = base + "%20[0-9]%1023s";
    while ((dp = readdir(dir)) != NULL)
    {
        if(sscanf(dp->d_name, pattern.c_str(), numbuf, extbuf) == 2)
        {
            if(strcmp(name_ext.c_str(), extbuf) == 0)
            {
                std::string num(numbuf);
                std::string name = name_base + num + name_ext;
                // skip matching files names that are not images
                if(isImage(name.c_str()))
                    result.push_back(num);
            }
        }
    }

    closedir(dir);

    vigra_precondition(errno == 0,
          "importVolume(): I/O error while searching for images.");

    std::sort(result.begin(), result.end(), detail::NumberCompare());
    numbers.swap(result);
}

#endif // _WIN32

// build a string from a sequence.
#if defined(_MSC_VER) && (_MSC_VER < 1300)
template <class iterator>
std::string stringify (const iterator &start, const iterator &end)
{
    return stringifyImpl(start, end, *start);
}

template <class iterator, class Value>
std::string stringifyImpl (const iterator &start, const iterator &end, Value const &)
{
    std::ostringstream out;
    // do not place a space character after the last sequence element.
    std::copy (start, end - 1,
               std::ostream_iterator <Value> (out, " "));
    out << *(end-1);
    return out.str ();
}

#else

template <class iterator>
std::string stringify (const iterator &start, const iterator &end)
{
    typedef typename std::iterator_traits<iterator>::value_type value_type;
    std::ostringstream out;
    // do not place a space character after the last sequence element.
    std::copy (start, end - 1,
               std::ostream_iterator <value_type> (out, " "));
    out << *(end-1);
    return out.str ();
}

#endif // _MSC_VER < 1300

void validate_filetype( std::string filetype )
{
    vigra_precondition( codecManager().fileTypeSupported(filetype),
                        "given file type is not supported" );
}

std::string impexListFormats()
{
    std::vector<std::string> ft = codecManager().supportedFileTypes();
    return stringify( ft.begin(), ft.end() );
}

std::string impexListExtensions()
{
    std::vector<std::string> ft = codecManager().supportedFileExtensions();
    return stringify( ft.begin(), ft.end() );
}

bool isImage(char const * filename)
{
    return CodecManager::manager().getFileTypeByMagicString(filename) != "";
}

// class ImageExportInfo

ImageExportInfo::ImageExportInfo( const char * filename )
    : m_filename(filename),
      m_x_res(0), m_y_res(0),
      fromMin_(0.0), fromMax_(0.0), toMin_(0.0), toMax_(0.0)
{}

ImageExportInfo::~ImageExportInfo()
{
}

ImageExportInfo & ImageExportInfo::setFileType( const char * filetype )
{
    m_filetype = filetype;
    return *this;
}

ImageExportInfo & ImageExportInfo::setForcedRangeMapping(double fromMin, double fromMax,
                                                     double toMin, double toMax)
{
    fromMin_ = fromMin;
    fromMax_ = fromMax;
    toMin_ = toMin;
    toMax_ = toMax;
    return *this;
}

bool ImageExportInfo::hasForcedRangeMapping() const
{
    return (fromMax_ > fromMin_) || (toMax_ > toMin_);
}

double ImageExportInfo::getFromMin() const
{
    return fromMin_;
}

double ImageExportInfo::getFromMax() const
{
    return fromMax_;
}

double ImageExportInfo::getToMin() const
{
    return toMin_;
}

double ImageExportInfo::getToMax() const
{
    return toMax_;
}

ImageExportInfo & ImageExportInfo::setCompression( const char * comp )
{
    m_comp = comp;
    return *this;
}

ImageExportInfo & ImageExportInfo::setFileName(const char * name)
{
    m_filename = name;
    return *this;
}

const char * ImageExportInfo::getFileName() const
{
    return m_filename.c_str();
}

const char * ImageExportInfo::getFileType() const
{
    return m_filetype.c_str();
}

ImageExportInfo & ImageExportInfo::setPixelType( const char * s )
{
    m_pixeltype = s;
    return *this;
}

const char * ImageExportInfo::getPixelType() const
{
    return m_pixeltype.c_str();
}

const char * ImageExportInfo::getCompression() const
{
    return m_comp.c_str();
}

float ImageExportInfo::getXResolution() const
{
    return m_x_res;
}

float ImageExportInfo::getYResolution() const
{
    return m_y_res;
}

ImageExportInfo & ImageExportInfo::setXResolution( float val )
{
    m_x_res = val;
    return *this;
}

ImageExportInfo & ImageExportInfo::setYResolution( float val )
{
    m_y_res = val;
    return *this;
}

ImageExportInfo & ImageExportInfo::setPosition(const vigra::Diff2D & pos)
{
    m_pos = pos;
    return *this;
}

vigra::Diff2D ImageExportInfo::getPosition() const
{
    return m_pos;
}

const ImageExportInfo::ICCProfile & ImageExportInfo::getICCProfile() const
{
    return m_icc_profile;
}

ImageExportInfo & ImageExportInfo::setICCProfile(
    const ImageExportInfo::ICCProfile &profile)
{
    m_icc_profile = profile;
    return *this;
}

// return an encoder for a given ImageExportInfo object
std::auto_ptr<Encoder> encoder( const ImageExportInfo & info )
{
    std::auto_ptr<Encoder> enc;

    std::string filetype = info.getFileType();
    if ( filetype != "" ) {
        validate_filetype(filetype);
        std::auto_ptr<Encoder> enc2
            = getEncoder( std::string( info.getFileName() ), filetype );
        enc = enc2;
    } else {
        std::auto_ptr<Encoder> enc2
            = getEncoder( std::string( info.getFileName() ) );
        enc = enc2;
    }

    std::string comp = info.getCompression();
    if ( comp != "" ) {

        // check for JPEG compression
        int quality = -1;
        std::istringstream compstream(comp.c_str());
        compstream >> quality;

        // FIXME: dangelo: This code might lead to strange effects (setting an invalid compression mode),
        // if other formats also support a numerical compression parameter.
        if ( quality != -1 ) {
            enc->setCompressionType( "JPEG", quality );
        } else {
            // leave any other compression type to the codec
            enc->setCompressionType(comp);
        }
    }

    std::string pixel_type = info.getPixelType();
    if ( pixel_type != "" ) {
        if(!isPixelTypeSupported( enc->getFileType(), pixel_type ))
        {
            std::string msg("exportImage(): file type ");
            msg += enc->getFileType() + " does not support requested pixel type "
                                      + pixel_type + ".";
            vigra_precondition(false, msg.c_str());
        }
        enc->setPixelType(pixel_type);
    }

    // set other properties
    enc->setXResolution(info.getXResolution());
    enc->setYResolution(info.getYResolution());
    enc->setPosition(info.getPosition());

    if ( info.getICCProfile().size() > 0 ) {
        enc->setICCProfile(info.getICCProfile());
    }

    return enc;
}

// class ImageImportInfo

ImageImportInfo::ImageImportInfo( const char * filename )
    : m_filename(filename)
{
    std::auto_ptr<Decoder> decoder = getDecoder(m_filename);

    m_filetype = decoder->getFileType();
    m_pixeltype = decoder->getPixelType();
    m_width = decoder->getWidth();
    m_height = decoder->getHeight();
    m_num_bands = decoder->getNumBands();
    m_num_extra_bands = decoder->getNumExtraBands();
    m_pos = decoder->getPosition();
    m_x_res = decoder->getXResolution();
    m_y_res = decoder->getYResolution();

    m_icc_profile = decoder->getICCProfile();

    decoder->abort(); // there probably is no better way than this
}

ImageImportInfo::~ImageImportInfo() {
}

const char * ImageImportInfo::getFileName() const
{
    return m_filename.c_str();
}

const char * ImageImportInfo::getFileType() const
{
    return m_filetype.c_str();
}

const char * ImageImportInfo::getPixelType() const
{
    return m_pixeltype.c_str();
}

ImageImportInfo::PixelType ImageImportInfo::pixelType() const
{
    const std::string pixeltype=ImageImportInfo::getPixelType();
   if (pixeltype == "UINT8")
     return UINT8;
   if (pixeltype == "INT16")
     return INT16;
   if (pixeltype == "UINT16")
     return UINT16;
   if (pixeltype == "INT32")
     return INT32;
   if (pixeltype == "UINT32")
     return UINT32;
   if (pixeltype == "FLOAT")
     return FLOAT;
   if (pixeltype == "DOUBLE")
     return DOUBLE;
   vigra_fail( "internal error: unknown pixel type" );
   return ImageImportInfo::PixelType();
}

int ImageImportInfo::width() const
{
    return m_width;
}

int ImageImportInfo::height() const
{
    return m_height;
}

int ImageImportInfo::numBands() const
{
    return m_num_bands;
}

int ImageImportInfo::numExtraBands() const
{
    return m_num_extra_bands;
}

Size2D ImageImportInfo::size() const
{
    return Size2D( m_width, m_height );
}

MultiArrayShape<2>::type ImageImportInfo::shape() const
{
    return MultiArrayShape<2>::type( m_width, m_height );
}

bool ImageImportInfo::isGrayscale() const
{
    return m_num_bands == 1;
}

bool ImageImportInfo::isColor() const
{
    return (m_num_bands - m_num_extra_bands) == 3;
}

bool ImageImportInfo::isByte() const
{
    return m_pixeltype == "UINT8";
}

Diff2D ImageImportInfo::getPosition() const
{
    return m_pos;
}

float ImageImportInfo::getXResolution() const
{
    return m_x_res;
}

float ImageImportInfo::getYResolution() const
{
    return m_y_res;
}

const ImageImportInfo::ICCProfile & ImageImportInfo::getICCProfile() const
{
    return m_icc_profile;
}

// return a decoder for a given ImageImportInfo object
std::auto_ptr<Decoder> decoder( const ImageImportInfo & info )
{
    std::string filetype = info.getFileType();
    validate_filetype(filetype);
    return getDecoder( std::string( info.getFileName() ), filetype );
}

// class VolumeExportInfo

VolumeExportInfo::VolumeExportInfo( const char * name_base, const char * name_ext ) : m_filename_base(name_base),
      m_filename_ext(name_ext), m_x_res(0), m_y_res(0), m_z_res(0), 
      fromMin_(0.0), fromMax_(0.0), toMin_(0.0), toMax_(0.0)
{
}

VolumeExportInfo::~VolumeExportInfo()
{
}

VolumeExportInfo & VolumeExportInfo::setFileType( const char * filetype )
{
    m_filetype = filetype;
    return *this;
}

VolumeExportInfo & VolumeExportInfo::setForcedRangeMapping(double fromMin, double fromMax,
                                                     double toMin, double toMax) 
{
    fromMin_ = fromMin;
    fromMax_ = fromMax;
    toMin_ = toMin;
    toMax_ = toMax;
    return *this;
}

bool VolumeExportInfo::hasForcedRangeMapping() const
{
    return (fromMax_ > fromMin_) || (toMax_ > toMin_);
}

double VolumeExportInfo::getFromMin() const
{
    return fromMin_;
}

double VolumeExportInfo::getFromMax() const
{
    return fromMax_;
}

double VolumeExportInfo::getToMin() const
{
    return toMin_;
}

double VolumeExportInfo::getToMax() const
{
    return toMax_;
}

VolumeExportInfo & VolumeExportInfo::setCompression( const char * comp )
{
    m_comp = comp;
    return *this;
}

VolumeExportInfo & VolumeExportInfo::setFileNameExt(const char * name_ext)
{
    m_filename_ext = name_ext;
    return *this;
}

const char * VolumeExportInfo::getFileNameExt() const
{
    return m_filename_ext.c_str();
}

VolumeExportInfo & VolumeExportInfo::setFileNameBase(const char * name_base)
{
    m_filename_base = name_base;
    return *this;
}

const char * VolumeExportInfo::getFileNameBase() const
{
    return m_filename_base.c_str();
}

const char * VolumeExportInfo::getFileType() const
{
    return m_filetype.c_str();
}

VolumeExportInfo & VolumeExportInfo::setPixelType( const char * s )
{
    m_pixeltype = s;
    return *this;
}

const char * VolumeExportInfo::getPixelType() const
{
    return m_pixeltype.c_str();
}

const char * VolumeExportInfo::getCompression() const
{
    return m_comp.c_str();
}

float VolumeExportInfo::getXResolution() const
{
    return m_x_res;
}

float VolumeExportInfo::getYResolution() const
{
    return m_y_res;
}

VolumeExportInfo & VolumeExportInfo::setXResolution( float val )
{
    m_x_res = val;
    return *this;
}

VolumeExportInfo & VolumeExportInfo::setYResolution( float val )
{
    m_y_res = val;
    return *this;
}

VolumeExportInfo & VolumeExportInfo::setZResolution( float val )
{
    m_z_res = val;
    return *this;
}

VolumeExportInfo & VolumeExportInfo::setPosition(const vigra::Diff2D & pos)
{
    m_pos = pos;
    return *this;
}

vigra::Diff2D VolumeExportInfo::getPosition() const
{
    return m_pos;
}

const VolumeExportInfo::ICCProfile & VolumeExportInfo::getICCProfile() const
{
    return m_icc_profile;
}

VolumeExportInfo & VolumeExportInfo::setICCProfile(
    const VolumeExportInfo::ICCProfile &profile)
{
    m_icc_profile = profile;
    return *this;
}

VolumeImportInfo::VolumeImportInfo(const std::string &filename)
: shape_(0, 0, 0),
  resolution_(1.f, 1.f, 1.f),
  numBands_(0)
{
    // first try image sequence loading
    std::string::const_reverse_iterator
        numBeginIt(filename.rbegin()), numEndIt(numBeginIt);

    do
    {
        numEndIt = std::find_if(numBeginIt, filename.rend(),(int (*)(int)) &isdigit);
        numBeginIt = std::find_if(numEndIt, filename.rend(), not1(std::ptr_fun((int (*)(int))&isdigit)));

        if(numEndIt != filename.rend())
        {
            std::string
                baseName(filename.begin(),
                         filename.begin() + (filename.rend()-numBeginIt)),
                extension(filename.begin() + (filename.rend()-numEndIt),
                          filename.end());

            std::vector<std::string> numbers;

            findImageSequence(baseName, extension, numbers);
            if(numbers.size() > 0)
            {
                getVolumeInfoFromFirstSlice(baseName + numbers[0] + extension);
                splitPathFromFilename(baseName, path_, name_);
                baseName_ = baseName;
                extension_ = extension;
                shape_[2] = numbers.size();
                std::swap(numbers, numbers_);

                break;
            }
        }
    }
    while(numEndIt != filename.rend());

    // no numbered images found, try .info file loading
    if(!numbers_.size())
    {
        std::ifstream stream(filename.c_str());

        while(stream.good())
        {
            char rawline[1024];
            stream.getline(rawline, 1024);

            // split off comments starting with '#':
            std::string line, comment;
            if(!detail::splitString(rawline, '#', line, comment))
                line = rawline;

            std::string key, value;
            if(detail::splitString(line, '=', key, value))
            {
                key = detail::trimString(key);
                value = detail::trimString(value);

                if(key == "width")
                    shape_[0] = atoi(value.c_str());
                else if(key == "height")
                    shape_[1] = atoi(value.c_str());
                else if(key == "depth")
                    shape_[2] = atoi(value.c_str());
                else if(key == "datatype")
                {
                    // FUTURE: store bit depth / signedness
                    if((value == "UNSIGNED_CHAR") || (value == "UNSIGNED_BYTE"))
                        numBands_ = 1;
                    else
                    {
                        std::cerr << "Unknown datatype '" << value << "'!\n";
                        break;
                    }
                }
                else if(key == "description")
                    description_ = value;
                else if(key == "name")
                    name_ = value;
                else if(key == "filename")
                    rawFilename_ = value;
                else
                {
                    std::cerr << "WARNING: Unknown key '" << key
                              << "' (value '" << value << "') in info file!\n";
                }
            }
            else
            {
                if(line[0]) // non-empty line?
                    std::cerr << "WARNING: could not parse line '" << line << "'!\n";
            }
        }

        if((shape_[0]*shape_[1]*shape_[2] > 0) && (rawFilename_.size() > 0))
        {
            if(!numBands_)
                numBands_ = 1; // default to UNSIGNED_CHAR datatype

            baseName_ = filename;
            if(name_.size() > 0)
            {
                std::string nameDummy;
                splitPathFromFilename(baseName_, path_, nameDummy);
            }
            else
            {
                splitPathFromFilename(baseName_, path_, name_);
            }
            return;
        }

        std::string message("VolumeImportInfo(): Unable to load volume '");
        message += filename + "'.";
        vigra_fail(message.c_str());
    }
}

VolumeImportInfo::VolumeImportInfo(const std::string &baseName, const std::string &extension)
: shape_(0, 0, 0),
  resolution_(1.f, 1.f, 1.f),
  numBands_(0)
{
    std::vector<std::string> numbers;
    findImageSequence(baseName, extension, numbers);

    std::string message("VolumeImportInfo(): No files matching '");
    message += baseName + "[0-9]+" + extension + "' found.";
    vigra_precondition(numbers.size() > 0, message.c_str());

    getVolumeInfoFromFirstSlice(baseName + numbers[0] + extension);

    splitPathFromFilename(baseName, path_, name_);
    baseName_ = baseName;
    extension_ = extension;
    shape_[2] = numbers.size();
    std::swap(numbers, numbers_);
}

void VolumeImportInfo::getVolumeInfoFromFirstSlice(const std::string &filename)
{
    ImageImportInfo info(filename.c_str());
    shape_[0] = info.width();
    shape_[1] = info.height();
    resolution_[1] = -1.f; // assume images to be right-handed
    pixelType_ = info.pixelType();
    numBands_ = info.numBands();
}

VolumeImportInfo::ShapeType VolumeImportInfo::shape() const { return shape_; }
VolumeImportInfo::Resolution VolumeImportInfo::resolution() const { return resolution_; }
VolumeImportInfo::PixelType VolumeImportInfo::pixelType() const
{
    const std::string pixeltype=VolumeImportInfo::getPixelType();
   if (pixeltype == "UINT8")
       return ImageImportInfo::UINT8;
   if (pixeltype == "INT16")
     return ImageImportInfo::INT16;
   if (pixeltype == "UINT16")
     return ImageImportInfo::UINT16;
   if (pixeltype == "INT32")
     return ImageImportInfo::INT32;
   if (pixeltype == "UINT32")
     return ImageImportInfo::UINT32;
   if (pixeltype == "FLOAT")
     return ImageImportInfo::FLOAT;
   if (pixeltype == "DOUBLE")
     return ImageImportInfo::DOUBLE;
   vigra_fail( "internal error: unknown pixel type" );
   return VolumeImportInfo::PixelType();
}
const char * VolumeImportInfo::getPixelType() const
{
    return pixelType_.c_str();
}
MultiArrayIndex VolumeImportInfo::numBands() const { return numBands_; }
bool VolumeImportInfo::isGrayscale() const { return numBands_ == 1; }
bool VolumeImportInfo::isColor() const { return numBands_ > 1; }
MultiArrayIndex VolumeImportInfo::width() const { return shape_[0]; }
MultiArrayIndex VolumeImportInfo::height() const { return shape_[1]; }
MultiArrayIndex VolumeImportInfo::depth() const { return shape_[2]; }
const std::string & VolumeImportInfo::name() const { return name_; }
const std::string & VolumeImportInfo::description() const { return description_; }

#ifdef HasHDF5

HDF5ImportInfo::HDF5ImportInfo(const std::string &filename, const std::string &path)
: m_file(0),
  m_dataset(0),
  m_filename(filename),
  m_datasetname(path)
{
	std::string path_name(path), group_name, data_set_name, message;
    hid_t group = 0, dspace = 0, dtype = 0, native_type = 0;

    m_file = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if(m_file < 0)
    {
        message = std::string("HDF5ImportInfo(): Unable to open file '") + filename + "'.";
        goto cleanup;
    }
    
	std::string::size_type delimiter = path_name.rfind('/');
	if(delimiter == std::string::npos)
	{
	    data_set_name = path_name;
	}
	else
	{
	    group_name = std::string(path_name.begin(), path_name.begin()+delimiter);
	    data_set_name = std::string(path_name.begin()+delimiter+1, path_name.end());
	}
	
	if(group_name != "")
	{
#if H5_VERS_MINOR <= 6
	    group = H5Gopen(m_file, group_name.c_str());
#else
	    group = H5Gopen(m_file, group_name.c_str(), H5P_DEFAULT);
#endif
	    if(group < 0)
	    {
	        message = std::string("HDF5ImportInfo(): Unable to open group '") + group_name + "'.";
	        goto cleanup;
	    }
	}
	else
	{
	    group = m_file;
	}
	
#if H5_VERS_MINOR <= 6
    m_dataset = H5Dopen(group, data_set_name.c_str());
#else
    m_dataset = H5Dopen(group, data_set_name.c_str(), H5P_DEFAULT);
#endif
    if(m_dataset < 0)
    {
        message = std::string("HDF5ImportInfo(): Unable to open data set '") + path + "'.";
        goto cleanup;
    }
    
    dspace = H5Dget_space(m_dataset);
    if(dspace < 0)
    {
        message = "HDF5ImportInfo(): Unable to open data space.";
        goto cleanup;
    }
    
    m_dimensions = H5Sget_simple_extent_ndims(dspace);
    if(m_dimensions < 2)
    {
        message = "HDF5ImportInfo(): Number of dimensions is lower than 2. Not an image!" ;
        goto cleanup;
    }
    
    m_dims.resize(m_dimensions);
    H5Sget_simple_extent_dims(dspace, m_dims.begin(), 0);
    
    dtype = H5Dget_type(m_dataset);
    if(dtype < 0)
    {
        message = "HDF5ImportInfo(): Unable to retrieve data type" ;
        goto cleanup;
    }
    native_type = H5Tget_native_type(dtype, H5T_DIR_ASCEND);
    if(native_type < 0)
    {
        message = "HDF5ImportInfo(): Unable to retrieve data type" ;
        goto cleanup;
    }

    if(native_type == H5T_NATIVE_CHAR)
        m_pixeltype = "INT8";
    if(native_type == H5T_NATIVE_SHORT)
        m_pixeltype = "INT16";
    if(native_type == H5T_NATIVE_INT && sizeof(int) == 4)
        m_pixeltype = "INT32";
    if(native_type == H5T_NATIVE_INT && sizeof(int) == 8)
        m_pixeltype = "INT64";
    if(native_type == H5T_NATIVE_LONG && sizeof(long) == 4)
        m_pixeltype = "INT32";
    if(native_type == H5T_NATIVE_LONG && sizeof(long) == 8)
        m_pixeltype = "INT64";
    if(native_type == H5T_NATIVE_LLONG && sizeof(long long) == 4)
        m_pixeltype = "INT32";
    if(native_type == H5T_NATIVE_LLONG && sizeof(long long) == 8)
        m_pixeltype = "INT64";

    if(native_type == H5T_NATIVE_UCHAR)
        m_pixeltype = "UINT8";
    if(native_type == H5T_NATIVE_USHORT)
        m_pixeltype = "UINT16";
    if(native_type == H5T_NATIVE_UINT && sizeof(unsigned int) == 4)
        m_pixeltype = "UINT32";
    if(native_type == H5T_NATIVE_UINT && sizeof(unsigned int) == 8)
        m_pixeltype = "UINT64";
    if(native_type == H5T_NATIVE_ULONG && sizeof(unsigned long) == 4)
        m_pixeltype = "UINT32";
    if(native_type == H5T_NATIVE_ULONG && sizeof(unsigned long) == 8)
        m_pixeltype = "UINT64";
    if(native_type == H5T_NATIVE_ULLONG && sizeof(unsigned long long) == 4)
        m_pixeltype = "UINT32";
    if(native_type == H5T_NATIVE_ULLONG && sizeof(unsigned long long) == 8)
        m_pixeltype = "UINT64";

    if(native_type == H5T_NATIVE_FLOAT)
        m_pixeltype = "FLOAT";
    if(native_type == H5T_NATIVE_DOUBLE)
        m_pixeltype = "DOUBLE";
    if(native_type == H5T_NATIVE_LDOUBLE)
        m_pixeltype = "LDOUBLE";
        
  cleanup:
    if(native_type)
        H5Tclose(native_type);
    if(dtype)
        H5Tclose(dtype);
    if(dspace)
        H5Sclose(dspace);
    if(group_name != "" && group)
        H5Gclose(group);
    if(message != "")
    {
        if(m_dataset)
            H5Dclose(m_dataset);
        if(m_file)
            H5Fclose(m_file);
        vigra_postcondition(false, message.c_str());
    }
}

HDF5ImportInfo::~HDF5ImportInfo()
{
    if(m_dataset)
        H5Dclose(m_dataset);
    if(m_file)
        H5Fclose(m_file);
}

HDF5ImportInfo::PixelType HDF5ImportInfo::pixelType() const
{
   const std::string pixeltype=HDF5ImportInfo::getPixelType();
   if (pixeltype == "UINT8")
       return HDF5ImportInfo::UINT8;
   if (pixeltype == "INT16")
     return HDF5ImportInfo::INT16;
   if (pixeltype == "UINT16")
     return HDF5ImportInfo::UINT16;
   if (pixeltype == "INT32")
     return HDF5ImportInfo::INT32;
   if (pixeltype == "UINT32")
     return HDF5ImportInfo::UINT32;
   if (pixeltype == "FLOAT")
     return HDF5ImportInfo::FLOAT;
   if (pixeltype == "DOUBLE")
     return HDF5ImportInfo::DOUBLE;
   vigra_fail( "internal error: unknown pixel type" );
   return HDF5ImportInfo::PixelType();
}

const char * HDF5ImportInfo::getPixelType() const
{
    return m_pixeltype.c_str();
}

MultiArrayIndex HDF5ImportInfo::shapeOfDimension(const int dim) const { return (MultiArrayIndex)m_dims[dim]; };
MultiArrayIndex HDF5ImportInfo::numDimensions() const { return m_dimensions; }
const std::string & HDF5ImportInfo::getDatasetName() const { return m_datasetname; }
const std::string & HDF5ImportInfo::getFileName() const { return m_filename; }
hid_t HDF5ImportInfo::getH5FileHandle() const { return m_file; }
hid_t HDF5ImportInfo::getDatasetHandle() const { return m_dataset; }

#endif // HasHDF5

} // namespace vigra
