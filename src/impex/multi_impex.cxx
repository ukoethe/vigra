#include "vigra/multi_impex.hxx"
#include "codecmanager.hxx"

#include <string>
#include <sys/types.h>
#include <dirent.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

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


#ifdef _MSC_VER
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

// find filenames matching the pattern "<path>/base[0-9]+ext" (Windows version)
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

#else // not _MSC_VER

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

// find filenames matching the pattern "<path>/base[0-9]+ext" (Unix version)
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

#endif // _MSC_VER

void validate_filetype( std::string filetype )
{
    vigra_precondition( codecManager().fileTypeSupported(filetype),
                        "given file type is not supported" );
}

// class ImageExportInfo

ImageExportInfo::ImageExportInfo( const char * filename, const char * mode )
    : m_filename(filename), m_mode(mode),
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

const char * ImageExportInfo::getMode() const
{
    return m_mode.c_str();
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

vigra::Size2D ImageExportInfo::getCanvasSize() const
{
    return m_canvas_size ;
}

ImageExportInfo & ImageExportInfo::setCanvasSize(const Size2D & size)
{
    m_canvas_size = size;
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
VIGRA_UNIQUE_PTR<Encoder> encoder( const ImageExportInfo & info )
{
    VIGRA_UNIQUE_PTR<Encoder> enc;

    std::string filetype = info.getFileType();
    if ( filetype != "" ) {
        validate_filetype(filetype);
        VIGRA_UNIQUE_PTR<Encoder> enc2 = getEncoder( std::string( info.getFileName() ), filetype, std::string( info.getMode() ) );
        std::swap(enc, enc2);
    } else {
        VIGRA_UNIQUE_PTR<Encoder> enc2 = getEncoder( std::string( info.getFileName() ), "undefined", std::string( info.getMode() ) );
        std::swap(enc, enc2);
    }

    std::string comp = info.getCompression();
    if ( comp != "" ) {

        // check for quality parameter of JPEG compression
        int quality = 0;

        // possibility 1: quality specified as "JPEG QUALITY=N" or "JPEG-ARITH QUALITY=N"
        // possibility 2 (deprecated): quality specified as just a number "10"
        std::string sq(" QUALITY="), parsed_comp;
        std::string::size_type pos = comp.rfind(sq), start = 0;

        if(pos != std::string::npos)
        {
            start = pos + sq.size();
            parsed_comp = comp.substr(0, pos);
        }

        std::istringstream compstream(comp.substr(start));
        compstream >> quality;
        if ( quality != 0 )
        {
            if(parsed_comp == "")
                parsed_comp = "JPEG";
             enc->setCompressionType( parsed_comp, quality );
        }
        else
        {
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
    enc->setCanvasSize(info.getCanvasSize());

    if ( info.getICCProfile().size() > 0 ) {
        enc->setICCProfile(info.getICCProfile());
    }

    return enc;
}

// class ImageImportInfo

ImageImportInfo::ImageImportInfo( const char * filename, unsigned int imageIndex )
    : m_filename(filename), m_image_index(imageIndex)
{
    readHeader_();
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

int ImageImportInfo::numImages() const
{
    return m_num_images;
}

void ImageImportInfo::setImageIndex(int index)
{
    m_image_index = index;
    readHeader_();
}

int ImageImportInfo::getImageIndex() const
{
    return m_image_index;
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
    return (m_num_bands - m_num_extra_bands) == 1;
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

Size2D ImageImportInfo::getCanvasSize() const
{
    return m_canvas_size;
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

void ImageImportInfo::readHeader_()
{
    VIGRA_UNIQUE_PTR<Decoder> decoder = getDecoder(m_filename, "undefined", m_image_index);
    m_num_images = decoder->getNumImages();

    m_filetype = decoder->getFileType();
    m_pixeltype = decoder->getPixelType();
    m_width = decoder->getWidth();
    m_height = decoder->getHeight();
    m_num_bands = decoder->getNumBands();
    m_num_extra_bands = decoder->getNumExtraBands();
    m_pos = decoder->getPosition();
    m_canvas_size = decoder->getCanvasSize();
    m_x_res = decoder->getXResolution();
    m_y_res = decoder->getYResolution();

    m_icc_profile = decoder->getICCProfile();

    decoder->abort(); // there probably is no better way than this
}

// return a decoder for a given ImageImportInfo object
VIGRA_UNIQUE_PTR<Decoder> decoder( const ImageImportInfo & info )
{
    std::string filetype = info.getFileType();
    validate_filetype(filetype);
    return getDecoder( std::string( info.getFileName() ), filetype, info.getImageIndex() );
}

// class VolumeExportInfo

VolumeExportInfo::VolumeExportInfo( const char * filename ) :
        m_x_res(0), m_y_res(0), m_z_res(0),
        m_filetype("MULTIPAGE"),
        m_filename_base(filename), m_filename_ext(".tif"),
        fromMin_(0.0), fromMax_(0.0), toMin_(0.0), toMax_(0.0)
{
}

VolumeExportInfo::VolumeExportInfo( const char * name_base, const char * name_ext ) :
        m_x_res(0), m_y_res(0), m_z_res(0),
        m_filename_base(name_base), m_filename_ext(name_ext),
        fromMin_(0.0), fromMax_(0.0), toMin_(0.0), toMax_(0.0)
{
    if(m_filename_ext == "")
    {
        m_filename_ext = ".tif";
        m_filetype = "MULTIPAGE";
    }
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
    std::string message;
    
#if 0 // defined(HasHDF5)
    // Deactivate this code because it effectively forces multi_impex.hxx to be compiled with HDF5 only.
    
    // try reading from HDF5
    // (do this first because it uses mangled 'filename/dataset_name' format)
    {
        // split the filename into external (true filename) and internal (dataset name) part
        // FIXME: at present, the delimiter must be the extension '.h5' or '.hdf5' - define a more robust rule
        std::string name, datasetName;
        std::size_t ext = filename.rfind(".h5");
        if(ext != std::string::npos)
        {
            name = filename.substr(0, ext+3);
            datasetName = filename.substr(ext+3); 
        }
        else
        {
            ext = filename.rfind(".hdf5");
            if(ext != std::string::npos)
            {
                name = filename.substr(0, ext+5);
                datasetName = filename.substr(ext+5); 
            }
            else
            {
                name = filename;
            }
        }
        
        if(H5Fis_hdf5(name.c_str()))
        {
            message = std::string("VolumeImportInfo(): File '");
            message += name + "' is HDF5, but filename doesn't contain an internal dataset path.";
            vigra_precondition(datasetName.size() > 0, message.c_str());

            HDF5File hdf5file(name, HDF5File::OpenReadOnly);
            
            message = std::string("VolumeImportInfo(): Dataset '");
            message += datasetName + "' not found in HDF5 file '" + name + "'.";
            vigra_precondition(hdf5file.existsDataset(datasetName), message.c_str());
            
            ArrayVector<hsize_t> shape(hdf5file.getDatasetShape(datasetName));
            message = std::string("VolumeImportInfo(): Dataset '");
            message += datasetName + "' in HDF5 file '" + name + "' is not a volume.";
            vigra_precondition(shape.size() == 3 || shape.size() == 4, message.c_str());
            
            shape_[0] = shape[0];
            shape_[1] = shape[1];
            shape_[2] = shape[2];
            pixelType_ = hdf5file.getDatasetType(datasetName);
            numBands_ = shape.size() == 4
                            ? shape[3]
                            : 1;
            baseName_ = name;
            extension_ = datasetName;
            fileType_ = "HDF5";
            return;
        }
    }
#endif // HasHDF5

    // check if file exists
    message = std::string("VolumeImportInfo(): File '");
    message += filename + "' not found.";
#ifdef _MSC_VER
    vigra_precondition(_access(filename.c_str(), 0) != -1, message.c_str());
#else
    vigra_precondition(access(filename.c_str(), F_OK) == 0, message.c_str());
#endif

    // try Andor SIF format
    {
        std::string magic_string;
        {
            std::ifstream siffile (filename.c_str());
            if( !siffile.is_open() )
            {
                message = std::string("VolumeImportInfo(): Unable to open file '");
                message += filename + "'.";
                vigra_precondition(false, message.c_str());
            }    
            
            getline(siffile, magic_string);
        }
        if(magic_string == "Andor Technology Multi-Channel File")
        {
            SIFImportInfo info(filename.c_str());
            shape_[0] = info.shapeOfDimension(0);
            shape_[1] = info.shapeOfDimension(1);
            shape_[2] = info.shapeOfDimension(2);
            pixelType_ = "FLOAT";
            numBands_ = 1;
            baseName_ = filename;
            fileType_ = "SIF";
            return;            
        }
    }

    // for TIFF files, check for single-file, 3D data, tiled TIFFs first:
    if(codecManager().getFileTypeByMagicString(filename) == "TIFF")
    {
        TIFFFile f(filename.c_str(), "r");

        shape_ = f.imageSize3D();

        pixelType_ = f.pixelType();
        numBands_ = 1; // FIXME

        path_ = filename;

        baseName_ = filename;
        fileType_ = "TILEDTIFF";
        return;
    }

    // try multi-page TIFF or image stack
    if(isImage(filename.c_str()))
    {
        ImageImportInfo info(filename.c_str());
        shape_[0] = info.width();
        shape_[1] = info.height();
        resolution_[1] = -1.f; // assume images to be right-handed
        pixelType_ = info.getPixelType();
        numBands_ = info.numBands();
        
        if(info.numImages() > 1)
        {
            // must be a multi-page TIFF
            splitPathFromFilename(filename, path_, name_);
            baseName_ = filename;
            fileType_ = "MULTIPAGE";
            shape_[2] = info.numImages();
            return;
        }
        else
        {
            // try image stack loading
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
                        splitPathFromFilename(baseName, path_, name_);
                        baseName_ = baseName;
                        extension_ = extension;
                        shape_[2] = numbers.size();
                        std::swap(numbers, numbers_);
                        fileType_ = "STACK";
                        return;
                    }
                }
            }
            while(numEndIt != filename.rend());
        }
        
        message = std::string("VolumeImportInfo(): File '");
        message += filename + "' is neither a multi-page TIFF nor a valid image stack.";
        vigra_precondition(false, message.c_str());
    }
    
    {
        static std::string pixelTypes[] = { std::string("UNSIGNED_CHAR"), 
                                            std::string("UNSIGNED_BYTE"), 
                                            std::string("UINT8"), 
                                            std::string("INT16"), 
                                            std::string("UINT16"), 
                                            std::string("INT32"), 
                                            std::string("UINT32"), 
                                            std::string("FLOAT"), 
                                            std::string("DOUBLE"), 
                                            std::string() };
        
        // try .info file loading
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
                    std::string * type = pixelTypes;
                    while(*type != "")
                    {
                        if(*type == value)
                        {
                            pixelType_ = value;
                            break;
                        }
                        ++type;
                    }
                    vigra_precondition(*type != "",
                        "VolumeImportInfo(): Invalid datatype '" + value +"' in .info file.");
                    if(pixelType_ == "UNSIGNED_CHAR" || pixelType_ == "UNSIGNED_BYTE")
                        pixelType_ = "UINT8";
                }
                else if(key == "description")
                    description_ = value;
                else if(key == "name")
                    name_ = value;
                else if(key == "filename")
                    rawFilename_ = value;
                else
                {
                    std::cerr << "VolumeImportInfo(): WARNING: Unknown key '" << key
                              << "' (value '" << value << "') in info file!\n";
                }
            }
            else
            {
                if(line[0]) // non-empty line?
                    std::cerr << "VolumeImportInfo(): WARNING: could not parse line '" << line << "'!\n";
            }
        }

        if((shape_[0]*shape_[1]*shape_[2] > 0) && (rawFilename_.size() > 0))
        {
            numBands_ = 1;

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
            fileType_ = "RAW";
            return;
        }
    }

    message = std::string("VolumeImportInfo(): Unable to load file '");
    message += filename + "' - not a recognized format.";
    vigra_precondition(false, message.c_str());
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
    fileType_ = "STACK";
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
const char * VolumeImportInfo::getFileType() const
{
    return fileType_.c_str();
}
MultiArrayIndex VolumeImportInfo::numBands() const { return numBands_; }
bool VolumeImportInfo::isGrayscale() const { return numBands_ == 1; }
bool VolumeImportInfo::isColor() const { return numBands_ > 1; }
MultiArrayIndex VolumeImportInfo::width() const { return shape_[0]; }
MultiArrayIndex VolumeImportInfo::height() const { return shape_[1]; }
MultiArrayIndex VolumeImportInfo::depth() const { return shape_[2]; }
const std::string & VolumeImportInfo::name() const { return name_; }
const std::string & VolumeImportInfo::description() const { return description_; }

} // namespace vigra
