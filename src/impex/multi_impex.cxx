#include "vigra/multi_impex.hxx"

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

// class VolumeExportInfo

VolumeExportInfo::VolumeExportInfo( const char * name_base, const char * name_ext ) :
        m_x_res(0), m_y_res(0), m_z_res(0),
        m_filename_base(name_base), m_filename_ext(name_ext),
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

} // namespace vigra
