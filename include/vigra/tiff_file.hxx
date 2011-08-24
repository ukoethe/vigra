#include <tiffio.h>
#include "vigra/tinyvector.hxx"

namespace vigra {

class TIFFFile
{
public:
	typedef TinyVector<uint32, 3> VolumeSize;

	TIFFFile(const char *filename, const char *mode);

	/**
	 * There are two possible TIFF file types that require different
	 * I/O APIs: striped files (which you know well) and tiled files
	 * (veeeery esoteric).
	 */
	enum FileType { INVALID = 0, STRIPED = 1, TILED = 2 };

	/**
	 * For usual TIFF files (no tile size fields), returns STRIPED.
	 * Returns TILED iff the TIFF file contains both tile width and
	 * height tags (and the tile depth tag for 3D files, i.e. if
	 * there's an image depth tag), INVALID otherwise (some, but not
	 * all extents specified).
	 */
	FileType fileType() const
	{
		return fileType_;
	}

	/**
	 * Returns true for 3D TIFFs, i.e. those with the SGI image/tile depth fields.
	 */
	bool hasDepth() const
	{
		return hasDepth_;
	}

	VolumeSize imageSize3D() const
	{
		return imageSize_;
	}

	VolumeSize tileSize3D() const
	{
		return tileSize_;
	}

private:
	TIFF *h_; // TIFF file handle (cf. libtiff library)

	TinyVector<uint32, 3> imageSize_, tileSize_;
	FileType fileType_;
	bool hasDepth_;
};

} // namespace vigra
