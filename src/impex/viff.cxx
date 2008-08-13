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

#include <iostream>
#include <fstream>
#include "vigra/config.hxx"
#include "vigra/sized_int.hxx"
#include "error.hxx"
#include "byteorder.hxx"
#include "void_vector.hxx"
#include "viff.hxx"

// VIFF definitions from the specification
// TODO: support complex types

/* definitions for version number,
   char release; */
#define XV_IMAGE_VER_NUM	3	/* Version 3 (3.1) */

/* definitions for release number,
   char version; */
#define XV_IMAGE_REL_NUM	1	/* Release 1   */

/* definitions for subimage information,
   long startx, starty; */
#define	VFF_NOTSUB   		~0	/* a negative number indicates that
					   the image is not a subimage	*/

/* definitions for machine dependencies,
   char machine_dep; */
#define	VFF_DEP_IEEEORDER   	0x2	/* IEEE byte ordering */
#define	VFF_DEP_DECORDER    	0x4	/* DEC (VAX) byte ordering */
#define VFF_DEP_NSORDER		0x8	/* NS32000 byte ordering */
#define VFF_DEP_CRAYORDER	0xA	/* Cray byte size and ordering */

#define	VFF_DEP_BIGENDIAN	VFF_DEP_IEEEORDER
#define	VFF_DEP_LITENDIAN	VFF_DEP_NSORDER


/* definitions for data storage type,
   unsigned long data_storage_type; */
#define	VFF_TYP_BIT		0	/* pixels are on or off (binary image)*/
                                        /* Note: This is an X11 XBitmap 
					   with bits packed into a byte and
					   padded to a byte */
#define	VFF_TYP_1_BYTE		1	/* pixels are byte (unsigned char) */
#define	VFF_TYP_2_BYTE		2	/* pixels are two byte (short int) */
#define	VFF_TYP_4_BYTE		4	/* pixels are four byte (integer) */
#define	VFF_TYP_FLOAT		5	/* pixels are float (single precision)*/
#define	VFF_TYP_COMPLEX		6	/* pixels are complex float */
#define VFF_TYP_DOUBLE		9	/* pixels are float (double precision)*/

#define VFF_TYP_DCOMPLEX	10	/* double complex */

/* definitions for data encoding scheme on disk - i.e. it may be
   compressed using RLE, or uncompressed (RAW).
   unsigned long data_encode_scheme; */
#define	VFF_DES_RAW		0	/* Raw - no compression */
#define VFF_DES_COMPRESS	1	/* Compressed using ALZ */
#define VFF_DES_RLE		2	/* Compressed using RLE */
#define VFF_DES_TRANSFORM	3	/* Transform based compression */
#define VFF_DES_CCITT		4	/* CCITT standard compression */
#define VFF_DES_ADPCM		5	/* ADPCM compression */
#define VFF_DES_GENERIC		6	/* User-specified compression */

/* definitions for map data or cells storage type,
   unsigned long map_storage_type; */
#define VFF_MAPTYP_NONE		0	/* No cell type is assigned  */
#define	VFF_MAPTYP_1_BYTE	1	/* cells are byte (unsigned char)    */
#define	VFF_MAPTYP_2_BYTE	2	/* cells are two byte (short int) */
#define	VFF_MAPTYP_4_BYTE	4	/* cells are four byte (integer) */
#define	VFF_MAPTYP_FLOAT	5	/* cells are float (single precision) */
#define	VFF_MAPTYP_COMPLEX	6	/* cells are complex FLOAT */
#define	VFF_MAPTYP_DOUBLE	7	/* cells are float (double precision) */

/* definitions for mapping schemes,
   unsigned long map_scheme; */
#define VFF_MS_NONE		0	/* No mapping is to be done, and no
					   maps are to be stored. */
#define	VFF_MS_ONEPERBAND	1	/* Each data band has its own map */
#define VFF_MS_CYCLE		2	/* An array of maps is selected in order
					   by groups of maps_per_cycle, allowing
					   "rotating the color map" */
#define	VFF_MS_SHARED		3	/* All data band share the same map */
#define VFF_MS_GROUP		4	/* All data bands are "grouped" 
					   together to point into one map */
/* definitions for enabling the map,
   unsigned long map_enable; */
#define VFF_MAP_OPTIONAL	1	/* The data is valid without being
					   sent thru the color map. If a
					   map is defined, the data may 
					   optionally be sent thru it. */
#define	VFF_MAP_FORCE		2	/* The data MUST be sent thru the map
					   to be interpreted */

/* definitions for color map models,
   unsigned long color_space_model; */

/*  the models use the following convention:
    ntsc: national television systems committee
    cie:  Commission Internationale de L'Eclairage
    ucs:  universal chromaticity scale
    RGB:  red band, green band, blue band
    CMY:  cyan band, magenta band, yellow band
    YIQ:  luminance, I and Q represent chrominance
    HSV:  hue, saturation, value
    IHS:  intensity, hue, saturation
    XYZ:
    UVW:  
    SOW:  
    Lab:
    Luv:

    For more information on how to interpret the combined meaning of the
    colorspace fields, see the document in $KHOROS_HOME/manual/viff_format,
    which contains detailed descriptions on the fields along with numerous
    examples of proper use.  */

#define VFF_CM_NONE	0
#define	VFF_CM_ntscRGB	1
#define	VFF_CM_ntscCMY	2
#define	VFF_CM_ntscYIQ	3
#define	VFF_CM_HSV	4
#define	VFF_CM_HLS	5
#define	VFF_CM_IHS	6
#define	VFF_CM_cieRGB	7
#define	VFF_CM_cieXYZ	8
#define	VFF_CM_cieUVW	9
#define	VFF_CM_cieucsUVW	10
#define	VFF_CM_cieucsSOW	11
#define	VFF_CM_cieucsLab	12
#define	VFF_CM_cieucsLuv	13
#define	VFF_CM_GENERIC		14	/* the color space is user defined */
#define VFF_CM_genericRGB	15	/* an RGB image but not conforming
					   to any standard */

/* definitions for location type,
   unsigned long location_type; */
#define	VFF_LOC_IMPLICIT	1	/*  The location of image pixels
					    or vector data is given by using
					    the implied 2D array given by
					    row_size and col_size.  */
#define	VFF_LOC_EXPLICIT	2	/*  The location of the image pixels
					    or the vectors is explicit */

namespace vigra {
    
    template< class T1, class T2 >
    class colormap
    {
    public:
        typedef T1 domain_type;
        typedef T2 value_type;
        typedef void_vector<value_type> vector_type;

    private:
        // size of lookup tables
        const unsigned int m_numTableElements;

        // number of lookup tables
        const unsigned int m_numTables;

        // number of bands in each lookup table
        const unsigned int m_numTableBands;

        // lookup tables
        vector_type m_tables;

    public:
        colormap( const unsigned int numTableElements,
                  const unsigned int numTables,
                  const unsigned int numTableBands )
            : m_numTableElements(numTableElements),
              m_numTables(numTables), m_numTableBands(numTableBands),
              m_tables( numTableElements * numTableBands )
        {
            vigra_precondition( ( numTables == 1 ) || ( numTableBands == 1 ),
                                "numTables or numTableBands must be 1" );
        }

        inline unsigned int getNumTables() const
        {
            return m_numTables;
        }

        inline unsigned int getNumElements() const
        {
            return m_numTableElements;
        }

        // initialize a single table, data array must have the correct length
        void initialize( const value_type * data, const unsigned int table )
        {
            vigra_precondition( table < m_numTables,
                                "table number out of range" );
            const unsigned int stride = m_numTableBands * m_numTableElements;
            std::copy( data, data + stride, m_tables.data() + stride * table );
        }

        // map a single value
        value_type operator()( const domain_type index,
                               const unsigned int band = 0 ) const
        {
            vigra_precondition( index < m_numTableElements,
                                "index out of range" );
            if ( m_numTables == 1 ) {
                // map bands with a single (interleaved or not) table
                vigra_precondition( band < m_numTableBands,
                                    "band out of range" );
                return m_tables[ m_numTableElements * band + index ];
            } else {
                // map bands with multiple, non-interleaved tables
                vigra_precondition( band < m_numTables, "band out of range" );
                const unsigned int stride
                    = m_numTableBands * m_numTableElements;
                return m_tables[ stride * band + index ];
            }
        }
    };


    // this function encapsulates the colormap functor
    template< class storage_type, class map_storage_type >
    void map_multiband( void_vector<map_storage_type> & dest,
                        unsigned int & dest_bands,
                        const void_vector<storage_type> & src,
                        unsigned int src_bands, unsigned int src_width, unsigned int src_height,
                        const void_vector<map_storage_type> & maps,
                        unsigned int map_bands, unsigned int map_width, unsigned int map_height )
    {
        typedef colormap< storage_type, map_storage_type > colormap_type;
        const unsigned int num_pixels = src_width * src_height;

        // build the color map
        const unsigned int map_band_size = map_width * map_height;
        colormap_type colormap( map_height, map_bands, map_width );
        for ( unsigned int i = 0; i < map_bands; ++i )
            colormap.initialize( maps.data() + map_band_size * i, i );

        // map each pixel
        const unsigned int band_size = src_width * src_height;
        dest_bands = map_bands * map_width;
        dest.resize( band_size * dest_bands );
        if ( map_width > 1 ) {
            // interleaved maps => there is only one band in the image
            for( unsigned int bandnum = 0; bandnum < dest_bands; ++bandnum )
                for( unsigned int i = 0; i < num_pixels; ++i )
                    dest[ bandnum * band_size + i ]
                        = colormap( src[i], bandnum );
        } else {
            // non-interleaved bands => map can be used per band
            for( unsigned int bandnum = 0; bandnum < dest_bands; ++bandnum )
                for( unsigned int i = 0; i < num_pixels; ++i )
                    dest[ bandnum * band_size + i ]
                        = colormap( src[ bandnum * band_size + i ], bandnum );
        }
    }

    CodecDesc ViffCodecFactory::getCodecDesc() const
    {
        CodecDesc desc;

        // init file type
        desc.fileType = "VIFF";

        // init pixel types
        desc.pixelTypes.resize(5);
        desc.pixelTypes[0] = "UINT8";
        desc.pixelTypes[1] = "INT16";
        desc.pixelTypes[2] = "INT32";
        desc.pixelTypes[3] = "FLOAT";
        desc.pixelTypes[4] = "DOUBLE";

        // init compression types
        desc.compressionTypes.resize(0);

        // init magic strings
        desc.magicStrings.resize(1);
        desc.magicStrings[0].resize(2);
        desc.magicStrings[0][0] = '\253'; // XV_FILE_MAGIC_NUM
        desc.magicStrings[0][1] = '\1';   // XV_FILE_TYPE_XVIFF

        // init file extensions
        desc.fileExtensions.resize(1);
        desc.fileExtensions[0] = "xv";

        desc.bandNumbers.resize(1);
        desc.bandNumbers[0] = 0;
        
        return desc;
    }

    std::auto_ptr<Decoder> ViffCodecFactory::getDecoder() const
    {
        return std::auto_ptr<Decoder>( new ViffDecoder() );
    }

    std::auto_ptr<Encoder> ViffCodecFactory::getEncoder() const
    {
        return std::auto_ptr<Encoder>( new ViffEncoder() );
    }

    class ViffHeader
    {
    public:
        typedef UInt32 field_type;

        field_type row_size, col_size, num_data_bands, data_storage_type,
            data_encode_scheme, map_scheme, map_storage_type, map_row_size,
            map_col_size;

        void from_stream( std::ifstream & stream, byteorder & bo );
        void to_stream( std::ofstream & stream, byteorder & bo ) const;
    };

    void ViffHeader::from_stream( std::ifstream & stream, byteorder & bo )
    {
        // scratch variables for values that do not need to be saved
        field_type scratch;

        // skip the magic number and the file type
        stream.seekg( 2, std::ios::cur );

        // check the release number
        vigra_precondition( stream.get() == XV_IMAGE_REL_NUM,
                            "file format release unsupported" );

        // check the version number
        vigra_precondition( stream.get() == XV_IMAGE_VER_NUM,
                            "file format version unsupported" );

        // check the endianess
        const char machine_dep = stream.get();
        if ( machine_dep == VFF_DEP_BIGENDIAN )
            bo.set("big endian");
        else if ( machine_dep == VFF_DEP_LITENDIAN )
            bo.set("little endian");
        else vigra_fail( "endianess unsupported" );

        // skip the comments
        stream.seekg( 0x208, std::ios::beg );

        // read the row size
        read_field( stream, bo, row_size );

        // read the col size
        read_field( stream, bo, col_size );

        // skip subrow size, subimage (x,y) and pixelsize (x,y)
        stream.seekg( 20, std::ios::cur );

        read_field( stream, bo, scratch );
        vigra_precondition( scratch != VFF_LOC_EXPLICIT,
                            "explicit locations are unsupported" );

        // skip location dim field (only used for explicit locations)
        stream.seekg( 4, std::ios::cur );

        // check number of images
        read_field( stream, bo, scratch );
        vigra_precondition( scratch < 2,
                            "multiple images are unsupported" );

        read_field( stream, bo, num_data_bands );

        read_field( stream, bo, data_storage_type );
        vigra_precondition( data_storage_type != VFF_TYP_BIT,
                            "bit storage type unsupported" );
        vigra_precondition( data_storage_type != VFF_TYP_COMPLEX,
                            "complex storage type unsupported" );
        vigra_precondition( data_storage_type != VFF_TYP_DCOMPLEX,
                            "double complex storage type unsupported" );

        read_field( stream, bo, data_encode_scheme );
        vigra_precondition( data_encode_scheme == VFF_DES_RAW,
                            "data compression unsupported" );

        read_field( stream, bo, map_scheme );
        vigra_precondition( map_scheme != VFF_MS_CYCLE,
                            "map cycling unsupported" );

        // only read the map fields if mapping is actually used
        if ( map_scheme != VFF_MS_NONE ) {

            read_field( stream, bo, map_storage_type );
            vigra_precondition( map_storage_type != VFF_MAPTYP_COMPLEX,
                                "complex storage type unsupported" );
            vigra_precondition( map_storage_type != VFF_MAPTYP_NONE,
                                "invalid maptype" );

            read_field( stream, bo, map_row_size );
            read_field( stream, bo, map_col_size );
        }
            
        // seek behind the header. (skip colorspace and pointers)
        stream.seekg( 1024, std::ios::beg );
    }

#if  __GNUC__ == 2
    #define VIGRA_STREAM_CHAR_TYPE char
#else
    #define VIGRA_STREAM_CHAR_TYPE std::ofstream::char_type
#endif

    void ViffHeader::to_stream( std::ofstream & stream, byteorder & bo ) const
    {
        field_type scratch, null = 0;

        // magic number
        stream.put((VIGRA_STREAM_CHAR_TYPE)0xAB);
            
        // file type
        stream.put((VIGRA_STREAM_CHAR_TYPE)0x01);

        // file format release number
        stream.put(XV_IMAGE_REL_NUM);

        // file format version number
        stream.put(XV_IMAGE_VER_NUM);

        // set the byte order
        if ( bo.get_host_byteorder() == "big endian" )
        {
            bo.set("big endian" );
            stream.put(VFF_DEP_BIGENDIAN);
        }
        else 
        {
            bo.set("little endian" );
            stream.put(VFF_DEP_LITENDIAN);
        }

        // pad, then zero out the comments
        unsigned int i;
        for(i = 0; i < 515; ++i )
            stream.put(0);

        // image size
        write_field( stream, bo, row_size );
        write_field( stream, bo, col_size );

        // zero out five fields
        for( i = 0; i < 20; ++i )
            stream.put(0);

        // using implicit locations
        write_field( stream, bo, scratch = VFF_LOC_IMPLICIT );

        // zero out one field
        write_field( stream, bo, null );

        // set number of images
        write_field( stream, bo, scratch = 1 );

        // image data format
        write_field( stream, bo, num_data_bands );
        write_field( stream, bo, data_storage_type );
        write_field( stream, bo, scratch = VFF_DES_RAW );

        // no color mapping
        write_field( stream, bo, scratch = VFF_MS_NONE );
        write_field( stream, bo, scratch = VFF_MAPTYP_NONE );

        // zero out five fields
        for( i = 0; i < 20; ++i )
            stream.put(0);

        // colorspace
        scratch = num_data_bands == 3 ? VFF_CM_genericRGB: VFF_CM_NONE;
        write_field( stream, bo, scratch );

        // zero out the last bytes of the header
        int offset = 1024 - stream.tellp();
        vigra_precondition( offset >= 0,
                            "machine is incapable to read viff" );
        for( int j = 0; j < offset; ++j )
            stream.put(0);
    }

    struct ViffDecoderImpl
    {
        unsigned int width, height, components, map_width,
            map_height, num_maps;
        std::string pixelType;
        int current_scanline;

        ViffHeader header;
        void_vector_base maps, bands;

        ViffDecoderImpl( const std::string & filename );

        void read_maps( std::ifstream & stream, byteorder & bo );
        void read_bands( std::ifstream & stream, byteorder & bo );
        void color_map();
    };

    ViffDecoderImpl::ViffDecoderImpl( const std::string & filename )
        : pixelType("undefined"), current_scanline(-1)
    {
#ifdef VIGRA_NEED_BIN_STREAMS
        std::ifstream stream( filename.c_str(), std::ios::binary );
#else
        std::ifstream stream( filename.c_str() );
#endif
        
        if(!stream.good())
        {
            std::string msg("Unable to open file '");
            msg += filename;
            msg += "'.";
            vigra_precondition(0, msg.c_str());
        }
        byteorder bo( "big endian" );

        // get header
        header.from_stream( stream, bo );
        width = header.row_size;
        height = header.col_size;
        components = header.num_data_bands;

        // read data and eventually map it
        if ( header.map_scheme != VFF_MS_NONE )
            read_maps( stream, bo );
        read_bands( stream, bo );
        if ( header.map_scheme != VFF_MS_NONE )
            color_map();
    }

    void ViffDecoderImpl::read_maps( std::ifstream & stream, byteorder & bo )
    {
        const bool shared_map = ( header.map_scheme == VFF_MS_SHARED );
        num_maps = shared_map ? 1 : header.num_data_bands;
        map_width = header.map_row_size;
        map_height = header.map_col_size;
        const unsigned int maps_size = map_width * map_height * num_maps;

        if ( header.map_storage_type == VFF_MAPTYP_1_BYTE ) {
            typedef void_vector<UInt8> maps_type;
            maps_type & castmaps = static_cast< maps_type & >(maps);
            castmaps.resize(maps_size);
            read_array( stream, bo, castmaps.data(), maps_size );
        } else if ( header.map_storage_type == VFF_MAPTYP_2_BYTE ) {
            typedef void_vector<Int16> maps_type;
            maps_type & castmaps = static_cast< maps_type & >(maps);
            castmaps.resize(maps_size);
            read_array( stream, bo, castmaps.data(), maps_size );
        } else if ( header.map_storage_type == VFF_MAPTYP_4_BYTE ) {
            typedef void_vector<Int32> maps_type;
            maps_type & castmaps = static_cast< maps_type & >(maps);
            castmaps.resize(maps_size);
            read_array( stream, bo, castmaps.data(), maps_size );
        } else if ( header.map_storage_type == VFF_MAPTYP_FLOAT ) {
            typedef void_vector<float> maps_type;
            maps_type & castmaps = static_cast< maps_type & >(maps);
            castmaps.resize(maps_size);
            read_array( stream, bo, castmaps.data(), maps_size );
        } else
            vigra_precondition( false, "map storage type unsupported" );
    }

    void ViffDecoderImpl::read_bands( std::ifstream & stream, byteorder & bo )
    {
        const unsigned int bands_size = width * height * components;

        if ( header.data_storage_type == VFF_TYP_1_BYTE ) {
            typedef void_vector<UInt8> bands_type;
            bands_type & castbands = static_cast< bands_type & >(bands);
            castbands.resize(bands_size);
            read_array( stream, bo, castbands.data(), bands_size );
            pixelType = "UINT8";
        } else if ( header.data_storage_type == VFF_TYP_2_BYTE ) {
            typedef void_vector<Int16> bands_type;
            bands_type & castbands = static_cast< bands_type & >(bands);
            castbands.resize(bands_size);
            read_array( stream, bo, castbands.data(), bands_size );
            pixelType = "INT16";
        } else if ( header.data_storage_type == VFF_TYP_4_BYTE ) {
            typedef void_vector<Int32> bands_type;
            bands_type & castbands = static_cast< bands_type & >(bands);
            castbands.resize(bands_size);
            read_array( stream, bo, castbands.data(), bands_size );
            pixelType = "INT32";
        } else if ( header.data_storage_type == VFF_TYP_FLOAT ) {
            typedef void_vector<float> bands_type;
            bands_type & castbands = static_cast< bands_type & >(bands);
            castbands.resize(bands_size);
            read_array( stream, bo, castbands.data(), bands_size );
            pixelType = "FLOAT";
        } else if ( header.data_storage_type == VFF_TYP_DOUBLE ) {
            typedef void_vector<double> bands_type;
            bands_type & castbands = static_cast< bands_type & >(bands);
            castbands.resize(bands_size);
            read_array( stream, bo, castbands.data(), bands_size );
            pixelType = "DOUBLE";
        } else
            vigra_precondition( false, "storage type unsupported" );
    }

    void ViffDecoderImpl::color_map()
    {
        void_vector_base temp_bands;
        unsigned int temp_num_bands;

        if ( header.map_storage_type == VFF_MAPTYP_1_BYTE ) {
            typedef UInt8 map_storage_type;

            if ( header.data_storage_type == VFF_TYP_1_BYTE ) {
                typedef UInt8 storage_type;
                typedef void_vector<storage_type> bands_type;
                typedef void_vector<map_storage_type> maps_type;
                map_multiband( static_cast< maps_type & >(temp_bands),
                               temp_num_bands,
                               static_cast< const bands_type & >(bands),
                               components, width, height,
                               static_cast< const maps_type & >(maps),
                               num_maps, map_width, map_height );

            } else if ( header.data_storage_type == VFF_TYP_2_BYTE ) {
                typedef UInt16 storage_type;
                typedef void_vector<storage_type> bands_type;
                typedef void_vector<map_storage_type> maps_type;
                map_multiband( static_cast< maps_type & >(temp_bands),
                               temp_num_bands,
                               static_cast< const bands_type & >(bands),
                               components, width, height,
                               static_cast< const maps_type & >(maps),
                               num_maps, map_width, map_height );

            } else if ( header.data_storage_type == VFF_TYP_4_BYTE ) {
                typedef UInt32 storage_type;
                typedef void_vector<storage_type> bands_type;
                typedef void_vector<map_storage_type> maps_type;
                map_multiband( static_cast< maps_type & >(temp_bands),
                               temp_num_bands,
                               static_cast< const bands_type & >(bands),
                               components, width, height,
                               static_cast< const maps_type & >(maps),
                               num_maps, map_width, map_height );

            } else
                vigra_precondition( false, "storage type unsupported" );

            pixelType = "UINT8";

        } else if ( header.map_storage_type == VFF_MAPTYP_2_BYTE ) {
            typedef UInt16 map_storage_type;

            if ( header.data_storage_type == VFF_TYP_1_BYTE ) {
                typedef UInt8 storage_type;
                typedef void_vector<storage_type> bands_type;
                typedef void_vector<map_storage_type> maps_type;
                map_multiband( static_cast< maps_type & >(temp_bands),
                               temp_num_bands,
                               static_cast< const bands_type & >(bands),
                               components, width, height,
                               static_cast< const maps_type & >(maps),
                               num_maps, map_width, map_height );

            } else if ( header.data_storage_type == VFF_TYP_2_BYTE ) {
                typedef UInt16 storage_type;
                typedef void_vector<storage_type> bands_type;
                typedef void_vector<map_storage_type> maps_type;
                map_multiband( static_cast< maps_type & >(temp_bands),
                               temp_num_bands,
                               static_cast< const bands_type & >(bands),
                               components, width, height,
                               static_cast< const maps_type & >(maps),
                               num_maps, map_width, map_height );

            } else if ( header.data_storage_type == VFF_TYP_4_BYTE ) {
                typedef UInt32 storage_type;
                typedef void_vector<storage_type> bands_type;
                typedef void_vector<map_storage_type> maps_type;
                map_multiband( static_cast< maps_type & >(temp_bands),
                               temp_num_bands,
                               static_cast< const bands_type & >(bands),
                               components, width, height,
                               static_cast< const maps_type & >(maps),
                               num_maps, map_width, map_height );

            } else
                vigra_precondition( false, "storage type unsupported" );

            pixelType = "INT16";

        } else if ( header.map_storage_type == VFF_MAPTYP_4_BYTE ) {
            typedef UInt32 map_storage_type;

            if ( header.data_storage_type == VFF_TYP_1_BYTE ) {
                typedef UInt8 storage_type;
                typedef void_vector<storage_type> bands_type;
                typedef void_vector<map_storage_type> maps_type;
                map_multiband( static_cast< maps_type & >(temp_bands),
                               temp_num_bands,
                               static_cast< const bands_type & >(bands),
                               components, width, height,
                               static_cast< const maps_type & >(maps),
                               num_maps, map_width, map_height );

            } else if ( header.data_storage_type == VFF_TYP_2_BYTE ) {
                typedef UInt16 storage_type;
                typedef void_vector<storage_type> bands_type;
                typedef void_vector<map_storage_type> maps_type;
                map_multiband( static_cast< maps_type & >(temp_bands),
                               temp_num_bands,
                               static_cast< const bands_type & >(bands),
                               components, width, height,
                               static_cast< const maps_type & >(maps),
                               num_maps, map_width, map_height );

            } else if ( header.data_storage_type == VFF_TYP_4_BYTE ) {
                typedef UInt32 storage_type;
                typedef void_vector<storage_type> bands_type;
                typedef void_vector<map_storage_type> maps_type;
                map_multiband( static_cast< maps_type & >(temp_bands),
                               temp_num_bands,
                               static_cast< const bands_type & >(bands),
                               components, width, height,
                               static_cast< const maps_type & >(maps),
                               num_maps, map_width, map_height );

            } else
                vigra_precondition( false, "storage type unsupported" );

            pixelType = "INT32";

        } else if ( header.map_storage_type == VFF_MAPTYP_FLOAT ) {
            typedef float map_storage_type;

            if ( header.data_storage_type == VFF_TYP_1_BYTE ) {
                typedef UInt8 storage_type;
                typedef void_vector<storage_type> bands_type;
                typedef void_vector<map_storage_type> maps_type;
                map_multiband( static_cast< maps_type & >(temp_bands),
                               temp_num_bands,
                               static_cast< const bands_type & >(bands),
                               components, width, height,
                               static_cast< const maps_type & >(maps),
                               num_maps, map_width, map_height );

            } else if ( header.data_storage_type == VFF_TYP_2_BYTE ) {
                typedef UInt16 storage_type;
                typedef void_vector<storage_type> bands_type;
                typedef void_vector<map_storage_type> maps_type;
                map_multiband( static_cast< maps_type & >(temp_bands),
                               temp_num_bands,
                               static_cast< const bands_type & >(bands),
                               components, width, height,
                               static_cast< const maps_type & >(maps),
                               num_maps, map_width, map_height );

            } else if ( header.data_storage_type == VFF_TYP_4_BYTE ) {
                typedef UInt32 storage_type;
                typedef void_vector<storage_type> bands_type;
                typedef void_vector<map_storage_type> maps_type;
                map_multiband( static_cast< maps_type & >(temp_bands),
                               temp_num_bands,
                               static_cast< const bands_type & >(bands),
                               components, width, height,
                               static_cast< const maps_type & >(maps),
                               num_maps, map_width, map_height );
 
            } else
                vigra_precondition( false, "storage type unsupported" );

            pixelType = "FLOAT";

        } else
            vigra_precondition( false, "map storage type unsupported" );

        swap_void_vector( bands, temp_bands );
        components = temp_num_bands;
    }

    void ViffDecoder::init( const std::string & filename )
    {
        pimpl = new ViffDecoderImpl(filename);
    }

    ViffDecoder::~ViffDecoder()
    {
        delete pimpl;
    }

    std::string ViffDecoder::getFileType() const
    {
        return "VIFF";
    }

    unsigned int ViffDecoder::getWidth() const
    {
        return pimpl->width;
    }

    unsigned int ViffDecoder::getHeight() const
    {
        return pimpl->height;
    }

    unsigned int ViffDecoder::getNumBands() const
    {
        return pimpl->components;
    }

    std::string ViffDecoder::getPixelType() const
    {
        return pimpl->pixelType;
    }

    unsigned int ViffDecoder::getOffset() const
    {
        return 1; // this is multiband data
    }

    const void * ViffDecoder::currentScanlineOfBand( unsigned int band ) const
    {
        const unsigned int index = pimpl->width
            * ( pimpl->height * band + pimpl->current_scanline );
        if ( pimpl->pixelType == "UINT8" ) {
            typedef void_vector<UInt8> bands_type;
            const bands_type & bands
                = static_cast< const bands_type & >(pimpl->bands);
            return bands.data() + index;
        } else if ( pimpl->pixelType == "INT16" ) {
            typedef void_vector<Int16> bands_type;
            const bands_type & bands
                = static_cast< const bands_type & >(pimpl->bands);
            return bands.data() + index;
        } else if ( pimpl->pixelType == "INT32" ) {
            typedef void_vector<Int32> bands_type;
            const bands_type & bands
                = static_cast< const bands_type & >(pimpl->bands);
            return bands.data() + index;
        } else if ( pimpl->pixelType == "FLOAT" ) {
            typedef void_vector<float> bands_type;
            const bands_type & bands
                = static_cast< const bands_type & >(pimpl->bands);
            return bands.data() + index;
        } else if ( pimpl->pixelType == "DOUBLE" ) {
            typedef void_vector<double> bands_type;
            const bands_type & bands
                = static_cast< const bands_type & >(pimpl->bands);
            return bands.data() + index;
        } else {
            vigra_fail( "PixelType was not set correctly" );
            return 0;
        }
    }

    void ViffDecoder::nextScanline()
    {
        ++(pimpl->current_scanline);
    }

    void ViffDecoder::close() {}
    void ViffDecoder::abort() {}

    // VIFF Encoder

    struct ViffEncoderImpl
    {
        std::ofstream stream;
        byteorder bo;
        std::string pixelType;
        int current_scanline;
        bool finalized;

        ViffHeader header;
        void_vector_base bands;

        ViffEncoderImpl( const std::string & filename )
#ifdef VIGRA_NEED_BIN_STREAMS
            : stream( filename.c_str(), std::ios::binary ), 
#else
            : stream( filename.c_str() ), 
#endif
              bo( "big endian" ),
              pixelType("undefined"), current_scanline(0), finalized(false)
        {
            if(!stream.good())
            {
                std::string msg("Unable to open file '");
                msg += filename;
                msg += "'.";
                vigra_precondition(false, msg);
            }
        }

    };

    void ViffEncoder::init( const std::string & filename )
    {
        pimpl = new ViffEncoderImpl(filename);
    }

    ViffEncoder::~ViffEncoder()
    {
        delete pimpl;
    }

    std::string ViffEncoder::getFileType() const
    {
        return "VIFF";
    }

    void ViffEncoder::setWidth( unsigned int width )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->header.row_size = static_cast<ViffHeader::field_type>(width);
    }

    void ViffEncoder::setHeight( unsigned int height )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->header.col_size = static_cast<ViffHeader::field_type>(height);
    }

    void ViffEncoder::setNumBands( unsigned int numBands )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->header.num_data_bands = static_cast<ViffHeader::field_type>(numBands);
    }

    void ViffEncoder::setCompressionType( const std::string & comp,
                                          int quality )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
    }

    void ViffEncoder::setPixelType( const std::string & pixelType )
    {
        VIGRA_IMPEX_FINALIZED(pimpl->finalized);
        pimpl->pixelType = pixelType;
        if ( pixelType == "UINT8" )
            pimpl->header.data_storage_type = VFF_TYP_1_BYTE;
        else if ( pixelType == "INT16" )
            pimpl->header.data_storage_type = VFF_TYP_2_BYTE;
        else if ( pixelType == "INT32" )
            pimpl->header.data_storage_type = VFF_TYP_4_BYTE;
        else if ( pixelType == "FLOAT" )
            pimpl->header.data_storage_type = VFF_TYP_FLOAT;
        else if ( pixelType == "DOUBLE" )
            pimpl->header.data_storage_type = VFF_TYP_DOUBLE;
    }

    unsigned int ViffEncoder::getOffset() const
    {
        return 1;
    }

    void ViffEncoder::finalizeSettings()
    {
        pimpl->header.to_stream( pimpl->stream, pimpl->bo );

        const unsigned int bands_size = pimpl->header.row_size
            * pimpl->header.col_size * pimpl->header.num_data_bands;

        if ( pimpl->header.data_storage_type == VFF_TYP_1_BYTE ) {
            typedef void_vector<UInt8> bands_type;
            bands_type & castbands = static_cast< bands_type & >(pimpl->bands);
            castbands.resize(bands_size);
        } else if ( pimpl->header.data_storage_type == VFF_TYP_2_BYTE ) {
            typedef void_vector<Int16> bands_type;
            bands_type & castbands = static_cast< bands_type & >(pimpl->bands);
            castbands.resize(bands_size);
        } else if ( pimpl->header.data_storage_type == VFF_TYP_4_BYTE ) {
            typedef void_vector<Int32> bands_type;
            bands_type & castbands = static_cast< bands_type & >(pimpl->bands);
            castbands.resize(bands_size);
        } else if ( pimpl->header.data_storage_type == VFF_TYP_FLOAT ) {
            typedef void_vector<float> bands_type;
            bands_type & castbands = static_cast< bands_type & >(pimpl->bands);
            castbands.resize(bands_size);
        } else if ( pimpl->header.data_storage_type == VFF_TYP_DOUBLE ) {
            typedef void_vector<double> bands_type;
            bands_type & castbands = static_cast< bands_type & >(pimpl->bands);
            castbands.resize(bands_size);
        } else
            vigra_precondition( false, "storage type unsupported" );

        pimpl->finalized = true;
    }

    void * ViffEncoder::currentScanlineOfBand( unsigned int band )
    {
        const unsigned int index = pimpl->header.row_size
            * ( pimpl->header.col_size * band + pimpl->current_scanline );
        if ( pimpl->pixelType == "UINT8" ) {
            typedef void_vector<UInt8> bands_type;
             bands_type & bands = static_cast<  bands_type & >(pimpl->bands);
            return bands.data() + index;
        } else if ( pimpl->pixelType == "INT16" ) {
            typedef void_vector<Int16> bands_type;
             bands_type & bands = static_cast<  bands_type & >(pimpl->bands);
            return bands.data() + index;
        } else if ( pimpl->pixelType == "INT32" ) {
            typedef void_vector<Int32> bands_type;
             bands_type & bands = static_cast<  bands_type & >(pimpl->bands);
            return bands.data() + index;
        } else if ( pimpl->pixelType == "FLOAT" ) {
            typedef void_vector<float> bands_type;
             bands_type & bands = static_cast<  bands_type & >(pimpl->bands);
            return bands.data() + index;
        } else if ( pimpl->pixelType == "DOUBLE" ) {
            typedef void_vector<double> bands_type;
             bands_type & bands = static_cast<  bands_type & >(pimpl->bands);
            return bands.data() + index;
        } else {
            vigra_postcondition( false, "PixelType was not set correctly" );
            return 0;
        }
    }

    void ViffEncoder::nextScanline()
    {
        ++(pimpl->current_scanline);
    }

    void ViffEncoder::close()
    {
        // write bands to disk
        const unsigned int bands_size = pimpl->header.row_size
            * pimpl->header.col_size * pimpl->header.num_data_bands;

        if ( pimpl->header.data_storage_type == VFF_TYP_1_BYTE ) {
            typedef void_vector<UInt8> bands_type;
            bands_type & castbands = static_cast< bands_type & >(pimpl->bands);
            write_array( pimpl->stream, pimpl->bo,
                         castbands.data(), bands_size );
        } else if ( pimpl->header.data_storage_type == VFF_TYP_2_BYTE ) {
            typedef void_vector<Int16> bands_type;
            bands_type & castbands = static_cast< bands_type & >(pimpl->bands);
            write_array( pimpl->stream, pimpl->bo,
                         castbands.data(), bands_size );
        } else if ( pimpl->header.data_storage_type == VFF_TYP_4_BYTE ) {
            typedef void_vector<Int32> bands_type;
            bands_type & castbands = static_cast< bands_type & >(pimpl->bands);
            write_array( pimpl->stream, pimpl->bo,
                         castbands.data(), bands_size );
        } else if ( pimpl->header.data_storage_type == VFF_TYP_FLOAT ) {
            typedef void_vector<float> bands_type;
            bands_type & castbands = static_cast< bands_type & >(pimpl->bands);
            write_array( pimpl->stream, pimpl->bo,
                         castbands.data(), bands_size );
        } else if ( pimpl->header.data_storage_type == VFF_TYP_DOUBLE ) {
            typedef void_vector<double> bands_type;
            bands_type & castbands = static_cast< bands_type & >(pimpl->bands);
            write_array( pimpl->stream, pimpl->bo,
                         castbands.data(), bands_size );
        } else
            vigra_precondition( false, "storage type unsupported" );
    }

    void ViffEncoder::abort() {}

} // namespace vigra
