/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#ifndef _viff_h_
#define _viff_h_


/*
******************************************************************
	Khoros Visualization/Image File Format		         
							         
	A Khoros data file is organized as a 1Kbyte header     
	followed by additional information.  The first two       
	bytes of the header tell what kind of stuff the          
	additional information is.  For Khoros image files,     
	the additional information consists of the maps,   
	the location data, and then the image or vector data.    
                                                                
	There is a supporting document for this file called      
	"The Khoros Visualization/Image File Format" that will help to        
	expalin the various fields, see $KHOROS_HOME/manual/viff_format.

	The header fields where carefully selected to            
	prevent contradictions between categories, i.e. they     
	were chosen to be orthogonal to each other.  However,
	in several situations this causes the fields to supply
	redundant information.
								 
	Note that the structure contains pointers to the various 
	chunks of data.  These will make sense ONLY when the 
	data is in memory.  

	*imagedata - points to a sequence of images, an image is
	made up of bands, and the bands are in a sequence; or
	it can point to vectors, where the vector dimension is
	the number of bands.

	*maps - points to a sequence of 2-dimensional maps, a map 
	is organized as stacked columns.  A data value indexes map rows.

	*location - points to bands of coordinate values, ie if
	two dimensional locations, then there would be a band
	of x's followed by a band of y's.
								 
	The Khoros convention for the image orientarion is	 
	with the image origin in the upper left hand corner.     
							         
******************************************************************
*/


/* 
 * The VIFF header size for a Cray is 4096 because of the
 * word size on a cray. So, set the headersize to be 4096
 * for the cray and 1024 of all other machines.
 */

#if defined(CRAY)
#   define VIFF_HEADERSIZE 4096
#else
#   define VIFF_HEADERSIZE 1024
#endif

#define VIFF_HEADER_RESERVESIZE (VIFF_HEADERSIZE - (21*sizeof(long))  \
			     - (520*sizeof(char)) - (2*sizeof(char *)) - (4*sizeof(float)) - (sizeof(float *)))


/* image structure definition */

struct xvimage	{

/*  Administrative or file management information */

	char	identifier;	/* a magic number that tells
				   the world that this is an
				   Khoros file */

	char	file_type;	/* tells if this file is a VIFF file */

	char 	release;	/* release number */

	char	version;	/* version number */

        char	machine_dep;	/* indicates peculiarities of */
				/* machine architecture */
	
	char    trash[3];	/* preserves word boundaries */
                                /* groups of 4 bytes */

	char	comment[512];	/* text for image commentary */

/* Things that specify the spatial properties of the image, pixel
   organization, and data storage arrangement.  */

	unsigned long	row_size;	/* length of row in pixels,
					   i.e. number of columns */

	unsigned long	col_size;	/* length of column in pixels,
				   	   i.e. number or rows */
					
	unsigned long	subrow_size;    /* Length of subrows. This is useful
                                           when one wants pixel vectors to
                                           represent 2D objects (images).
                                           The size of each pixel "image"
                                           would be subrow_size (columns)
                                           by num_data_bands/subrow_size (rows).
                                           This field may be ignored except
                                           by routines that need the 2D
                                           interpretation. */

/* The product of row_size and col_size is used to indicate 
   the number of locations when the location type is explicit,
   the product also indicates the number of pixels in a band,
   or the number of vectors.  */

	long	startx, starty;	/* subimage starting position (upper
				   left hand corner), negative indicates
				   that it is not a subimage */

	float 	pixsizx, pixsizy;	/* Actual size of pixel at time of
					   digitization in meters */

	unsigned long	location_type;	/* implied or explicit location
					   data	(implied locations are
					   derived from row_size and
					   col_size */

	unsigned long	location_dim;	/* explicit locations can be of
					   any dimension  */

	unsigned long	num_of_images;     /* number of images
					   pointed to by *imagedata,
					   do not confuse with number of
					   bands */

        unsigned long	num_data_bands;	/* Number of bands per data pixel,
					   or number of bands per image, or
					   dimension of vector data, or
					   number of elements in a vector */

	unsigned long	data_storage_type;  /* storage type for disk data */

	unsigned long	data_encode_scheme; /* encoding scheme of disk data */

/* Things that determine how the mapping (if any) of data bands is 
   to be done to obtain the actual "image" or data.  */

        unsigned long	map_scheme;	 /* How mapping (if any) is to occur */

	unsigned long	map_storage_type;/* Storage type of cells in the maps */

	unsigned long	map_row_size; /* number of columns in map array */

	unsigned long	map_col_size;	/* number of entries in map (rows) */

	unsigned long	map_subrow_size; /* Length of subrows. This is useful
				            when using the output vector from
                                            the map as a 2-D image, rather
					    than just a vector. The size of
					    the 2-D image would be:
					    map_subrow_size (columns) by
					    map_row_size/map_subrow_size
                                            (rows). This field may be ignored
                                            except by routines that need the 2D
                                            interpretation */

	unsigned long	map_enable;	/* Tells if the disk data is valid
					   with or without being sent thru the
					   map. Some data MUST be mapped to be
					   valid. */

	unsigned long	maps_per_cycle;	/* number of maps to constitue a "cycle"
 					   for VFF_MS_CYCLE */

/* Specification of the particular color model in use when working with a
   color image. This just tells what the coordinate system and axis orientation
   of the color space is.  */

        unsigned long	color_space_model;

/* Extra fields for use by the user as needed. These are NOT SUPPORTED
   in any way, except for being read and written correctly with respect
   to machine dependencies.  */

	unsigned long 	ispare1,ispare2;	/* Spare long ints */

	float		fspare1,fspare2;	/* Spare floats */

/* Pointers to the actual data - these are valid only when in memory! */

	char	reserve[VIFF_HEADER_RESERVESIZE];
				/* maximum header information is 
				   1024 bytes, what is not currently
				   used is saved in reserve */

	char *maps;		/* a pointer to the maps, must be cast into 
				   the proper type */

	float *location;	/* a pointer to the location data (for
				   explicit locations, each location is
				   paired with data pointed to by 
				   *imagedata,  all locations are 
				   in float  */

	char *imagedata;	/* a pointer to the input data (straight off
				   of disk, must be cast into the proper type */

	} ;
        
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


/* The following are a couple of defines that make it easier
   to use the image size data */
#define number_of_rows col_size
#define number_of_cols row_size


/********************************************************/
/*                                                      */
/*              read/write/create xvimages              */
/*                                                      */
/********************************************************/

#ifdef __cplusplus

extern "C" {

typedef xvimage ViffImage;

/** \defgroup VIFFFunctions Read/write/create/delete VIFF images 
    See the <a href="documents/viff.ps">VIFF format decription</a> for details.
   @{
*/
    /** \brief Create ViffImage and read contents from given file.
    
	\#include: "<a href="viff_8h-source.html">vigra/viff.h</a>"
	
    */
ViffImage * readViffImage(char * filename);


    /** \brief Delete ViffImage and free memory.
    
	\#include: "<a href="viff_8h-source.html">vigra/viff.h</a>"
	
    */
void freeViffImage(ViffImage * img);

    /** Write ViffImage to file.
    
        <b> Declaration: </b>
	\code
	void writeViffImage(char * filename, ViffImage * img);
	\endcode

	\#include: "<a href="viff_8h-source.html">vigra/viff.h</a>"
	
    */
int writeViffImage(char * filename, ViffImage * img);

    /** \brief Create ViffImage with given characteristics. 
    
    This is normally 
    only called within the function \ref createViffImage. See
    the <a href="documents/viff.ps">VIFF format decription</a> for more details.
    
            
	\#include: "<a href="viff_8h-source.html">vigra/viff.h</a>"
	
    */
ViffImage * createimage(int height, int width, int type,
                             int number_of_images, int number_of_bands,
			     char * comment,
			     int map_row_size, int map_col_size, 
			     int map_scheme, int map_storage_type,
			     int location_type, int location_dim);
		     
    /** \brief Create ViffImage with identical characteristics.
    
        (Data are not copied.) See
        the <a href="documents/viff.ps">VIFF format decription</a> for more details.
    
	\#include: "<a href="viff_8h-source.html">vigra/viff.h</a>"
	
    */
ViffImage * createsameimage(ViffImage * img);

    /** \brief Create single-band xvimage with given size and type. 
    
       This is normally 
       only called within the function \ref createViffImage. See
       the <a href="documents/viff.ps">VIFF format decription</a> for more details.
    
	\#include: "<a href="viff_8h-source.html">vigra/viff.h</a>"
	
    */
ViffImage * createsimpleimage(int height, int width, int type);

    /** \brief Create multi-band (multi-spectral) ViffImage with given 
        characteristics. 
        
        This is normally 
        only called within the function \ref createViffImage. See
        the <a href="documents/viff.ps">VIFF format decription</a> for more details.
    
        'type' may be one of the following:
        <DL>
        <DT>VFF_TYP_1_BYTE:<DD> pixels are byte (unsigned char) 
        <DT>VFF_TYP_2_BYTE:<DD> pixels are two byte (short int) 
        <DT>VFF_TYP_4_BYTE:<DD> pixels are four byte (integer) 
        <DT>VFF_TYP_FLOAT:<DD> pixels are float (single precision)
        <DT>VFF_TYP_COMPLEX:<DD> pixels are complex float 
        <DT>VFF_TYP_DOUBLE:<DD> pixels are float (double precision)
        </DL>
        
	\#include: "<a href="viff_8h-source.html">vigra/viff.h</a>"
	
    */
ViffImage * 
createmultibandimage(int height, int width, int type, int number_of_bands);

/**@}*/

} /* extern "C" */

#else /* ! __cplusplus */

struct xvimage * readViffImage(/* char * filename*/);

void freeViffImage(/* struct xvimage * img */);

int writeViffImage(/* char * filename, struct xvimage * img */);

struct xvimage * createimage(/* int height, int width, int type,
                                 int number_of_images, int number_of_bands,
			         char * comment,
			         int map_row_size, int map_col_size, 
			         int map_scheme, int map_storage_type,
			         int location_type, int location_dim */);
		     
struct xvimage * createsameimage(/* xvimage * img */);

struct xvimage * createsimpleimage(/* int height, int width, int type */);

struct xvimage * 
createmultibandimage(/* int height, int width, int type, int number_of_bands */);

#endif /* __cplusplus */

#endif /* _viff_h_ */
/* Don't add after the endif */
