/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

/*INCLUDE*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vigra/viff.h"
#include "vdefines.h"

#define LENGTH 512

/**************************************************************
*
* MODULE NAME: createimage
*
*     PURPOSE: Create a generic image
*
*       INPUT: 	col_size -- the size of a column
*		row_size -- the size of a row
*		data_storage_type -- the VFF_TYP_* define of image
*		num_of_images -- # of images pointed to by imagedata
*		num_data_bands -- # of bands/pixel, /image or dim vec data
*		comment -- description of image
*		map_row_size -- # of columns in map array
*		map_col_size -- # of rows in map array
*		map_storage_type -- Storage type of cells in the maps
*		location_type -- implied or explicit location data
*		location_dim -- explicit locations can be of any dimension
*
*      OUTPUT: 	1.  returns a pointer to an xvimage with image defined
*
* CALLED FROM: 
*
* ROUTINES CALLED: 
*
**************************************************************/
struct xvimage *
createimage(col_size, row_size, data_storage_type, num_of_images,
            num_data_bands, comment, map_row_size, map_col_size,
	    map_scheme, map_storage_type, location_type, location_dim)
unsigned
long	col_size,
	row_size,
	data_storage_type,
	num_of_images,
	num_data_bands,
	map_row_size,
	map_col_size,
	map_scheme,
	map_storage_type,
	location_type,
	location_dim;
char	*comment;
{
struct
xvimage *image;
char	*maps, 
	*imagedata, 
	tmp_comment[LENGTH];
long 	machtype();
float	*location;
int	cstrlen,
	image_data_size_bytes,		/* # data bytes */
	image_data_count_pixels,	/* # data pixels */
	map_size_bytes,			/* # map bytes */
	map_count_cells,		/* # map cells */
	location_size_bytes,		/* # location bytes */
	location_count_objects;		/* # location objs */
int imagesize();

/* malloc room for the xvimage structure */

    if ((image=(struct xvimage *)calloc(1,sizeof(struct xvimage)))==NULL)
    {
        (void) fprintf(stderr,"createimage: No space for image \
- malloc failed!\n");
        return(0);
    }

/* setup the comment (can only be 511 chars) */

    cstrlen = VStrlen(comment);
    if (cstrlen > 0)
    {
       if (cstrlen < 512)
          (void) strcpy(tmp_comment, comment);
       else
       {
          (void) strncpy(tmp_comment, comment, LENGTH - 1);
          (void) strcat(tmp_comment, "");
       }
    }
    else
       (void) strcpy(tmp_comment, "");

/* Load the image header with the values. These can be over-ridden by
   giving them a different value after returning to the calling routine.
 */
    image->identifier = (char)XV_FILE_MAGIC_NUM;
    image->file_type = XV_FILE_TYPE_XVIFF;
    image->release = XV_IMAGE_REL_NUM;
    image->version = XV_IMAGE_VER_NUM;
    image->machine_dep = (char)machtype(NULL);
    (void) strcpy(image->comment, tmp_comment);
    image->row_size = row_size;
    image->col_size = col_size;
    image->startx = VFF_NOTSUB;
    image->starty = VFF_NOTSUB;
    image->pixsizx = 1.0;
    image->pixsizy = 1.0;
    image->location_type = location_type;
    image->location_dim = location_dim;
    image->num_of_images = num_of_images;
    image->num_data_bands = num_data_bands;
    image->data_storage_type = data_storage_type;
    image->data_encode_scheme = VFF_DES_RAW;
    image->map_scheme = map_scheme;
    image->map_storage_type = map_storage_type;
    image->map_row_size = map_row_size;
    image->map_col_size = map_col_size;
    image->map_subrow_size = 0;
    image->map_enable = VFF_MAP_OPTIONAL;
    image->maps_per_cycle = 0;      /* Don't care */
    image->color_space_model = VFF_CM_NONE;
    image->ispare1 = 0;
    image->ispare2 = 0;
    image->fspare1 = 0;
    image->fspare2 = 0;

/* get the sizes for the image data, map data, and location data */

    if (! imagesize(image, 			/* xvimage */
		    &image_data_size_bytes,	/* # data bytes */
		    &image_data_count_pixels,	/* # data pixels */
		    &map_size_bytes,		/* # map bytes */
		    &map_count_cells,		/* # map cells */
		    &location_size_bytes,	/* # location bytes */
		    &location_count_objects	/* # location objs */
		   ))
    {
	fprintf(stderr, "createimage: Uninterpretable image \
specificationa\n");
	return(0);
    }

/* malloc room for the image data */

    if (image_data_size_bytes > 0)
    {
       if ((imagedata = (char *)
	  malloc((unsigned int)image_data_size_bytes)) == NULL)
       {
	  (void) fprintf(stderr,"createimage: Not enough memory for image\
 data!\n");
	  return(0);
       }
    }
    else
    {
       imagedata = NULL;
    }

/* malloc room for the color map data */

    if (map_size_bytes > 0)
    {
       if ((maps = (char *)malloc((unsigned int)map_size_bytes)) == NULL)
       {
	   (void) fprintf(stderr,"createimage: Not enough memory for maps\
 data!\n");
	   return(0);
       }
    }
    else
    {
       maps = NULL;
    }


/* malloc room for the location */

    if (location_size_bytes)
    {
       if ((location = (float *)
	   malloc((unsigned int)location_size_bytes))==NULL)
       {
           (void) fprintf(stderr,"createimage: Not enough memory \
 for location data!\n");
           return(0);
       }
    }
    else
    {
       location = NULL;
    }


/* Load the image data, color map data, and location data */

    image->maps = maps;
    image->location = location;
    image->imagedata = imagedata;

    return(image);
}

/**************************************************************
*
* MODULE NAME: createsameimage
*
*     PURPOSE: Create a generic image using the header information 
*              of a given image
*
*       INPUT: 	pointer to an xvimage struct
*
*      OUTPUT: 	1.  returns a pointer to an xvimage with image defined
*
* CALLED FROM: 
*
* ROUTINES CALLED: 
*
**************************************************************/
struct xvimage *
createsameimage(template)
struct xvimage * template;
{
    return createimage(template->col_size, template->row_size,
                       template->data_storage_type,
                       template->num_of_images,template->num_data_bands,
                       NULL,    /* no comment */
                       template->map_row_size, template->map_col_size,
	                   template->map_scheme, template->map_storage_type,
	                   template->location_type, template->location_dim);
}

/**************************************************************
*
* MODULE NAME: createsimpleimage
*
*     PURPOSE: Create a plain single band image of specified type and size
*
*       INPUT: 	height, width, and type
*
*      OUTPUT: 	1.  returns a pointer to an xvimage with image defined
*
* CALLED FROM: 
*
* ROUTINES CALLED: 
*
**************************************************************/
struct xvimage *
createsimpleimage(height, width, type)
int height, width, type;
{
    return createimage(height, width, type,
                       1, 1,
                       NULL,    /* no comment */
                       0, 0,
	               VFF_MS_NONE, VFF_MAPTYP_NONE,
	               VFF_LOC_IMPLICIT, 0);
}

/**************************************************************
*
* MODULE NAME: createmultibandimage
*
*     PURPOSE: Create a plain multiband image of specified type and size
*
*       INPUT: 	height, width, type, and number of bands
*
*      OUTPUT: 	1.  returns a pointer to an xvimage with image defined
*
* CALLED FROM: 
*
* ROUTINES CALLED: 
*
**************************************************************/
struct xvimage *
createmultibandimage(height, width, type, bands)
int height, width, type, bands;
{
    return createimage(height, width, type,
                       1, bands,
                       NULL,    /* no comment */
                       0, 0,
	               VFF_MS_NONE, VFF_MAPTYP_NONE,
	               VFF_LOC_IMPLICIT, 0);
}




