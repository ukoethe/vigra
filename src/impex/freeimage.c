/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "vigra/viff.h"	
#include "vdefines.h"	


/************************************************************
*
*  MODULE NAME: freeimage
*
*      PURPOSE: This routine frees an khoros xvimage structure.
*
*        INPUT: image --  a pointer to an khoros xvimage structure that 
*			  contains the image structure to be freed.
*
*       OUTPUT: (none)  since all we are doing is freeing as much of
*		a image structure as we can.
*
*    CALLED BY: any routine that wishes to free an xvimage structure
*
*
*************************************************************/

void freeViffImage(image)
struct xvimage *image;
{
	unsigned char id1,id2;

	/*
	 *  Free as much of the xvimage structure as we can.  But first check to
	 *  make sure the image pointer is not NULL.
	 */
	if (image != NULL)
	{
	   /* Now see of the image itself is legal. This catches accidental
	      attempts to free an image already turned loose by freeimage(). */
           id1 = image->identifier;
	   id2 = XV_FILE_MAGIC_NUM;
	   if (id1 != id2)
	     {
	       (void)fprintf(stderr,
       "freeimage: Attempt to free an object that is not a VIFF image.\n");
	       (void)fprintf(stderr,
       "freeimage: Object may be a VIFF image that has already been free'd.\n");
	       (void)fprintf(stderr,"freeimage: Attempt aborted.\n");
	       return;
	     }

	   /*  free image data */
	   if (image->imagedata != NULL && (image->row_size * 
	       image->col_size) > 0)
	      free ((char *) image->imagedata);

	   /*  free map data */
	   if (image->maps != NULL && image->map_row_size > 0)
	      free ((char *) image->maps);

	   /*  free location data */
	   if (image->location != NULL && image->row_size *
               image->col_size > 0 &&
               image->location_type == VFF_LOC_EXPLICIT)
	      free ((char *) image->location);

	   /*
	    *  Get rid of the image structure itself.  We know it is not
	    *  NULL, since we checked it above. BUT: before we do this,
            *  fill the header with zeros so that should an unspecting
            *  programmer try to use an already free'd image it will not
            *  just appear that all is well!
	    */
           memset(image, 0, sizeof(struct xvimage));
	   free((char *) image);

	}
}

void freeimage(image)
struct xvimage *image;
{
    freeViffImage(image);
}

