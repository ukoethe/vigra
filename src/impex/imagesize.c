/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/


/*

    IMAGESIZE - Compute the size in bytes and objects of the
                components of an KHOROS 3 image.

    Usage:
       i = imagesize(&image,&dsize,&dcount,&msize,&mcount,&lsize,&lcount);
       if (i == 0) printf("Invalid image configuration\n");
       else {
         printf("Number of data "pixels" of specified data type: %d\n,dcount);
         printf("Number of map cells of specified map data type: %d\n,mcount);
         printf("Number of location objects (floats): %d\n,lcount);
         printf("Number of data bytes: %d\n,dsize);
         printf("Number of map bytes: %d\n,msize);
         printf("Number of location bytes: %d\n,lsize);
       }
*/

#include "vigra/viff.h"	
#include <stdio.h>

int
imagesize(image,dsize,dcount,msize,mcount,lsize,lcount)
struct xvimage *image;
int	*dsize,*dcount,*msize,*mcount,*lsize,*lcount;
{
    long rows,cols,getmachorder();
    unsigned long getmachsize(),mach;
    int datasize,datacount;
    int mapsize,mapcount;
    int locsize,loccount;

    cols = image->row_size;
    rows = image->col_size;

    mach = (long)image->machine_dep;

    /* 
    ** Compute total size of DATA in bytes 
    */
    if(image->data_storage_type==VFF_TYP_BIT){
       datasize = ((cols+7)/8)*rows;
       datacount = datasize;
    }else{
       datasize = cols*rows * getmachsize(mach,
                  (unsigned long)image->data_storage_type);
       datacount = cols*rows;
    }

    datasize *= image->num_of_images*image->num_data_bands;
    datacount *= image->num_of_images*image->num_data_bands;
 
    /* 
    ** Compute number of MAP data objects 
    */
    switch(image->map_scheme)
      {
        case VFF_MS_NONE:
          mapcount = 0;
          break;
        case VFF_MS_ONEPERBAND:
        case VFF_MS_CYCLE:
          mapcount = 
              image->num_data_bands*image->map_row_size*image->map_col_size;
          break;
        case VFF_MS_SHARED:
        case VFF_MS_GROUP:
          mapcount = image->map_row_size*image->map_col_size;
          break;
        default:
          (void)fprintf(stderr,"\nimagesize: Unknown mapping scheme:");
          (void)fprintf(stderr," %d\n",image->map_scheme);
          return(0);
          /* break; */
      }

    /* 
    ** mapcount now contains the number of CELLS, so convert to bytes 
    */ 
    if(image->map_storage_type==VFF_MAPTYP_NONE){
       mapsize = 0;
    }else{
       mapsize = mapcount*getmachsize(mach,
                 (unsigned long)image->map_storage_type);
    }

    /* 
    ** Compute size of LOCATION data in bytes and floats
    */
    loccount = rows*cols*image->location_dim;
    locsize  = loccount*getmachsize(mach,(long)VFF_TYP_FLOAT);

    *dsize = datasize;
    *dcount = datacount;
    *msize = mapsize;
    *mcount = mapcount;
    *lsize = locsize;
    *lcount = loccount;
    return(1);
}

