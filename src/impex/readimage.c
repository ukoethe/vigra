/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/


/*
 *
 *   READIMAGE - Read an KHOROS image from the specified file. Special
 *              filenames are "-" which is STDIN and "#" which is
 *              STDERR.
 *
 *              Machine dependencies are
 *              taken care of by converting the incoming data into
 *              the data format of the host.
 *
 *  NOTE:       THIS CODE IS MACHINE DEPENDENT!
 *
 *             
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifdef _MSC_VER
#    include <io.h>
#else
#    include <sys/uio.h>
#endif

#include "vigra/viff.h"
#include "vdefines.h"

struct xvimage *readViffImage(filename)
char *filename;
{
    char           buf[512];
    long           machtype(),getmachorder();
    long           src_order,dest_order;
    unsigned long  src_machine,dest_machine,dimension,getmachsize();
    int            imgsize,mapsize,locsize,imgcnt,mapcnt,loccnt,storage_type;
    struct xvimage *image, *readheader();
    int            block_read();
    int            read_compressed(),elem_size;
    int            file;
    int            cast_data(),i;
	int imagesize();
	void freeimage();

    /* 
    ** if file cannot be found, return NULL indicating failure.
    */    
    if (!vfullpath(filename,NULL,buf)) return(NULL);

#ifdef _MSC_VER
     if((file = open(buf,O_RDONLY | O_BINARY)) < 0){
        (void)fprintf(stderr,"\nreadimage: unable to access file %s\n",
                      filename);
        return(NULL);
     }
#else
     if((file = open(buf,O_RDONLY)) < 0){
        (void)fprintf(stderr,"\nreadimage: unable to access file %s\n",
                      filename);
        return(NULL);
     }
#endif
    
    /* 
     * Read the image header 
     */
    if((image = readheader(file)) == NULL){
       close(file);
       return(NULL);
    }

    /* 
     * Get size of image components in bytes 
     */
    if(!imagesize(image,&imgsize,&imgcnt,&mapsize,&mapcnt,&locsize,&loccnt)){
       free((void *)image);
       close(file);
       return(NULL);
    }

/*
 * this is where things start to get hairy.
 */

    src_machine = (unsigned long)image->machine_dep;
    dest_machine = machtype(NULL);
    
    src_order = getmachorder((long int)src_machine);
    dest_order = getmachorder((long int)dest_machine);

    /* 
     * Read  maps 
     */
    if(mapsize!=0){
       if((image->maps =(char*)malloc((unsigned)mapsize*sizeof(char)))==NULL){
          fprintf(stderr,"\nreadimage:  Memory allocation failure\n");
          freeimage(image);
          close(file);
          return(NULL);
       }

       if(image->data_encode_scheme == VFF_DES_RAW) {
           if (block_read(file,(char *)image->maps, mapsize) != mapsize) {
               fprintf(stderr,"\nreadimage: ");
    	       fprintf(stderr,"Incorrect byte count reading maps!\n");
               freeimage(image);
               close(file);
               return(NULL);
           }
       }else{
          if(read_compressed(file,image->maps,mapsize,image) != mapsize) {
              fprintf(stderr,"\nreadimage: ");
              fprintf(stderr,"Unable to interpret compressed maps!\n");
              freeimage(image);
              close(file);
              return(NULL);
          }
       } /* fi */

       /*
        * this obnoxious if statement handles the fact that double has
        * a different value for maps and data in vigra/viff.h!
        */
       if (image->map_storage_type == VFF_MAPTYP_DOUBLE)
          storage_type = VFF_TYP_DOUBLE;
       else
          storage_type = image->map_storage_type;

       elem_size = getmachsize(src_machine,(unsigned long)storage_type);
       dimension = mapsize/elem_size;
       if(storage_type!=VFF_MAPTYP_NONE){
          i = cast_data((unsigned char **)&(image->maps),
                         (unsigned int)dimension,
                         (unsigned int)storage_type,
                         (unsigned int)src_machine,
                         (unsigned int)dest_machine);
          if(!i){
            fprintf(stderr,"\nData format conversion failed in readimage().\n");
            fprintf(stderr,"Error occured while reading map data.\n");
            return(0);
          }
       } /* end if storage_type */
    } /* end if maptype!=0 */

    /* 
     * Read  locations 
     */
    if(locsize!=0){
       if((image->location = (float *)malloc(locsize*sizeof(char)))==NULL){
           fprintf(stderr,"\nreadimage:  Memory allocation failure\n");
           freeimage(image);
           close(file);
           return(NULL);
       }

       if(image->data_encode_scheme == VFF_DES_RAW) {
           if (block_read(file,(char *)(image->location), locsize) != locsize) {
               fprintf(stderr,"\nreadimage: ");
    	       fprintf(stderr,"Incorrect byte count reading locations!\n");
               freeimage(image);
               close(file);
               return(NULL);
           }
       } else {
           if (read_compressed(file,(char *)image->location,locsize,image) 
                != locsize) {
               fprintf(stderr,"\nreadimage: ");
    	       fprintf(stderr,"Unable to interpret compressed locations!\n");
               freeimage(image);
               close(file);
               return(NULL);
            }
       } /* fi */

       elem_size = getmachsize(src_machine,(unsigned long)VFF_TYP_FLOAT);
       dimension = locsize/elem_size;
       i = cast_data((unsigned char **)&(image->location),
                     (unsigned int)dimension, 
                     (unsigned int)VFF_TYP_FLOAT, 
                     (unsigned int)src_machine, 
                     (unsigned int)dest_machine);
       if(!i){
          fprintf(stderr,"\nData format conversion failed in readimage().\n");
          fprintf(stderr,"Error occured while reading location data.\n");
          return(0);
       }

    } /* end if locsize!=0 */

    /* 
     * Read  image data 
     */
    if(imgsize!=0){
       if((image->imagedata = (char *)malloc((unsigned)imgsize*sizeof(char)))
          ==NULL){
          fprintf(stderr,"\nreadimage:  Memory allocation failure\n");
          freeimage(image);
          close(file);
          return(NULL);
       }

       if(image->data_encode_scheme == VFF_DES_RAW) {
           if (block_read(file,image->imagedata, imgsize) != imgsize) {
               fprintf(stderr,"\nreadimage: ");
       	       fprintf(stderr,"Incorrect byte count reading image data!\n");
               freeimage(image);
               close(file);
               return(NULL);
           }
       } else {
          if (read_compressed(file,image->imagedata,imgsize,image) != imgsize) {
             fprintf(stderr,"\nreadimage: ");
             fprintf(stderr,"Unable to interpret compressed image data!\n");
             freeimage(image);
             close(file);
             return(NULL);
           }
       } /* fi */

       if(image->data_storage_type!=VFF_TYP_BIT){
          elem_size = getmachsize(src_machine,image->data_storage_type);
          dimension = imgsize/elem_size;
          i = cast_data((unsigned char **) &(image->imagedata),
                            (unsigned int)dimension,
                            (unsigned int)image->data_storage_type, 
                            (unsigned int)src_machine, 
                            (unsigned int)dest_machine); 
          if(!i){
             fprintf(stderr,
                     "\nData format conversion failed in readimage().\n");
             fprintf(stderr,"Error occurred while reading image data.\n");
             return(0);
          }

       }
    } /* end if imgsize!=0 */

  close(file);

  image->machine_dep = (char)machtype(NULL);

  /* 
   * Return the whole mess to the caller 
   */
  return(image);
}

int read_compressed()
{
    return 0;
}


