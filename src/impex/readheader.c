/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/


/*
 *
 *  READHEADER- Read an KHOROS image header from the file pointed 
 *              to by the supplied file descriptor.
 *
 *  NOTE:       THIS CODE IS MACHINE DEPENDENT!
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "vigra/viff.h"
#include "vdefines.h"

#define BYTE_FIELD_SIZE 520
#define OTHER_FIELD_CNT 25 
#define VIFF_HEADER_SIZE 1024

struct xvimage
*readheader(file)
int file;
{
    struct xvimage *imageptr;
    unsigned char *temp,*temp2,*temp3,*temp4;
    int convert_order(), convert_long(), convert_float();
    long machtype(),current_machine,source_machine;
    long current_order,source_order,word_size,getmachorder();
    unsigned char id1,id2;
    int numread,headersize;  /* Added by Ron. E. Neher to fix pipe problem */
    int i;
	int block_read();

    /* 
     * Grab space for the image structure and read it 
     */

    if((imageptr = (struct xvimage *)calloc(1,sizeof(struct xvimage)))==NULL){
       (void)fprintf(stderr,"\nreadheader: Not enough memory ");
       (void)fprintf(stderr,"for image header\n");
       return(NULL);
    }

    /*
     * malloc some temp space that we'll be needing later
     */
    headersize = VIFF_HEADER_SIZE;
    if((temp = (unsigned char *)calloc(1,(unsigned)headersize))==NULL){
       (void)fprintf(stderr,"\nreadheader: Not enough memory to read ");
       (void)fprintf(stderr,"image header\n");
       return(NULL);
    }

    if((temp2 = (unsigned char *)calloc(1,5*sizeof(long)))==NULL){
       (void)fprintf(stderr,"\nreadheader: Not enough memory to read ");
       (void)fprintf(stderr,"image header\n");
       return(NULL);
    }
    if((temp3 = (unsigned char *)calloc(1,2*sizeof(float)))==NULL){
       (void)fprintf(stderr,"\nreadheader: Not enough memory to read ");
       (void)fprintf(stderr,"image header\n");
       return(NULL);
    }
    if((temp4 = (unsigned char *)calloc(1,18*sizeof(long)))==NULL){
       (void)fprintf(stderr,"\nreadheader: Not enough memory to read ");
       (void)fprintf(stderr,"image header\n");
       return(NULL);
    }

    /*
     * attempt to read the whole header in.
     */

    if((numread = block_read(file,(char *)temp, headersize))!=headersize){
        (void) fprintf(stderr,"\nreadheader: Incorrect header byte count: ");
	(void) fprintf(stderr,"found %d, should be %d\n",numread,headersize);
        (void) fprintf(stderr,"Possible invalid input VIFF data or an empty ");
        (void) fprintf(stderr,"file.\n");
        free(imageptr);
        return(NULL);
    }

    /*
     * fill in the portions of the header that we never need conversion
     */

    (void)memmove((char *)imageptr,(char *)temp,BYTE_FIELD_SIZE);

    /* 
     * check to see if we are an VIFF file 
     */

    id1 = imageptr->identifier;
    id2 = XV_FILE_MAGIC_NUM;
    if (id1 != id2 || imageptr->file_type != XV_FILE_TYPE_XVIFF) {
        (void)fprintf(stderr,"\nreadheader: Cannot read non-VIFF file!\n");
        (void)free(imageptr);
        return(NULL);
    }
    
    
    /* 
     * Check the release and version numbers 
     */

    if (imageptr->release != XV_IMAGE_REL_NUM ||
        imageptr->version != XV_IMAGE_VER_NUM)
    {
        (void) fprintf(stderr,"\nreadheader: Release or version number ");
	(void) fprintf(stderr,"mismatch!\n");
        (void) fprintf(stderr,"readheader: Release is %d, should be %d\n",
                imageptr->release,XV_IMAGE_REL_NUM);
        (void) fprintf(stderr,"readheader: Version is %d, should be %d\n",
                imageptr->version,XV_IMAGE_VER_NUM);
        free(imageptr);
        return(NULL);
    }

    /* 
     * See if we need to convert the header to the current machine type 
     */

    current_machine = machtype(NULL);
    current_order =   getmachorder(current_machine);
    source_machine  = imageptr->machine_dep;
    source_order =    getmachorder(source_machine);

    word_size = 4;

    /*
     * we use IEEE format on crays for the header since their wordsize
     * is different
     */
    if(source_machine==VFF_DEP_CRAYORDER){
       source_machine = VFF_DEP_IEEEORDER;
       source_order = VFF_DEP_BIGENDIAN;
    }

    if(current_machine==source_machine){
        memmove((char *)&imageptr->row_size,
              (char *)(temp + BYTE_FIELD_SIZE),
              (VIFF_HEADER_SIZE - BYTE_FIELD_SIZE)); 
    }else{
        i = convert_order((unsigned char *)(temp+BYTE_FIELD_SIZE),
                      (unsigned int)source_order,
                      (unsigned int)VFF_DEP_BIGENDIAN,
                      (unsigned int)OTHER_FIELD_CNT,
                      (unsigned int)word_size);
        if(!i){
           fprintf(stderr,"\nOrder conversion failed in readheader.\n");
           free(temp);
           free(temp2);
           free(temp3);
           free(temp4);
           return(0);
        }

        /*
         * this call converts row_size through starty to new encoding
         */

        memmove((char *)temp2,(char *)(temp+BYTE_FIELD_SIZE),20);
        i = convert_long((unsigned char **)&temp2,
                             (unsigned int)source_machine,
                             (unsigned int)current_machine,
                             (unsigned int)5);
        if(!i){
           fprintf(stderr,"\nData format conversion failed in readheader.\n");
           free(temp);
           free(temp2);
           free(temp3);
           free(temp4);
           return(0);
        }

        /* 
         * this call converts pixsizx and pixsizy to new encoding 
         */

        memmove((char *)temp3,(char *)(temp+BYTE_FIELD_SIZE+20),8);
        i = convert_float((unsigned char **)&temp3,
                             (unsigned int)source_machine,
                             (unsigned int)current_machine,
                             (unsigned int)2);
        if(!i){
           fprintf(stderr,"\nData format conversion failed in readheader.\n");
           free(temp);
           free(temp2);
           free(temp3);
           free(temp4);
           return(0);
        }
 
        /* 
         * this call converts location_type through fspare2 to new encoding
         */

        memmove((char *)temp4,(char *)(temp+BYTE_FIELD_SIZE+28),72);
        i = convert_long((unsigned char **)&temp4,
                             (unsigned int)source_machine,
                             (unsigned int)current_machine,
                             (unsigned int)18);
        if(!i){
           fprintf(stderr,"\nData format conversion failed in readheader.\n");
           free(temp);
           free(temp2);
           free(temp3);
           free(temp4);
           return(0);
        }


        memmove((char *)&imageptr->row_size,(char *)temp2,5*sizeof(long));
        memmove((char *)&imageptr->pixsizx,(char *)temp3,2*sizeof(float));
        memmove((char *)&imageptr->location_type,(char *)temp4,18*sizeof(long));

        i = convert_order((unsigned char *)(&imageptr->row_size),
                      (unsigned int)VFF_DEP_BIGENDIAN,
                      (unsigned int)current_order,
                      (unsigned int)OTHER_FIELD_CNT,
                      (unsigned int)word_size);
        if(!i){
           fprintf(stderr,"\nData format conversion failed in readheader.\n");
           free(temp);
           free(temp2);
           free(temp3);
           free(temp4);
           return(0);
        }

        /*
         * this last call converts the resulting header to the
         * correct byte order
         */

        if(current_machine==VFF_DEP_CRAYORDER)
           word_size = 8;
        else
           word_size = 4;


    } /* fi big if */

    free(temp);
    free(temp2);
    free(temp3);
    free(temp4);

    /* 
     * Return the whole mess to the caller 
     */
    return(imageptr);
}
