/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

/*

    WRITEIMAGE - Write an KHOROS image to the file pointed to by the
                 supplied file descriptor.

    Returns:  1 if successful
	      0 if unsuccessful

*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifdef _MSC_VER
#    include <io.h>
#else
#    include <unistd.h>
#    include <sys/uio.h>
#endif

#include <string.h>

char * vfullpath();

#include "vigra/viff.h"
#include "vencode.h"

#define VWRITE_SIZE 65536		/* Write 64KB at a time */
#define BYTE_FIELD_SIZE 520
#define OTHER_FIELD_CNT 25
#define VIFF_HEADER_SIZE 1024

int writeViffImage(filename,imageptr)
struct xvimage *imageptr;
char *filename;
  {
    char buf[512];
    int imgsize,locsize,mapsize,imgcnt,loccnt,mapcnt,index,n;
    char *temp1,*temp3;
    float *temp2;
    unsigned char *temp,*ptemp,*ptemp2,*ptemp3,*ptemp4;
    int headersize,i;
    int file;
	int imagesize();
	int convert_long(), convert_float();
	int write_compressed();

    /* 
     * Check the image pointer 
     */
    if(imageptr == NULL){
       fprintf(stderr,"writeimage: NULL image pointer, no image written!\n");
       return(0);
    }

    /* 
     * Open the output file 
     */
    if(!vfullpath(filename,NULL,buf))return(0);

#ifdef _MSC_VER
    if((file = open(buf,O_WRONLY|O_TRUNC|O_CREAT|O_BINARY,0664)) < 0){
       fprintf(stderr,"writeimage: unable to access %s for writing\n",
                           filename);
       return(0);
    }
#else
    if((file = open(buf,O_WRONLY|O_TRUNC|O_CREAT,0664)) < 0){
       fprintf(stderr,"writeimage: unable to access %s for writing\n",
                           filename);
       return(0);
    }
#endif

    if(!imagesize(imageptr,&imgsize,&imgcnt,&mapsize,&mapcnt,&locsize, 
       &loccnt)) {
       fprintf(stderr,"writeimage: Uninterpretable image specification\n");
       close(file);
       return(0);
    }

    /* 
     * Sanity check 
     */
    if(imgsize != 0 && imageptr->imagedata == NULL){
       fprintf(stderr,
          "writeimage: Bogus image - no image data but nonzero image size!\n");
       fprintf(stderr,"writeimage: Image not written!\n");
       close(file);
       return(0);
    }

    if(mapsize != 0 && imageptr->maps == NULL){
       fprintf(stderr,
           "writeimage: Bogus image - no map data but nonzero map size!\n");
       fprintf(stderr,"writeimage: Image not written!\n");
       close(file);
       return(0);
    }

    if(locsize != 0 && imageptr->location == NULL){
       fprintf(stderr,"writeimage: Bogus image - no location ");
       fprintf(stderr,"data but nonzero location size!\n");
       fprintf(stderr,"writeimage: Image not written!\n");
       close(file);
       return(0);
    }

    /* 
     * Write the image header 
     * Note that the pointer fields in the header must be written out
     * as NULLs in order that an image file on the disk can be
     * compared using diff or cmp! 
     */
    temp1 = imageptr->maps;      /* Save pointers to data areas */
    temp2 = imageptr->location;
    temp3 = imageptr->imagedata;

    imageptr->maps = NULL;       /* Set 'em to NULL */
    imageptr->location = NULL;
    imageptr->imagedata = NULL;

    if(imageptr->machine_dep!=VFF_DEP_CRAYORDER){
	   int w = write(file,(char *)imageptr,sizeof(struct xvimage));
       if(w != sizeof(struct xvimage)) {
          fprintf(stderr,"writeimage: Unable to write image header\n");
          close(file);
          return(0);
       }
    }else{
        temp = (unsigned char *)imageptr;

        headersize = VIFF_HEADER_SIZE;
        if((ptemp = (unsigned char *)malloc((unsigned)headersize))==NULL){
           (void)fprintf(stderr,"\nwriteimage: Not enough memory to write ");
           (void)fprintf(stderr,"image header\n");
           return(0);
        }

        if((ptemp2 = (unsigned char *)malloc(5*sizeof(long)))==NULL){
           (void)fprintf(stderr,"\nwriteimage: Not enough memory to write ");
           (void)fprintf(stderr,"image header\n");
           return(0);
        }
        if((ptemp3 = (unsigned char *)malloc(2*sizeof(float)))==NULL){
           (void)fprintf(stderr,"\nwriteimage: Not enough memory to write ");
           (void)fprintf(stderr,"image header\n");
           return(0);
        }
        if((ptemp4 = (unsigned char *)malloc(18*sizeof(long)))==NULL){
           (void)fprintf(stderr,"\nwriteimage: Not enough memory to write ");
           (void)fprintf(stderr,"image header\n");
           return(0);
        }

        memmove((char *)ptemp2,(char *)(temp+BYTE_FIELD_SIZE),40);
        memmove((char *)ptemp3,(char *)(temp+BYTE_FIELD_SIZE+40),16);
        memmove((char *)ptemp4,(char *)(temp+BYTE_FIELD_SIZE+56),144);

        /*
         * first convert long stuff
         */
        i = convert_long((unsigned char **)&ptemp2,
                             (unsigned int)VFF_DEP_CRAYORDER,
                             (unsigned int)VFF_DEP_IEEEORDER,
                             (unsigned int)5);
        if(!i){
           fprintf(stderr,"\nData format conversion failed in readheader.\n");
           free(ptemp);
           free(ptemp2);
           free(ptemp3);
           free(ptemp4);
           return(0);
        }

        /*
         * then convert float stuff
         */
        i = convert_float((unsigned char **)&ptemp3,
                             (unsigned int)VFF_DEP_CRAYORDER,
                             (unsigned int)VFF_DEP_IEEEORDER,
                             (unsigned int)2);
        if(!i){
           fprintf(stderr,"\nData format conversion failed in readheader.\n");
           free(ptemp);
           free(ptemp2);
           free(ptemp3);
           free(ptemp4);
           return(0);
        }

        /*
         * and then convert last long stuff 
         */
        i = convert_long((unsigned char **)&ptemp4,
                             (unsigned int)VFF_DEP_CRAYORDER,
                             (unsigned int)VFF_DEP_IEEEORDER,
                             (unsigned int)18);
        if(!i){
           fprintf(stderr,"\nData format conversion failed in readheader.\n");
           free(ptemp);
           free(ptemp2);
           free(ptemp3);
           free(ptemp4);
           return(0);
        }

        /*
         * now construct a new header and write it
         */
        memmove(ptemp,(char *)imageptr,BYTE_FIELD_SIZE);
        memmove(ptemp+BYTE_FIELD_SIZE, (char *)ptemp2,20);
        memmove(ptemp+BYTE_FIELD_SIZE+20,(char *)ptemp3,8);
        memmove(ptemp+BYTE_FIELD_SIZE+28,(char *)ptemp4,72);

        if(write(file,(char *)ptemp,headersize)!=headersize) {
           fprintf(stderr,"writeimage: Unable to write image header\n");
           close(file);
           return(0);
        }
    } /* end cray special case */

    imageptr->maps = temp1;      /* Restore the previous values */
    imageptr->location = temp2;
    imageptr->imagedata = temp3;

    /* 
     * Write maps 
     */

    index = imageptr->data_encode_scheme;

    /* 
     * Write to normal file descriptor 
     */
    switch((int)imageptr->data_encode_scheme){
       case VFF_DES_RAW:
            n = write(file,imageptr->maps,mapsize);
	    if(n != mapsize){
	       fprintf(stderr,"writeimage: Unable to write maps\n");
               close(file);
               return(0);
            }
            break;
       case VFF_DES_COMPRESS:
       case VFF_DES_RLE:
       case VFF_DES_TRANSFORM:
       case VFF_DES_CCITT:
       case VFF_DES_ADPCM:
       case VFF_DES_GENERIC:
            n = write_compressed(file,(char *)imageptr->maps,mapsize,
                 compr_cmd[index]);
            if (n != mapsize){
                fprintf(stderr,"writeimage: Unable to write maps\n");
                close(file);
                return(0);
            }
	    break;
       default:
            fprintf(stderr,
               "writeimage: Unknown encoding scheme: %d maps not written!\n",
               (int)imageptr->data_encode_scheme);
            close(file);
            return(0);
            /* break; */
    } /* end switch */
                 
    /* 
     * Write the location data 
     */

    switch((int)imageptr->data_encode_scheme){
       case VFF_DES_RAW: 
            n = write(file,(char *)imageptr->location,locsize);
            if(n != locsize){
               fprintf(stderr,"writeimage: Unable to write location data\n");
               close(file);
               return(0);
            }
            break;
       case VFF_DES_COMPRESS:
       case VFF_DES_RLE:
       case VFF_DES_TRANSFORM:
       case VFF_DES_CCITT:
       case VFF_DES_ADPCM:
       case VFF_DES_GENERIC:
            n = write_compressed(file,(char *)imageptr->location,locsize,
                                 compr_cmd[index]);
            if(n != locsize){
               fprintf(stderr,"writeimage: Unable to write location data\n");
               close(file);
               return(0);
            }
	    break;
       default:
            fprintf(stderr,
            "writeimage: Unknown encoding scheme: %d locations not written!\n",
            (int)imageptr->data_encode_scheme);
            close(file);
            return(0);
            /* break; */
    } /* end switch */
                 
    /* 
     * write image data 
     */

    switch((int)imageptr->data_encode_scheme){
       case VFF_DES_RAW:
            n = write(file,imageptr->imagedata,imgsize);
            if(n != imgsize){
               fprintf(stderr,"writeimage: Unable to write image data\n");
               close(file);
               return(0);
            }
            break;
       case VFF_DES_COMPRESS:
       case VFF_DES_RLE:
       case VFF_DES_TRANSFORM:
       case VFF_DES_CCITT:
       case VFF_DES_ADPCM:
       case VFF_DES_GENERIC:
            n = write_compressed(file,imageptr->imagedata,imgsize,                                               compr_cmd[index]);
            if(n != imgsize){
               fprintf(stderr,"writeimage: Unable to write image data\n");
               close(file);
               return(0);
            }
	    break;
       default:
            fprintf(stderr,
              "writeimage: Unknown encoding scheme: %d imagedata not written!\n",
              index);
            close(file);
            return(0);
            /* break; */
    }

    close(file);
    return(1);
}

int write_compressed()
{
    fprintf(stderr,"Khoros::write_compressed(): Sorry, not implemented.\n");
    return -1;
}
