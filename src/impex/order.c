/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/


/*

      ORDER.C - A library of routines for converting values between
                machine storage formats.


      NOTE:    THIS CODE IS MACHINE ARCHITECTURE DEPENDENT!

*/

#include <stdlib.h>
#include <stdio.h>
#include "vigra/viff.h"

/***********************************************************************
*
*  Routine Name: convert_order()
*
*          Date: Sat Sep 21 10:26:06 MDT 1991
*        
*       Purpose: converts a string of MSB or LSB words of any length to
*                another order.  i.e. this can be used to convert MSB
*                words to LSB or LSB to MSB, etc.
*
*         Input: data:           input data cast to unsigned char
*		 src_mach:       source machine type
*		 dest_mach:      destination machine type 
*                num_elem:       number of nbyte elements in the string
*                bytes_per_elem: number of bytes per element.
*
*        Output: data:           reordered data.
*
*    Written By: Mark Young, Jeremy Worley 
*
* Modifications:
*
***********************************************************************/

int convert_order(data,src_mach,dest_mach,num_elem,bytes_per_elem)
   register unsigned char *data;
   unsigned int src_mach,dest_mach,num_elem,bytes_per_elem;
{
   
   register unsigned int i,j,indx1,indx2,num_swaps;
   unsigned char tmp;
   long src_order,dest_order,getmachorder();


   src_order=getmachorder((long int)src_mach);
   dest_order=getmachorder((long int)dest_mach);
   if (src_order != dest_order)
   {
      num_swaps = bytes_per_elem >> 1;
      if(num_swaps==0)return(1);

      for(i=0,indx1=0;i<num_elem;i++){
         indx2 = indx1 + (bytes_per_elem-1);
         for(j=0;j<num_swaps;j++){
	    tmp = data[indx1];
            data[indx1++] = data[indx2];
            data[indx2--] = tmp;
         }
         indx1 += num_swaps;
      }
   } /* end switch */
   return(1);
}
 
/********************************************************************/
/************************ C A S T   N U M B E R *********************/
/********************************************************************/

/***********************************************************************
*
*  Routine Name: cast_data()
*
*          Date: Mon Sep 23 10:44:30 MDT 1991
*        
*       Purpose: casts an array of data from one machine's representation
*                to another's.
*
*                This function assumes that data is in MSB representation.
*
*                This function also assumes that the variable "data" is
*                already malloced.  
*
*         Input: data
*
*        Output: returns a 1 for success and a 0 for failure.
*
*    Written By: Jeremy Worley 
*
* Modifications:
*
***********************************************************************/

int cast_data(data,dimension,type,src_mach,dest_mach)
   unsigned char **data;
   unsigned int dimension,type,src_mach,dest_mach;
{
   int i;
   unsigned int elem_size,ord_size,ord_dim;
   unsigned long  getmachsize();
   int convert_short(), convert_long(), convert_float(), convert_double();

   if(dimension<=0 || src_mach==dest_mach)return(1);

   elem_size = getmachsize(src_mach,(unsigned long)type);

   /*
   ** complex and double complex are treated as floats and doubles
   ** respectively.  This is because convert_order will rearrage 
   ** the real and imaginary components if we pass in the size of
   ** the whole complex structure.
   */
   if(type==VFF_TYP_COMPLEX || type==VFF_TYP_DCOMPLEX){
      ord_dim  = dimension*2;
      ord_size = elem_size/2;
   }else{
      ord_dim  = dimension;
      ord_size = elem_size;
   }

   i = convert_order((unsigned char *)(*data),
                     (unsigned int)src_mach,
                     (unsigned int)VFF_DEP_BIGENDIAN, 
                     (unsigned int)ord_dim,
                     (unsigned int)ord_size);
   if(!i) return(0);

   switch(type){
      /* 
      ** first case is 1_byte.  return same data under all circumstances
      */
      case VFF_TYP_1_BYTE:
          return(1);
      /* 
      ** 
      */
      case VFF_TYP_2_BYTE:
          i = convert_short(data,src_mach,dest_mach,dimension);
          break; 
      /*
      ** 
      */
      case VFF_TYP_4_BYTE:
	   i = convert_long(data,src_mach,dest_mach,dimension);
           break; 
      /*
      ** float: 
      */
      case VFF_TYP_FLOAT:
          i = convert_float(data,src_mach,dest_mach,dimension);
          break;
      /*
      **
      */
      case VFF_TYP_COMPLEX:
          i = convert_float(data,src_mach,dest_mach,2*dimension);
          break;
      /*
      **
      */
      case VFF_TYP_DOUBLE:
          i = convert_double(data,src_mach,dest_mach,dimension);
          break;
      /*
      **
      */
      case VFF_TYP_DCOMPLEX:
          i = convert_double(data,src_mach,dest_mach,2*dimension);
          break;
   } /* end switch */

   if (!i) return(0);

   elem_size = getmachsize(dest_mach,(unsigned long)type);
   i = convert_order((unsigned char *)(*data),
                         (unsigned int)VFF_DEP_BIGENDIAN,
                         (unsigned int)dest_mach, 
                         (unsigned int)ord_dim, 
                         (unsigned int)ord_size);
   return(i);
}

/********************************************************************/
/************************ S H O R T    I N T  ***********************/
/********************************************************************/

/***********************************************************************
*
*  Routine Name: convert_short()
*
*          Date: Fri Oct  4 16:14:20 MDT 1991
*        
*       Purpose: in general, this routine is designed to convert one
*                short integer representation to another.  There
*                are special cases involving 64 bit machines, however. 
*                If this routine gets a src=32bit and a dest=64bit, 
*                then it is expected to "stretch" the data.  In other
*                words, in that case, the input is in 16 bit components
*                and the output is in 32 bit components.
*
*                Huge and unjustifiable assumption:  ALL shorts are signed.
*                (sign extension is performed in all cases)
*
*         Input: data:  address of the input data cast to (unsigned char **)  
*                src:   unsigned int specifying source architecture.
*                dest:  another unsigned int specifying the destination 
*                       (usually current) architecture.
*                num:   unsigned int specifying number of short elements
*                       int *data (NOT the number of bytes in the array)
*
*        Output: This function returns a 1 for success and a 0 for
*                failure.  Failure indicates (1) that the programmer
*                didn't know what the heck he/she was doing; (2) that
*                this function wasn't able to realloc the data if
*                that was necessary; or (3) that the src or dest
*                arguments were invalid. 
*
*    Written By: Jeremy Worley 
*
* Modifications:
*
***********************************************************************/

int convert_short(data,src,dest,num)
   unsigned char **data;
   unsigned int src,dest,num;
{
   unsigned long src_size,dest_size,i,j,idxd,idxs;
   unsigned long getmachsize();
   unsigned char *temp,*iput;

/*
** test to see if the goofy programmer called this routine with
** src == dest 
*/

  if(src==dest)return(1);

/*
** go ahead and return saying that everything is cool if the 
** calling function passed in num=0.  
*/

  if(num==0)return(1);

/*
** test to see if src or dest are valid
*/

  if(src>=255 || dest>=255)return(0);   

/*
** test to see if the calling routine "knew" to malloc data before calling
** this one.
*/

   if(data == NULL)return(0);
   if(*data == NULL)return(0);

/*
** because we deal with byte ordering as a separate issue, the only
** real cases we gotta deal with are cases where the source and 
** destination architectures are of different word lengths and 
** those cases where the "short" is defined differently on the
** two machines.  (They are really the same problem.) 
** NOTE:  Because the way the system is implemented, VFF_TYP_2_BYTE
** really means "short".
*/

   src_size = getmachsize((unsigned long)src,(unsigned long)VFF_TYP_2_BYTE);
   dest_size = getmachsize((unsigned long)dest,(unsigned long)VFF_TYP_2_BYTE);

   if(src_size==dest_size)return(1);

   /* malloc tempspace to do the conversion */
   if((temp = (unsigned char *)malloc((unsigned)(num*dest_size)))==NULL) 
      return(0);

   iput = *data;

   for(i=0;i<num;i++){
       idxd = i*dest_size;   
       idxs = i*src_size;
       
       if(src_size>dest_size){
          temp[idxd] = (iput[idxs+src_size-dest_size] & (unsigned char)0x7f) + 
                       (iput[idxs] & (unsigned char)0x80);
          for(j=1;j<dest_size;j++)
              temp[idxd+j] = iput[idxs+src_size-dest_size+j];
       }else{
          if(iput[idxs] & (unsigned char)0x80){
             for(j=0;j<dest_size-src_size;j++)
                 temp[idxd+j] = (unsigned char)0xff;
          }else{
             for(j=0;j<dest_size-src_size;j++)
                 temp[idxd+j] = (unsigned char)0x0;
          }
          for(j=dest_size-src_size;j<dest_size;j++)
              temp[idxd+j] = iput[idxs-src_size+j];
       }
   }

   free(*data);
   *data = temp;

/*
** well, it all looks cool.
*/
 
   return(1);  
}
 
/********************************************************************/
/************************** L O N G    I N T  ***********************/
/********************************************************************/

/***********************************************************************
*
*  Routine Name: convert_long()
*
*          Date:
*        
*       Purpose:  
*
*         Input: data:  address of the input data cast to (unsigned char **)  
*                src:   unsigned int specifying source architecture.
*                dest:  another unsigned int specifying the destination 
*                       (usually current) architecture.
*                num:   unsigned int specifying number of short elements
*                       int *data (NOT the number of bytes in the array)
*
*        Output: This function returns a 1 for success and a 0 for
*                failure.  Failure indicates (1) that the programmer
*                didn't know what the heck he/she was doing; (2) that
*                this function wasn't able to realloc the data if
*                that was necessary; or (3) that the src or dest
*                arguments were invalid. 
*
*    Written By: Jeremy Worley 
*
* Modifications:
*
***********************************************************************/

int convert_long(data,src,dest,num)
   unsigned char **data;
   unsigned int src,dest,num;
{
   unsigned long src_size,dest_size,i,j,idxd,idxs;
   unsigned long getmachsize();
   unsigned char *temp,*iput;

/*
** test to see if the goofy programmer called this routine with
** src == dest 
*/

  if(src==dest)return(1);

/*
** go ahead and return saying that everything is cool if the 
** calling function passed in num=0.  
*/

  if(num==0)return(1);

/*
** test to see if src or dest are valid
*/

  if(src>=255 || dest>=255)return(0);   

/*
** test to see if the calling routine "knew" to malloc data before calling
** this one.
*/

   if(data == NULL)return(0);
   if(*data == NULL)return(0);

/*
** because we deal with byte ordering as a separate issue, the only
** real cases we gotta deal with are cases where the source and 
** destination architectures are of different word lengths and 
** those cases where the "short" is defined differently on the
** two machines.  (They are really the same problem.) 
** NOTE:  Because the way the system is implemented, VFF_TYP_4_BYTE
** really means "long".
*/

   src_size = getmachsize((unsigned long)src,(unsigned long)VFF_TYP_4_BYTE);
   dest_size = getmachsize((unsigned long)dest,(unsigned long)VFF_TYP_4_BYTE);

   if(src_size==dest_size)return(1);

   /* malloc tempspace to do the conversion */
   if((temp = (unsigned char *)malloc((unsigned)(num*dest_size)))==NULL) 
      return(0);

   iput = *data;

   for(i=0;i<num;i++){
       idxd = i*dest_size;   
       idxs = i*src_size;
       
       if(src_size>dest_size){
          temp[idxd] = (iput[idxs+src_size-dest_size] & (unsigned char)0x7f) + 
                       (iput[idxs] & (unsigned char)0x80);
          for(j=1;j<dest_size;j++)
              temp[idxd+j] = iput[idxs+src_size-dest_size+j];
       }else{
          if(iput[idxs] & (unsigned char)0x80){
             for(j=0;j<dest_size-src_size;j++)
                 temp[idxd+j] = (unsigned char)0xff;
          }else{
             for(j=0;j<dest_size-src_size;j++)
                 temp[idxd+j] = (unsigned char)0x0;
          }
          for(j=dest_size-src_size;j<dest_size;j++)
              temp[idxd+j] = iput[idxs-src_size+j];
       }
   }

   free(*data);
   *data = temp;

/*
** well, it all looks cool.
*/
 
   return(1);  
}

/***********************************************************************
*
*  Routine Name: convert_float()
*
*          Date:
*
*       Purpose:  
*
*         Input: data:  address of the input data cast to (unsigned char **)  
*                src:   unsigned int specifying source architecture.
*                dest:  another unsigned int specifying the destination 
*                       (usually current) architecture.
*                num:   unsigned int specifying number of short elements
*                       int *data (NOT the number of bytes in the array)
*
*        Output: This function returns a 1 for success and a 0 for
*                failure.  Failure indicates (1) that the programmer
*                didn't know what the heck he/she was doing; (2) that
*                this function wasn't able to realloc the data if
*                that was necessary; or (3) that the src or dest
*                arguments were invalid. 
*
*    Written By: Jeremy Worley 
*
* Modifications:
*
***********************************************************************/

int convert_float(data,src,dest,num)
   unsigned char **data;
   unsigned int src,dest,num;
{
   unsigned char c,*kern;
   unsigned int i,j,idx;
   long l;

/*
** test to see if the goofy programmer called this routine with
** src == dest or if the src and dest equals either VFF_DEP_NSORDER
** and VFF_DEP_IEEEORDER (which are the same).
*/

 if((src==dest) ||
     (src==VFF_DEP_IEEEORDER && dest==VFF_DEP_NSORDER) ||
     (src==VFF_DEP_NSORDER && dest==VFF_DEP_IEEEORDER))
  {
     return(1);
  }
 
/*
** go ahead and return saying that everything is cool if the 
** calling function passed in num=0.  
*/

  if(num==0)return(1);

/*
** test to see if src or dest are valid
*/

  if(src>=255 || dest>=255)return(0);   

/*
** test to see if the calling routine "knew" to malloc data before calling
** this one.
*/

   if(data == NULL)return(0);
   if(*data == NULL)return(0);

/*
** normal cases first:
*/
   kern = (*data);
 
   if((src==VFF_DEP_IEEEORDER || src==VFF_DEP_NSORDER) &&
      dest==VFF_DEP_DECORDER){
      for(i=0;i<num*4;i+=4){
	  /* check for negative zero */
          if(kern[i]==0x80 && (kern[i+1]|kern[i+2]|kern[i+3])==0x00)
             kern[i]=kern[i]&0x7f;
        
          c = (kern[i]<<1) + (kern[i+1]>>7) + 2;

          kern[i] = (kern[i]&0x80) + (c>>1);
          if(c&0x01)kern[i+1] = kern[i+1] | 0x80;

          
          c = kern[i];
          kern[i] = kern[i+2];
          kern[i+2] = c;

          c = kern[i+1];
          kern[i+1] = kern[i+3];
          kern[i+3] = c;
      }
   }else if(src==VFF_DEP_DECORDER && 
      (dest==VFF_DEP_IEEEORDER || dest==VFF_DEP_NSORDER)){
      for(i=0;i<num*4;i+=4){ 
          c = kern[i];
          kern[i] = kern[i+2];
          kern[i+2] = c;

          c = kern[i+1];
          kern[i+1] = kern[i+3];
          kern[i+3] = c;

          c = (kern[i]<<1) + (kern[i+1]>>7); /* get exponent */
          if(c==0)
             kern[i] = kern[i+1] = kern[i+2] = kern[i+3] = 0x0;
          else{
             c -= 2;
             kern[i] = (kern[i]&0x80) + (c>>1);
             if(c&0x01)kern[i+1] = kern[i+1] | 0x80;
          }
      }
   }else if(src==VFF_DEP_CRAYORDER && dest==VFF_DEP_DECORDER){
      for(i=0,j=0;i<num*8;i+=8,j+=4){
          l = ((long)(kern[i] & 0x7f))*256 + (long)kern[i+1] - 0x3f80;
          c = l & 0xff;
   
          kern[j] = (kern[i]&0x80) + (c>>1);
          if(c&0x01)
             kern[j+1] = kern[i+2] | 0x80; /* funny business here */
          else
             kern[j+1] = kern[i+2] & 0x7f;
          kern[j+2] = kern[i+3];
          kern[j+3] = kern[i+4];

          c = kern[j];
          kern[j] = kern[j+2];
          kern[j+2] = c;

          c = kern[j+1];
          kern[j+1] = kern[j+3];
          kern[j+3] = c;
       }
       if((kern = (unsigned char *)realloc(kern,num*4))==NULL){
          return(0);
       }
       *data = kern;
   }else if(src==VFF_DEP_CRAYORDER && 
     (dest==VFF_DEP_IEEEORDER || dest==VFF_DEP_NSORDER)){
      for(i=0,j=0;i<num*8;i+=8,j+=4){
          l = ((long)(kern[i] & 0x7f))*256 + (long)kern[i+1] - 0x3f82;
          c = l & 0xff;
   
          if(l==0)
             kern[j] = kern[j+1] = kern[j+2] = kern[j+3] = 0x0;
          else{
             kern[j] = (kern[i]&0x80) + (c>>1);
             if(c&0x01)
                kern[j+1] = kern[i+2] | 0x80; /* funny business here */
             else
                kern[j+1] = kern[i+2] & 0x7f;
             kern[j+2] = kern[i+3];
             kern[j+3] = kern[i+4];
          }
       }
       if((kern = (unsigned char *)realloc(kern,num*4))==NULL){
          return(0);
       }
       *data = kern;
   }else if(src==VFF_DEP_DECORDER && dest==VFF_DEP_CRAYORDER){
       for(i=0;i<4*num;i+=4){
          c = kern[i];
          kern[i] = kern[i+2];
          kern[i+2] = c;

          c = kern[i+1];
          kern[i+1] = kern[i+3];
          kern[i+3] = c;
       }
       if((kern = (unsigned char *)realloc(kern,num*8))==NULL){
          return(0);
       }
       *data = kern;
       for(idx=0,i=4*num-4,j=8*num-8;idx<num;idx++,i-=4,j-=8){
          l = ((kern[i] & 0x7f)<<1) + ((kern[i+1]&0x80)>>7) + 0x3f80;
         
          kern[j+5] = kern[j+6] = kern[j+7] = 0x0;
          kern[j+4] = kern[i+3];
          kern[j+3] = kern[i+2];
          kern[j+2] = (kern[i+1]&0x7f) + 0x80;
          kern[j+1] = (l & 0xff);
          kern[j] = (kern[i] & 0x80) + ((l>>8)&0x7f);
       
       }
   }else{ /* src=ieee/ns dest=cray */
       if((kern = (unsigned char *)realloc(kern,num*8))==NULL){
          return(0);
       }
       *data = kern;
       for(idx=0,i=4*num-4,j=8*num-8;idx<num;idx++,i-=4,j-=8){
          l = ((kern[i] & 0x7f)<<1) + ((kern[i+1]&0x80)>>7) + 0x3f82;
         
          kern[j+5] = kern[j+6] = kern[j+7] = 0x0;
          kern[j+4] = kern[i+3];
          kern[j+3] = kern[i+2];
          kern[j+2] = (kern[i+1]&0x7f) + 0x80;
          kern[j+1] = (l & 0xff);
          kern[j] = (kern[i] & 0x80) + ((l>>8)&0x7f);
       
       }
   }
   return(1);
}

/***********************************************************************
*
*  Routine Name: convert_double()
*
*          Date:
*        
*       Purpose:  
*
*         Input: data:  address of the input data cast to (unsigned char **)  
*                src:   unsigned int specifying source architecture.
*                dest:  another unsigned int specifying the destination 
*                       (usually current) architecture.
*                num:   unsigned int specifying number of short elements
*                       int *data (NOT the number of bytes in the array)
*
*        Output: This function returns a 1 for success and a 0 for
*                failure.  Failure indicates (1) that the programmer
*                didn't know what the heck he/she was doing; (2) that
*                this function wasn't able to realloc the data if
*                that was necessary; or (3) that the src or dest
*                arguments were invalid. 
*
*    Written By: Jeremy Worley 
*
* Modifications:
*
***********************************************************************/

int convert_double(data,src,dest,num)
   unsigned char **data;
   unsigned int src,dest,num;
{

/*
** test to see if the goofy programmer called this routine with
** src == dest 
*/

  if(src==dest)return(1);

/*
** go ahead and return saying that everything is cool if the 
** calling function passed in num=0.  
*/

  if(num==0)return(1);

/*
** test to see if src or dest are valid
*/

  if(src>=255 || dest>=255)return(0);   

/*
** test to see if the calling routine "knew" to malloc data before calling
** this one.
*/

   if(data == NULL)return(0);
   if(*data == NULL)return(0);

/*
** insert stuff here
*/

   if((src==VFF_DEP_IEEEORDER && dest==VFF_DEP_NSORDER) ||
	   (dest==VFF_DEP_IEEEORDER && src==VFF_DEP_NSORDER))
   {
	   unsigned char * cp = (unsigned char *)*data;
	   int i,j;

	   for(i=0; i<(int)num; ++i, cp+=8)
	   {
		   for(j=0; j<4; ++j)
		   {
			   unsigned char c = cp[j];
			   cp[j] = cp[7-j];
			   cp[7-j] = c;
		   }
	   }

	   return 1;
   }
/*
** case where nothing worked.
*/

   fprintf(stderr,"Sorry, byte order conversion for VFF_TYP_DOUBLE not implemented for this source machine.\n");
   return(0);
}

/***********************************************************************
*
* What follows are a set of obsolete routines that are kept in place
* just in case someone actually had the nerve and patience to use them.
*
***********************************************************************/

void ieeetonss(s)  /* From IEEE short int to NS  short int on a NS  machine */
short *s;
  {
    unsigned char *c,ctmp;
    c = (unsigned char *)s;
    ctmp = *c;
    *c = *(c+1);
    *(c+1) = ctmp;
    return;
  }

void nstoieees(s)  /* From NS short int to IEEE short int on a IEEE machine */
short *s;
  {
    unsigned char *c,ctmp;
    c = (unsigned char *)s;
    ctmp = *c;
    *c = *(c+1);
    *(c+1) = ctmp;
    return;
  }

void ieeetodecs(s) /* From IEEE short int to DEC short int on a DEC machine */
short *s;
  {
    unsigned char *c,ctmp;
    c = (unsigned char *)s;
    ctmp = *c;
    *c = *(c+1);
    *(c+1) = ctmp;
    return;
  }

void dectoieees(s)  /* From DEC short int to IEEE short int on an IEEE machine */
short *s;
  {
    unsigned char *c,ctmp;
    c = (unsigned char *)s;
    ctmp = *c;
    *c = *(c+1);
    *(c+1) = ctmp;
    return;
  }

void dectonss(s) /* From DEC short int to NS short int on a NS machine */
short *s;
  {
    return; /* No byte swap necessary */
  }

void nstodecs(s) /* From NS short int to DEC short int on a DEC machine */
short *s;
  {
    return; /* No byte swap necessary */
  }

void ieeetonsl(l) /* From IEEE long int to NS long int on a NS machine */
long *l;
  {
    unsigned char *c,ctmp;
    c = (unsigned char *)l;
    ctmp = *c;
    *c = *(c+3);
    *(c+3) = ctmp;
    ctmp = *(c+1);
    *(c+1) = *(c+2);
    *(c+2) = ctmp;
    return;
  } 

void nstoieeel(l) /* From NS long int to IEEE long int on a NS machine */
long *l;
  {
    unsigned char *c,ctmp;
    c = (unsigned char *)l;
    ctmp = *c;
    *c = *(c+3);
    *(c+3) = ctmp;
    ctmp = *(c+1);
    *(c+1) = *(c+2);
    *(c+2) = ctmp;
    return;
  } 

void ieeetodecl(l) /* From IEEE long int to DEC long int on a DEC machine */
long *l;
  {
    unsigned char *c,ctmp;
    c = (unsigned char *)l;
    ctmp = *c;
    *c = *(c+3);
    *(c+3) = ctmp;
    ctmp = *(c+1);
    *(c+1) = *(c+2);
    *(c+2) = ctmp;
    return;
  } 

void dectoieeel(l)  /* From DEC long int to IEEE long int on an IEEE machine */
long *l;
  {
    unsigned char *c,ctmp;
    c = (unsigned char *)l;
    ctmp = *c;
    *c = *(c+3);
    *(c+3) = ctmp;
    ctmp = *(c+1);
    *(c+1) = *(c+2);
    *(c+2) = ctmp;
    return;
  }

void dectonsl(l) /* From DEC long int to NS long int on a NS machine */
long *l;
  {
    return; /* No byte swap necessary */
  }

void nstodecl(l) /* From NS long int to IEEE long int on a NS machine */
long *l;
  {
    return; /* No byte swap necessary */
  }

/*************************************************************************/
/******************************* F L O A T *******************************/

void ieeetonsf(f) /* From IEEE float to NS float on a NS machine */
float *f;
  {
    unsigned char *c,ctmp;
    c = (unsigned char *)f;
    ctmp = *c;
    *c = *(c+3);
    *(c+3) = ctmp;
    ctmp = *(c+1);
    *(c+1) = *(c+2);
    *(c+2) = ctmp;
    return;
  } 

void nstoieeef(f) /* From NS float to IEEE float on a IEEE machine */
float *f;
  {
    unsigned char *c,ctmp;
    c = (unsigned char *)f;
    ctmp = *c;
    *c = *(c+3);
    *(c+3) = ctmp;
    ctmp = *(c+1);
    *(c+1) = *(c+2);
    *(c+2) = ctmp;
    return;
  } 

void ieeetodecf(f)  /* From IEEE float to DEC float on a DEC machine */ 
float *f;
  {
    unsigned char *c,ctmp;
    unsigned short *s;

    c = (unsigned char *)f; /* Do byte swaps */
    s = (unsigned short *)c;
    ctmp = *c;
    *c = *(c+1);
    *(c+1) = ctmp;
    ctmp = *(c+2);
    *(c+2) = *(c+3);
    *(c+3) = ctmp;
    if (*s == 0x8000 && *(s+1) == 0x0000)
      {
        /* Have detected an IEEE negative zero. Make this a positive
	   zero to keep the VAX frow giving a Reserved Operand trap. */
        *s = 0x0000;
      }
    *f *= 4;                /* Do left shift */
  }

void dectoieeef(f) /* From DEC float to IEEE float on an IEEE machine */
float *f;
  {
    unsigned char *c,ctmp;
    c = (unsigned char *)f;   /* Do byte swaps */
    ctmp = *c;
    *c = *(c+1);
    *(c+1) = ctmp;
    ctmp = *(c+2);
    *(c+2) = *(c+3);
    *(c+3) = ctmp;
    *f /= 4;                  /* Do right shift */
  }

void nstodecf(f)  /* From NS float to DEC float on a DEC machine */ 
float *f;
  {
    unsigned short *s,stmp;

    s = (unsigned short *)f;   /* Do byte swaps */
    stmp = *s;
    *s = *(s+1);
    *(s+1) = stmp;
    if (*s == 0x8000 && *(s+1) == 0x0000)
      {
        /* Have detected an IEEE negative zero. Make this a positive
	   zero to keep the VAX frow giving a Reserved Operand trap. */
        *s = 0x0000;
      }
    *f *= 4;                /* Do left shift */
  }

void dectonsf(f)  /* From DEC float to NS float on a NS machine */ 
float *f;
  {
    unsigned short *s,stmp;

    s = (unsigned short *)f;   /* Do byte swaps */
    stmp = *s;
    *s = *(s+1);
    *(s+1) = stmp;
    *f /= 4;                /* Do left shift */
  }

