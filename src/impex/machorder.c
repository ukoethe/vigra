/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   >>>>                                                       <<<<
   >>>>	    file name: machorder                              <<<<
   >>>>                                                       <<<<
   >>>>   description: Machine Order utilities                <<<<
   >>>>                                                       <<<<
   >>>>      routines: machorder()                            <<<<
   >>>>                                                       <<<<
   >>>> modifications:					      <<<<
   >>>>                                                       <<<<
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "machdefs.h"
#include "vdefines.h"

/***********************************************************************
*
*  Routine Name: machorder()
*
*          Date: 6/29/90
*        
*       Purpose:  This routine will return the machine order
*		  given a machine name. The Machine name must
*		  be contained in machdefs.h in include.
*
*         Input: mach_name - name of the machine
*
*        Output: machine order.
*
*    Written By:  Tom Sauer
*
* Modifications:
*
***********************************************************************/

long
machorder(mach_name)

char *mach_name;
{
int i;

  for (i = 0; i < NUMBER_MACHINES; i++)
  {
     if (strcmp(mach_name, machine_defs[i].hosttype) == 0)
       return(machine_defs[i].order);
  }
  return(255);  /* an error has occured. This keeps consistent with
                 * machtype() 
                 */
 
}

/***********************************************************************
*
*  Routine Name: getmachorder()
*
*          Date: Tue Sep 24 16:20:33 MDT 1991
*        
*       Purpose: returns the word ordering for a given machine.
*                this routine probably ought to be called machorder(),
*                but that name was stolen by the above routine to 
*                return the machtype. Hmmm.
*
*         Input: mtype - type of machine...returned by machorder().
*
*        Output: returns a number corresponding to VFF_DEP_BIGENDIAN or
*                VFF_DEP_LITENDIAN from the table machorder[] in
*                machdefs.h
*
*    Written By: Jeremy Worley 
*
* Modifications:
*
***********************************************************************/

long getmachorder(mtype)
  long int mtype;
{

 switch((int)mtype){
    case VFF_DEP_CRAYORDER:
    case VFF_DEP_IEEEORDER: return(VFF_DEP_BIGENDIAN);
    case VFF_DEP_NSORDER:    
    case VFF_DEP_DECORDER:  return(VFF_DEP_LITENDIAN);
 }

 return((long)255);
   
}


/***********************************************************************
*
*  Routine Name: 
*
*          Date: Thu Oct  3 15:34:48 MDT 1991
*        
*       Purpose: Returns the length in bytes of a given data type
*                that is supported in khoros on a particular
*                architecture. 
*
*         Input: mtype - machine type as defined by VFF_DEP's
*                dtype - data type as defined by VFF_TYP's
*
*        Output: size of that data type in bytes 
*
*    Written By:  Jeremy Worley
*
* Modifications: J.W. Thu Oct 31 07:37:35 MST 1991
*			Changed complex and double complex cases to 
*		 	return the same values as float.  This is because
*		 	we want to know the size of each component during 
*		 	conversion.
*
***********************************************************************/

unsigned long getmachsize(mtype,dtype)
   unsigned long mtype,dtype;
{
   unsigned long tmp = (mtype==VFF_DEP_CRAYORDER) + 1;
   switch((int)dtype){
      case VFF_TYP_BIT     : return((unsigned long)0);
      case VFF_TYP_1_BYTE  : return((unsigned long)1);
      case VFF_TYP_2_BYTE  : if(mtype==VFF_DEP_CRAYORDER)
                                return((unsigned long)8);
                             else
                                return((unsigned long)2);
      case VFF_TYP_4_BYTE  : return((unsigned long)4*tmp);
      case VFF_TYP_FLOAT   : return((unsigned long)4*tmp);
      case VFF_TYP_DOUBLE  : return((unsigned long)8);
      case VFF_TYP_COMPLEX : return((unsigned long)8*tmp);
      case VFF_TYP_DCOMPLEX: return((unsigned long)16);
      default: return((unsigned long)255);
   }
}

/***********************************************************************
*
*  Routine Name: getmachlist()
*
*          Date: Sat Nov  9 16:35:39 MST 1991
*        
*       Purpose: getmachlist() returns all the host type names to
*
*         Input: None
*
*        Output: returns the number of corresponding machine entries
*		 and there number.
*
*    Written By: Jeremy Worley 
*
* Modifications:
*
***********************************************************************/


#define ALLOCSIZE 100
char **getmachlist(num_entries)

int	*num_entries;
{
	int	i, j, num = ALLOCSIZE;
	char	**list;


	if ((list = (char **) malloc(num * sizeof(char *))) == NULL)
	{
	   fprintf(stderr,"getmachlist:  Not enough memory....\n\n");
	   fprintf(stderr,"  Unable to malloc (%d) bytes for the khoros \
machine defintion list.\n", num*sizeof(char *));
	   return(NULL);
	}

	for (i = 0, j = 0; i < NUMBER_MACHINES; i++)
	{
	    if (machine_defs[i].hosttype != NULL)
	    {
	       if (j >= num)
	       {
	          num += ALLOCSIZE;
		  if ((list = (char **) realloc(list, num * sizeof(char *))) ==
				NULL)
		  {
		     fprintf(stderr,"getmachlist:  Not enough memory..\n\n");
		     fprintf(stderr,"  Unable to malloc (%d) bytes for the \
khoros transport list.\n", num*sizeof(char *));
		     return(NULL);
		  }
	       }
	       list[j++] = VStrcpy(machine_defs[i].hosttype);
	    }
	}
	*num_entries = j;
	return(list);
}
