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
#include "vdefines.h"
#include "vgparm.h"	


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   >>>>                                                       <<<<
   >>>>	    file name: fullpath.c                             <<<<
   >>>>                                                       <<<<
   >>>>   description: File utility                           <<<<
   >>>>                                                       <<<<
   >>>>      routines: vfullpath()			      <<<<
   >>>>                                                       <<<<
   >>>> modifications:					      <<<<
   >>>>                                                       <<<<
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

#define MAX_VARIABLES 1
static char *environ_variables[] =
{
		"KHOROS_HOME",
};

/**************************************************************
*
* MODULE NAME: vfullpath
*
*     PURPOSE: This function returns the full path to the
*	       user supplied file.  "vfullpath" expands
*	       path and returns a full path back the user.
*	       If the file is specified by user name then
*	       the path will be expanded.  The routine will
*	       also expand the KHOROS_HOME environmental
*	       variable if defined by the user.
*
*	       The user can can pass in a global directory "global_dir"
*	       which will be prepended to the front of the filename.
*	       The global directory will be prepended to the filename if
*	       
*	       
*
*       INPUT: filename    - the file to be expanded
*	       global_dir  - the global directory is used as a prefix
*			     to the filename.  If global_dir is
*	       return_file - if the return file is not NULL then
*			     this will be where we place the vfullpath
*			     to the expanded file.  Otherwise we malloc
*			     space and return the fille.
*	       
*        
*      OUTPUT: the full path to the expanded file.
*
* CALLED FROM: 
*
*  WRITTEN BY:  Mark Young
*
**************************************************************/


char *vfullpath(filename, global_dir, return_file)

char	*filename, *global_dir, *return_file;
{
	int	i;
	char	file[LENGTH], directory[LENGTH], *dir, *buffer,
		*vreplace_keyword();
	char    *_cleanup_string(), *_expand_variable(), *_expand_tilda();


	/*
	 *  If the filename is null then error and return NULL
	 */
	if (filename == NULL)
	{
	   (void) fprintf(stderr,"\nvfullpath:\n");
	   (void) fprintf(stderr,"   Error!  NULL input file encountered.\n");
	   return(NULL);
	}

	/*
	 *  If the filename is  empty then error and return NULL
	 */
	strcpy(file, filename);
	buffer = _cleanup_string(file);
	if (VStrlen(buffer) == 0)
	{
	   (void) fprintf(stderr,"\nvfullpath:\n");
	   (void) fprintf(stderr,"   Error!  Empty input file encountered.\n");
	   return(NULL);
	}


       if (global_dir == NULL)
	  dir = NULL;
       else if (VStrlen(global_dir) > 0)
       {
	  strcpy(directory, global_dir);
	  dir = _cleanup_string(directory);
       }
       else
	  dir = NULL;

       /*
	*  Now check to see if the filename or global_dir contains
	*  KHOROS_HOME, in which case we will need 
	*  to expand the string first.
	*/
       if ((buffer = _expand_variable(buffer, NULL)) == NULL)
	  return(NULL);

       if (dir != NULL)
       {
	  dir = _expand_variable(dir, NULL);
	  if (dir == NULL) return(NULL);
       }

       for (i = 0; i < MAX_VARIABLES; i++)
       {
	   buffer = _expand_variable(buffer, environ_variables[i]);
	   if (buffer == NULL) return(NULL);

	   if (dir != NULL)
	   {
	      dir = _expand_variable(dir, environ_variables[i]);
	      if (dir == NULL) return(NULL);
	   }
       }

       /*
	*  Expand the username if the string begins with a '~'.
	*/
       if (buffer[0] == '~' || buffer[0] == '/' || buffer[0] == '.')
       {
	  if (!(buffer = _expand_tilda(buffer)))
	     return(NULL);

	  dir = NULL;
       }
       else if (dir != NULL)
       {
	  dir = _expand_tilda(dir);
       }

       /*
	*  Return the fullpath back to the user.
	*/
       if (dir != NULL)
       {
	  i = VStrlen(dir);
	  if (i > 0 && dir[i -1] != '/')
	     strcat(dir, "/");

	  strcat(dir, buffer);
       }
       else
	  dir = buffer;

	if (return_file == NULL)
	   return(VStrcpy(dir));
	else
	{
	   strcpy(return_file, dir);
	   return(return_file);
	}
}
