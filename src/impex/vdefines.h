/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   >>>>                                                          <<<<
   >>>>      file: vdefines.h                                    <<<<
   >>>>                                                          <<<<
   >>>>      contains:  common definitions needed by             <<<<
   >>>>			Khoros programs				 <<<<
   >>>>                                                          <<<<
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

#ifndef _vdefines_h_
#define _vdefines_h_

#define	XV_FILE_MAGIC_NUM	0xab	/* Khoros file identifier */
#define	XV_FILE_TYPE_XVIFF	1	/* indicates an image file */

#undef  TRUE
#define TRUE 1

#undef  FALSE
#define FALSE 0

#undef LENGTH
#define LENGTH 512

#undef  COMMENTLEN
#define COMMENTLEN 512


/*
 *  Need to see if we are on a system V machine.  If so then we need to define
 *  bcopy(), bzero(), bcmp() to use memcpy(), memset(), memcmp(), LOCK_EX,
 *  LOCK_UN, lockf
 */
#ifdef SYSV

#ifndef bcopy
#define bcopy(s1, s2, n)	memcpy(s2, s1, n)
#endif

#ifndef bzero
#define bzero(s1, n)		memset(s1, 0, n)
#endif

#ifndef bcmp
#define bcmp(s1, s2, n)		memcmp(s2, s1, n)
#endif

#if !defined(_IBMR2) && !defined(sgi)
#ifndef LOCK_EX
#define LOCK_EX F_TLOCK
#endif


#ifndef LOCK_UN
#define LOCK_UN F_ULOCK
#endif

#ifndef LOCK_SH
#define LOCK_SH F_TLOCK
#endif

#ifndef flock
#define flock(fids, function) lockf(fids, function, 0) 
#endif
#endif

#endif

#ifdef BSDTYPES

#ifndef u_char
typedef unsigned char u_char;
#endif
#ifndef u_short
typedef unsigned short u_short;
#endif
#ifndef u_int
typedef unsigned int u_int;
#endif
#ifndef u_long
typedef unsigned long u_long;
#endif

#endif

/*
 *  Some machines like the luna88k do not define errno
 *  if errno is not defined, then we define it
 */

#ifndef errno
extern int errno;
#endif


/*
 *  Defines used for vlistdir for it's list mode
 */
#define XV_PATH	 (1 << 0)
#define XV_FILE	 (1 << 1)
#define XV_DIR	 (1 << 2)
#define XV_DOT	 (1 << 3)
#define XV_LINK	 (1 << 4)
#define XV_SOCK	 (1 << 5)
#define XV_READ	 (1 << 6)
#define XV_WRITE (1 << 7)
#define XV_EXEC	 (1 << 8)


/*
 *  Need to make sure that S_IFSOCK & S_IFLNK are defined.  Normally these
 *  defines are defined in <sys/stat.h>, but on certain POSIX machines these
 *  are not defined.  If they are not defined, we then define them to something
 *  meaningless (ie. -1, -2).  Also, make sure that MAXPATHLEN is defined.
 */
#ifndef S_IFSOCK
#define S_IFSOCK -1
#endif

#ifndef S_IFLNK
#define S_IFLNK -2
#endif

#ifdef UNICOS
#define MAXPATHLEN    PATHSIZE
#endif

#ifndef MAXPATHLEN
#define MAXPATHLEN    1024
#endif


/*
 *  Defines for the dsp routines
 */
#undef  DSP_VECTOR 
#define DSP_VECTOR  0

#undef  DSP_BAND 
#define DSP_BAND 1


/*
 *  Declare some of the basic basic string manipulation routines
 */
#if 0
#ifndef AUX
#  ifndef malloc
       char *malloc();
#  endif

#  ifndef calloc
       char *calloc();
#  endif

#  ifndef strchr
       char *strchr();
#  endif

#  ifndef strrchr
       char *strrchr();
#  endif
#endif /* AUX */
#endif /* 0 */

/*  This is a KLUDGE to take care of a new read_raw for patch 3 */

#define read_raw(fid, offset, points, type, mach_type) \
         new_read_raw(fid, offset, points, type, mach_type, TRUE)

#define MAXBUFF 	20      /* max # of command line string args */

/**************************************************************
*
* MODULE NAME: VStrlen
*     PURPOSE: a better strlen (will check for null pointers)
*       INPUT: a character pointer
*      OUTPUT: the length of the string
* CALLED FROM: EVERY place that uses strlen
*
**************************************************************/

#define	   VStrlen(strng)	((strng != NULL) ? strlen(strng) : 0)

/**************************************************************
*
* MODULE NAME: VStrcat
*     PURPOSE: a memory allocation  strcat (will malloc for me)
*       INPUT: s1 - a character pointer to the first string
*              s2 - a character pointer to the second string
*      OUTPUT: a pointer to the concatenated string
* CALLED FROM: EVERY place that needs to malloc inside of strcpy
*
**************************************************************/

#define VStrcat(s1, s2) (strcat(strcpy(malloc(VStrlen(s1) + VStrlen(s2) +1), \
				(s1)), (s2)))

/**************************************************************
*
* MODULE NAME: VStr3cat
*     PURPOSE: a memory allocation  strcat (will malloc for me)
*	       of three strings
*       INPUT: s1 - a character pointer to the first string
*              s2 - a character pointer to the second string
*	       s3 - a character pointer to the third string
*      OUTPUT: a pointer to the concatenated string
* CALLED FROM: EVERY place that needs to malloc inside of strcpy
*
**************************************************************/

#define VStr3cat(s1, s2, s3) (strcat(strcat(strcpy(malloc(VStrlen(s1) + \
			    VStrlen(s2) + VStrlen(s3)+1), (s1)), (s2)),(s3)))

/**************************************************************
*
* MODULE NAME: VStrcpy
*     PURPOSE: a memory allocation  strcpy (will calloc for me)
*       INPUT: a character pointer
*      OUTPUT: a pointer to a calloc'ed string
*
*
**************************************************************/

#define	VStrcpy(strng) ((strng != NULL) ? (strcpy(calloc((unsigned)1, \
				(unsigned)(VStrlen(strng)+1)),strng)) : NULL)

/**************************************************************
*
* MODULE NAME: VStrncpy
*     PURPOSE: a memory allocation  strncpy (will calloc for me)
*       INPUT: strng - a character pointer to string to copy
*	       n - integer containing number of chars to copy
*      OUTPUT: a pointer to a calloc'ed string
*
**************************************************************/

#define	VStrncpy(strng,n) ((strng != NULL && n > 0) ? \
				(strncpy(calloc((unsigned)1, \
				(unsigned)(n+1)),strng,n)) : NULL)


/**************************************************************
*
* MODULE NAME: VStrcmp
*     PURPOSE: VStrcmp is a frontend for the system's strcmp().  The
*	       only difference is that VStrcmp() protects against
*	       NULL strings.  The routine compares the two input
*	       strings and returns whether str1 is less than str2 "-1",
*	       str1 equal to str2 "0", str1 greater than str2 "1".  If
*	       both or either of the strings are NULL then the define
*	       protects against it by saying that the NULL string
*	       is less than the other (unless both strings are NULL
*	       in which case they are returned as being equal "0").
*
*			if str1 is less str2     - return -1
*			if str1 is equal str2    - return  0
*			if str1 is greater str2  - return  1
*
*       INPUT: str1 - a character pointer to the first string
*              str2 - a character pointer to the second string
*      OUTPUT: returns 
*
**************************************************************/

#define	VStrcmp(str1, str2) ((str1 != NULL && str2 != NULL) ? strcmp(str1,str2)\
	 : ((str1 == NULL && str2 == NULL) ? 0 : ((str1 == NULL) ? -1 : 1)))

/**************************************************************
*
* MODULE NAME: VStrncmp
*     PURPOSE: VStrncmp is a frontend for the system's strncmp().  The
*	       only difference is that VStrncmp() protects against
*	       NULL strings.  The routine compares the two input
*	       strings and returns whether str1 is less than str2 "-1",
*	       str1 equal to str2 "0", str1 greater than str2 "1".  If
*	       both or either of the strings are NULL then the define
*	       protects against it by saying that the NULL string
*	       is less than the other (unless both strings are NULL
*	       in which case they are returned as being equal "0").
*
*			if str1 is less str2     - return -1
*			if str1 is equal str2    - return  0
*			if str1 is greater str2  - return  1
*
*       INPUT: 1. str1 - a character pointer to the first string
*       	  str2 - a character pointer to the second string
*		  n - number of chars to compare
*
*      OUTPUT: returns 
*
*
**************************************************************/

#define	VStrncmp(str1, str2, n) ((str1 != NULL && str2 != NULL) ? \
	  strncmp(str1,str2,((n <= 0) ? 1 : n)) \
	 : ((str1 == NULL && str2 == NULL) ? 0 : ((str1 == NULL) ? -1 : 1)))


/************************************************************
*
*  MODULE NAME: kmalloc
*
*      PURPOSE: Calls malloc() to allocate the requested data.  Currently
*		the only check is to make sure that the user doesn't alloc
*		0 bytes, since this will fail on certain architectures.
*	        If the user requests 0 bytes the request will be changed
*		to alloc a single byte.
*
*        INPUT: size --  the number of desired bytes
*       OUTPUT: Returns the number of bytes alloc'ed or NULL upon failure.
*
*************************************************************/

#define kmalloc(size) malloc(((size) == 0) ? 1 : (size))


/************************************************************
*
*  MODULE NAME: kalloca
*
*      PURPOSE: Calls alloca() to allocate the requested data.  Currently
*		the only check is to make sure that the user doesn't alloc
*		0 bytes, since this will fail on certain architectures.
*	        If the user requests 0 bytes the request will be changed
*		to alloc a single byte.
*
*        INPUT: size --  the number of desired bytes
*       OUTPUT: Returns the number of bytes alloc'ed or NULL upon failure.
*
*************************************************************/

#define kalloca(size) alloca(((size) == 0) ? 1 : (size))


/************************************************************
*
*  MODULE NAME: kcalloc
*
*      PURPOSE: Calls calloc() to allocate the requested data.  Currently
*		the only check is to make sure that the user doesn't alloc
*		0 bytes, since this will fail on certain architectures.
*	        If the user requests 0 bytes the request will be changed
*		to alloc a single byte.
*
*        INPUT: nelem  --  the number of desired elements
*        	elsize --  the number of bytes per element
*       OUTPUT: Returns the number of bytes alloc'ed or NULL upon failure.
*
*************************************************************/

#define kcalloc(nelem, elsize) calloc((((nelem) == 0) ? 1 : (nelem)), \
				      (((elsize) == 0) ? 1 : (elsize)))

/************************************************************
*
*  MODULE NAME: krealloc
*
*      PURPOSE: Calls realloc() to re-allocate the requested data.  Currently
*		the only check is to make sure that the user doesn't re-alloc
*		0 bytes, since this will fail on certain architectures.
*	        If the user requests 0 bytes the request will be changed
*		to alloc a single byte.
*
*        INPUT: ptr -- the previously allocated pointer
*		size --  the number of desired bytes
*       OUTPUT: Returns the number of bytes alloc'ed or NULL upon failure.
*
*************************************************************/

#define krealloc(ptr, size) (((ptr) == NULL) ? \
				(malloc(((size) == 0) ? 1 : (size))) : \
				((char *)realloc(ptr,((size)==0) ? 1 : (size))))


/************************************************************
*
*  MODULE NAME: kfree
*
*      PURPOSE: Calls free() to free previously allocated data.  Currently
*		the only check is to make sure that the user doesn't free
*		a NULL pointer.
*
*        INPUT: 1.  ptr --  the previously allocated pointer
*
*       OUTPUT: 
*
*
*************************************************************/

#define kfree(ptr) if ((ptr) != NULL) free(ptr)


/*
 * library routine declarations
 */
char *vbasename();
char *vfullpath();
char *vtempnam();
struct xvimage *readimage();
struct xvimage *createimage();
struct xvimage *createsameimage();
struct xvimage *copyimage();
struct xvimage **create_image_list();
float *new_read_raw();
float *read_ascii();
char **load_vector();
char **dload_vector();
char *unload_vector();

char *vlower_string();
char *vupper_string();
char *vreplace_string();
char *vreplace_char();
char *vcleanup_string();
char *vstrstr();
char *vstrtok();
char *vstrpbrk();

char **vlistdir();
char **vlistfile();
char **vmergelist();
char **vsortlist();
char **vcopylist();
int  vfreelist();

/*
 * Define for use with all the Khoros tools
 */

#define MACH_FILE_PATH "repos/config/src_conf"

/*
 * Handy defines for getting around in image data arrays
 */

/* PIXEL - For indexing into SINGLE BAND 2D arrays, incurs less computation
           then BPIXEL. */
#define PIXEL(x,y,rows,cols) (y)*(cols)+(x) 

/* BPIXEL - For indexing into MULTI-BAND 2D arrays. */
#define BPIXEL(b,x,y,rows,cols) (b)*(rows)*(cols)+(y)*(cols)+(x)

#endif /* _vdefines_h_ */
/* Don't add after this point. */
