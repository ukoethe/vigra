/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   >>>>                                                       <<<<
   >>>>	    file name: vgparm.h                               <<<<
   >>>>                                                       <<<<
   >>>>   description: indescribable                          <<<<
   >>>>                                                       <<<<
   >>>> modifications:					      <<<<
   >>>>                                                       <<<<
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


#ifndef _VGPARMS_H_
#define _VGPARMS_H_


/*
 * Varargs parameter list handling...YECH!!!!
 */
#if !defined(__STDC__)
#include <varargs.h>
#define	VA_START(arg, param)	va_start(arg)
#else
#include <stdarg.h>
#define	VA_START(arg, param)	va_start(arg, param)
#endif

#define VA_END(arg) va_end(arg)

#define MAX_ARGS 200
#define MAX_ARG_LEN 512
#define MAX_LINE_LEN 512
#define MAX_KEY_LEN 40
#define BUFFERSIZE 512
#define MODE 00644

typedef struct arg_entry {
		char key[MAX_KEY_LEN];
		char sarg[MAX_ARG_LEN];
		struct arg_entry *next;
		} ARG_ENTRY;

char *_cleanup_string();
char *_expand_tilda();
char *_expand_variable();

#endif /* _VGPARMS_H_ */
/* don`t add after the endif */
