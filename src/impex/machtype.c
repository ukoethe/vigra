/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/


#include <stdlib.h>
#include <string.h>
#include "machdefs.h"	

long
machtype(name)
char *name;
{
	if (name != NULL)
	   (void) strcpy(name, LocalDef.hosttype);

	return(LocalDef.order);
}
