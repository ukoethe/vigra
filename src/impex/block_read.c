/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/

#include <sys/types.h>

#ifdef _MSC_VER
#    include <io.h>
#else
#    include <unistd.h>
#    include <sys/uio.h>
#endif

#include "readbuffer.h"

/************************************************************
*
*  MODULE NAME: block_read
*
*      PURPOSE: block_read is a stupid little routine that performs
*		a block read on the number of desired blocks.  This
*		useful when the input is coming from standard in, since
*		read can return without reading all the desired blocks.
*
*        INPUT: 1.  fid   --  the file descriptor in which to read from
*		2.  input --  a character pointer where the bytes will be 
*			      read and stored to
*		3.  size --   an integer specifying the number of bytes to 
*			      to be read.
*
*       OUTPUT: Returns the number of bytes read.
*
*    CALLED BY: any routine that wishes to block read
*
*
*************************************************************/


int block_read(fid, input, size)

int	fid;
char	*input;
int	size;
  {
    int	numread = 0, nbytes = 1, num = 0;

    if (fid == 0 && num_buffer > 0)
      {
	while (num_buffer > 0 && num < size)
          {
	    *(input + num) = read_buffer[num];
	    num_buffer--; num++; numread++;
	  }
      }

        while ((numread < size) && (nbytes >0))
          {
            nbytes = read(fid, (char *) (input + numread), size - numread);
	    numread += nbytes;
          }
    return(numread);
  }
