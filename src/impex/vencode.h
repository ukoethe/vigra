/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/



#ifndef	_vencode_h_
#define _vencode_h_

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   >>>>                                                              <<<<
   >>>>      file: vencode.h					     <<<<
   >>>>                                                              <<<<
   >>>>      contains: 
   >>>>                This file contains the tables of commands     <<<<
   >>>>      used to implement the automatic encoding/decoding       <<<<
   >>>>      of viff images.					     <<<<
   >>>>                                                              <<<<
   >>>>      The table entitled "compr_cmd" associates the actual    <<<<
   >>>>      command required to encode the data with the encoding   <<<<
   >>>>      scheme declared in KHOROS_HOME/include/viff.h.          <<<<
   >>>>                                                              <<<<
   >>>>      The table entitled "uncompr_cmd" associates the actual  <<<<
   >>>>      command required to decode the data with the encoding   <<<<
   >>>>      scheme declared in  KHOROS_HOME/include/viff.h.         <<<<
   >>>>                                                              <<<<
   >>>>      In both cases, the actual command should be written to  <<<<
   >>>>      accept stdin for the input, and stdout for the output   <<<<
   >>>>      (this is the only requirement as far as viff support is <<<<
   >>>>      required)  For an example of this, see the UNIX         <<<<
   >>>>      compress(1) command.				     <<<<
   >>>>                                                              <<<<
   >>>>      NOTE: The command given will be implanted in a command  <<<<
   >>>>            command string similar to:			     <<<<
   >>>>            % cat input | your-command-here | > output        <<<<
   >>>>            so you can put whatever switches, etc, you want!  <<<<
   >>>>		   Just remember that >> sh <<, not >> csh << will   <<<<
   >>>>		   execute your command!  If the command does not    <<<<
   >>>>		   exist, the support routines will try to detect    <<<<
   >>>>		   this fact and exit cleanly.			     <<<<
   >>>>								     <<<<
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
  
#ifndef	_compr_cmd_
#define _compr_cmd_
static char *compr_cmd[] =
	{"This_should_never_happen",      /* DES-RAW */
         "compress",                      /* DES_COMPRESS */
         "rle",                           /* DES_RLE */
         "transform",                     /* DES_TRANSFORM */
         "ccitt",                         /* DES_CCITT */
         "adpcm",                         /* DES_ADPCM */
         "gencomp"};                      /* DES_GENERIC */
#endif

#ifndef	_uncompr_cmd_
#define _uncompr_cmd_
static char *uncompr_cmd[] =
	{"This_should_never_happen",      /* DES-RAW */
         "uncompress",                    /* DES_COMPRESS */
         "unrle",                         /* DES_RLE */
         "untransform",                   /* DES_TRANSFORM */
         "unccitt",                       /* DES_CCITT */
         "unadpcm",                       /* DES_ADPCM */
         "ungencomp"};                    /* DES_GENERIC */
#endif

/*
  In case you're wondering why This_should_never_happen, you'll
  find that readimage() and writeimage() will trap for VFF_DES_RAW
  and call write() or block_read() directly for efficiency!

*/
  
#endif /* _vencode_h_ */
/* Don't add after the endif */
