/************************************************************************/
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/


#ifndef _machdefs_h_
#define _machdefs_h_

#include "vigra/viff.h"

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   >>>>                                                       <<<<
   >>>>	    file name: machdefs.h                             <<<<
   >>>>                                                       <<<<
   >>>>      contains: machine definitions for correct	      <<<<
   >>>>                byte ordering and known		      <<<<
   >>>>                                                       <<<<
   >>>>      written by:                                      <<<<
   >>>>                                                       <<<<
   >>>>      date: 3/5/90                                     <<<<
   >>>>                                                       <<<<
   >>>>      modifications:				      <<<<
   >>>>                                                       <<<<
   >>>>>>>>>>>>>>>>>>>>>>>>>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


typedef struct machine_definitions
{
	char	*hosttype;
	long	order;
	long	byte_order;
	long	float_order;
} MachineDefs;


static MachineDefs machine_defs[] =
{
#define APOLLO_DEF	0
   { "Apollo",			  VFF_DEP_IEEEORDER,	0,	0 },

#define DECSTATION_DEF	1
   { "DECstation",		  VFF_DEP_NSORDER,	0,	0 },

#define NS32000_DEF	2
   { "Encore, Sequent, NS32000",  VFF_DEP_NSORDER,	0,	0 },

#define IBMRT_DEF	3
   { "IBM RT",			  VFF_DEP_IEEEORDER,	0,	0 },

#define IBM032_DEF	4
   { "IBM032",			  VFF_DEP_IEEEORDER,	0,	0 },

#define MC680x0_DEF	5
   { "MC680x0", 		  VFF_DEP_IEEEORDER,	0,	0 },

#define MC68k32_DEF	6
   { "MC68k32",			  VFF_DEP_IEEEORDER,	0,	0 },

#define SONY_DEF	7
   { "Sony News",		  VFF_DEP_IEEEORDER,	0,	0 },

#define SUN_DEF		8
   { "Sun",			  VFF_DEP_IEEEORDER,	0,	0 },

#define TITAN_DEF	9
   { "Titan",			  VFF_DEP_IEEEORDER,	0,	0 },

#define VAX_DEF		10
   { "Vax",			  VFF_DEP_DECORDER,	0,	0 },

#define MIPSEB_DEF	11
   { "Mips cpu, big endian mode", VFF_DEP_IEEEORDER,	0,	0 },

#define MIPS_DEF	12
   { "MIPS",			  VFF_DEP_IEEEORDER,	0,	0 },

#define SGI_DEF		13
   { "Silicon Graphics",	  VFF_DEP_IEEEORDER,	0,	0 },

#define NEXT_DEF	14
   { "NeXT",			  VFF_DEP_IEEEORDER,	0,	0 },

#define HP9000_DEF	15
   { "hp9000",			  VFF_DEP_IEEEORDER,	0,	0 },

#define RISC6000_DEF	16
   { "RISC 6000",		  VFF_DEP_IEEEORDER,	0,	0 },

#define CRAY_DEF	17
   { "Cray",			  VFF_DEP_CRAYORDER,	0,	0 },

#define CONVEX_DEF	18
   { "Convex",			  VFF_DEP_IEEEORDER,	0,	0 },

#define MACII_DEF	19
   { "MacII",			  VFF_DEP_IEEEORDER,	0,	0 },

#define DGUX_DEF	20
   { "Data General",		  VFF_DEP_IEEEORDER,	0,	0 },

#define I386_DEF	21
   { "386/486",			  VFF_DEP_NSORDER, 	0, 	0 },
 
#define M88K_DEF	22	
   { "Motorola 88000",		  VFF_DEP_IEEEORDER,	0,	0 },

#define LUNA88K_DEF	23
   { "luna",			  VFF_DEP_IEEEORDER,	0,	0 },

#define UNKNOWN_DEF	24
   { "Unknown",			  255,			-1,	-1 },
         /* this should always be the last entry in the list */
};

#if defined(_MSC_VER)
#    define i386
#endif

                     /* This one for the Vax */
#if defined(vax) && !defined(LocalDef)
#define LocalDef	machine_defs[VAX_DEF]
#endif

                     /* This one for the Sun */
#if defined(sun) && !defined(LocalDef)
#define LocalDef	machine_defs[SUN_DEF]
#endif

                     /* This one for Sony */
#if defined(sony_news) && !defined(LocalDef)
#define LocalDef	machine_defs[SONY_DEF]
#endif

                     /* This one for Silicon Graphics */ 
#if defined(sgi) && !defined(LocalDef) 
#define LocalDef	machine_defs[SGI_DEF]
#endif

                    /* This one for NeXT */
#if defined(NeXT) && !defined(LocalDef) 
#define LocalDef	machine_defs[NEXT_DEF]
#endif

#if defined(mc68000) && !defined(LocalDef)
#define LocalDef	machine_defs[MC680x0_DEF]
#endif

#if defined(mc68k32) && !defined(LocalDef)
#define LocalDef	machine_defs[MC68k32_DEF]
#endif

                   /* This one is for the ENCORE and SEQUENT */
#if defined(ns32000) && !defined(LocalDef)
#define LocalDef	machine_defs[NS32000_DEF]
#endif

                   /* This one is for the Ardent Titan   */
#if defined(titan) && !defined(LocalDef)
#define LocalDef	machine_defs[TITAN_DEF]
#endif

			 /* This one is for the DECstation 3100's */
#if defined(MIPSEL) && !defined(LocalDef)
#define LocalDef	machine_defs[DECSTATION_DEF]
#endif

			 /* This one is for the DECstation 3100's */
#if defined(MIPSEB) && !defined(LocalDef)
#define LocalDef	machine_defs[MIPSEB_DEF]
#endif

			 /* This one is for generic MIPS machines */
#if defined(mips) && !defined(LocalDef)
#define LocalDef	machine_defs[MIPS_DEF]
#endif

                     /* This one for IBM-RT's running AIX */
#if defined(aiws) && !defined(LocalDef)
#define LocalDef	machine_defs[IBMRT_DEF]
#endif

                    /* This one for IBM-RT's running IBM 4.3 */
#if defined(ibm032) && !defined(LocalDef)
#define LocalDef	machine_defs[IBM032_DEF]
#endif

                   /* This works for APOLLO DN10000 */
#if defined(apollo)  && !defined(LocalDef)
#define LocalDef	machine_defs[APOLLO_DEF]
#endif

			/* this works for hp9000s{300,400,700,800} */
#if defined(hp9000) && !defined(LocalDef)
#define LocalDef	machine_defs[HP9000_DEF]
#endif

		/* this works for hp9000s{300,400,700,800} */
#if defined(hp9000s300) && !defined(LocalDef)
#define LocalDef	machine_defs[HP9000_DEF]
#endif

		/* this works for hp9000s{300,400,700,800} */
#if defined(hp9000s400) && !defined(LocalDef)
#define LocalDef	machine_defs[HP9000_DEF]
#endif

		/* this works for hp9000s{300,400,700,800} */
#if defined(hp9000s700) && !defined(LocalDef)
#define LocalDef	machine_defs[HP9000_DEF]
#endif

		/* this works for hp9000s{300,400,700,800} */
#if defined(hp9000s800) && !defined(LocalDef)
#define LocalDef	machine_defs[HP9000_DEF]
#endif

                    /* This works for ibm 6000 */
#if defined(ibm) && !defined(LocalDef)
#define LocalDef	machine_defs[RISC6000_DEF]
#endif

                    /* This works for the cray */
#if defined(cray) && !defined(LocalDef)
#define LocalDef	machine_defs[CRAY_DEF]
#endif

               /* This works for the Convex */
#if defined(__convex__) && !defined(LocalDef)
#define LocalDef	machine_defs[CONVEX_DEF]
#endif

               /* This works for the Mac */
#if defined(AUX) && !defined(LocalDef)
#define LocalDef	machine_defs[MACII_DEF]
#endif

               /* This works for the Data General */
#if defined(DGUX) && !defined(LocalDef)
#define LocalDef	machine_defs[DGUX_DEF]
#endif

               /* This works for the Intel 386/486*/
#if defined(i386) && !defined(LocalDef)
#define LocalDef	machine_defs[I386_DEF]
#endif
               /* This works for the Motorola 88000 */
#if defined(m88k) && !defined(LocalDef)
#define LocalDef	machine_defs[M88K_DEF]
#endif

               /* This works for the Luna 88000 */
#if defined(luna) && !defined(LocalDef)
#define LocalDef	machine_defs[LUNA88K_DEF]
#endif

               /* unknown type */
#if !defined(LocalDef)
#define LocalDef	machine_defs[UNKNOWN]
#endif


/*
 *  Compute the number of current Machine definitions
 */
#define NUMBER_MACHINES ((int) sizeof(machine_defs)/sizeof(machine_defs[0]))


/*
 *  routine declarations
 */
char	**getmachlist();

#endif /* _machdefs_h_ */
/* do not add after this line */
