/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2000 by Ullrich Koethe                  */
/*       Cognitive Systems Group, University of Hamburg, Germany        */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    You may use, modify, and distribute this software according       */
/*    to the terms stated in the LICENSE file included in               */
/*    the VIGRA distribution.                                           */
/*                                                                      */
/*    The VIGRA Website is                                              */
/*        http://kogs-www.informatik.uni-hamburg.de/~koethe/vigra/      */
/*    Please direct questions, bug reports, and contributions to        */
/*        koethe@informatik.uni-hamburg.de                              */
/*                                                                      */
/*  THIS SOFTWARE IS PROVIDED AS IS AND WITHOUT ANY EXPRESS OR          */
/*  IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED      */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. */
/*                                                                      */
/************************************************************************/
 
 
#ifndef VIGRA_BORDERTREATMENT_HXX
#define VIGRA_BORDERTREATMENT_HXX

/********************************************************/
/*                                                      */
/*                      BorderTreatmentMode             */
/*                                                      */
/********************************************************/

/** Choose between different border treatment modes.
    In the convolution algorithms, these modes apply to 
    all image pixels where the kernel does not completely fit inside 
    the image.
    
    Include-File:
    \URL[vigra/bordertreatment.hxx]{../include/vigra/bordertreatment.hxx}
*/   
enum BorderTreatmentMode 
{
      /** do not operate on a pixel where the kernel does not fit in the
          image
          @memo
      */
   BORDER_TREATMENT_AVOID, 
   
      /** clip kernel at image border (this is only useful if the
          kernel is >= 0 everywhere)
          @memo
      */
   BORDER_TREATMENT_CLIP, 
   
      /** repeat the last valid pixel
          @memo
      */
   BORDER_TREATMENT_REPEAT,
   
      /** reflect image at last line 
          @memo
      */
   BORDER_TREATMENT_REFLECT, 
   
      /** wrap image around (periodic boundary conditions)
          @memo
      */
   BORDER_TREATMENT_WRAP
};

#endif // VIGRA_BORDERTREATMENT_HXX
