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

namespace vigra {


/*! \page BorderTreatmentMode BorderTreatmentMode

    Choose between different border treatment modes. In the convolution 
    algorithms, these modes apply to 
    all image pixels where the kernel does not completely fit inside 
    the image.
    
    <b>\#include</b> "<a href="bordertreatment_8hxx-source.html">vigra/bordertreatment.hxx</a>"<br>
    Namespace: vigra
    
    \code
    enum BorderTreatmentMode 
    {
          // do not operate on a pixel where the kernel does 
          // not fit in the image
       BORDER_TREATMENT_AVOID, 

          // clip kernel at image border (this is only useful if the
          //  kernel is >= 0 everywhere)
       BORDER_TREATMENT_CLIP, 

          // repeat the nearest valid pixel
       BORDER_TREATMENT_REPEAT,

          // reflect image at last row/column 
       BORDER_TREATMENT_REFLECT, 

          // wrap image around (periodic boundary conditions)
       BORDER_TREATMENT_WRAP
    };
    \endcode
*/   
enum BorderTreatmentMode 
{
   BORDER_TREATMENT_AVOID, 
   BORDER_TREATMENT_CLIP, 
   BORDER_TREATMENT_REPEAT,
   BORDER_TREATMENT_REFLECT, 
   BORDER_TREATMENT_WRAP
};

} // namespace vigra

#endif // VIGRA_BORDERTREATMENT_HXX
