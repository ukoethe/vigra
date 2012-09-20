/************************************************************************/
/*                                                                      */
/*               Copyright 1998-2002 by Ullrich Koethe                  */
/*                                                                      */
/*    This file is part of the VIGRA computer vision library.           */
/*    The VIGRA Website is                                              */
/*        http://hci.iwr.uni-heidelberg.de/vigra/                       */
/*    Please direct questions, bug reports, and contributions to        */
/*        ullrich.koethe@iwr.uni-heidelberg.de    or                    */
/*        vigra@informatik.uni-hamburg.de                               */
/*                                                                      */
/*    Permission is hereby granted, free of charge, to any person       */
/*    obtaining a copy of this software and associated documentation    */
/*    files (the "Software"), to deal in the Software without           */
/*    restriction, including without limitation the rights to use,      */
/*    copy, modify, merge, publish, distribute, sublicense, and/or      */
/*    sell copies of the Software, and to permit persons to whom the    */
/*    Software is furnished to do so, subject to the following          */
/*    conditions:                                                       */
/*                                                                      */
/*    The above copyright notice and this permission notice shall be    */
/*    included in all copies or substantial portions of the             */
/*    Software.                                                         */
/*                                                                      */
/*    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    */
/*    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   */
/*    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          */
/*    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       */
/*    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      */
/*    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      */
/*    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     */
/*    OTHER DEALINGS IN THE SOFTWARE.                                   */                
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
    
    <b>\#include</b> \<vigra/bordertreatment.hxx\><br>
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

          // assume that all outside points have value zero
       BORDER_TREATMENT_ZEROPAD
    };
    \endcode
*/   
enum BorderTreatmentMode 
{
   BORDER_TREATMENT_AVOID, 
   BORDER_TREATMENT_CLIP, 
   BORDER_TREATMENT_REPEAT,
   BORDER_TREATMENT_REFLECT, 
   BORDER_TREATMENT_WRAP,
   BORDER_TREATMENT_ZEROPAD
};

} // namespace vigra

#endif // VIGRA_BORDERTREATMENT_HXX
