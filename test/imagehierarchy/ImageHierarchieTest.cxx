#include <iostream>
#include <unittest.h>
#include <functional>
#include <algorithm>
#include "vigra/imagehierarchie.hxx"



/*
*Die Konstruktoren von SingleBandImage sind protected, wie werden sie geprueft
*Was macht man mit virtual. 
*SelectBandImage beinhaltet Unterklassen sollte man sie auch testen, ich denke
*schon.  
*627 , 716 , 722 wo kommt das copyImage her?, 1185, 1192, 1208 
*ImageHandle wird wie Object in Java verwendet?
*Was sind Proxy
*
*
*/
class ImageHierarchieTest
{

}



class VariableBandsImageTest
{
    /*
    alle unsere Initialisierungen
    */
    void testWidth(){}
    void testHeight(){}
    void testSize(){}
    void testActualWidth(){}
    void testActualHeight(){}
    void testActualSize(){}
    void testBands(){}
    void testActualBands(){}
    void testRoiUpperLeft(){}
    
    void testOperatorIndex2D(){}    
    void testOperatorIndex2DConst(){}
    void testOperatorFunctionCallInt(){}
    void testOperatorFunctionCallIntConst(){}
    void testLowerRight(){}
    void testConstLowerRight(){}
    void testUpperLeft(){}
    void testConstUpperLeft(){}
    void testAccessor(){}
    void testAccessorConst(){}
    
    /*
    und alles andere
    */
    
}


class FixedBandsImageTest
: VariableBandsImageTest
{
    /*
    *die Funktionen muessen ueberladen werden, da sie andere
    *Rueckgabewerte liefern oder auch nicht??
    */
    
    
    /*
    *hier noch die ganzen virtuellen sachen
    */
    void testCopyConstructor()
    void testImage2DConstuctor()
    void testImageIntConstuctor()
    void testImageInnerConstuctor()
    
    void testOperatorAssignmentImage(){}
    
    void testOperatorIndex2D(){}    
    void testOperatorIndex2DConst(){}
    void testOperatorFunctionCallInt(){}
    void testOperatorFunctionCallIntConst(){}
    
    void testLowerRight(){}
    void testConstLowerRight(){}
    void testUpperLeft(){}
    void testConstUpperLeft(){}
    void testAccessor(){}
    void testAccessorConst(){}

}

class FixedRGBImageTest                                   
: public FixedBandsImageTest
{    
    /*
    *hier noch die ganzen virtuellen sachen
    */
    void testCopyConstructor(){}
    void testImage2DConstuctor(){}
    void testImageIntConstuctor(){}
    void testImageInnerConstuctor(){}

    void testOperatorAssignmentImage(){}

    void testAccessor(){}
    void testAccessorConst(){}  
}








