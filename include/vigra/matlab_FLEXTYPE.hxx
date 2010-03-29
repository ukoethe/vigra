/************************************************************************/
/*                                                                      */
/*                    Copyright 2009 by Rahul Nair                      */
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

#include <string>
#ifndef MATLAB_FLEXTYPE
#define MATLAB_FLEXTYPE



#define NOFLEXIBLETYPE vigraMain(outputs, inputs);


#define FLEX_TYPE(inClass, pos, name)\
    std::string name_str(#name);\
    if ((1##pos##1) != 0 && ((pos < inputs.size()-1) || pos == (inputs.size()-1) && !inputs.options_.isValid()))\
    {\
        inClass = mxGetClassID(inputs[pos]);\
    }\
    else\
    {\
        if(!inputs.options_.isValid())\
            mexErrMsgTxt("Input at (Position: pos, Name in Struct: name) not set");\
        inClass = mxGetClassID(inputs.options_[#name]);\
    }\

#define DEFAULT_ERROR\
    default:\
    std::string msg = "Invalid Inputtype for data element '" + name_str + "' - see documentation for supported Types";\
    mexErrMsgTxt(msg.c_str());\


#define FLEXIBLE_TYPE_START(pos, name)\
    mxClassID inClass;\
    std::string name_str(#name);\
    if ((1##pos##1) != 0 && ((pos < inputs.size()-1) || pos == (inputs.size()-1) && !inputs.options_.isValid()))\
    {\
        inClass = mxGetClassID(inputs[pos]);\
    }\
    else\
    {\
        if(!inputs.options_.isValid())\
            mexErrMsgTxt("Input at (Position: pos, Name in Struct: name) not set");\
        inClass = mxGetClassID(inputs.options_[#name]);\
    }\
    switch(inClass){\



#define FLEXIBLE_TYPE_END\
    default:\
    std::string msg = "Invalid Inputtype for data element " + name_str + " - see documentation for supported Types";\
    mexErrMsgTxt(msg.c_str());\
    }


/*UINT*/

#define ALLOW_UINT\
        case mxUINT8_CLASS:\
            vigraMain<UInt8>(outputs, inputs);     break;\
        case mxUINT16_CLASS:\
            vigraMain<UInt16>(outputs, inputs);    break;\
        case mxUINT32_CLASS:\
            vigraMain<UInt32>(outputs, inputs);    break;\
        case mxUINT64_CLASS:\
            vigraMain<UInt64>(outputs, inputs);    break;


#define ALLOW_UINT_8\
        case mxUINT8_CLASS:\
            vigraMain<UInt8>(outputs, inputs);     break;

#define ALLOW_UINT_16\
        case mxUINT16_CLASS:\
            vigraMain<UInt16>(outputs, inputs);     break;

#define ALLOW_UINT_32\
        case mxUINT32_CLASS:\
            vigraMain<UInt32>(outputs, inputs);     break;

#define ALLOW_UINT_64\
        case mxUINT64_CLASS:\
            vigraMain<UInt64>(outputs, inputs);     break;

#define ALLOW_UINT_8_16\
        case mxUINT8_CLASS:\
            vigraMain<UInt8>(outputs, inputs);     break;\
        case mxUINT16_CLASS:\
            vigraMain<UInt16>(outputs, inputs);    break;

#define ALLOW_UINT_16_32\
        case mxUINT16_CLASS:\
            vigraMain<UInt16>(outputs, inputs);     break;\
        case mxUINT32_CLASS:\
            vigraMain<UInt32>(outputs, inputs);    break;

#define ALLOW_UINT_32_64\
        case mxUINT16_CLASS:\
            vigraMain<UInt32>(outputs, inputs);     break;\
        case mxUINT32_CLASS:\
            vigraMain<UInt64>(outputs, inputs);    break;

#define ALLOW_UINT_8_32\
        case mxUINT8_CLASS:\
            vigraMain<UInt8>(outputs, inputs);     break;\
        case mxUINT16_CLASS:\
            vigraMain<UInt16>(outputs, inputs);    break;\
        case mxUINT32_CLASS:\
            vigraMain<UInt32>(outputs, inputs);    break;

#define ALLOW_UINT_16_64\
        case mxUINT16_CLASS:\
            vigraMain<UInt16>(outputs, inputs);    break;\
        case mxUINT32_CLASS:\
            vigraMain<UInt32>(outputs, inputs);    break;\
        case mxUINT64_CLASS:\
            vigraMain<UInt64>(outputs, inputs);     break;


/*INT*/
#define ALLOW_INT\
        case mxINT8_CLASS:\
            vigraMain<Int8>(outputs, inputs);     break;\
        case mxINT16_CLASS:\
            vigraMain<Int16>(outputs, inputs);    break;\
        case mxINT32_CLASS:\
            vigraMain<Int32>(outputs, inputs);    break;\
        case mxINT64_CLASS:\
            vigraMain<Int64>(outputs, inputs);    break;

#define ALLOW_INT_8_64 ALLOW_INT;
#define ALLOW_UINT_8_64 ALLOW_UINT;

#define ALLOW_INT_8\
        case mxINT8_CLASS:\
            vigraMain<Int8>(outputs, inputs);     break;

#define ALLOW_INT_16\
        case mxINT16_CLASS:\
            vigraMain<Int16>(outputs, inputs);     break;

#define ALLOW_INT_32\
        case mxINT32_CLASS:\
            vigraMain<Int32>(outputs, inputs);     break;

#define ALLOW_INT_64\
        case mxINT64_CLASS:\
            vigraMain<Int64>(outputs, inputs);     break;

#define ALLOW_INT_8_16\
        case mxINT8_CLASS:\
            vigraMain<Int8>(outputs, inputs);     break;\
        case mxINT16_CLASS:\
            vigraMain<Int16>(outputs, inputs);    break;

#define ALLOW_INT_16_32\
        case mxINT16_CLASS:\
            vigraMain<Int16>(outputs, inputs);     break;\
        case mxINT32_CLASS:\
            vigraMain<Int32>(outputs, inputs);    break;

#define ALLOW_INT_32_64\
        case mxINT16_CLASS:\
            vigraMain<Int32>(outputs, inputs);     break;\
        case mxINT32_CLASS:\
            vigraMain<Int64>(outputs, inputs);    break;

#define ALLOW_INT_8_32\
        case mxINT8_CLASS:\
            vigraMain<Int8>(outputs, inputs);     break;\
        case mxINT16_CLASS:\
            vigraMain<Int16>(outputs, inputs);    break;\
        case mxINT32_CLASS:\
            vigraMain<Int32>(outputs, inputs);    break;

#define ALLOW_INT_16_64\
        case mxINT16_CLASS:\
            vigraMain<Int16>(outputs, inputs);    break;\
        case mxINT32_CLASS:\
            vigraMain<Int32>(outputs, inputs);    break;\
        case mxINT64_CLASS:\
            vigraMain<Int64>(outputs, inputs);     break;

/*Float double*/

#define ALLOW_FD\
        case mxDOUBLE_CLASS:\
            vigraMain<double>(outputs, inputs);    break;\
        case mxSINGLE_CLASS:\
            vigraMain<float>(outputs, inputs);     break;

#define ALLOW_F\
        case mxSINGLE_CLASS:\
            vigraMain<float>(outputs, inputs);     break;

#define ALLOW_D\
        case mxDOUBLE_CLASS:\
            vigraMain<double>(outputs, inputs);    break;

#endif
