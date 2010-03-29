#######################################################################
#                                                                      
#         Copyright 1998-2003 by Ullrich Koethe, Hans Meine            
#                                                                      
#    This file is part of the VIGRA computer vision library.           
#    The VIGRA Website is                                              
#        http://hci.iwr.uni-heidelberg.de/vigra/                       
#    Please direct questions, bug reports, and contributions to        
#        ullrich.koethe@iwr.uni-heidelberg.de    or                    
#        vigra@informatik.uni-hamburg.de                               
#                                                                      
#    Permission is hereby granted, free of charge, to any person       
#    obtaining a copy of this software and associated documentation    
#    files (the "Software"), to deal in the Software without           
#    restriction, including without limitation the rights to use,      
#    copy, modify, merge, publish, distribute, sublicense, and/or      
#    sell copies of the Software, and to permit persons to whom the    
#    Software is furnished to do so, subject to the following          
#    conditions:                                                       
#                                                                      
#    The above copyright notice and this permission notice shall be    
#    included in all copies or substantial portions of the             
#    Software.                                                         
#                                                                      
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND    
#    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   
#    OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          
#    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT       
#    HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,      
#    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      
#    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR     
#    OTHER DEALINGS IN THE SOFTWARE.                                   
#                                                                      
#######################################################################

import numpy

class Function(object):
    test_types = numpy.typecodes['AllInteger'][:-2] + numpy.typecodes['AllFloat']+'O'
    len_test_types = len(test_types)
    kindToNumber = {'b': 1, 'u': 2, 'i': 2, 'f': 3, 'c': 4}
    boolFunctions = ['equal', 'greater', 'greater_equal', 'less', 'less_equal', 'not_equal',
                     'logical_and', 'logical_not', 'logical_or', 'logical_xor']
    
    def __init__(self, function):
        self.function = function
        self.is_bool = function.__name__ in self.boolFunctions
        self.is_abs  = function.__name__ == "absolute"
        self.__doc__ = function.__doc__
        
    def __getattr__(self, name):
        return getattr(self.function, name)
        
    def __repr__(self):
        return "<vigra.ufunc '%s'>" % self.__name__
    
    def priorities(self, *args):
        '''Among the inputs with largest size, find the one with highest 
           __array_priority__. Return this input, or None if there is no 
           inputs with 'size' and '__array_priority__' defined.'''
        maxSize = max([getattr(x, 'size', 0) for x in args])
        if maxSize == 0:
            return None
        priorities = [(getattr(x, '__array_priority__', -1.0), x) for x in args if getattr(x, 'size', 0) == maxSize]
        priorities.sort(key = lambda (p, x): p)
        if priorities[-1][0] == -1.0:
            return None
        else:
            return priorities[-1][1]
    
    def common_type_numpy(self, *args):
        '''Find a common type for the given inputs.
           This function will become obsolete when numpy.find_common_type() will be fixed.
           Code taken from the partial fix in numpy changeset 7133.
        '''
        arrayTypes = [x.dtype for x in args if hasattr(x, 'dtype')]
        N = len(arrayTypes)
        if N == 1:
            highestArrayType = arrayTypes[0]
        else:
            k = 0
            while k < self.len_test_types:
                highestArrayType = numpy.dtype(self.test_types[k])
                numcoerce = len([x for x in arrayTypes if highestArrayType >= x])
                if numcoerce == N:
                    break
                k += 1
        scalarTypes = [numpy.dtype(type(x)) for x in args if numpy.isscalar(x)]
        N = len(scalarTypes)
        if N == 0 or highestArrayType >= scalarTypes[0]:
            return (highestArrayType, highestArrayType)
        else:
            h, s = highestArrayType.kind, scalarTypes[0].kind
            if (h in ['i', 'u'] and s in ['i', 'u']) or \
               (h == 'f' and s == 'f') or \
               (h == 'c' and s == 'c'):
                return (highestArrayType, highestArrayType)
            return (highestArrayType, scalarTypes[0])
        
    def common_type(self, *args):
        '''Find the appropriate pair (in_dtype, out_dtype) according to 
           vigranumpy typecasting rules. in_dtype is the type into which 
           the arguments will be casted before performing the operation
           (to prevent possible overflow), out_type is the type the output
           array will have (unless an explicit out-argument is provided). 
           
           The ideas behind the vigranumpy typcasting rules are (i) to represent
           data with at most 32 bit, when possible, (ii) to reduce the number of
           types that occur as results of mixed expressions, and (iii) to minimize 
           the chance of bad surprises. Default output types are thus determined 
           according to the following rules:
           
           1. The output type does not depend on the order of the arguments::
           
                 a + b results in the same type as b + a
           
           2. With exception of logical functions and abs(), the output type 
              does not depend on the function to be executed. The output type 
              of logical functions is bool. The output type of abs() follows
              general rules unless the input is complex, in which case the
              output type is the corresponding float type::
              
                 a + b results in the same type as a / b
                 a == b => bool
                 abs(complex128) => float64
                 
           3. If the inputs have the same type, the type is preserved::
           
                 uint8 + uint8 => uint8
           
           4. If (and only if) one of the inputs has at least 64 bits, the output 
              will also have at least 64 bits::
              
                 int64 + uint32 => int64
                 int64 + 1.0    => float64
                 
           5. If an array is combined with a scalar of the same kind (integer,
              float, or complex), the array type is preserved. If an integer 
              array with at most 32 bits is combined with a float scalar, the 
              result is float32 (and rule 4 kicks in if the array has 64 bits)::
              
                 uint8   + 1   => uint8
                 uint8   + 1.0 => float32
                 float32 + 1.0 => float32
                 float64 + 1.0 => float64
                 
           6. Integer expressions with mixed types always produce signed results.
              If the arguments have at most 32 bits, the result will be int32, 
              otherwise it will be int64 (cf. rule 4)::
              
                 int8  + uint8  => int32
                 int32 + uint8  => int32
                 int32 + uint32 => int32
                 int32 + int64  => int64
                 int64 + uint64 => int64
                 
           7. In all other cases, the output type is equal to the highest input 
              type::
              
                 int32   + float32    => float32
                 float32 + complex128 => complex128
                 
           8. All defaults can be overridden by providing an explicit output array::
           
                 ufunc.add(uint8, uint8, uint16) => uint16
                 
              In order to prevent overflow, necessary upcasting is performed before 
              the function is executed.
        '''
        if self.is_abs and args[0].dtype.kind == "c" and args[1] is None:
            dtype = args[0].dtype
            if dtype == numpy.complex64:
                return dtype, numpy.float32
            if dtype == numpy.complex128:
                return dtype, numpy.float64
            if dtype == numpy.clongdouble:
                return dtype, numpy.longdouble
        arrayTypes = [(self.kindToNumber[x.dtype.kind], x.dtype.itemsize, x.dtype) for x in args if hasattr(x, 'dtype')]
        arrayTypes.sort()
        if arrayTypes[0] != arrayTypes[-1] and arrayTypes[-1][0] == 2:
            if arrayTypes[-1][1] <= 4: 
                highestArrayType = (2, 4, numpy.int32)
            else: 
                highestArrayType = (2, 8, numpy.int64)
        else:
            highestArrayType = arrayTypes[-1]
            
        if self.is_bool:
            return (highestArrayType[-1], numpy.bool)

        scalarType = [numpy.dtype(type(x)) for x in args if numpy.isscalar(x)]
        if not scalarType:
            return (highestArrayType[-1], highestArrayType[-1])
        scalarType = (self.kindToNumber[scalarType[0].kind], scalarType[0].itemsize, scalarType[0])
        if highestArrayType[0] >= scalarType[0]:
            return (highestArrayType[-1], highestArrayType[-1])
        elif scalarType[0] == 3 and highestArrayType[1] <= 4:
            return (highestArrayType[-1], numpy.float32)        
        else:
            return (highestArrayType[-1], scalarType[-1])        
        
    def permutation(self, p):
        '''Find the axis permutation that makes p as close to C-order
           as possible. Return the permutation, its inverse, and the 
           permuted shape object.'''
        permutation  = [i for i, s in sorted(enumerate(p.strides), key = lambda (i, s): -s)]
        inversePermutation = [permutation.index(i) for i in range(p.ndim)]
        pshape       = [p.shape[j] for j in permutation]
        return permutation, inversePermutation, pshape

class UnaryFunction(Function):
    def __call__(self, a, out = None):
        if not isinstance(a, numpy.ndarray):
            return self.function(a, out)
        p = self.priorities(a, out)
        if p is None:
            return self.function(a, out)

        permutation, inversePermutation, pshape = self.permutation(p)
        dtype, out_dtype = self.common_type(a, out)

        a = a.transpose(permutation)
        a = numpy.require(a, dtype)
        o = out.transpose(permutation) if hasattr(out, 'transpose') else \
            numpy.ndarray(pshape, out_dtype, order='C').view(p.__class__)
        self.function(a, o)

        if out is None:
            out = o.transpose(inversePermutation) 
        return out

class UnaryFunctionOut2(Function):
    def __call__(self, a, out1 = None, out2 = None):
        if not isinstance(a, numpy.ndarray):
            return self.function(a, out1, out2)
        p = self.priorities(a, out1, out2)
        if p is None:
            return self.function(a, out1, out2)

        permutation, inversePermutation, pshape = self.permutation(p)
        dtype, out_dtype = self.common_type(a, out1, out2)

        a = a.transpose(permutation)
        a = numpy.require(a, dtype)
        o1 = out1.transpose(permutation) if hasattr(out1, 'transpose') else \
             numpy.ndarray(pshape, out_dtype, order='C').view(p.__class__)
        o2 = out2.transpose(permutation) if hasattr(out2, 'transpose') else \
             numpy.ndarray(pshape, out_dtype, order='C').view(p.__class__)
        self.function(a, o1, o2)

        if out1 is None:
            out1 = o1.transpose(inversePermutation) 
        if out2 is None:
            out2 = o2.transpose(inversePermutation) 
        return out1, out2

class BinaryFunction(Function):
    def __call__(self, a, b, out = None):
        a_isarray, b_isarray = isinstance(a, numpy.ndarray), isinstance(b, numpy.ndarray)
        if not a_isarray and not b_isarray:
            return self.function(a, b, out)
        p = self.priorities(a, b, out)
        if p is None:
            return self.function(a, b, out)

        ndim = p.ndim
        permutation, inversePermutation, pshape = self.permutation(p)
        dtype, out_dtype = self.common_type(a, b, out)

        if a_isarray:
            if a.ndim < ndim:
                a = a.reshape(((1,)*ndim + a.shape)[-ndim:])
            a = a.transpose(permutation)
        if b_isarray:
            if b.ndim < ndim:
                b = b.reshape(((1,)*ndim + b.shape)[-ndim:])
            b = b.transpose(permutation)

        # make sure that at least one input array has type dtype
        if a_isarray and b_isarray:
            if a.dtype != dtype and b.dtype != dtype:
                if b.size < a.size:
                    b = numpy.require(b, dtype)
                else:
                    a = numpy.require(a, dtype)
        elif a_isarray and a.dtype != dtype:
            a = numpy.require(a, dtype)
        elif b_isarray and b.dtype != dtype:
            b = numpy.require(b, dtype)

        o = out.transpose(permutation) if hasattr(out, 'transpose') else \
            numpy.ndarray(pshape, out_dtype, order='C').view(p.__class__)
        
        self.function(a, b, o)

        if out is None:
            out = o.transpose(inversePermutation) 
        return out

        
__all__ = []

for k in numpy.__dict__.itervalues():
     if type(k) == numpy.ufunc:
        if k.nin == 1 and k.nout == 1:
            exec k.__name__ + " = UnaryFunction(k)"
        if k.nin == 1 and k.nout == 2:
            exec k.__name__ + " = UnaryFunctionOut2(k)"
        if k.nin == 2:
            exec k.__name__ + " = BinaryFunction(k)"
        __all__.append(k.__name__)

__all__.sort()

__doc__ = '\nThe following mathematical functions are available in this module::\n\n'

for k in range(0, len(__all__), 7):
    __doc__ += '        ' + '   '.join(__all__[k:k+7]) + '\n'

__doc__ += '''
Some of these functions are also provided as member functions of the vigra array types::

        __abs__   __add__   __and__   __div__   __divmod__   __eq__   __floordiv__
        __ge__   __gt__   __invert__   __le__   __lshift__   __lt__   __mod__
        __mul__   __ne__   __neg__   __or__   __pos__   __pow__   __radd__
        __radd__   __rand__   __rdiv__   __rdivmod__   __rfloordiv__   __rlshift__
        __rmod__   __rmul__   __ror__   __rpow__   __rrshift__   __rshift__
        __rsub__   __rtruediv__   __rxor__   __sub__   __truediv__   __xor__

'''
