﻿#######################################################################
#
#         Copyright 2009-2010 by Ullrich Koethe
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
import copy

vigraTypecastingRules = '''
Default output types are thus determined according to the following rules:

   1. The output type does not depend on the order of the arguments::

         a + b results in the same type as b + a

   2.a With exception of logical functions and abs(), the output type
       does not depend on the function to be executed.

   2.b The output type of logical functions is bool.

   2.c The output type of abs() follows general rules unless the
       input contains complex numbers, in which case the output type
       is the corresponding float number type::

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
        self.nin = function.nin
        self.nout = function.nout

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
        priorities = sorted(priorities,key = lambda tuplepx: tuplepx[0])
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

           See ufunc.vigraTypecastingRules for detailed information on coercion rules.
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
            return (highestArrayType[-1], numpy.bool8)

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

class UnaryFunction(Function):
    def __call__(self, arg, out=None):
        a = arg.squeeze().transposeToNumpyOrder()
        dtype, out_dtype = self.common_type(a, out)

        if out is None:
            out = arg.__class__(arg, dtype=out_dtype, order='A', init=False)
            o = out.squeeze().transposeToNumpyOrder()
        else:
            o = out.squeeze().transposeToNumpyOrder()
            if not a.axistags.compatible(o.axistags):
                raise RuntimeError("%s(): axistag mismatch" % self.function.__name__)

        a = numpy.require(a, dtype).view(numpy.ndarray) # view(ndarray) prevents infinite recursion
        self.function(a, o)
        return out

class UnaryFunctionOut2(Function):
    def __call__(self, arg, out1=None, out2=None):
        a = arg.squeeze().transposeToNumpyOrder()
        dtype, out_dtype = self.common_type(a, out1, out2)

        if out1 is None:
            out1 = arg.__class__(arg, dtype=out_dtype, order='A', init=False)
            o1 = out1.squeeze().transposeToNumpyOrder()
        else:
            o1 = out1.squeeze().transposeToNumpyOrder()
            if not a.axistags.compatible(o1.axistags):
                raise RuntimeError("%s(): axistag mismatch" % self.function.__name__)

        if out2 is None:
            out2 = arg.__class__(arg, dtype=out_dtype, order='A', init=False)
            o2 = out2.squeeze().transposeToNumpyOrder()
        else:
            o2 = out2.squeeze().transposeToNumpyOrder()
            if not a.axistags.compatible(o2.axistags):
                raise RuntimeError("%s(): axistag mismatch" % self.function.__name__)

        a = numpy.require(a, dtype).view(numpy.ndarray) # view(ndarray) prevents infinite recursion
        self.function(a, o1, o2)
        return out1, out2

class BinaryFunction(Function):
    def __call__(self, arg1, arg2, out=None):
        if arg1.__class__ is numpy.ndarray or arg2.__class__ is numpy.ndarray:
            return self.function(arg1, arg2, out)

        dtype, out_dtype = self.common_type(arg1, arg2, out)

        if isinstance(arg1, numpy.ndarray):
            a1 = arg1.transposeToNumpyOrder()
            if isinstance(arg2, numpy.ndarray):
                a2 = arg2.transposeToNumpyOrder()

                if arg1.__array_priority__ == arg2.__array_priority__:
                    priorityArg = arg2 if arg1.ndim < arg2.ndim else arg1
                else:
                    priorityArg = arg2 if arg1.__array_priority__ < arg2.__array_priority__ else arg1

                if a1.ndim < a2.ndim:
                    a1 = a1.insertChannelAxis(order='C')
                elif a1.ndim > a2.ndim:
                    a2 = a2.insertChannelAxis(order='C')

                axistags = a1.axistags

                if not axistags.compatible(a2.axistags):
                    raise RuntimeError("%s(): input axistag mismatch %r vs. %r" %
                                         (self.function.__name__, axistags, a2.axistags))
                shape = tuple(max(k) for k in zip(a1.shape, a2.shape))
                a2 = numpy.require(a2, dtype).view(numpy.ndarray)
            else:
                priorityArg = arg1
                axistags = a1.axistags
                shape = a1.shape
                a2 = arg2
            a1 = numpy.require(a1, dtype).view(numpy.ndarray) # view(ndarray) prevents infinite recursion
        else:
            a1 = arg1
            a2 = arg2.transposeToNumpyOrder()
            axistags = a2.axistags
            shape = a2.shape
            priorityArg = arg2
            a2 = numpy.require(a2, dtype).view(numpy.ndarray)

        if out is None:
            outClass = priorityArg.__class__
            inversePermutation = priorityArg.permutationFromNumpyOrder()
            o = outClass(shape, dtype=out_dtype, order='C', axistags=axistags, init=False)
            if priorityArg.ndim < o.ndim:
                out = o.dropChannelAxis().transpose(inversePermutation)
            else:
                out = o.transpose(inversePermutation)
        else:
            o = out.transposeToNumpyOrder()
            if o.ndim < len(shape):
                o = o.insertChannelAxis(order='C')
            if not axistags.compatible(o.axistags):
                raise RuntimeError("%s(): output axistag mismatch %r vs. %r" %
                                         (self.function.__name__, axistags, o.axistags))
        self.function(a1, a2, o)
        return out

__all__ = []

for _k in numpy.__dict__.values():
     if type(_k) == numpy.ufunc:
        if _k.nin == 1 and _k.nout == 1:
            exec(_k.__name__ + " = UnaryFunction(_k)")
        if _k.nin == 1 and _k.nout == 2:
            exec(_k.__name__ + " = UnaryFunctionOut2(_k)")
        if _k.nin == 2:
            exec(_k.__name__ + " = BinaryFunction(_k)")
        __all__.append(_k.__name__)

__all__ = sorted(__all__)

def _prepareDoc():
    doc = '''
The following mathematical functions are available in this module
(refer to numpy for detailed documentation)::

'''

    k = 0
    while k < len(__all__):
        t = 8
        while True:
            d = '    ' + '   '.join(__all__[k:k+t]) + '\n'
            if len(d) <= 80:
                break
            t -= 1
        doc += d
        k += t

    return doc + '''
Some of these functions are also provided as member functions of
VigraArray::

    __abs__   __add__   __and__   __div__   __divmod__   __eq__
    __floordiv__   __ge__   __gt__   __invert__   __le__   __lshift__
    __lt__   __mod__   __mul__   __ne__   __neg__   __or__   __pos__
    __pow__   __radd__   __radd__   __rand__   __rdiv__   __rdivmod__
    __rfloordiv__   __rlshift__   __rmod__   __rmul__   __ror__   __rpow__
    __rrshift__   __rshift__   __rsub__   __rtruediv__   __rxor__   __sub__
    __truediv__   __xor__

As usual, these functions are applied independently at each pixel.

Vigranumpy overloads the numpy-versions of these functions in order to make their
behavior more suitable for image analysis. In particular, we changed two aspects:

* Axistag consistency is checked, and the order of axes and strides is
  preserved in the result array. (In contrast, plain numpy functions
  always create C-order arrays, disregarding the stride order of the
  inputs.)
* Typecasting rules are changed such that (i) data are represented with
  at most 32 bits, when possible, (ii) the number of types that occur as
  results of mixed expressions is reduced, and (iii) the chance of bad
  surprises is minimized.

''' + vigraTypecastingRules

__doc__ = _prepareDoc()