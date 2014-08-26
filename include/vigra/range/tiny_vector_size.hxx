#ifndef TINY_VECTOR_SIZE_HXX_
#define TINY_VECTOR_SIZE_HXX_

namespace vigra
{

template <class T>
struct TinyVectorSize;
template<template <class, int> class TinyVectorType, class T, int N>  
struct TinyVectorSize<TinyVectorType<T, N> >
{
    enum{ value = N };
};

}

#endif

