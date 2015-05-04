/*
    (c) Sergey Ryazanov (http://home.onego.ru/~ryazanov)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

For details, see "The Impossibly Fast C++ Delegates" at
http://www.codeproject.com/Articles/11015/The-Impossibly-Fast-C-Delegates
*/

#ifndef VIGRA_DELEGATE_INCLUDED
#define VIGRA_DELEGATE_INCLUDED

namespace vigra
{
#ifdef VIGRA_DELEGATE_PREFERRED_SYNTAX
    template <typename TSignature> class delegate;
    template <typename TSignature> class delegate_invoker;
#endif
}

#ifdef _MSC_VER
#define VIGRA_DELEGATE_CALLTYPE __fastcall
#else
#define VIGRA_DELEGATE_CALLTYPE
#endif

#include "detail/delegate_list.hxx"

#undef VIGRA_DELEGATE_CALLTYPE

#endif//VIGRA_DELEGATE_INCLUDED
