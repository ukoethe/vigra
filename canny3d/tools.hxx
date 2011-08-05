#ifndef TOOLS_HXX
#define TOOLS_HXX

#include "vigra/tinyvector.hxx"
#include "vigra/mathutil.hxx"

namespace vigra {
  
  /// component-wise max of one vector
  template <class V, int SIZE, class D1, class D2>
  inline
  V max(TinyVectorBase<V, SIZE, D1, D2> const & v);

  template <class V, class D1, class D2>
  inline
  V max(TinyVectorBase<V, 3, D1, D2> const & v) {
    return VIGRA_CSTD::max(v[2],VIGRA_CSTD::max(v[1],v[0]));
  }

  /// component-wise max of two vectors
  template <class V1, int SIZE, class D1, class D2, class V2, class D3, class D4>
  inline
  TinyVector<typename PromoteTraits<V1, V2>::Promote, SIZE>
  max(TinyVectorBase<V1, SIZE, D1, D2> const & r1,
      TinyVectorBase<V2, SIZE, D3, D4> const & r2);

  template <class V1, class D1, class D2, class V2, class D3, class D4>
  inline
  TinyVector<typename PromoteTraits<V1, V2>::Promote, 3>
  max(TinyVectorBase<V1, 3, D1, D2> const & r1,
      TinyVectorBase<V2, 3, D3, D4> const & r2)
  {
    typedef TinyVector<typename PromoteTraits<V1, V2>::Promote, 3> Ret;
    return Ret(VIGRA_CSTD::max(r1[0],r2[0]),
           VIGRA_CSTD::max(r1[1],r2[1]),
           VIGRA_CSTD::max(r1[2],r2[2]));
  }
} // namespace vigra
#endif
