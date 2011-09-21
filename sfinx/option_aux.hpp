#pragma once
#include <cmath>

namespace sfinx {
namespace aux {

template <typename Num>
inline Num d1(Num S, Num K, Num T, Num r, Num v, Num q = 0.0)
{
  return (log(S / K) + (r - q + v * v / 2) * T) / (v * sqrt(T));
}

template <typename Num>
inline Num d2(Num d1, Num T, Num v)
{
  return d1 - v * sqrt(T);
}

template <typename Num>
inline Num d2(Num S, Num K, Num T, Num r, Num v, Num q = 0.0)
{
  return d2(d1(S, K, T, r, v, q), T, v);
}

template <typename Num>
inline std::pair<Num, Num> d(Num S, Num K, Num T, Num r, Num v, Num q = 0.0)
{
  Num d_1 = d1(S, K, T, r, v, q);
  return std::make_pair(d_1, d2(d_1, T, v));
}

} // namespace aux
} // namespace sfinx

