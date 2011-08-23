#pragma once
#include <cstddef>
#include <utility>
#include <cmath>

namespace sfinx {

template <typename F, typename Range, typename T>
auto bisect(F f, Range const& r, T eps, size_t maxIter = 1000)
  -> std::pair<T, T>
{
  T a = std::get<0>(r), b = std::get<1>(r);
  T fa = f(a), fb = f(b);
  if (std::fabs(fa) < eps) return std::make_pair(a, fa);
  if (std::fabs(fb) < eps) return std::make_pair(b, fb);
  if (fa * fb > 0)
    return std::fabs(fa) < std::fabs(fb) ? std::make_pair(a, fa) : std::make_pair(b, fb);
  else if (fb < 0)
    std::swap(a, b);

  while (maxIter--) {
    T c = (a + b) / 2, fc = f(c);
    if (std::fabs(fc) < eps) return std::make_pair(c, fc);
    if (f(a) * fc < 0.0)
      b = c;
    else
      a = c;
  }
  return std::make_pair((a + b) / 2, f((a + b) / 2));
}

} // namespace sfinx

