#pragma once
#include <cstddef>
#include <utility>
#include <cmath>

namespace sfinx { namespace solver {

template <typename F, typename Range, typename Num>
auto bisection(F f, Range const& r, Num eps, size_t maxIter = 1000)
  -> std::pair<Num, Num>
{
  Num a = std::get<0>(r), b = std::get<1>(r);
  Num fa = f(a), fb = f(b);
  if (std::abs(fa) < eps) return std::make_pair(a, fa);
  if (std::abs(fb) < eps) return std::make_pair(b, fb);
  if (fa * fb > 0)
    return abs(fa) < abs(fb) ? std::make_pair(a, fa) : std::make_pair(b, fb);
  else if (fb < 0)
    std::swap(a, b);

  while (maxIter--) {
    Num c = (a + b) / 2, fc = f(c);
    if (std::abs(fc) < eps) return std::make_pair(c, fc);
    if (f(a) * fc < 0.0)
      b = c;
    else
      a = c;
  }
  return std::make_pair((a + b) / 2, f((a + b) / 2));
}

template <typename F, typename DF, typename Num>
auto newton(F f, DF df, Num x0, Num eps, size_t maxIter = 1000)
  -> std::pair<Num, Num>
{
  while (maxIter--) {
    Num fx0 = f(x0);
    if (std::abs(fx0) < eps)
      return std::make_pair(x0, fx0);
    x0 -= fx0 / df(x0);
  }
  // Fail to converge within maxIter
  return std::make_pair(x0, f(x0));
}

/**************** can also add *****************
 * Secant, Broyden, Brent
 ***********************************************/

} } // namespace sfinx::solver


