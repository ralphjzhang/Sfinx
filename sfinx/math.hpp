#pragma once
#include <cmath>
#include <utility>
#include <numeric>

namespace sfinx {

template <typename Ret, typename T, typename U, typename Prod, typename Accum>
Ret inner_product(Ret init, T const& t, U const& u, Prod prod, Accum accum)
{
  return std::inner_product(std::begin(t), std::end(t), std::begin(u), init, accum, prod);
}

template <typename Ret, typename T, typename U, typename Prod>
Ret inner_product(Ret init, T const& t, U const& u, Prod prod)
{
  return inner_product(init, t, u, prod, std::plus<Ret>());
}

template <typename Ret, typename X, typename Y>
auto linear_interpolate(Ret x, X const& xs, Y const& ys)
  -> std::pair<Ret, bool>
{
  using namespace std;

  auto x0 = begin(xs), y0 = begin(ys);
  if (x0 == end(xs) || y0 == end(ys))
    return make_pair(Ret(), false);
  if (x == *x0)
    return make_pair(*y0, true);
  auto x1 = x0, y1 = y0;
  while (++x1 != end(xs) && ++y1 != end(ys))
  {
    if (*x0 < x && x <= *x1)
      return make_pair(*y0 + (x - *x0) * (*y1 - *y0) / (*x1 - *x0), true);
    ++x0; ++y0;
  }
  return make_pair(Ret(), false);
}

} // namespace sfinx

