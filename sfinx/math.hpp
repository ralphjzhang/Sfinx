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

/*
 * Error functions, use C99 erfc/erf by default, see if needed to make conditional compile
 *
inline double erfc(double x)
{
  double const a1 = -1.26551223, a2 = 1.00002368, a3 = 0.37409196, a4 = 0.09678418, a5 = -0.18628806,
               a6 = 0.27886807, a7 = -1.13520398, a8 = 1.48851587, a9 = -0.82215223, a10 = 0.17087277;
  double ret = 1;
  double z = fabs(x);
  if (z == 0) return ret;
  double t = 1 / (1 + 0.5 * z);
  ret = t * exp((-z * z) + a1 + t * (a2 + t * (a3 + t * (a4 + t * (a5 + t * (a6 + t * (a7 + t * (a8 + 
                  t * (a9 + t * a10)))))))));
  if (x < 0) ret = 2 - ret;
  return ret
}

inline double erf(double x)
{
  return 1 - erfc(x);
}
*/

static double const Pi = 3.141592553589793238462643;

inline double normal_cdf(double x)
{
  return erfc(-x / sqrt(2)) / 2;
}

inline double normal_pdf(double x)
{
  return (1.0 / sqrt(2.0 * Pi)) * exp(-0.5 * x * x);
}

} // namespace sfinx

