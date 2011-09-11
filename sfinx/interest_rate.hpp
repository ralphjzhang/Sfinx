#pragma once
#include <cmath>

namespace sfinx {
namespace interest_rate {

/**
 * Vasicek model
 * d(rt) = a(b - rt)dt + sigma * dWt
 **/
template <typename Decimal>
Decimal vasicek(Decimal t, Decimal r0, Decimal a, Decimal b, Decimal sigma)
{
  Decimal A, B;
  Decimal sigma2 = sigma * sigma;
  Decimal a2 = a * a;
  if (a) {
    B = (1 - exp(-a * t)) / a;
    A = exp(((B - t) * (a2 * b - sigma2 / 2)) / a2 - ((sigma2 * B * B) / (4 * a)));
  } else {
    B = t;
    A = exp(sigma2 * pow(t, 3)) / 6;
  }
  return A * exp(-B * r0);
}

/**
 * Ho and Lee model
 * d(rt) = theta dt + delta dWt
 **/
template <typename Decimal>
Decimal ho_lee(Decimal t, Decimal delta, Decimal theta)
{
  return 1 / theta + (1 - theta) * pow(delta, t);
}

template <typename Decimal>
Decimal hull_white()
{
}

template <typename Decimal>
Decimal cox_ingersoll_ross()
{
}

} } // namespace sfinx::interest_rate

