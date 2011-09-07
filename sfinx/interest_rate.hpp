#pragma once
#include <cmath>

namespace sfinx {
namespace interest_rate {

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

template <typename Decimal>
Decimal ho_lee()
{
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

