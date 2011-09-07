#pragma once
#include <cmath>

namespace sfinx {
namespace term_structure {

/**
 * Nelson-Siegel Model
 *                         /1 - exp(-t/lambda)\
 * r(t) = b0 + (b1 + b2) *(--------------------) + b2 * exp(-t/lambda)
 *                         \     t/lambda     /
 **/
template <typename Decimal>
Decimal nelson_siegel(Decimal t, Decimal b0, Decimal b1, Decimal b2, Decimal lambda)
{
  if (t == 0) return b0;
  Decimal t1 = t / lambda;
  return b0 + (b1 + b2) * ((1 - exp(-t1)) / t1) + b2 * exp(-t1);
}

/**
 * Nelson-Siegel-Svensson Model
 *                /1 - exp(-t/t1)\       /1 - exp(-t/t1)            \        /1 - exp(-t/t2)            \
 * r(t) = b0 + b*(----------------)+ b2*(--------------- - exp(-t/t1)) + b3*(--------------- - exp(-t/t2))
 *                \     t/t1     /       \     t/t1                 /        \     t/t2                 /
 **/
template <typename Decimal>
Decimal svensson(Decimal t, Decimal b0, Decimal b1, Decimal b2, Decimal b3, Decimal t1, Decimal t2)
{
  if (t == 0) return b0;
  Decimal x1 = t / t1, x2 = t / t2;
  return b0 + b1 * (1 - exp(-x1)) / x1 
            + b2 * (((1 - exp(-x1)) / x1) - exp(-x1))
            + b3 * (((1 - exp(-x2)) / x2) - exp(-x2));
}

} // namespace term_structure
} // namespace sfinx

