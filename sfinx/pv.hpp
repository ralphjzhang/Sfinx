#pragma once
#include <cmath>
#include <type_traits>
#include <utility>
#include "math.hpp"
#include "discount_factor.hpp"

namespace sfinx {

enum class Flow : unsigned
{
  Discrete,
  Continuous
};

template <Flow F, typename T, typename U>
auto pv(T const& times, U const& amounts, double r)
  -> typename std::enable_if<F == Flow::Discrete, double>::type
{
  return inner_product(0.0, times, amounts,
      [r](double t, double c) { return c * discount_factor(r, t, 1); }
    );
}

template <Flow F, typename T, typename U, typename Ret>
auto pv(T const& times, U const& amounts, Ret r)
  -> typename std::enable_if<F == Flow::Continuous, Ret>::type
{
  return inner_product(0.0, times, amounts,
      [r](Ret t, Ret c) { return c * discount_factor(r, t); }
    );
}

} // namespace sfinx

