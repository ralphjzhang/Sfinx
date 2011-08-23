#pragma once
#include <cmath>
#include <type_traits>
#include <utility>
#include "math.hpp"

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
  return sigma_func(0.0, times, amounts,
      [r](double t, double c) { return c / std::pow(1.0 + r, t); }
    );
}

template <Flow F, typename T, typename U>
auto pv(T const& times, U const& amounts, double r)
  -> typename std::enable_if<F == Flow::Continuous, double>::type
{
  return sigma_func(0.0, times, amounts,
      [r](double t, double c) { return c * std::exp(-r * t); }
    );
}

} // namespace sfinx
