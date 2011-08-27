#pragma once
#include <cmath>
#include <type_traits>

namespace sfinx {

template <typename T>
T discount_factor(T rate, T t, size_t n = 0)
{
  return (n == 0) ? exp(-rate * t) : pow(1 + rate / n, -(n * t));
}

double yield(double discount_factor, double t)
{
  return -log(discount_factor) / t;
}

/*
template <typename T>
T yield(T discount_factor, T t)
{
  return -log(discount_factor) / t;
}
*/


} // namespace sfinx

