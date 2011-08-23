#pragma once
#include <cmath>

namespace sfinx {

double discount_factor(double rate, double t)
{
  return exp(-rate * t);
}

double yield(double discount_factor, double t)
{
  return -log(discount_factor) / t;
}

/*
template <typename T>
T discount_factor(T r, T t)
{
  return exp(-r * t);
}

template <typename T>
T yield(T discount_factor, T t)
{
  return -log(discount_factor) / t;
}
*/


} // namespace sfinx

