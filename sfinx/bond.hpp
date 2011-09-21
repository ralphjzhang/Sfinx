#pragma once
#include <functional>
#include "discount_factor.hpp"
#include "pv.hpp"
#include "solver.hpp"
#include "math.hpp"

namespace sfinx {

template <Flow F, typename T, typename U, typename Ret>
Ret bond_price(T const& times, U const& amounts, Ret r)
{
  return pv<F>(times, amounts, r);
}

template <Flow F, typename T, typename U>
double ytm(T const& times, U const& amounts, double price)
{
  using namespace std::placeholders;
  auto range = std::make_pair(0.0, 1.0);
  double eps = 1.0e-5;
  auto f = [&](double r){ return price - pv<F>(times, amounts, r); };
  auto res = solver::bisection(f, range, eps);
  return res.second <= eps ? res.first : std::numeric_limits<double>::quiet_NaN();
}

template <Flow F, typename T, typename U>
auto bond_duration(T const& times, U const& amounts, double r)
  -> typename std::enable_if<F == Flow::Discrete, double>::type
{
  auto f = [r](double t, double c) { return t * c * discount_factor(r, t, 1); };
  double pv_time = inner_product(0.0, times, amounts, f);
  return pv_time / pv<F>(times, amounts, r);
}

template <Flow F, typename T, typename U, typename Ret>
auto bond_duration(T const& times, U const& amounts, Ret r)
  -> typename std::enable_if<F == Flow::Continuous, Ret>::type
{
  auto f = [r](Ret t, Ret c) { return t * c * discount_factor(r, t); };
  double pv_time = inner_product(0.0, times, amounts, f);
  return pv_time / pv<F>(times, amounts, r);
}

template <Flow F, typename T, typename U>
double bond_macaulay_duration(T const& times, U const& amounts, double price)
{
  double r = ytm<F>(times, amounts, price);
  return bond_duration<F>(times, amounts, r);
}

template <Flow F, typename T, typename U>
double bond_modified_duration(T const& times, U const& amounts, double price)
{
  double r = ytm<F>(times, amounts, price);
  return bond_duration<F>(times, amounts, r) / (1 + r);
}

template <Flow F, typename T, typename U>
auto bond_convexity(T const& times, U const& amounts, double r)
  -> typename std::enable_if<F == Flow::Discrete, double>::type
{
  auto f = [r](double t, double c) { return c * t * (t + 1) * discount_factor(r, t, 1); };
  double Cx = inner_product(0.0, times, amounts, f);
  double B = bond_price<F>(times, amounts, r);
  return (Cx / std::pow(1 + r, 2)) / B;
}

template <Flow F, typename T, typename U, typename Ret>
auto bond_convexity(T const& times, U const& amounts, Ret r)
  -> typename std::enable_if<F == Flow::Continuous, Ret>::type
{
  auto f = [r](Ret t, Ret c) { return c * t * t * discount_factor(r, t); };
  double C = inner_product(0.0, times, amounts, f);
  double B = bond_price<F>(times, amounts, r);
  return C / B;
}

} // namespace sfinx

