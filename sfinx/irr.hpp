#pragma once
#include <functional>
#include <limits>
#include "pv.hpp"
#include "solver.hpp"

namespace sfinx {

template <Flow F, typename T, typename U>
double irr(T const& times, U const& amounts)
{
  using namespace std::placeholders;
  auto range = std::make_pair(0.0, 1.0);
  double eps = 1.0e-5;
  auto f = [&](double r){ return pv<F>(times, amounts, r); };
  auto res = solver::bisection(f, range, eps);
  return res.second <= eps ? res.first : std::numeric_limits<double>::quiet_NaN();
}


} // namespace sfinx

