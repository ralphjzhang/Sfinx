#pragma once
#include "solver.hpp"
#include "black_scholes.hpp"

namespace sfinx { namespace baw {

namespace aux {

template <typename Num>
inline Num q2(Num N, Num M, Num K)
{
  return (-(N - 1) + sqrt((N - 1) * (N - 1) + 4 * M / K)) / 2;
}

template <typename Num>
inline Num q1(Num N, Num M, Num K)
{
  return (-(N - 1) - sqrt((N - 1) * (N - 1) + 4 * M / K)) / 2;
}

template <typename Num>
inline Num rhs(Num d1, Num Si, Num X, Num T, Num r, Num b, Num v, Num q2)
{
  auto d = std::make_pair(d1, d1 - v * sqrt(T));
  auto C = bsm_general::aux::call(d, Si, X, T, r, b);
  return C + (1 - exp((b - r) * T) * normal_cdf(d1)) * Si / q2;
}

template <typename Num>
inline Num hs(Num d1, Num Sj, Num X, Num T, Num r, Num b, Num v, Num q1)
{
  auto d = std::make_pair(d1, d1 - v * sqrt(T));
  auto P = bsm_general::aux::put(d, Sj, X, T, r, b);
  return P - (1 - exp((b - r) * T) * normal_cdf(-d1)) * Sj / q1;
}

template <typename Num>
inline Num bi(Num d1, Num T, Num r, Num b, Num v, Num q2)
{
  auto N = normal_cdf(d1), n = normal_pdf(d1);
  return exp((b - r) * T) * N * (1 - 1 / q2) + (1 - exp((b - r) * T) * n / (v * sqrt(T))) / q2;
}

template <typename Num>
inline Num bj(Num d1, Num T, Num r, Num b, Num v, Num q1)
{
  auto N = normal_cdf(-d1), n = normal_pdf(-d1);
  return -exp((b - r) * T) * N * (1 - 1 / q1) - (1 + exp((b - r) * T) * n / (v * sqrt(T))) / q1;
}

template <typename Num>
Num seed_Ss(Num N, Num M, Num X, Num T, Num b, Num v)
{
  Num q2u = q2(N, M, 1.0);
  Num su = X / (1 - 1 / q2u);
  Num h2 = -(b * T + 2 * v * sqrt(T)) * X / (su - X);
  return X + (su - X) * (1 - exp(h2));
}

template <typename Num>
Num seed_Sss(Num N, Num M, Num X, Num T, Num b, Num v)
{
  Num q1u = q1(N, M, 1.0);
  Num su = X / (1 - 1 / q1u);
  Num h1 = (b * T - 2 * v * sqrt(T)) * X / (X - su);
  return su + (X - su) * exp(h1); 
}

template <typename Num>
Num solve_Ss(Num X, Num T, Num r, Num b, Num v)
{
  Num N = 2 * b / (v * v);
  Num M = 2 * r / (v * v);
  Num Si = seed_Ss(N, M, X, T, b, v);

  Num K = 1 - exp(-r * T);
  Num d1 = bsm_general::d1(Si, X, T, b, v);
  Num q2_ = q2(N, M, K);
  Num LHS = Si - X;
  Num RHS = rhs(d1, Si, X, T, r, b, v, q2_);
  Num bi_ = bi(d1, T, r, b, v, q2_);
  Num eps = 1e-6;
  while (std::abs(LHS - RHS) / X > eps){
    Si = (X + RHS - bi_ * Si) / (1 - bi_);
    d1 = bsm_general::d1(Si, X, T, b, v);
    LHS = Si - X;
    RHS = rhs(d1, Si, X, T, r, b, v, q2_);
    bi_ = bi(d1, T, r, b, v, q2_);
  }
  return Si;
}

template <typename Num>
Num solve_Sss(Num X, Num T, Num r, Num b, Num v)
{
  Num N = 2 * b / (v * v);
  Num M = 2 * r / (v * v);
  Num Sj = seed_Sss(N, M, X, T, b, v);

  Num K = 1 - exp(-r * T);
  Num d1 = bsm_general::d1(Sj, X, T, b, v);
  Num q1_ = q1(N, M, K);
  Num VS = X - Sj;
  Num HS = hs(d1, Sj, X, T, r, b, v, q1_);
  Num bj_ = bj(d1, T, r, b, v, q1_);
  Num eps = 1e-6;
  while (std::abs(VS - HS) / X > eps) {
    Sj = (X - HS + bj_ * Sj) / (1 + bj_);
    d1 = bsm_general::d1(Sj, X, T, b, v);
    VS = X - Sj;
    HS = hs(d1, Sj, X, T, r, b, v, q1_);
    bj_ = bj(d1, T, r, b, v, q1_);
  }
  return Sj;
}

} // namespace sfinx::baw::aux

template <typename Num>
Num call(Num S, Num X, Num T, Num r, Num b, Num v)
{
  if (b >= r)
    return bsm_general::call(S, X, T, r, b, v);
  Num Ss = aux::solve_Ss(X, T, r, b, v);
  Num N = 2 * b / (v * v);
  Num M = 2 * r / (v * v);
  Num K = 1 - exp(-r * T);
  Num q2 = aux::q2(N, M, K);
  Num d1 = bsm_general::d1(Ss, X, T, b, v);
  Num A2 = (Ss / q2) * (1 - exp((b - r) * T) * normal_cdf(d1));
  if (S < Ss)
    return bsm_general::call(S, X, T, r, b, v) + A2 * pow(S / Ss, q2);
  else
    return S - X;
}

template <typename Num>
Num put(Num S, Num X, Num T, Num r, Num b, Num v)
{
  Num Sss = aux::solve_Sss(X, T, r, b, v);
  Num N = 2 * b / (v * v);
  Num M = 2 * r / (v * v);
  Num K = 1 - exp(-r * T);
  Num q1 = aux::q1(N, M, K);
  Num d1 = bsm_general::d1(Sss, X, T, b, v);
  Num A1 = -(Sss / q1) * (1 - exp((b - r) * T) * normal_cdf(-d1));
  if (S > Sss)
    return bsm_general::put(S, X, T, r, b, v) + A1 * pow(S / Sss, q1);
  else
    return X - S;
}

} } // namespace sfinx::baw

