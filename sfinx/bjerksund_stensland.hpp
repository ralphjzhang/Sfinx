#pragma once
#include "math.hpp"
#include "black_scholes.hpp"

namespace sfinx { namespace bs93 {

namespace aux {

template <typename Num>
Num phi(Num S, Num T, Num gamma, Num H, Num I, Num r, Num b, Num v)
{
  Num lambda = (-r + gamma * b + 0.5 * gamma * (gamma - 1) * v * v) * T;
  Num d = -(log(S / H) + (b + (gamma - 0.5) * v * v) * T) / (v * sqrt(T));
  Num kappa = 2 * b / (v * v) + (2 * gamma - 1);
  return exp(lambda) * pow(S, gamma)
      * (normal_cdf(d) - pow(I / S, kappa) * normal_cdf(d - 2 * log(I / S) / (v * sqrt(T))));
}

} // namespace sfinx::bs93::aux

template <typename Num>
Num call(Num S, Num X, Num T, Num r, Num b, Num v)
{
  if (b >= r)
    return bsm_general::call(S, X, T, r, b, v); // Not optimal to exercise early
  Num bv = b / (v * v) - 0.5;
  Num beta = -bv + sqrt(bv * bv + 2 * r / (v * v));
  Num B_inf = beta / (beta - 1) * X;
  Num B0 = std::max(X, r / (r - b) * X);
  Num hT = -(b * T + 2 * v * sqrt(T)) * B0 / (B_inf - B0);
  Num I = B0 + (B_inf - B0) * (1 - exp(hT));
  Num alpha = (I - X) * pow(I, -beta);
  if (S >= I)
    return S - X;
  else
    return alpha * pow(S, beta) - alpha * aux::phi(S, T, beta, I, I, r, b, v) 
        + aux::phi(S, T, 1.0, I, I, r, b, v) - aux::phi(S, T, 1.0, X, I, r, b, v)
        - X * aux::phi(S, T, 0.0, I, I, r, b, v) + X * aux::phi(S, T, 0.0, X, I, r, b, v);
}

} } // namespace sfinx::bs93

