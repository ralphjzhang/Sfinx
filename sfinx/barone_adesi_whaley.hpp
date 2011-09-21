#pragma once
#include "solver.hpp"
#include "option_aux.hpp"
#include "black_scholes.hpp"

namespace sfinx {

namespace baw {

template <typename Num>
Num q2(Num K, Num r, Num v, Num b)
{
  Num v2 = v * v;
  Num M = 2 * r / v2, n = 2 * b / v2 - 1;
  return (-n + sqrt((n * n) + 4 * M / K)) / 2;
}

template <typename Num>
Num A2(Num Ss, Num T, Num b, Num r, Num q2, Num d1)
{
  return Ss * (1 - exp((b - r) * T) * normal_cdf(d1)) / q2;
}

template <typename Num>
Num call(Num S, Num X, Num T, Num r, Num v, Num b)
{
  Num K = 1 - exp(-r * T);
  auto g = [&](Num Si) { 
    Num d1 = aux::d1(Si, X, T, b, v);
    Num c = black_scholes::call(Si, X, T, r, v, b);
    Num q = q2(K, r, v, b);
    Num ret = Si - X - c - Si * (1 - exp((b - r) * T) * normal_cdf(d1)) / q;
    std::cout<<"Si="<<Si<<", c="<<c<<", d1="<<d1<<" ,q="<<q<<", g="<<ret;
    return ret;
  };
  auto dg = [&](Num Si) {
    Num d1 = aux::d1(Si, X, T, b, v);
    Num q = q2(K, r, v, b);
    Num ret = (1 - 1 / q) * (1 - exp((b - r) * T) * normal_cdf(d1)) + 
           (exp((b - r) * T) * normal_pdf(d1)) / (q * v * sqrt(T));
    std::cout << ", dg=" << ret<<std::endl;
    return ret;
  };
  Num seed = 114.482; // TODO use the seed from paper
  return solver::newton(g, dg, seed, 1e-5, 100).first;
}

template <typename Num>
Num put()
{
}

}

} // namespace sfinx

