#pragma once
#include <type_traits>
#include <cmath>
#include "math.hpp"

namespace sfinx { namespace option {

namespace aux {

inline double d1(double S, double K, double T, double r, double v)
{
  return (log(S / K) + (r + v * v / 2) * T) / (v * sqrt(T));
}

inline double d2(double d1, double T, double v)
{
  return d1 - v * sqrt(T);
}

inline double d2(double S, double K, double T, double r, double v)
{
  return d2(d1(S, K, T, r, v), T, v);
}

inline std::pair<double, double> d(double S, double K, double T, double r, double v)
{
  double d_1 = d1(S, K, T, r, v);
  return std::make_pair(d_1, d2(d_1, T, v));
}

inline double call(double d1, double d2, double S, double K, double T, double r)
{
  return S * normal_cdf(d1) - K * exp(-r * T) * normal_cdf(d2);
}

inline double put(double d1, double d2, double S, double K, double T, double r)
{
  return K * exp(-r * T) * normal_cdf(-d2) - S * normal_cdf(-d1);
}


} // namespace aux

enum class Type 
{
  Call, Put, Both
};

enum class Exercise
{
  European, American
};

template <Type type = Type::Both>
auto black_scholes(double S, double K, double T, double r, double v)
  -> typename std::enable_if<type == Type::Both, std::pair<double, double>>::type
{
  auto d = aux::d(S, K, T, r, v);
  double call = aux::call(d.first, d.second, S, K, T, r), put = aux::put(d.first, d.second, S, K, T, r);
  return std::make_pair(call, put);
}

template <Type type>
auto black_scholes(double S, double K, double T, double r, double v)
  -> typename std::enable_if<type == Type::Call, double>::type
{
  auto d = aux::d(S, K, T, r, v);
  return aux::call(d.first, d.second, S, K, T, r);
}

template <Type type>
auto black_scholes(double S, double K, double T, double r, double v)
  -> typename std::enable_if<type == Type::Put, double>::type
{
  auto d = aux::d(S, K, T, r, v);
  return aux::put(d.first, d.second, S, K, T, r);
}

template <Type type>
auto delta(double S, double K, double T, double r, double v)
  -> typename std::enable_if<type == Type::Call, double>::type
{
  return normal_cdf(aux::d1(S, K, T, r, v));
}

template <Type type>
auto delta(double S, double K, double T, double r, double v)
  -> typename std::enable_if<type == Type::Put, double>::type
{
  return normal_cdf(aux::d1(S, K, T, r, v)) - 1;
}

template <Type type>
auto delta(double S, double K, double T, double r, double v)
  -> typename std::enable_if<type == Type::Both, std::pair<double, double>>::type
{
  double d = delta<Type::Call>(S, K, T, r, v);
  return std::make_pair(d, d - 1);
}

inline double gamma(double S, double K, double T, double r, double v)
{
  return normal_pdf(aux::d1(S, K, T, r, v)) / (S * v * sqrt(T));
}

inline double vega(double S, double K, double T, double r, double v)
{
  return S * normal_pdf(aux::d1(S, K, T, r, v)) * sqrt(T);
}

template <Type type>
auto theta(double S, double K, double T, double r, double v)
  -> typename std::enable_if<type == Type::Call, double>::type
{
  double d1 = aux::d1(S, K, T, r, v);
  double d2 = aux::d2(d1, T, v);
  return -S * normal_pdf(d1) * v / (2 * sqrt(T)) - r * K * exp(-r * T) * normal_cdf(d2);
}

template <Type type>
auto theta(double S, double K, double T, double r, double v)
  -> typename std::enable_if<type == Type::Put, double>::type
{
  double d1 = aux::d1(S, K, T, r, v);
  double d2 = aux::d2(d1, T, v);
  return -S * normal_pdf(d1) * v / (2 * sqrt(T)) + r * K * exp(-r * T) * normal_cdf(-d2);
}

template <Type type>
auto theta(double S, double K, double T, double r, double v)
  -> typename std::enable_if<type == Type::Both, std::pair<double, double>>::type
{
  double d1 = aux::d1(S, K, T, r, v);
  double d2 = aux::d2(d1, T, v);
  double t1 = -S * normal_pdf(d1) * v / (2 * sqrt(T)), t2 = r * K * exp(-r * T);
  return std::make_pair(t1 - t2 * normal_cdf(d2), t1 + t2 * normal_cdf(-d2));
}

template <Type type>
auto rho(double S, double K, double T, double r, double v)
  -> typename std::enable_if<type == Type::Call, double>::type
{
  return K * T * exp(-r * T) * normal_cdf(aux::d2(S, K, T, r, v));
}

template <Type type>
auto rho(double S, double K, double T, double r, double v)
  -> typename std::enable_if<type == Type::Put, double>::type
{
  return -K * T * exp(-r * T) * normal_cdf(-aux::d2(S, K, T, r, v));
}

template <Type type>
auto rho(double S, double K, double T, double r, double v)
  -> typename std::enable_if<type == Type::Both, std::pair<double, double>>::type
{
  double d2 = aux::d2(S, K, T, r, v);
  double t = K * T * exp(-r * T);
  return std::make_pair(t * normal_cdf(d2), -t * normal_cdf(-d2));
}

} // namespace option
} // namespace sfinx

