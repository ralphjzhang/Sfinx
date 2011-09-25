#pragma once
#include <type_traits>
#include <cmath>
#include "math.hpp"

namespace sfinx { namespace option {

enum class Type 
{
  Call, Put, Both
};

enum class Exercise
{
  European, American
};

} // namespace sfinx::option

/// Generalized Black Scholes Merton
//
namespace bsm_general {

namespace aux {
template <typename Num>
inline Num call(std::pair<Num, Num> d, Num S, Num X, Num T, Num r, Num b)
{
  return S * exp((b - r) * T) * normal_cdf(d.first) - X * exp(-r * T) * normal_cdf(d.second);
}

template <typename Num>
inline Num put(std::pair<Num, Num> d, Num S, Num X, Num T, Num r, Num b)
{
  return X * exp(-r * T) * normal_cdf(-d.second) - S * exp((b - r) * T) * normal_cdf(-d.first);
}
} // namespace sfinx::bsm_general::aux

template <typename Num>
inline Num d1(Num S, Num X, Num T, Num b, Num v)
{
  return (log(S / X) + (b + v * v / 2) * T) / (v * sqrt(T));
}

template <typename Num>
inline Num d2(Num d1, Num T, Num v)
{
  return d1 - v * sqrt(T);
}

template <typename Num>
inline std::pair<Num, Num> d(Num S, Num X, Num T, Num r, Num b, Num v)
{
  auto d1_ = d1(S, X, T, b, v);
  return std::make_pair(d1_, d2(d1_, T, v));
}

/**
 * Call option 
 **/
template <typename Num>
inline Num call(Num S, Num X, Num T, Num r, Num b, Num v)
{
  return aux::call(d(S, X, T, r, b, v), S, X, T, r, b);
}

/**
 * Put option
 **/
template <typename Num>
inline Num put(Num S, Num X, Num T, Num r, Num b, Num v)
{
  return aux::put(d(S, X, T, r, b, v), S, X, T, r, b);
}

} // namespace bsm_general


/// Standard Black-Scholes
//
namespace bs {

namespace aux {
template <typename Num>
inline Num call(std::pair<Num, Num> d, Num S, Num X, Num T, Num r)
{
  return S * normal_cdf(d.first) - X * exp(-r * T) * normal_cdf(d.second);
}

template <typename Num>
inline Num put(std::pair<Num, Num> d, Num S, Num X, Num T, Num r)
{
  return X * exp(-r * T) * normal_cdf(-d.second) - S * normal_cdf(-d.first);
}

} // namespace sfinx::bs::aux

template <typename Num>
inline Num d1(Num S, Num X, Num T, Num r, Num v)
{
  return bsm_general::d1(S, X, T, r, v);
}

template <typename Num>
inline Num d2(Num d1, Num T, Num v)
{
  return bsm_general::d2(d1, T, v);
}

template <typename Num>
inline Num d2(Num S, Num X, Num T, Num r, Num v)
{
  return d1(S, X, T, r, v) - v * sqrt(T);
}


template <typename Num>
inline std::pair<Num, Num> d(Num S, Num X, Num T, Num r, Num v)
{
  auto d1_ = d1(S, X, T, r, v);
  return std::make_pair(d1_, d2(d1_, T, v));
}

/**
 * Call option value 
 **/
template <typename Num>
Num call(Num S, Num X, Num T, Num r, Num v)
{
  return aux::call(d(S, X, T, r, v), S, X, T, r);
}

/**
 * Put option value
 **/
template <typename Num>
Num put(Num S, Num X, Num T, Num r, Num v)
{
  return aux::put(d(S, X, T, r, v), S, X, T, r);
}

template <option::Type type, typename Num>
auto value(Num S, Num X, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Both, std::pair<Num, Num>>::type
{
  auto d = bs::d(S, X, T, r, v);
  return std::make_pair(aux::call(d, S, X, T, r), aux::put(d, S, X, T, r));
}

template <option::Type type, typename Num>
auto value(Num S, Num X, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Call, Num>::type
{
  return call(S, X, T, r, v);
}

template <option::Type type, typename Num>
auto value(Num S, Num X, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Put, Num>::type
{
  return put(S, X, T, r, v);
}

/**
 * Greeks
 * */
template <option::Type type, typename Num>
auto delta(Num S, Num X, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Call, Num>::type
{
  return normal_cdf(bs::d1(S, X, T, r, v));
}

template <option::Type type, typename Num>
auto delta(Num S, Num X, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Put, Num>::type
{
  return normal_cdf(bs::d1(S, X, T, r, v)) - 1;
}

template <option::Type type, typename Num>
auto delta(Num S, Num X, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Both, std::pair<Num, Num>>::type
{
  Num d = delta<option::Type::Call>(S, X, T, r, v);
  return std::make_pair(d, d - 1);
}

template <typename Num>
inline Num gamma(Num S, Num X, Num T, Num r, Num v)
{
  return normal_pdf(bs::d1(S, X, T, r, v)) / (S * v * sqrt(T));
}

template <typename Num>
inline Num vega(Num S, Num X, Num T, Num r, Num v)
{
  return S * normal_pdf(bs::d1(S, X, T, r, v)) * sqrt(T);
}

template <option::Type type, typename Num>
auto theta(Num S, Num X, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Call, Num>::type
{
  Num d1 = bs::d1(S, X, T, r, v);
  Num d2 = bs::d2(d1, T, v);
  return -S * normal_pdf(d1) * v / (2 * sqrt(T)) - r * X * exp(-r * T) * normal_cdf(d2);
}

template <option::Type type, typename Num>
auto theta(Num S, Num X, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Put, Num>::type
{
  Num d1 = bs::d1(S, X, T, r, v);
  Num d2 = bs::d2(d1, T, v);
  return -S * normal_pdf(d1) * v / (2 * sqrt(T)) + r * X * exp(-r * T) * normal_cdf(-d2);
}

template <option::Type type, typename Num>
auto theta(Num S, Num X, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Both, std::pair<Num, Num>>::type
{
  Num d1 = bs::d1(S, X, T, r, v);
  Num d2 = bs::d2(d1, T, v);
  Num t1 = -S * normal_pdf(d1) * v / (2 * sqrt(T)), t2 = r * X * exp(-r * T);
  return std::make_pair(t1 - t2 * normal_cdf(d2), t1 + t2 * normal_cdf(-d2));
}

template <option::Type type, typename Num>
auto rho(Num S, Num X, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Call, Num>::type
{
  return X * T * exp(-r * T) * normal_cdf(bs::d2(S, X, T, r, v));
}

template <option::Type type, typename Num>
auto rho(Num S, Num X, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Put, Num>::type
{
  return -X * T * exp(-r * T) * normal_cdf(-bs::d2(S, X, T, r, v));
}

template <option::Type type, typename Num>
auto rho(Num S, Num X, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Both, std::pair<Num, Num>>::type
{
  Num d2 = bs::d2(S, X, T, r, v);
  Num t = X * T * exp(-r * T);
  return std::make_pair(t * normal_cdf(d2), -t * normal_cdf(-d2));
}

} // namespace bs

} // namespace sfinx

