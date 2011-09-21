#pragma once
#include <type_traits>
#include <cmath>
#include "math.hpp"
#include "option_aux.hpp"

namespace sfinx { namespace option {

enum class Type 
{
  Call, Put, Both
};

enum class Exercise
{
  European, American
};

} // namespace option

namespace black_scholes {

/**
 * Call option value 
 **/
template <typename Num>
inline Num call(Num d1, Num d2, Num S, Num K, Num T, Num r)
{
  return S * normal_cdf(d1) - K * exp(-r * T) * normal_cdf(d2);
}

/**
 * Call option value with payout q
 **/
template <typename Num>
inline Num call(Num d1, Num d2, Num S, Num K, Num T, Num r, Num q)
{
  return call(d1, d2, S * exp(-q * T), K, T, r);
}

/**
 * Put option value
 **/
template <typename Num>
inline Num put(Num d1, Num d2, Num S, Num K, Num T, Num r)
{
  return K * exp(-r * T) * normal_cdf(-d2) - S * normal_cdf(-d1);
}

/**
 * Put option value with payout q
 **/
template <typename Num>
inline Num put(Num d1, Num d2, Num S, Num K, Num T, Num r, Num q)
{
  return put(d1, d2, S * exp(-q * T), K, T, r);
}

template <option::Type type, typename Num>
auto value(Num S, Num K, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Both, std::pair<Num, Num>>::type
{
  auto d = aux::d(S, K, T, r, v);
  Num c = call(d.first, d.second, S, K, T, r), p = put(d.first, d.second, S, K, T, r);
  return std::make_pair(c, p);
}

template <option::Type type, typename Num>
auto value(Num S, Num K, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Call, Num>::type
{
  auto d = aux::d(S, K, T, r, v);
  return call(d.first, d.second, S, K, T, r);
}

template <option::Type type, typename Num>
auto value(Num S, Num K, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Put, Num>::type
{
  auto d = aux::d(S, K, T, r, v);
  return put(d.first, d.second, S, K, T, r);
}

template <option::Type type, typename Num>
auto delta(Num S, Num K, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Call, Num>::type
{
  return normal_cdf(aux::d1(S, K, T, r, v));
}

template <option::Type type, typename Num>
auto delta(Num S, Num K, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Put, Num>::type
{
  return normal_cdf(aux::d1(S, K, T, r, v)) - 1;
}

template <option::Type type, typename Num>
auto delta(Num S, Num K, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Both, std::pair<Num, Num>>::type
{
  Num d = delta<option::Type::Call>(S, K, T, r, v);
  return std::make_pair(d, d - 1);
}

template <typename Num>
inline Num gamma(Num S, Num K, Num T, Num r, Num v)
{
  return normal_pdf(aux::d1(S, K, T, r, v)) / (S * v * sqrt(T));
}

template <typename Num>
inline Num vega(Num S, Num K, Num T, Num r, Num v)
{
  return S * normal_pdf(aux::d1(S, K, T, r, v)) * sqrt(T);
}

template <option::Type type, typename Num>
auto theta(Num S, Num K, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Call, Num>::type
{
  Num d1 = aux::d1(S, K, T, r, v);
  Num d2 = aux::d2(d1, T, v);
  return -S * normal_pdf(d1) * v / (2 * sqrt(T)) - r * K * exp(-r * T) * normal_cdf(d2);
}

template <option::Type type, typename Num>
auto theta(Num S, Num K, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Put, Num>::type
{
  Num d1 = aux::d1(S, K, T, r, v);
  Num d2 = aux::d2(d1, T, v);
  return -S * normal_pdf(d1) * v / (2 * sqrt(T)) + r * K * exp(-r * T) * normal_cdf(-d2);
}

template <option::Type type, typename Num>
auto theta(Num S, Num K, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Both, std::pair<Num, Num>>::type
{
  Num d1 = aux::d1(S, K, T, r, v);
  Num d2 = aux::d2(d1, T, v);
  Num t1 = -S * normal_pdf(d1) * v / (2 * sqrt(T)), t2 = r * K * exp(-r * T);
  return std::make_pair(t1 - t2 * normal_cdf(d2), t1 + t2 * normal_cdf(-d2));
}

template <option::Type type, typename Num>
auto rho(Num S, Num K, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Call, Num>::type
{
  return K * T * exp(-r * T) * normal_cdf(aux::d2(S, K, T, r, v));
}

template <option::Type type, typename Num>
auto rho(Num S, Num K, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Put, Num>::type
{
  return -K * T * exp(-r * T) * normal_cdf(-aux::d2(S, K, T, r, v));
}

template <option::Type type, typename Num>
auto rho(Num S, Num K, Num T, Num r, Num v)
  -> typename std::enable_if<type == option::Type::Both, std::pair<Num, Num>>::type
{
  Num d2 = aux::d2(S, K, T, r, v);
  Num t = K * T * exp(-r * T);
  return std::make_pair(t * normal_cdf(d2), -t * normal_cdf(-d2));
}

} // namespace black_scholes
} // namespace sfinx

