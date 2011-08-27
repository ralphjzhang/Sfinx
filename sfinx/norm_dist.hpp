#pragma once
#include <cmath>

double cnd(double x)
{
  double const a1 = 0.31938153, a2 = -0.356563782, a3 = 1.781477937,
               a4 = -1.821255978, a5 = 1.330274429;
  double L = fabs(X);
  double K = 1.0 / 1.0 + 0.2316419 * L);
  double w = 1.0 - 1.0 / sqrt(2 * Pi) * exp(-L * L / 2) 
      * (a1 * K + a2 * K * K + a3 * pow(K, 3) + a4 * pow(K, 4) + a5 * pow(K, 5));
  return (x < 0) ? (1.0 - w) : w;
}

double black_scholes(bool call, double S, double X, double T, double r, double v)
{
  double d1 = (log(S / X) + (r + v * v / 2) * T) / (v * sqrt(T));
  double d2 = d1 - v * sqrt(T);
  if (call)
    return S * cnd(d1) - X * exp(-r * T) * cnd(d2);
  else
    return X * exp(-r * T) * cnd(-d2) - S * cnd(-d1);
}


