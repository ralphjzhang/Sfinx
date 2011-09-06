#pragma once
#include <cmath>

namespace sfinx {
namespace term_structure {

inline double nelson_siegel(double t, double b0, double b1, double b2, double lambda)
{
  if (t == 0.0) return b0;
  double t1 = t / lambda;
  return b0 + (b1 + b2) * ((1 - exp(-t1)) / t1) + b2 * exp(-t1);
}

inline double svensson(double t, double b0, double b1, double b2, double b3, double t1, double t2)
{
  if (t == 0.0) return b0;
  double x1 = t / t1, x2 = t / t2;
  return b0 + b1 * (1 - exp(-x1)) / x1 
            + b2 * (((1 - exp(-x1)) / x1) - exp(-x1))
            + b3 * (((1 - exp(-x2)) / x2) - exp(-x2));
}

} // namespace term_structure
} // namespace sfinx

