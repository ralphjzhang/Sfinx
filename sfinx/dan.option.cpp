#include <gtest/gtest.h>
#include "black_scholes.hpp"
#include "barone_adesi_whaley.hpp"
#include "bjerksund_stensland.hpp"


TEST(option, black_scholes)
{
  using namespace sfinx::option;
  using namespace sfinx::bs;

  double S = 45, X = 50, T = 0.50, r = 0.01, v = 0.20;
  double eps = 1.0e-4;
  EXPECT_NEAR(value<Type::Call>(S, X, T, r, v), 0.9392, eps);
  EXPECT_NEAR(value<Type::Put>(S, X, T, r, v), 5.6898, eps);
  EXPECT_NEAR(delta<Type::Call>(S, X, T, r, v), 0.2615, eps);
  EXPECT_NEAR(delta<Type::Put>(S, X, T, r, v), -0.7385, eps);
  EXPECT_NEAR(gamma(S, X, T, r, v), 0.0511, eps);

  EXPECT_NEAR(vega(S, X, T, r, v), 10.3504, eps);
  EXPECT_NEAR(theta<Type::Call>(S, X, T, r, v), -2.1783, eps);
  EXPECT_NEAR(theta<Type::Put>(S, X, T, r, v), -1.6809, eps);
  EXPECT_NEAR(rho<Type::Call>(S, X, T, r, v), 5.4126, eps);
  EXPECT_NEAR(rho<Type::Put>(S, X, T, r, v), -19.4628, eps);
}

TEST(option, black_scholes_merton_generalized)
{
  using namespace sfinx::bsm_general;
  double S = 75, X = 70, T = 0.50, r = 0.1, b = 0.05, v = 0.35;
  double eps = 1.0e-4;
  EXPECT_NEAR(put(S, X, T, r, b, v), 4.0870, eps);
}

TEST(option, barone_adesi_whaley)
{
  using namespace sfinx;
  double S = 90, X = 100, r = 0.1, b = 0;
  double eps = 1.0e-4;
  {
    double T = 0.1, v = 0.15;
    EXPECT_NEAR(baw::call(S, X, T, r, b, v), 0.0206, eps);
    v = 0.25;
    EXPECT_NEAR(baw::call(S, X, T, r, b, v), 0.3159, eps);
    v = 0.35;
    EXPECT_NEAR(baw::call(S, X, T, r, b, v), 0.9495, eps);
    T = 0.5, v = 0.15;
    EXPECT_NEAR(baw::call(S, X, T, r, b, v), 0.8208, eps);
    v = 0.25;
    EXPECT_NEAR(baw::call(S, X, T, r, b, v), 2.7437, eps);
    v = 0.35;
    EXPECT_NEAR(baw::call(S, X, T, r, b, v), 5.0062, eps);
  }
  {
    double T = 0.1, v = 0.15;
    EXPECT_NEAR(baw::put(S, X, T, r, b, v), 10.0000, eps);
    v = 0.25;
    EXPECT_NEAR(baw::put(S, X, T, r, b, v), 10.2531, eps); //10.2533
    v = 0.35; 
    EXPECT_NEAR(baw::put(S, X, T, r, b, v), 10.8786, eps);
    T = 0.5, v = 0.15;
    EXPECT_NEAR(baw::put(S, X, T, r, b, v), 10.5593, eps); //10.5595
    v = 0.25;
    EXPECT_NEAR(baw::put(S, X, T, r, b, v), 12.4417, eps); //10.4419
    v = 0.35;
    EXPECT_NEAR(baw::put(S, X, T, r, b, v), 14.6944, eps); //12.6945
  }
  S = 110;
  {
    double T = 0.1, v = 0.15;
    EXPECT_NEAR(baw::call(S, X, T, r, b, v), 10.0061, eps); //10.0089
    v = 0.25;
    EXPECT_NEAR(baw::call(S, X, T, r, b, v), 10.3902, eps); //10.3919 ???
    v = 0.35;
    EXPECT_NEAR(baw::call(S, X, T, r, b, v), 11.1679, eps);
    T = 0.5, v = 0.15;
    EXPECT_NEAR(baw::call(S, X, T, r, b, v), 10.8086, eps); //10.8087
    v = 0.25;
    EXPECT_NEAR(baw::call(S, X, T, r, b, v), 13.0168, eps); //13.0170
    v = 0.35;
    EXPECT_NEAR(baw::call(S, X, T, r, b, v), 15.5685, eps); //15.5689
  }
  {
    double T = 0.1, v = 0.15;
    EXPECT_NEAR(baw::put(S, X, T, r, b, v), 0.0410, eps);
    v = 0.25;
    EXPECT_NEAR(baw::put(S, X, T, r, b, v), 0.4562, eps);
    v = 0.35;
    EXPECT_NEAR(baw::put(S, X, T, r, b, v), 1.2402, eps);
    T = 0.5, v = 0.15;
    EXPECT_NEAR(baw::put(S, X, T, r, b, v), 1.0822, eps);
    v = 0.25;
    EXPECT_NEAR(baw::put(S, X, T, r, b, v), 3.3226, eps);
    v = 0.35;
    EXPECT_NEAR(baw::put(S, X, T, r, b, v), 5.8822, eps); //5.8823
  }
}

TEST(option, bjerksund_stensland)
{
  using namespace sfinx;
  double eps = 1.0e-4;
  EXPECT_NEAR(bs93::call(42., 40., 0.75, 0.04, -0.04, 0.35), 5.2704, eps);

  double S = 90, X = 100, r = 0.1, b = 0;
  {
    double T = 0.1, v = 0.15;
    EXPECT_NEAR(bs93::call(S, X, T, r, b, v), 0.0206, eps);
    v = 0.25;
    EXPECT_NEAR(bs93::call(S, X, T, r, b, v), 0.3159, eps);
    v = 0.35;
    EXPECT_NEAR(bs93::call(S, X, T, r, b, v), 0.9495, eps);
    T = 0.5, v = 0.15;
    EXPECT_NEAR(bs93::call(S, X, T, r, b, v), 0.8208, eps);
    v = 0.25;
    EXPECT_NEAR(bs93::call(S, X, T, r, b, v), 2.7436, eps);
    v = 0.35;
    EXPECT_NEAR(bs93::call(S, X, T, r, b, v), 5.0062, eps);
  }
}

