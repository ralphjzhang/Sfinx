#include <gtest/gtest.h>
#include "black_scholes.hpp"
#include "barone_adesi_whaley.hpp"
#include "bjerksund_stensland.hpp"


TEST(option, black_scholes)
{
  using namespace sfinx::option;
  using namespace sfinx::black_scholes;

  double S = 45, K = 50, T = 0.50, r = 0.01, v = 0.20;
  double eps = 1.0e-4;
  EXPECT_NEAR(value<Type::Call>(S, K, T, r, v), 0.9392, eps);
  EXPECT_NEAR(value<Type::Put>(S, K, T, r, v), 5.6898, eps);
  EXPECT_NEAR(delta<Type::Call>(S, K, T, r, v), 0.2615, eps);
  EXPECT_NEAR(delta<Type::Put>(S, K, T, r, v), -0.7385, eps);
  EXPECT_NEAR(gamma(S, K, T, r, v), 0.0511, eps);

  EXPECT_NEAR(vega(S, K, T, r, v), 10.3504, eps);
  EXPECT_NEAR(theta<Type::Call>(S, K, T, r, v), -2.1783, eps);
  EXPECT_NEAR(theta<Type::Put>(S, K, T, r, v), -1.6809, eps);
  EXPECT_NEAR(rho<Type::Call>(S, K, T, r, v), 5.4126, eps);
  EXPECT_NEAR(rho<Type::Put>(S, K, T, r, v), -19.4628, eps);
}

TEST(option, barone_adesi_whaley)
{
  using namespace sfinx::baw;
  double S = 100, K = 100, r = 0.08, v = 0.20, b = -0.04, T = 0.25;
  double eps = 1.0e-5;
  EXPECT_NEAR(call(S, K, T, r, v, b), 5.74339, eps);
}

