#include <array>
#include <tuple>
#include <gtest/gtest.h>
#include "solver.hpp"

using namespace sfinx::solver;

TEST(solver, bisection)
{
  double eps = 0.00001;
  auto f = [](double x) { return x * x - 4; };
  double root = bisection(f, std::array<double, 2>{{ 0.0, 10.0 }}, eps).first;
  EXPECT_NEAR(root, 2.0, eps);
  root = bisection(f, std::make_pair(1.0, 3.0), eps).first;
  EXPECT_NEAR(root, 2.0, eps);
  root = bisection(f, std::make_tuple(3.0, 1.0), eps).first;
  EXPECT_NEAR(root, 2.0, eps);

  // This range contains no root
  auto res = bisection(f, std::make_tuple(3.0, 10.0), eps);
  EXPECT_EQ(res.first, 3.0);
  EXPECT_EQ(res.second, 5.0);
  // This range contains no root
  res = bisection(f, std::make_tuple(10.0, 3.0), eps);
  EXPECT_EQ(res.first, 3.0);
  EXPECT_EQ(res.second, 5.0);
}

TEST(solver, newton)
{
  double eps = 0.00001;
  auto f = [](double x) { return x * x - 4; }; // f(x) = x^2 - 4
  auto df = [](double x) { return 2 * x; };    // f'(x) = 2x
  double root = newton(f, df, 5.0, eps).first;
  EXPECT_NEAR(root, 2.0, eps);
}

