#include <array>
#include <tuple>
#include <gtest/gtest.h>
#include "pv.hpp"
#include "bisect.hpp"
#include "irr.hpp"
#include "bond.hpp"

using namespace sfinx;

TEST(sfinx, linear_interpolate)
{
  double eps = 1.0e-5;
  double xs[] = { 0, 1.0, 2.0 };
  double ys[] = { 0, 2.0, 4.0 };
  EXPECT_FALSE(linear_interpolate(-1.0, xs, ys).second);
  EXPECT_LT(std::fabs(linear_interpolate(0.0, xs, ys).first - 0.0), eps); 
  EXPECT_LT(std::fabs(linear_interpolate(0.5, xs, ys).first - 1.0), eps); 
  EXPECT_LT(std::fabs(linear_interpolate(1.7, xs, ys).first - 3.4), eps); 
  EXPECT_LT(std::fabs(linear_interpolate(2.0, xs, ys).first - 4.0), eps); 
  EXPECT_FALSE(linear_interpolate(2.1, xs, ys).second);
}

TEST(sfinx, pv)
{
  double eps = 1.0;
  {
    double times[] = { 0, 1.0, 2.0 };
    double amounts[] = { -100.0, 75.0, 75.0 };
    EXPECT_LT(std::fabs(pv<Flow::Discrete>(times, amounts, 0.1) - 30), eps);
    EXPECT_LT(std::fabs(pv<Flow::Continuous>(times, amounts, 0.1) - 29.3), eps) ;
  }
  {
    double times[] = { 1, 2, 3 };
    double amounts[] = { 10, 10, 110 };
    EXPECT_LT(std::fabs(pv<Flow::Discrete>(times, amounts, 0.09) - 102), eps);
  }
}

TEST(sfinx, bisect)
{
  double eps = 0.00001;
  auto f = [](double x) { return x * x - 4; };
  double root = bisect(f, std::array<double, 2>{ {0.0, 10.0} }, eps).first;
  EXPECT_LT(std::fabs(root - 2.0), eps);
  root = bisect(f, std::make_pair(1.0, 3.0), eps).first;
  EXPECT_LT(std::fabs(root - 2.0), eps);
  root = bisect(f, std::make_tuple(3.0, 1.0), eps).first;
  EXPECT_LT(std::fabs(root - 2.0), eps);

  // This range contains no root
  auto res = bisect(f, std::make_tuple(3.0, 10.0), eps);
  EXPECT_EQ(res.first, 3.0);
  EXPECT_EQ(res.second, 5.0);
  // This range contains no root
  res = bisect(f, std::make_tuple(10.0, 3.0), eps);
  EXPECT_EQ(res.first, 3.0);
  EXPECT_EQ(res.second, 5.0);
}

TEST(sfinx, irr)
{
  double eps = 1.0e-5;
  double times[] = { 0, 1.0, 2.0 };
  double amounts[] = { -100.0, 10.0, 110.0 };
  EXPECT_LT(std::fabs(irr<Flow::Discrete>(times, amounts) - 0.1), eps);
}

TEST(sfinx, bond_price)
{
  double eps = 1.0e-3;
  double times[] = { 1.0, 2.0, 3.0 };
  double amounts[] = { 10.0, 10.0, 110.0 };
  EXPECT_LT(std::fabs(bond_price<Flow::Discrete>(times, amounts, 0.09) - 102.531), eps);
  EXPECT_LT(std::fabs(bond_price<Flow::Continuous>(times, amounts, 0.09) - 101.464), eps);
}

TEST(sfinx, ytm)
{
  double eps = 0.01;
  double times[] = { 1.0, 2.0, 3.0 };
  double amounts[] = { 10.0, 10.0, 110.0 };
  EXPECT_LT(std::fabs(ytm<Flow::Discrete>(times, amounts, 102.531) - 0.09), eps);
}

TEST(sfinx, bond_duration)
{
  double eps = 1.0e-5;
  double times[] = { 1.0, 2.0, 3.0 };
  double amounts[] = { 10.0, 10.0, 110.0 };
  EXPECT_LT(std::fabs(bond_duration<Flow::Discrete>(times, amounts, 0.09) - 2.73895), eps);
  EXPECT_LT(std::fabs(bond_duration<Flow::Continuous>(times, amounts, 0.09) - 2.73753), eps);
}

TEST(sfinx, bond_macaulay_duration)
{
  double eps = 1.0e-5;
  double times[] = { 1.0, 2.0, 3.0 };
  double amounts[] = { 10.0, 10.0, 110.0 };
  EXPECT_LT(std::fabs(bond_macaulay_duration<Flow::Discrete>(times, amounts, 102.531) - 2.73895), eps);
}

TEST(sfinx, bond_modified_duration)
{
  double eps = 1.0e-4;
  double times[] = { 1.0, 2.0, 3.0 };
  double amounts[] = { 10.0, 10.0, 110.0 };
  EXPECT_LT(std::fabs(bond_modified_duration<Flow::Discrete>(times, amounts, 102.531) - 2.5128), eps);
}

TEST(sfinx, bond_convexity)
{
  double eps = 1.0e-5;
  double times[] = { 1.0, 2.0, 3.0 };
  double amounts[] = { 10.0, 10.0, 110.0 };
  EXPECT_LT(std::fabs(bond_convexity<Flow::Discrete>(times, amounts, 0.09) - 8.93248), eps);
  EXPECT_LT(std::fabs(bond_convexity<Flow::Continuous>(times, amounts, 0.09) - 7.86779), eps);
}


