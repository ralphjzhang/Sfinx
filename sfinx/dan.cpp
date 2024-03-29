#include <cmath>
#include <array>
#include <tuple>
#include <gtest/gtest.h>
#include "discount_factor.hpp"
#include "math.hpp"
#include "pv.hpp"
#include "irr.hpp"
#include "bond.hpp"
#include "black_scholes.hpp"
#include "term_structure.hpp"
#include "interest_rate.hpp"


using namespace sfinx;

TEST(sfinx, linear_interpolate)
{
  double eps = 1.0e-5;
  double xs[] = { 0, 1.0, 2.0 };
  double ys[] = { 0, 2.0, 4.0 };
  EXPECT_FALSE(linear_interpolate(-1.0, xs, ys).second);
  EXPECT_LT(fabs(linear_interpolate(0.0, xs, ys).first - 0.0), eps); 
  EXPECT_LT(fabs(linear_interpolate(0.5, xs, ys).first - 1.0), eps); 
  EXPECT_LT(fabs(linear_interpolate(1.7, xs, ys).first - 3.4), eps); 
  EXPECT_LT(fabs(linear_interpolate(2.0, xs, ys).first - 4.0), eps); 
  EXPECT_FALSE(linear_interpolate(2.1, xs, ys).second);
}

TEST(sfinx, erfc)
{
  double eps = 1.0e-15;
  EXPECT_LT(fabs(erfc(-1.0) - 1.8427007929497148), eps);
  EXPECT_LT(fabs(erfc(0) - 1), eps);
  EXPECT_LT(fabs(erfc(2) - 0.004677734981047265), eps);
}

TEST(sfinx, normal_pdf)
{
  double eps = 1.0e-6;
  EXPECT_LT(fabs(normal_pdf(0) - 0.398942), eps);
  EXPECT_LT(fabs(normal_pdf(1) - 0.241971), eps);
  EXPECT_LT(fabs(normal_pdf(2) - 0.053991), eps);
  EXPECT_LT(fabs(normal_pdf(-0.5) - 0.352065), eps);
  EXPECT_LT(fabs(normal_pdf(-0.1) - 0.396953), eps);
}

TEST(sfinx, normal_cdf)
{
  double eps = 1.0e-4;
  EXPECT_LT(fabs(normal_cdf(-1) - 0.1587), eps);
  EXPECT_LT(fabs(normal_cdf(0) - 0.5), eps);
  EXPECT_LT(fabs(normal_cdf(2) - 0.9772), eps);
}

TEST(sfinx, discount_factor)
{
  double eps = 1.0e-6;
  EXPECT_LT(fabs(discount_factor(0.05, 1.0) - 0.951229), eps);
  EXPECT_LT(fabs(discount_factor(0.05, 1.0, 1) - 0.952381), eps);
  EXPECT_LT(fabs(yield(0.9, 2.0) - 0.0526803), eps);
}

TEST(sfinx, pv)
{
  double eps = 1.0;
  {
    double times[] = { 0, 1.0, 2.0 };
    double amounts[] = { -100.0, 75.0, 75.0 };
    EXPECT_LT(fabs(pv<Flow::Discrete>(times, amounts, 0.1) - 30), eps);
    EXPECT_LT(fabs(pv<Flow::Continuous>(times, amounts, 0.1) - 29.3), eps) ;
  }
  {
    double times[] = { 1, 2, 3 };
    double amounts[] = { 10, 10, 110 };
    EXPECT_LT(fabs(pv<Flow::Discrete>(times, amounts, 0.09) - 102), eps);
  }
}

TEST(sfinx, irr)
{
  double eps = 1.0e-5;
  double times[] = { 0, 1.0, 2.0 };
  double amounts[] = { -100.0, 10.0, 110.0 };
  EXPECT_LT(fabs(irr<Flow::Discrete>(times, amounts) - 0.1), eps);
}

TEST(sfinx, bond_price)
{
  double eps = 1.0e-3;
  double times[] = { 1.0, 2.0, 3.0 };
  double amounts[] = { 10.0, 10.0, 110.0 };
  EXPECT_LT(fabs(bond_price<Flow::Discrete>(times, amounts, 0.09) - 102.531), eps);
  EXPECT_LT(fabs(bond_price<Flow::Continuous>(times, amounts, 0.09) - 101.464), eps);
}

TEST(sfinx, ytm)
{
  double eps = 0.01;
  double times[] = { 1.0, 2.0, 3.0 };
  double amounts[] = { 10.0, 10.0, 110.0 };
  EXPECT_LT(fabs(ytm<Flow::Discrete>(times, amounts, 102.531) - 0.09), eps);
}

TEST(sfinx, bond_duration)
{
  double eps = 1.0e-5;
  double times[] = { 1.0, 2.0, 3.0 };
  double amounts[] = { 10.0, 10.0, 110.0 };
  EXPECT_LT(fabs(bond_duration<Flow::Discrete>(times, amounts, 0.09) - 2.73895), eps);
  EXPECT_LT(fabs(bond_duration<Flow::Continuous>(times, amounts, 0.09) - 2.73753), eps);
}

TEST(sfinx, bond_macaulay_duration)
{
  double eps = 1.0e-5;
  double times[] = { 1.0, 2.0, 3.0 };
  double amounts[] = { 10.0, 10.0, 110.0 };
  EXPECT_LT(fabs(bond_macaulay_duration<Flow::Discrete>(times, amounts, 102.531) - 2.73895), eps);
}

TEST(sfinx, bond_modified_duration)
{
  double eps = 1.0e-4;
  double times[] = { 1.0, 2.0, 3.0 };
  double amounts[] = { 10.0, 10.0, 110.0 };
  EXPECT_LT(fabs(bond_modified_duration<Flow::Discrete>(times, amounts, 102.531) - 2.5128), eps);
}

TEST(sfinx, bond_convexity)
{
  double eps = 1.0e-5;
  double times[] = { 1.0, 2.0, 3.0 };
  double amounts[] = { 10.0, 10.0, 110.0 };
  EXPECT_LT(fabs(bond_convexity<Flow::Discrete>(times, amounts, 0.09) - 8.93248), eps);
  EXPECT_LT(fabs(bond_convexity<Flow::Continuous>(times, amounts, 0.09) - 7.86779), eps);
}

TEST(sfinx, nelson_siegel)
{
  using namespace sfinx::term_structure;
  double eps = 1.0e-7;
  double b0 = 0.01, b1 = 0.01, b2 = 0.01, lambda = 5.0, t = 1;
  EXPECT_LT(fabs(nelson_siegel(t, b0, b1, b2, lambda) - 0.0363142), eps);
}

TEST(sfinx, vasicek)
{
  using namespace sfinx::interest_rate;
  double eps = 1.0e-6;
  double t = 1, r0 = 0.05, a = -0.1, b = 0.1, sigma = 0.1;
  EXPECT_LT(fabs(vasicek(t, r0, a, b, sigma) - 0.955408), eps);
}

