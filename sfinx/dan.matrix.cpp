#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "matrix.hpp"

using namespace sfinx::matrix;
using namespace Eigen;

TEST(matrix, size)
{
  Matrix<int, 2, 3> m;
  m << 1, 2, 3, 4, 5, 6;
  EXPECT_EQ(rows(m), 2u);
  EXPECT_EQ(cols(m), 3u);
  EXPECT_EQ(size(m), 6u);
}

TEST(matrix, col_row)
{
  Matrix<int, 2, 3> m;
  m << 1, 2, 3, 4, 5, 6;
  decltype(col(m, 0)) c0, c2;
  c0 << 1, 4; c2 << 3, 6;
  EXPECT_TRUE(c0 == col(m, 0));
  EXPECT_TRUE(c2 == col(m, 2));

  decltype(row(m, 0)) r0, r1;
  r0 << 1, 2, 3; r1 << 4, 5, 6;
  EXPECT_TRUE(r0 == row(m, 0));
  EXPECT_TRUE(r1 == row(m, 1));
}

TEST(matrix, transpose)
{
  {
    Matrix<int, 2, 3> m;
    m << 1, 2, 3, 4, 5, 6;
    decltype(transpose(m)) mt;
    mt << 1, 4, 2, 5, 3, 6;
    EXPECT_TRUE(mt == transpose(m));
  }
  {
    Matrix<int, 2, Dynamic> m;
    m.resize(2, 4);
    m << 1, 2, 3, 4, 5, 6, 7, 8;
    decltype(transpose(m)) mt;
    mt.resize(4, 2);
    mt << 1, 5, 2, 6, 3, 7, 4, 8;
    EXPECT_TRUE(mt == transpose(m));
  }
}

TEST(matrix, inverse)
{
  Matrix<float, 2, 2> m, mr;
  m << 5, 4, 6, 5;
  mr << 5, -4, -6, 5;
  EXPECT_TRUE(mr == inverse(m));
}

TEST(matrix, eigenvalues)
{
  Matrix<double, 2, 2> m;
  m << 1, 2, 2, 1;
  decltype(eigenvalues(m)) m1;
  m1 << 3.0, -1.0;
  auto m2 = eigenvalues(m);
  EXPECT_LT(fabs(m1(0) - m2(0)), std::numeric_limits<double>::epsilon() * 10);
  EXPECT_LT(fabs(m1(1) - m2(1)), std::numeric_limits<double>::epsilon() * 10);
}

TEST(matrix, eigenvectors)
{
  Matrix<double, 2, 2> m;
  m << 1, 2, 2, 1;
  auto m1 = eigenvectors(m);
  double eps = 1e-6;
  EXPECT_LT(fabs(m1(0, 0) + 0.707107), eps);
  EXPECT_LT(fabs(m1(0, 1) + 0.707107), eps);
  EXPECT_LT(fabs(m1(1, 0) - 0.707107), eps);
  EXPECT_LT(fabs(m1(1, 1) + 0.707107), eps);
}

