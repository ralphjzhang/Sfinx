/// Common testcases for matrix
//

TEST(matrix, size)
{
  auto m = matrix<int, 2, 3>();
  m << 1, 2, 3, 4, 5, 6;
  EXPECT_EQ(rows(m), 2u);
  EXPECT_EQ(cols(m), 3u);
  EXPECT_EQ(size(m), 6u);
}

TEST(matrix, col_row)
{
  auto m = matrix<int, 2, 3>();
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
    auto m = matrix<int, 2, 3>();
    m << 1, 2, 3, 4, 5, 6;
    decltype(transpose(m)) mt;
    mt << 1, 4, 2, 5, 3, 6;
    EXPECT_TRUE(mt == transpose(m));
  }
  {
    auto m = matrix<float, 2, dynamic>();
    resize(m, 2, 4);
    m << 1, 2, 3, 4, 5, 6, 7, 8;
    decltype(transpose(m)) mt;
    resize(mt, 4, 2);
    mt << 1, 5, 2, 6, 3, 7, 4, 8;
    EXPECT_TRUE(mt == transpose(m));
  }
}

TEST(matrix, inverse)
{
  auto m = matrix<float, 2, 2>(), mr = matrix<float, 2, 2>();
  m << 5, 4, 6, 5;
  mr << 5, -4, -6, 5;
  EXPECT_TRUE(mr == inverse(m));
}

TEST(matrix, eigenvalues)
{
  auto m = matrix<double, 2, 2>();
  m << 1, 2, 2, 1;
  decltype(eigenvalues(m)) m1;
  m1 << 3.0, -1.0;
  auto m2 = eigenvalues(m);
  EXPECT_DOUBLE_EQ(m1(0), m2(0));
  EXPECT_DOUBLE_EQ(m1(1), m2(1));
}

TEST(matrix, eigenvectors)
{
  auto m = matrix<double, 2, 2>();
  m << 1, 2, 2, 1;
  auto m1 = eigenvectors(m);
  double eps = 1e-6;
  EXPECT_NEAR(m1(0, 0), -0.707107, eps);
  EXPECT_NEAR(m1(0, 1), -0.707107, eps);
  EXPECT_NEAR(m1(1, 0), 0.707107, eps);
  EXPECT_NEAR(m1(1, 1), -0.707107, eps);
}

