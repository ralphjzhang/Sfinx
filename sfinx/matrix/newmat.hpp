#pragma once
#include <newmat.h>

template <typename T>
MatrixInput& operator , (MatrixInput& m, T t)
{
  m << t;
  return m;
}

namespace sfinx { namespace matrix {

static const auto dyna = 0;

template <typename T, size_t Rows, size_t Cols>
auto matrix() -> Matrix
{
  return Matrix(Rows, Cols);
}

inline size_t rows(Matrix const& m)
{
  return m.Nrows();
}

inline size_t cols(Matrix const& m)
{
  return m.Ncols();
}

inline size_t size(Matrix const& m)
{
  return m.Nrows() * m.Ncols();
}

inline void resize(Matrix& m, size_t row, size_t col)
{
  m.ReSize(row, col);
}

inline RowVector row(Matrix const& m, size_t n)
{
  return m.Row(n);
}

inline ColumnVector col(Matrix const& m, size_t n)
{
  return m.Column(n);
}

inline Matrix transpose(Matrix const& m)
{
  return m.t();
}

inline Matrix inverse(Matrix const& m)
{
  return m.i();
}

inline RowVector eigenvalues(Matrix const& m)
{
  //DiagonalMatrix d;
  //EigenValues(m, d);
  RowVector d;
  return d;
}

inline Matrix eigenvectors(Matrix const& m)
{
  DiagonalMatrix d;
  Matrix v;
  //EigenValues(m, d, v);
  return v;
}

} } // namespace sfinx::matrix

