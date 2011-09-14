#pragma once
#include <vector>
#include "matrix_meta.hpp"

namespace sfinx { namespace matrix {

/// Functions
//
template <typename Mx>
inline size_t rows(Mx const& m)
{
  return m.rows();
}

template <typename Mx>
inline size_t cols(Mx const& m)
{
  return m.cols();
}

template <typename Mx>
inline size_t size(Mx const& m)
{
  return m.size();
}

template <typename Mx>
inline auto row(Mx const& m, size_t n) -> decltype(row_type(m))
{
  return m.row(n);
}

template <typename Mx>
inline auto col(Mx const& m, size_t n) -> decltype(col_type(m))
{
  return m.col(n);
}

template <typename Mx>
inline auto transpose(Mx const& m) -> decltype(transpose_type(m))
{
  return m.transpose();
}

template <typename Mx>
inline Mx inverse(Mx const& m)
{
  return m.inverse();
}

template <typename Mx>
inline auto eigenvalues(Mx const& m) -> decltype(col_type(m))
{
  auto ev = m.eigenvalues();
  decltype(col_type(m)) ret;
  for (int i = 0; i < ev.rows(); ++i)
    ret(i) = ev(i).real();
  return ret;
}

template <typename Mx>
inline auto eigenvectors(Mx const& m) -> Mx
{
  Eigen::SelfAdjointEigenSolver<Mx> solver(m);
  return solver.eigenvectors();
}

} // namespace matrix 
} // namespace sfinx

