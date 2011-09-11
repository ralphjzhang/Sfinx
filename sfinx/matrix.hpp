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
inline auto row(Mx const& m, size_t n)
  -> typename row_type<Mx>::type
{
  return m.row(n);
}

template <typename Mx>
inline auto col(Mx const& m, size_t n)
  -> typename col_type<Mx>::type
{
  return m.col(n);
}

template <typename Mx>
inline auto transpose(Mx const& m)
  -> typename transpose_type<Mx>::type
{
  return m.transpose();
}

template <typename Mx>
inline size_t inverse(Mx const& m)
{
  return m.inverse();
}

template <typename Mx>
inline auto eigenvalues(Mx const& m)
  -> std::vector<typename scalar_type<Mx>::type>
{
  return m.eigenvalues();
}

template <typename Mx>
inline auto eigenvectors(Mx const& m)
  -> std::vector<typename vector_type<Mx>::type>
{
  return m.eigenvectors();
}

} // namespace matrix 
} // namespace sfinx

