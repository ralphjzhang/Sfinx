#pragma once

namespace sfinx {

template <typename Matrix>
size_t rows(Matrix const& m)
{
  return m.rows();
}

template <typename Matrix>
size_t cols(Matrix const& m)
{
  return m.cols();
}

template <typename Matrix>
size_t size(Matrix const& m)
{
  return m.size();
}

template <typename Matrix>
size_t transpose(Matrix const& m)
{
  return m.transpose();
}

template <typename Matrix>
size_t inverse(Matrix const& m)
{
  return m.inverse();
}

} // namespace sfinx

