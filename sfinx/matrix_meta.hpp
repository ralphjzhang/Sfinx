#pragma once

namespace sfinx { namespace matrix {

/**
 * Meta functions, for type extraction, etc.
 */

template <typename Mx>
struct scalar_type;

/// vector_type<Mx>
//
template <typename Mx>
struct vector_type;

/// row_type<Mx>
//
template <typename Mx>
struct row_type;

template <template <typename T, int rows, int cols, int, int, int> class Mx,
          typename T, int rows, int cols, int options, int maxrows, int maxcols>
struct row_type<Mx<T, rows, cols, options, maxrows, maxcols>>
{
  typedef Mx<T, 1, cols, !options, 1, maxcols> type;
};

/// col_type<Mx>
//
template <typename Mx>
struct col_type;

template <template <typename T, int rows, int cols, int, int, int> class Mx,
          typename T, int rows, int cols, int options, int maxrows, int maxcols>
struct col_type<Mx<T, rows, cols, options, maxrows, maxcols>>
{
  typedef Mx<T, rows, 1, options, maxrows, 1> type;
};

/// transpose_type<Mx>
//
template <typename Mx>
struct transpose_type;

template <template <typename T, int rows, int cols, int, int, int> class Mx,
          typename T, int rows, int cols, int options, int maxrows, int maxcols>
struct transpose_type<Mx<T, rows, cols, options, maxrows, maxcols>>
{
  typedef Mx<T, cols, rows, options, maxcols, maxrows> type;
};


} } // namespace sfinx::matrix

