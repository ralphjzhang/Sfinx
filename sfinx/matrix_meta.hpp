#pragma once

namespace sfinx { namespace matrix {

/**
 * Meta functions, for type computation, etc.
 */

template <typename Mx>
struct scalar_type;

/// vector_type<Mx>
//
template <typename Mx>
struct vector_type;

/// row_type<Mx>
//
template <template <typename T, int rows, int cols, int, int, int> class Mx,
          typename T, int rows, int cols, int options, int maxrows, int maxcols>
auto row_type(Mx<T, rows, cols, options, maxrows, maxcols>)
  -> Mx<T, 1, cols, !options, 1, maxcols>;
                  //TODO !options changes col-major to row-major, but should be a better way

/// col_type(Mx)
//
template <template <typename T, int rows, int cols, int, int, int> class Mx,
          typename T, int rows, int cols, int options, int maxrows, int maxcols>
auto col_type(Mx<T, rows, cols, options, maxrows, maxcols>)
  -> Mx<T, rows, 1, options, maxrows, 1>;

/// transpose_type<Mx>
//
template <template <typename T, int rows, int cols, int, int, int> class Mx,
          typename T, int rows, int cols, int options, int maxrows, int maxcols>
auto transpose_type(Mx<T, rows, cols, options, maxrows, maxcols>)
  -> Mx<T, cols, rows, options, maxcols, maxrows>;


} } // namespace sfinx::matrix

