#pragma once

namespace sfinx { namespace matrix {

/**
 * Meta functions, for type computation, etc.
 */

auto row_type(Matrix const&) -> RowVector;

auto col_type(Matrix const&) -> ColumnVector;

auto transpose_type(Matrix const&) -> Matrix;

} } // namespace sfinx::matrix

