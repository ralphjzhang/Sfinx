#pragma once
#include <iterator>
#include <utility>

namespace std {
  template <typename T> typename T::const_iterator begin(T const& t) { return t.begin(); }
}

namespace sfinx {

template <typename Ret, typename Comb, typename Accum, typename... Ranges>
Ret combine(Comb comb, Accum accum, Ret init, Ranges const&... ranges)
{
  comb((*std::begin(ranges))...);
  comb((*std::advance(ranges, 1))...);
  return init;
}

} // namespace sfinx


