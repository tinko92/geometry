// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
#ifndef BOOST_GEOMETRY_TEST_EXACT_TEST_CONTEXT_HPP
#define BOOST_GEOMETRY_TEST_EXACT_TEST_CONTEXT_HPP

#include <boost/geometry/strategies/cartesian/exact/operation_context.hpp>

template <typename C = double>
inline void prepare_exact_cache_repeat()
{
    namespace op = boost::geometry::strategies::detail::exact_cartesian::operation_detail;
    auto* const owner = op::current<C>();
    if (!owner) return;
    // Cache growth is permitted. Reserve outside the allocation census for a
    // repeat of the priming call; unrelated numerical allocations still count.
    auto& cache = owner->storage();
    cache.native.reserve(2 * cache.native.size() + 8);
    cache.auxiliary.reserve(2 * cache.auxiliary.size() + 8);
}
#endif
