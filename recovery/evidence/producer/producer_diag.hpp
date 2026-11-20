// External, single-thread counter replay. Not a production Geometry header.
#ifndef PRODUCER_BOUND_DIAGNOSTIC_HPP
#define PRODUCER_BOUND_DIAGNOSTIC_HPP

#ifndef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#error Compile the diagnostic with -DBOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#endif
#include <boost/geometry/strategies/cartesian/exact/filter_statistics.hpp>
#include <array>
#include <cstdint>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>

namespace producer_diag {
namespace ec = boost::geometry::strategies::detail::exact_cartesian;
enum kind { orientation, direction, signed_measure, unsigned_measure, mixed_measure,
            coordinate_order, coordinate_equal, kinds };
inline constexpr char const* names[] = {"orientation", "direction", "signed_measure",
    "unsigned_measure", "mixed_measure", "coordinate_order", "coordinate_equal"};
struct stages {
    std::uint64_t first_calls{}, first_misses{}, refinements{}, refined_misses{};
};
struct counts : stages {
    stages axes{};
    std::uint64_t exact_after_retry{};
};
inline thread_local std::array<counts, kinds> rows{};
inline void reset() {
    rows = {};
    ec::filtering_statistics() = {};
}
inline bool observe(stages& s, bool refined, bool hit) {
    ++(refined ? s.refinements : s.first_calls);
    (refined ? s.refined_misses : s.first_misses) += !hit;
    return hit;
}
inline void write_header(std::ostream& o) {
    o << "workload,operation,predicate,first_calls,first_misses,refinements,refined_misses,"
         "shortcut_after_miss,exact_fallbacks,small_integer_exact,post_refine_exact,tagged_exact,"
         "axis_first_calls,axis_first_misses,axis_refined_calls,axis_refined_misses\n";
}
inline void dump(std::ostream& o, std::string const& workload, std::string const& operation) {
    auto const& s = ec::filtering_statistics();
    if (rows[signed_measure].exact_after_retry + rows[unsigned_measure].exact_after_retry
        + rows[mixed_measure].exact_after_retry != s.compact_measure_comparisons + s.full_measure_comparisons)
        throw std::logic_error("measure exact exits disagree with Geometry counters");
    for (int k = 0; k < kinds; ++k) {
        auto const& r = rows[k];
        auto const small = k == unsigned_measure ? s.small_integer_measure_comparisons : 0;
        auto const tagged = k == orientation ? s.tagged_orientation_calls : 0;
        if (r.first_misses > r.first_calls || small + r.refinements > r.first_misses
            || r.refined_misses > r.refinements || r.refined_misses != r.exact_after_retry
            || tagged > r.exact_after_retry)
            throw std::logic_error(std::string("bad stage accounting: ") + names[k]);
        o << std::quoted(workload, '"', '"') << ',' << std::quoted(operation, '"', '"') << ','
          << names[k] << ',' << r.first_calls << ',' << r.first_misses << ',' << r.refinements
          << ',' << r.refined_misses << ',' << r.first_misses - small - r.refinements
          << ',' << r.exact_after_retry + small << ',' << small << ',' << r.exact_after_retry
          << ',' << tagged << ',' << r.axes.first_calls << ',' << r.axes.first_misses
          << ',' << r.axes.refinements << ',' << r.axes.refined_misses << '\n';
    }
}
} // namespace producer_diag
#endif
