// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <boost/geometry/util/exact_coordinate.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/strategies/cartesian/exact/coordinate_predicates.hpp>

#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using boost::multiprecision::cpp_int;
using rational = boost::multiprecision::cpp_rational;
using coordinate = bg::exact_coordinate<rational>;
using decimal = bg::exact_coordinate<std::string>;
using limits = std::numeric_limits<double>;

static_assert(!std::is_convertible<coordinate, double>::value, "The binary64 view is explicit");
static_assert(std::is_same<decltype(std::declval<decimal&>().scalar()), std::string const&>::value,
              "The authoritative scalar is read-only");

rational power_of_two(int exponent)
{
    if (exponent >= 0) { return rational(cpp_int(1) << exponent); }
    return rational(1, cpp_int(1) << -exponent);
}

void check(rational const& value, double expected, bool exact = false)
{
    coordinate const c(value);
    BOOST_TEST_EQ(c.scalar(), value);
    BOOST_TEST_EQ(c.binary64(), expected);
    BOOST_TEST_EQ(std::signbit(c.binary64()), std::signbit(expected));
    BOOST_TEST_EQ(c.binary64_is_exact(), exact);
}

void ties_and_extremes()
{
    rational const half_ulp = power_of_two(-53);
    rational const epsilon = power_of_two(-150);
    double const odd = std::nextafter(1.0, 2.0);
    double const even = std::nextafter(odd, 2.0);
    check(1 + half_ulp - epsilon, 1);
    check(1 + half_ulp, 1);
    check(1 + half_ulp + epsilon, odd);
    check(1 + 3 * half_ulp, even);
    check(-1 - 3 * half_ulp, -even);
    check(2 - half_ulp, 2);
    check(rational(odd), odd, true);

    double const min = limits::denorm_min();
    rational const half_min = power_of_two(-1075);
    check(rational(0), 0, true);
    check(half_min, 0);
    check(-half_min, -0.0);
    check(half_min + power_of_two(-1200), min);
    check(-half_min - power_of_two(-1200), -min);
    check(3 * half_min, 2 * min);
    check(rational(min), min, true);
    check(power_of_two(-2000), 0);
    check(power_of_two(-1022) - half_min, (limits::min)());

    double const max = (limits::max)();
    rational const overflow = rational(max) + power_of_two(970);
    check(rational(max), max, true);
    check(overflow - 1, max);
    check(overflow, limits::infinity());
    check(-overflow, -limits::infinity());
    check(power_of_two(4000), limits::infinity());
}

void rounding_modes()
{
    int const saved = std::fegetround();
    for (int mode : {FE_TONEAREST, FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO})
    {
        BOOST_TEST_EQ(std::fesetround(mode), 0);
        check(1 + power_of_two(-53), 1);
        check(1 + power_of_two(-53) + power_of_two(-90), std::nextafter(1.0, 2.0));
        check(-power_of_two(-1075), -0.0);
        check(power_of_two(-1074), limits::denorm_min(), true);
    }
    BOOST_TEST_EQ(std::fesetround(saved), 0);
}

void decimal_intent_and_invariants()
{
    decimal d;
    BOOST_TEST_EQ(d.scalar(), "0");
    BOOST_TEST_EQ(d.binary64(), 0);
    BOOST_TEST(d.binary64_is_exact());

    d.assign("+001.00e-1");
    rational const tenth(1, 10);
    BOOST_TEST_EQ(d.rational(), tenth);
    BOOST_TEST_NE(d.rational(), rational(0.1));
    BOOST_TEST_EQ(d.binary64(), 0.1);
    BOOST_TEST(!d.binary64_is_exact());

    decimal const copy(std::move(d));
    BOOST_TEST_EQ(copy.scalar(), "+001.00e-1");
    BOOST_TEST_EQ(d.rational(), tenth);
    d = copy;
    decimal const& self = d;
    d = self;
    BOOST_TEST_EQ(d.scalar(), copy.scalar());

    for (std::string const text : {"", "+", ".", "1e", "1e-", "1..2", "1e2e3", " 1", "1 ",
                                  "NaN", "inf", "0x1p0", "2/3", "1,2"})
    {
        BOOST_TEST_THROWS(d.assign(text), std::invalid_argument);
        BOOST_TEST_EQ(d.scalar(), copy.scalar());
        BOOST_TEST_EQ(d.binary64(), copy.binary64());
        BOOST_TEST(!d.binary64_is_exact());
    }
    BOOST_TEST_THROWS(d.assign("1e10001"), std::out_of_range);
    BOOST_TEST_THROWS(d.assign("1e-9999999999999999999999"), std::out_of_range);
    BOOST_TEST_THROWS(d.assign("." + std::string(10001, '0')), std::out_of_range);
    BOOST_TEST_THROWS(d.assign(std::string("1\0suffix", 8)), std::invalid_argument);
    BOOST_TEST_EQ(d.rational(), tenth);
    BOOST_TEST_EQ(d.binary64(), 0.1);

    decimal target("3.");
    target = std::move(d);
    BOOST_TEST_EQ(target.scalar(), d.scalar());
    BOOST_TEST_EQ(target.binary64(), d.binary64());
    BOOST_TEST_EQ(decimal("-.125E+1").rational(), rational(-5, 4));
    BOOST_TEST(decimal(".125").binary64_is_exact());
    BOOST_TEST_EQ(decimal("-0.0").scalar(), "-0.0");
    BOOST_TEST(!std::signbit(decimal("-0.0").binary64()));
    BOOST_TEST(std::signbit(decimal("-1e-10000").binary64()));
    BOOST_TEST(std::isinf(decimal("1e10000").binary64()));

    decimal const tie("1.00000000000000011102230246251565404236316680908203125");
    BOOST_TEST_EQ(tie.rational(), rational(1) + power_of_two(-53));
    BOOST_TEST_EQ(tie.binary64(), 1);
    BOOST_TEST_EQ(decimal(tie.scalar() + "1").binary64(), std::nextafter(1.0, 2.0));

    bg::model::point<decimal, 2, bg::cs::cartesian> point;
    bg::set<1>(point, target);
    BOOST_TEST_EQ(bg::get<0>(point).rational(), rational(0));
    BOOST_TEST_EQ(bg::get<1>(point).scalar(), target.scalar());
    BOOST_TEST_EQ(bg::get<1>(point).binary64(), target.binary64());
}

struct counted_conversion : bg::traits::exact_coordinate_conversion<rational>
{
    static int conversions;

    static rational const& to_rational(rational const& value)
    {
        ++conversions;
        return value;
    }
};

int counted_conversion::conversions = 0;

void predicate_routing()
{
    using counted = bg::exact_coordinate<rational, counted_conversion>;
    counted const half(rational(1, 2));
    counted const near_half(rational(1, 2) + power_of_two(-56));
    counted_conversion::conversions = 0;
    BOOST_TEST_EQ(ec::compare_coordinates(decimal("0.1"), half), -1);
    BOOST_TEST_EQ(ec::compare_coordinates(decimal("0.5"), half), 0);
    BOOST_TEST_EQ(counted_conversion::conversions, 0);
    BOOST_TEST(!ec::compare(half, near_half).certain);
    BOOST_TEST_EQ(ec::compare_coordinates(half, near_half), -1);
    BOOST_TEST_EQ(counted_conversion::conversions, 2);
    BOOST_TEST_EQ(ec::compare_coordinates(decimal("0.1"), coordinate(rational(0.1))), -1);
    BOOST_TEST_EQ(ec::compare_coordinates(decimal("1e500"), decimal("2e500")), -1);
    BOOST_TEST_EQ(ec::compare_coordinates(decimal("-1e-500"), decimal("0")), -1);

    std::array<decimal, 2> const a{{decimal("0"), decimal("0")}};
    std::array<decimal, 2> const b{{decimal("0.1"), decimal("0.1")}};
    std::array<counted, 2> const left{{counted(rational(1, 10)), counted(rational(1, 5))}};
    counted_conversion::conversions = 0;
    BOOST_TEST(ec::orientation_filter(a, b, left).certain);
    BOOST_TEST_EQ(ec::orientation(a, b, left), 1);
    BOOST_TEST_EQ(ec::orientation(a, left, b), -1);
    BOOST_TEST_EQ(counted_conversion::conversions, 0);

    std::array<decimal, 2> const collinear{{decimal(".2"), decimal("0.20")}};
    BOOST_TEST(!ec::orientation_filter(a, b, collinear).certain);
    BOOST_TEST_EQ(ec::orientation(a, b, collinear), 0);
    std::array<coordinate, 2> const off{{coordinate(rational(0.1)), coordinate(rational(1, 10))}};
    BOOST_TEST(!ec::orientation_filter(a, b, off).certain);
    BOOST_TEST_EQ(ec::orientation(a, b, off), -1);

    std::array<decimal, 2> const huge{{decimal("1e500"), decimal("1e501")}};
    BOOST_TEST(!ec::orientation_filter(a, b, huge).certain);
    BOOST_TEST_EQ(ec::orientation(a, b, huge), 1);
    std::array<decimal, 2> const tiny{{decimal("1e-500"), decimal("2e-500")}};
    BOOST_TEST(!ec::orientation_filter(a, b, tiny).certain);
    BOOST_TEST_EQ(ec::orientation(a, b, tiny), 1);
}

void test_forward_signs()
{
    unsigned accepted = 0, uncertain = 0;
    for (int exponent : {-160, -128, -22, 0, 44, 120, 160})
    {
        rational const step = power_of_two(exponent) / 10;
        auto point = [&](int x, int y)
        { return std::array<coordinate, 2>{{coordinate(step*x), coordinate(step*y)}}; };
        for (int i = -9; i <= 9; ++i)
        {
            auto const a = point(i, 2*i), b = point(i + 3, 2*i - 5);
            for (int k : {-1, 0, 1})
            {
                auto const c = point(i + 6, 2*i - 10 + k);
                auto const f = ec::orientation_filter(a, b, c);
                if (f.certain) { ++accepted; BOOST_TEST_EQ(f.value, k); }
                else { ++uncertain; }
                BOOST_TEST_EQ(ec::orientation(a, b, c), k);
            }
        }
    }
    BOOST_TEST(accepted != 0 && uncertain != 0);
}

int main()
{
    ties_and_extremes();
    rounding_modes();
    decimal_intent_and_invariants();
    predicate_routing();
    test_forward_signs();
    return boost::report_errors();
}
