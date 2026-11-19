// Boost.Geometry
// Unit Test
// Copyright (c) 2026
// Distributed under the Boost Software License, Version 1.0.

#include <boost/geometry/strategies/cartesian/exact/bounded_radical_field.hpp>
#include <boost/core/lightweight_test.hpp>

namespace ec = boost::geometry::strategies::detail::exact_cartesian;
using field = ec::round_buffer::bounded_radical_field;
using number = ec::exact_number;
using boost::multiprecision::cpp_int;

void test_square_classes()
{
    field f;
    field::value a, b, c, d;
    BOOST_TEST(f.try_sqrt(2, a));
    BOOST_TEST(f.try_sqrt(5, b));
    BOOST_TEST(f.try_sqrt(10, c));
    BOOST_TEST(f.try_sqrt(number(8) / 9, d));
    BOOST_TEST_EQ(f.rank(), 2u);
    BOOST_TEST_EQ(f.sign(f.subtract(c, f.multiply(a, b))), 0);
    BOOST_TEST_EQ(f.sign(f.subtract(d, f.scale(a, number(2) / 3))), 0);
    BOOST_TEST_EQ(f.sign(f.subtract(f.multiply(a, a), f.rational(2))), 0);
    auto const s = f.add(a, b);
    BOOST_TEST_EQ(f.sign(f.subtract(f.multiply(s, s), f.add(f.rational(7), f.scale(c, 2)))), 0);

    number rational = -7;
    BOOST_TEST(f.try_rational_quotient(f.scale(s, number(3) / 7), s, rational));
    BOOST_TEST_EQ(rational, number(3) / 7);
    BOOST_TEST_NOT(f.try_rational_quotient(a, b, rational));
    BOOST_TEST_EQ(rational, number(3) / 7);
    BOOST_TEST_THROWS(f.try_rational_quotient(a, field::value{}, rational), std::domain_error);
    BOOST_TEST_THROWS(f.try_sqrt(-1, c), std::domain_error);
    field other;
    BOOST_TEST_THROWS(other.sign(a), std::invalid_argument);
}

void test_close_and_high_rank_signs()
{
    field f;
    field::value a, b, c, d, e;
    BOOST_TEST(f.try_sqrt(2, a));
    BOOST_TEST(f.try_sqrt(3, b));
    BOOST_TEST(f.try_sqrt(5, c));
    BOOST_TEST(f.try_sqrt(7, d));
    // A positive factor from each independent generator. Their product
    // expands to all 16 terms with mixed signs and a small net value.
    auto const product = f.multiply(f.multiply(f.subtract(a, f.rational(1)), f.subtract(b, f.rational(1))),
                                    f.multiply(f.subtract(c, f.rational(2)), f.subtract(d, f.rational(2))));
    BOOST_TEST_EQ(f.sign(product), 1);
    BOOST_TEST_EQ(f.sign(f.scale(product, -1)), -1);
    BOOST_TEST_EQ(f.sign(f.subtract(product, product)), 0);

    cpp_int const denominator = cpp_int(1) << 300;
    cpp_int const lower = boost::multiprecision::sqrt(cpp_int(2) * denominator * denominator);
    BOOST_TEST_EQ(f.sign(f.subtract(a, f.rational(number(lower) / denominator))), 1);
    BOOST_TEST_EQ(f.sign(f.subtract(a, f.rational(number(lower + 1) / denominator))), -1);

    e = f.rational(42);
    BOOST_TEST_NOT(f.try_sqrt(11, e));
    BOOST_TEST_EQ(f.rank(), field::max_rank);
    BOOST_TEST_EQ(f.sign(f.subtract(e, f.rational(42))), 0);
    BOOST_TEST(f.try_sqrt(210, e)); // dependent even at the rank cap
    BOOST_TEST_EQ(f.sign(f.subtract(e, f.multiply(f.multiply(a, b), f.multiply(c, d)))), 0);
}

int main()
{
    test_square_classes();
    test_close_and_high_rank_signs();
    return boost::report_errors();
}
