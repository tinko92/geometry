// Boost.Geometry
// Distributed under the Boost Software License, Version 1.0.

#ifndef BOOST_GEOMETRY_TEST_PROJECTED_BENCHMARK_FIXTURE_HPP
#define BOOST_GEOMETRY_TEST_PROJECTED_BENCHMARK_FIXTURE_HPP

#include <random>

double constexpr easting = 450000, northing = 5200000;

template <typename Function>
void for_each_input_point(input& g, Function transform)
{
    bg::for_each_point(g.a, transform);
    bg::for_each_point(g.b, transform);
    bg::for_each_point(g.la, transform);
    bg::for_each_point(g.lb, transform);
    transform(g.query);
}

point surveyed_point(double x, double y)
{
    // Snap once to the binary64 ULP of each projected axis (2^(18-52)
    // and 2^(22-52) metres). The same vertices then fit both frames exactly.
    return point(coordinate(std::ldexp(std::round(std::ldexp(x, 34)), -34)),
                 coordinate(std::ldexp(std::round(std::ldexp(y, 30)), -30)));
}

template <typename Ring, typename Center>
void corridor(Ring& ring, std::size_t n, double west, double east, Center center)
{
    double const tau = 2 * std::acos(-1.0);
    for (bool upper : {true, false})
    {
        auto const count = upper ? n / 2 : n - n / 2;
        for (std::size_t i = 0; i < count; ++i)
        {
            double const t = double(upper ? i : count - i - 1) / double(count - 1);
            double const x = west + (east - west) * (t + 0.12 * std::sin(tau * t) / tau);
            double const half_width = 64 + 13 * std::cos(2 * tau * t);
            ring.push_back(surveyed_point(x, center(x) + (upper ? half_width : -half_width)));
        }
    }
    ring.push_back(ring.front());
}

input projected_input(std::size_t vertices, std::size_t index, std::mt19937_64& random,
                      bool local)
{
    input g;
    double const dx = double(std::int64_t(random() % 401) - 200);
    double const dy = double(std::int64_t(random() % 201) - 100);
    double const phase = double(random() % 1000) / 1000;
    auto const center = [=](double x) { return dy + 84 * std::sin((x - dx) / 380 + phase); };
    corridor(g.a.outer(), vertices, dx - 1500, dx + 1250, center);
    corridor(g.b.outer(), vertices, dx - 1130, dx + 1550,
        [=](double x) { return center(x) + 210 * std::sin((x - dx) / 295 + phase + 0.41); });

    // Crossing routes vs. roughly parallel routes, both with bends and
    // staggered vertices. These exercise linear relate as well as clipping.
    auto const count = (std::max)(std::size_t(2), vertices / 8);
    for (std::size_t i = 0; i < count; ++i)
    {
        double const t = (double(i) + 0.37) / double(count);
        double const x = dx - 1330 + 2380 * t;
        double const y = dy - 240 + 510 * t;
        g.la.push_back(line{surveyed_point(x - 137, dy - 570), surveyed_point(x - 39, dy - 171),
                            surveyed_point(x + 41, dy + 138), surveyed_point(x + 106, dy + 543)});
        g.lb.push_back(line{surveyed_point(dx - 1610, y - 57), surveyed_point(dx - 581, y + 21),
                            surveyed_point(dx + 472, y - 29), surveyed_point(dx + 1710, y + 53)});
    }
    g.query = surveyed_point(dx, center(dx) + (index % 3 == 0 ? 230 : 0));

    for_each_input_point(g, [=](point& p)
    {
        double const x = bg::get<0>(p), y = bg::get<1>(p);
        double const east = x + easting, north = y + northing;
        // Each axis stays in one positive binade, so any same-axis
        // subtraction of projected input coordinates satisfies Sterbenz.
        if (east < (1 << 18) || east >= (1 << 19) || north < (1 << 22) || north >= (1 << 23)
            || east - easting != x || north - northing != y)
        {
            throw std::logic_error("Projected fixture left its Sterbenz box or lost precision");
        }
        if (!local) { p = point(coordinate(east), coordinate(north)); }
    });
    return g;
}

#endif
