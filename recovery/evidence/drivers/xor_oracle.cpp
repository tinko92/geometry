// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <vector>

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>



namespace bg = boost::geometry;
namespace ex = bg::strategies::detail::exact_cartesian;
using point = bg::model::d2::point_xy<double>;
using polygon = bg::model::polygon<point>;
using multi = bg::model::multi_polygon<polygon>;
using rational = ex::exact_number;
using exact_point = bg::model::d2::point_xy<rational>;
using exact_polygon = bg::model::polygon<exact_point>;
using exact_multi = bg::model::multi_polygon<exact_polygon>;
using exact_open = bg::model::multi_polygon<bg::model::polygon<exact_point, false, false>>;
using strategy = bg::strategies::relate::cartesian_exact<double, bg::strategies::relate::nearest_output>;

template <typename G = polygon>
G wkt(char const* text)
{
    G result;
    bg::read_wkt(text, result);
    bg::correct(result);
    return result;
}

template <typename G>
std::string signature(G const& result)
{
    std::ostringstream out;
    out << std::hexfloat << bg::wkt(result);
    return out.str();
}

// Independent oracle: sweep open strips between vertex/crossing abscissae.
// At a strip midpoint, compare the complete vertical occupancy profile and
// boundary multiplicities. A second sweep with the axes swapped checks vertical
// boundaries too. This does not use the noder, its atom labels, or BG predicates.
struct sweep
{
    struct line { rational low, high, slope, intercept; std::size_t owner; };
    std::set<rational> breaks;
    std::vector<line> lines;

    template <bool Swap, typename G>
    void add(G const& g, std::size_t owner)
    {
        auto add_ring = [&](auto const& r)
        {
            for (std::size_t i = 0; i < r.size(); ++i)
            {
                auto const& a = r[i]; auto const& b = r[(i + 1) % r.size()];
                rational const x0 = bg::get<Swap ? 1 : 0>(a), x1 = bg::get<Swap ? 1 : 0>(b);
                breaks.insert(x0);
                if (x0 == x1) continue;
                rational const y0 = bg::get<Swap ? 0 : 1>(a), y1 = bg::get<Swap ? 0 : 1>(b);
                rational const slope = (y1 - y0) / (x1 - x0);
                lines.push_back({(std::min)(x0, x1), (std::max)(x0, x1), slope, y0 - slope * x0, owner});
            }
        };
        if constexpr (std::is_same<bg::tag_t<G>, bg::multi_polygon_tag>::value)
            for (auto const& p : g) add<Swap>(p, owner);
        else
        {
            add_ring(bg::exterior_ring(g));
            for (auto const& hole : bg::interior_rings(g)) add_ring(hole);
        }
    }

    void check(std::size_t shape_count)
    {
        for (std::size_t i = 0; i < lines.size(); ++i)
        for (std::size_t j = i + 1; j < lines.size(); ++j)
        {
            auto const& a = lines[i]; auto const& b = lines[j];
            if (a.slope == b.slope) continue;
            rational const x = (b.intercept - a.intercept) / (a.slope - b.slope);
            if (a.low < x && x < a.high && b.low < x && x < b.high) breaks.insert(x);
        }
        for (auto it = breaks.begin(); it != breaks.end() && std::next(it) != breaks.end(); ++it)
        {
            rational const x = (*it + *std::next(it)) / 2;
            std::map<rational, std::vector<std::size_t>> events;
            for (auto const& l : lines)
                if (l.low < x && x < l.high) events[l.slope * x + l.intercept].push_back(l.owner);
            std::vector<bool> inside(shape_count, false);
            for (auto const& e : events)
            {
                bool const before = inside[0] != inside[1];
                std::size_t emitted = 0;
                for (auto id : e.second) { inside[id] = !inside[id]; emitted += id >= 2; }
                bool const after = inside[0] != inside[1];
                BOOST_TEST_EQ(emitted, unsigned(before != after));
                BOOST_TEST_EQ(std::count(inside.begin() + 2, inside.end(), true), int(after));
            }
            BOOST_TEST(std::none_of(inside.begin(), inside.end(), [](bool x) { return x; }));
        }
    }
};

template <typename A, typename B, typename Result>
void verify(A const& a, B const& b, Result const& result)
{
    sweep x, y;
    x.add<false>(a, 0); x.add<false>(b, 1);
    y.add<true>(a, 0); y.add<true>(b, 1);
    auto check_ring = [](auto const& r, bool shell)
    {
        using ring = std::decay_t<decltype(r)>;
        constexpr bool closed = bg::closure<ring>::value == bg::closed;
        if (r.size() < 3u + closed) { BOOST_ERROR("short output ring"); return; }
        if constexpr (closed)
        {
            BOOST_TEST(bg::get<0>(r.front()) == bg::get<0>(r.back()));
            BOOST_TEST(bg::get<1>(r.front()) == bg::get<1>(r.back()));
        }
        std::set<std::pair<rational, rational>> seen;
        rational signed_area = 0;
        auto const n = r.size() - closed;
        for (std::size_t i = 0; i < n; ++i)
        {
            auto const& a = r[i]; auto const& b = r[(i + 1) % n];
            BOOST_TEST(seen.emplace(bg::get<0>(a), bg::get<1>(a)).second);
            signed_area += rational(bg::get<0>(a)) * bg::get<1>(b) - rational(bg::get<1>(a)) * bg::get<0>(b);
        }
        bool const positive = (bg::point_order<ring>::value == bg::counterclockwise) == shell;
        BOOST_TEST(positive ? signed_area > 0 : signed_area < 0);
    };
    std::size_t id = 2;
    for (auto const& p : result)
    {
        check_ring(bg::exterior_ring(p), true);
        for (auto const& h : bg::interior_rings(p)) check_ring(h, false);
        x.add<false>(p, id); y.add<true>(p, id++);
    }
    x.check(id); y.check(id);
}

template<class A,class B> void run(A const& a,B const& b, unsigned& admitted){
    strategy s;
    bg::strategies::detail::with_operation_sources(s,a,b,[&](auto const& st){admitted+=st.no_common_boundary_edge(a,b);});
    exact_multi out; bg::sym_difference(a,b,out,s); verify(a,b,out);
}
int main(){
    unsigned admitted=0, n=0;
    std::vector<polygon> shapes;
    auto vertex=[](int x,int y){return point(2*x+y,x+2*y);};
    for(int l=0;l<4;++l)for(int r=l+1;r<4;++r)
    for(int b=0;b<4;++b)for(int t=b+1;t<4;++t)
        shapes.emplace_back().outer()={vertex(l,b),vertex(l,t),vertex(r,t),vertex(r,b),vertex(l,b)};
    for(auto const& a:shapes)for(auto const& b:shapes){run(a,b,admitted);++n;}
    std::mt19937 rng(1517);
    auto triangle=[&]{polygon p;do{
        auto& r=p.outer(); r.clear();
        for(int i=0;i<3;++i)r.emplace_back(int(rng()%17)-8,int(rng()%17)-8);
        r.push_back(r.front());bg::correct(p);
    }while(bg::area(p)==0);return p;};
    auto const holes=wkt("POLYGON((-7 -7,-7 7,7 7,7 -7,-7 -7),(-3 -3,3 -3,3 3,-3 3,-3 -3))");
    for(int i=0;i<2048;++i){auto a=triangle(),b=triangle();run(a,b,admitted);++n;if(i<256){run(holes,b,admitted);++n;}}
    std::cout<<"independent edge certificate "<<admitted<<"/"<<n<<" polygon pairs\n";
    return boost::report_errors();
}
