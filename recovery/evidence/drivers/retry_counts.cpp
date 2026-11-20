// Diagnostic only, not part of a timed call. Include the exact benchmark corpus.
#define main benchmark_main
#include "simple_cartesian_benchmark.cpp"
#undef main
#include <boost/geometry/strategies/cartesian/exact/construction_rounding.hpp>

namespace ec = bg::strategies::detail::exact_cartesian;
using edge = std::pair<point,point>;
using pair_input = std::array<double,8>;
template <typename R> void edges(std::vector<edge>& e, R const& r)
{
    for (std::size_t i=1; i<r.size(); ++i) e.emplace_back(r[i-1],r[i]);
}
std::vector<edge> edges(polygon const& p)
{
    std::vector<edge> e; edges(e,p.outer());
    for (auto const& r:p.inners()) edges(e,r);
    return e;
}
std::vector<edge> edges(lines const& m)
{
    std::vector<edge> e;
    for (auto const& r:m) edges(e,r);
    return e;
}
template <unsigned N> struct totals {
    unsigned axes = 0, fast_off = 0, fast_on = 0, fallback_off = 0, fallback_on = 0;
    void add(pair_input const& s)
    {
        ec::forward_error::point_approximation p, q;
        auto const a=ec::construction_rounding::detail::bounded_intersection<N,false,true,true,false,false>(s,p);
        auto const b=ec::construction_rounding::detail::bounded_intersection<N,false,true,true,false,true>(s,q);
        for (int i=0; i<2; ++i) {
            if (a[i] && (!b[i] || p.value[i]!=q.value[i])) throw std::logic_error("Lost first certificate");
            ++axes; fast_off += a[i]; fast_on += b[i];
        }
        fallback_off += !a[0] || !a[1]; fallback_on += !b[0] || !b[1];
    }
    void print(char const* w) const
    { std::cout<<w<<','<<N<<','<<axes<<','<<fast_off<<','<<fast_on<<','<<fallback_off<<','<<fallback_on<<'\n'; }
};
int main()
{
    using side=bg::strategy::side::side_robust<double,bg::strategy::side::fp_equals_policy,3>;
    std::cout<<"workload,N,axes,fast_off,fast_on,points_requiring_refinement_off,points_requiring_refinement_on\n";
    for (char const* w:{"degenerate","general","decimal","projected","projected-local"}) {
        totals<8> eight; totals<100> hundred;
        auto count=[&](auto const& first,auto const& second) {
            for(auto const& [a,b]:first) for(auto const& [c,d]:second) {
                if (side::apply(a,b,c)*side::apply(a,b,d)>=0 || side::apply(c,d,a)*side::apply(c,d,b)>=0) continue;
                pair_input s{{a.x(),a.y(),b.x(),b.y(),c.x(),c.y(),d.x(),d.y()}};
                eight.add(s); hundred.add(s);
            }
        };
        for(auto const& g:make_inputs(w)) {
            count(edges(g.a),edges(g.b)); count(edges(g.la),edges(g.a)); count(edges(g.la),edges(g.lb));
        }
        eight.print(w); hundred.print(w);
    }
}
