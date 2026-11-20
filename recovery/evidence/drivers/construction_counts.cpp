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
template <unsigned N> unsigned accepted(pair_input const& s)
{
    ec::forward_error::point_approximation p;
    auto const yes=ec::construction_rounding::detail::bounded_intersection<N,false>(s,p);
    return unsigned(yes[0])+unsigned(yes[1]);
}
int main()
{
    using side=bg::strategy::side::side_robust<double,bg::strategy::side::fp_equals_policy,3>;
    std::cout<<"workload,axes,N1,N2,N4,N8,N100,pass4fail8_axes\n";
    for (auto w:{"degenerate","general","decimal","projected","projected-local"})
    {
        std::array<unsigned,7> n{};
        auto count=[&](std::vector<edge> const& first,std::vector<edge> const& second){
            for(auto const& x:first) for(auto const& y:second) {
                auto const& [a,b]=x; auto const& [c,d]=y;
                if(side::apply(a,b,c)*side::apply(a,b,d)>=0 || side::apply(c,d,a)*side::apply(c,d,b)>=0)
                    continue;
                pair_input s{{a.x(),a.y(),b.x(),b.y(),c.x(),c.y(),d.x(),d.y()}};
                n[0]+=2; n[1]+=accepted<1>(s); n[2]+=accepted<2>(s);
                n[3]+=accepted<4>(s); n[4]+=accepted<8>(s); n[5]+=accepted<100>(s);
                ec::forward_error::point_approximation p;
                auto const four=ec::construction_rounding::detail::bounded_intersection<4,false>(s,p);
                auto const eight=ec::construction_rounding::detail::bounded_intersection<8,false>(s,p);
                for(int i=0;i<2;++i) n[6]+=four[i]&&!eight[i];
            }
        };
        for(auto const& g:make_inputs(w)) {
            count(edges(g.a),edges(g.b)); count(edges(g.la),edges(g.a)); count(edges(g.la),edges(g.lb));
        }
        std::cout<<w; for(auto x:n) std::cout<<','<<x; std::cout<<'\n';
    }
}
