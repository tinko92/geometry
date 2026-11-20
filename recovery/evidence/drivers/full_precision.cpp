// Replay the public benchmark inputs, comparing full outputs without quantization.
#define main frozen_benchmark_main
#include "simple_cartesian_benchmark.cpp"
#undef main
#ifdef EXPECT_CANONICAL_NEAREST
static_assert(std::is_same<typename strategy::output_policy, bg::strategies::relate::nearest_output>::value);
#endif

int main()
{
    auto path = [](auto const& r) {
        std::cout << r.size() << ':';
        for (auto const& p : r)
            std::cout << std::hexfloat << bg::get<0>(p) << ',' << bg::get<1>(p) << ';';
    };
    for (auto w : {"degenerate","general","decimal","projected","projected-local"})
    for (auto const& in : make_inputs(w))
    for (auto op : operation_names)
    {
        auto id = operation_code(op);
        polygons p; lines l; std::cout << w << ' ' << op << ' ';
        if (id == operation_id::relate_aa) std::cout << bg::relation(in.a,in.b,strategy{}).str();
        else if (id == operation_id::relate_la) std::cout << bg::relation(in.la,in.a,strategy{}).str();
        else if (id == operation_id::relate_ll) std::cout << bg::relation(in.la,in.lb,strategy{}).str();
        else {
            std::cout << execute(id,in,p,l) << ' ';
            for (auto const& g : p) {
                std::cout << 'P'; path(g.outer());
                for (auto const& h : g.inners()) { std::cout << 'H'; path(h); }
            }
            for (auto const& line : l) { std::cout << 'L'; path(line); }
        }
        std::cout << '\n';
    }
}
