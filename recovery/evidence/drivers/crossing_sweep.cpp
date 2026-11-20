// A deliberately near-degenerate corpus, not part of the public timing suite.
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include "producer_diag.hpp"
#endif
#include <array>
#include <chrono>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
namespace sr = bg::strategies::relate;
using point = bg::model::d2::point_xy<double>;
using line = bg::model::linestring<point>;

template <unsigned N> void sweep()
{
    using root = ec::borrowed_root<double, sr::lean_bounded_ulp_output<N>, line, line>;
    point const o(0,0), d(1,1), u(1,0), v(0,2);
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
    char const* const operation[]={"side","dot","compare_signed","compare_distance"};
    producer_diag::write_header(std::cout);
#endif
    for (unsigned op=0; op!=4; ++op)
    {
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
        producer_diag::reset();
#endif
        for (int m = -768; m <= 768; ++m)
        {
        // P=(1, 1 + (m+1/3)*2^-52), exactly. Each endpoint is binary64.
        line const first{{0, 1 + m*0x1p-52}, {3, 1 + (m+1)*0x1p-52}};
        line const second{{1,0}, {1,2}};
        int const want = m < 0 ? -1 : 1;
            ec::point_cache_scope<double> scope;
            ec::borrowed_context<root> context;
            ec::borrowed_stage<root> stage(context,first,second);
            std::array<typename root::operation_type,2> ops;
            ops[0].seg_id = bg::segment_identifier(0,-1,-1,0);
            ops[1].seg_id = bg::segment_identifier(1,-1,-1,0);
            ec::borrowed_record<root> record;
            record.stage = &stage;
            record.operations = {{&ops[0],&ops[1]}};
            ec::borrowed_point<root> p(record);
            record.numeric.init(p.producers());
            int got;
            if (op==0) got=ec::root_orientation<root>(o,d,p);
            else if (op==1) got=-ec::root_direction<root>(v,d,p);
            else if (op==2) got=ec::root_compare_distance(ec::root_signed_determinant<root>(o,u,p),
                                                          ec::root_signed_determinant<root>(o,u,d));
            else got=ec::root_compare_distance(ec::root_squared_distance<root>(o,p),
                                               ec::root_squared_distance<root>(o,d));
            if (got!=want) {
                std::fprintf(stderr,"N=%u m=%d probe=%u got=%d wanted=%d\n",N,m,op,got,want);
                throw std::logic_error("crossing sweep sign");
            }
        }
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
        producer_diag::dump(std::cout,"targeted",operation[op]);
#endif
    }
}

int main(int argc, char** argv)
{
    if (argc!=2) return 2;
    unsigned n=std::atoi(argv[1]);
    void (*run)();
    switch (n) {
        case 1: run=sweep<1>; break; case 2: run=sweep<2>; break;
        case 4: run=sweep<4>; break; case 8: run=sweep<8>; break;
        case 100: run=sweep<100>; break; default: return 2;
    }
#ifdef BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
    run();
#else
    using clock=std::chrono::steady_clock;
    auto const start=clock::now();
    unsigned reps=0;
    double secs;
    do { run(); ++reps; secs=std::chrono::duration<double>(clock::now()-start).count(); } while(secs<.15);
    std::printf("%.9f\n", secs*1e9/(reps*6148.));
#endif
}
