// External public-call microbenchmarks. Preparation and output fingerprinting are untimed.
#include <boost/geometry.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(CAND) && defined(DECIMAL)
#include <boost/geometry/algorithms/prepare_filtered_decimal_polygon.hpp>
#include <boost/geometry/strategies/relate/cartesian_decimal.hpp>
#endif
namespace bg = boost::geometry;
using q = boost::multiprecision::cpp_rational;
template<class T> using point = bg::model::d2::point_xy<T>;
template<class P> using poly = bg::model::polygon<P>;
template<class P> using multi = bg::model::multi_polygon<poly<P>>;
template<class P> using line = bg::model::linestring<P>;
std::string dataset, op;
int rep;
volatile std::size_t sink = 0;
std::uint64_t hash(std::string const& str) {
    std::uint64_t h=14695981039346656037ULL;
    for(unsigned char c:str) h=(h^c)*1099511628211ULL;
    return h;
}
template<class G> G read(std::string const& wkt) { G g; bg::read_wkt(wkt,g); return g; }
template<class G> std::string wkt(G const& g) {
    std::ostringstream s; s<<std::setprecision(17)<<bg::wkt(g); return s.str();
}
std::string packed(std::vector<std::string> x) {
    std::sort(x.begin(),x.end()); std::string s;
    for(auto& v:x) s+=std::to_string(v.size())+":"+v;
    return s;
}
template<class Ring> std::string ring_key(Ring const& r) {
    std::vector<std::string> edges;
    auto key=[](auto const& p) {
        using P=std::decay_t<decltype(p)>;
        if constexpr(std::is_same_v<bg::coordinate_type_t<P>,double>)
            return std::to_string(std::llround(bg::get<0>(p)*65536))+":"
                 + std::to_string(std::llround(bg::get<1>(p)*65536));
        else return wkt(p);
    };
    for(std::size_t i=1;i<r.size();++i) {
        auto a=key(r[i-1]), b=key(r[i]); if(b<a) std::swap(a,b);
        edges.push_back(a+":"+b);
    }
    return packed(edges);
}
template<class Multi> std::string poly_key(Multi const& m) {
    std::vector<std::string> polys;
    for(auto& p:m) {
        std::vector<std::string> holes;
        for(auto& r:p.inners()) holes.push_back(ring_key(r));
        polys.push_back(ring_key(p.outer())+"|"+packed(holes));
    }
    return packed(polys);
}
template<class Call, class Observe>
void sample(int size, std::string const& inputs, Call call, Observe observe) {
    sink=call(); auto initial=observe();
    std::uint64_t calls=0; auto c=std::clock(); auto start=std::chrono::steady_clock::now();
    double s=0;
    do {
        sink=call(); calls+=size;
        s=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    } while(s<.13);
    double cpu=double(std::clock()-c)/CLOCKS_PER_SEC;
    if(observe()!=initial) throw std::logic_error("unstable observation");
    std::cout<<"workload,operation,repeat,corpus_size,calls,seconds,cpu_seconds,ns_per_call,input_fingerprint,output_control\n"
        <<std::setprecision(16)<<dataset<<','<<op<<','<<rep<<','<<size<<','<<calls<<','<<s<<','<<cpu<<','
        <<s*1e9/calls<<','<<hash(inputs)<<','<<initial<<'\n';
}
#if defined(DECIMAL)
#include "decimal-fixtures.hpp"
#ifdef CAND
auto prepare(std::string const& s) {
    return bg::prepare_filtered_decimal_polygon(read<poly<point<bg::exact_coordinate<std::string>>>>(s));
}
using strategy=bg::strategies::relate::cartesian_decimal;
#else
auto prepare(std::string const& s) {
    // Multiprecision rational extraction accepts integer fractions, not decimal tokens.
    std::regex token("(-?)([0-9]+)\\.([0-9]+)");
    std::string expanded; std::size_t pos=0;
    for(std::sregex_iterator it(s.begin(),s.end(),token), end; it!=end; ++it) {
        auto m=*it; expanded+=s.substr(pos,m.position()-pos);
        auto n=m[2].str()+m[3].str(); auto first=n.find_first_not_of('0');
        expanded+=m[1].str()+(first==std::string::npos?"0":n.substr(first))+"/1"+std::string(m[3].length(),'0');
        pos=m.position()+m.length();
    }
    return read<poly<point<q>>>(expanded+s.substr(pos));
}
using strategy=bg::strategies::relate::cartesian<q>;
#endif
void run() {
    using g=decltype(prepare("")); std::vector<std::pair<g,g>> cases; std::string inputs;
    for(auto const& r:input_wkts) if(r.category==dataset) {
        cases.emplace_back(prepare(r.a),prepare(r.b));
        inputs+=wkt(cases.back().first)+"|"+wkt(cases.back().second)+";";
    }
    if(cases.size()!=3) throw std::logic_error("expected literal three-input corpus");
    std::array<multi<point<q>>,3> out;
    std::array<std::string,3> matrices;
    strategy s;
    auto call=[&]() {
        std::size_t count=0;
        for(unsigned i=0;i<3;++i) {
            auto const& a=cases[i].first; auto const& b=cases[i].second; auto& m=out[i]; m.clear();
            if(op=="relate_aa") { matrices[i]=bg::relation(a,b,s).str(); count+=matrices[i][0]; }
            else if(op=="union") bg::union_(a,b,m,s);
            else if(op=="intersection") bg::intersection(a,b,m,s);
            else if(op=="difference") bg::difference(a,b,m,s);
            else if(op=="sym_difference") bg::sym_difference(a,b,m,s);
            else throw std::logic_error("unknown op");
            count+=m.size();
        }
        return count;
    };
    sample(3, inputs, call, [&] {
        if(op=="relate_aa") return matrices[0]+":"+matrices[1]+":"+matrices[2];
        std::string sig;
        for(auto const& m:out) sig+=poly_key(m)+";";
        return std::to_string(hash(sig));
    });
}
#endif
int main(int argc,char**argv) {
    try { if(argc!=4) return 2; dataset=argv[1]; op=argv[2]; rep=std::stoi(argv[3]); run(); return 0; }
    catch(std::exception const& e){ std::cerr<<e.what()<<'\n'; return 1; }
}
