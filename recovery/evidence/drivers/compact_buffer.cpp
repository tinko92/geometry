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
#ifdef CAND
#include <boost/geometry/algorithms/buffer_point_certified.hpp>
#include <boost/geometry/algorithms/buffer_linestring_certified.hpp>
#include <boost/geometry/algorithms/buffer_polygon_certified.hpp>
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
template<class G> std::string input_wkt(G const& g) {
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
        else return input_wkt(p);
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
template<class G> void buffer(G const& g) {
#ifdef CAND
    poly<point<q>> one; multi<point<q>> many;
#else
    multi<point<double>> many;
#endif
    auto call=[&]()->std::size_t{
#ifdef CAND
        one={}; many.clear(); bg::certified_buffer_status status;
        if constexpr(std::is_same_v<G,point<double>>) status=bg::buffer_point_certified(g,one,2,.125);
        else if constexpr(std::is_same_v<G,line<point<double>>>) status=bg::buffer_linestring_certified(g,one,2,.125);
        else status=bg::buffer_polygon_certified(g,many,2,.125);
        if(status!=bg::certified_buffer_status::success) throw std::runtime_error("buffer not certified");
        return one.outer().size()+many.size();
#else
        many.clear();
        bg::buffer(g,many,bg::strategy::buffer::distance_symmetric<double>(2),bg::strategy::buffer::side_straight{},
                   bg::strategy::buffer::join_round(16),bg::strategy::buffer::end_round(16),bg::strategy::buffer::point_circle(16));
        return many.size();
#endif
    };
    sample(1,input_wkt(g),call,[&]{
        std::string key=poly_key(many);
#ifdef CAND
        key+=ring_key(one.outer());
#endif
        return std::to_string(hash(key));
    });
}
void run(){
    if(dataset=="point") buffer(point<double>{3,4});
    else if(dataset=="bent-line") buffer(line<point<double>>{{-6,-8},{0,0},{9,0},{12,4}});
    else buffer(read<poly<point<double>>>("POLYGON((0 0,0 30,40 30,40 0,0 0),(8 8,26 8,26 20,8 20,8 8))"));
}
int main(int argc,char**argv) {
    try { if(argc!=4) return 2; dataset=argv[1]; op=argv[2]; rep=std::stoi(argv[3]); run(); return 0; }
    catch(std::exception const& e){ std::cerr<<e.what()<<'\n'; return 1; }
}
