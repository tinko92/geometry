// Larger, file-backed public-call workload. All parsing, validation and result hashing are untimed.
#include <boost/geometry.hpp>
#ifdef SIMPLE
#include <boost/geometry/strategies/relate/simple_cartesian.hpp>
#elif defined(EXACT)
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#elif defined(RATIONAL)
#include <boost/multiprecision/cpp_int.hpp>
#ifdef FILTERED
#include <boost/geometry/algorithms/prepare_filtered_decimal_polygon.hpp>
#include <boost/geometry/strategies/relate/cartesian_decimal.hpp>
#endif
#endif
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace bg=boost::geometry;
template<class C> using pt=bg::model::d2::point_xy<C>;
template<class C> using poly=bg::model::polygon<pt<C>>;
#ifdef RATIONAL
using Q=boost::multiprecision::cpp_rational;
using C=Q;
#ifdef FILTERED
using strategy=bg::strategies::relate::cartesian_decimal;
#else
using strategy=bg::strategies::relate::cartesian<Q>;
#endif
#else
using C=double;
#ifdef EXACT
#if defined(NEAREST)
using policy=bg::strategies::relate::nearest_output;
#elif defined(LEAN)
using policy=bg::strategies::relate::lean_bounded_ulp_output<STEPS>;
#elif defined(STEPS)
using policy=bg::strategies::relate::bounded_ulp_output<STEPS>;
#else
using policy=bg::strategies::relate::approximate_output;
#endif
using strategy=bg::strategies::relate::cartesian_exact<double,policy>;
#elif defined(SIMPLE)
using strategy=bg::strategies::relate::simple_cartesian;
#else
using strategy=bg::strategies::relate::cartesian<>;
#endif
#endif
using polygons=bg::model::multi_polygon<poly<C>>;
using lines=bg::model::multi_linestring<bg::model::linestring<pt<C>>>;
using words=std::vector<std::uint64_t>;
std::string op,phase,folder;
int grid;
volatile std::size_t sink=0;
template<class T> void emit(std::string const& key,T const& v){std::cout<<key<<'\t'<<v<<'\n';}
std::uint64_t mix(std::uint64_t h,std::uint64_t x){return (h^x)*1099511628211ULL;}
std::uint64_t text_hash(std::string const& s){std::uint64_t h=1;for(unsigned char c:s)h=mix(h,c);return h;}
std::uint64_t group(words v){std::sort(v.begin(),v.end());std::uint64_t h=v.size();for(auto x:v)h=mix(h,x);return h;}
template<class T> std::uint64_t coord(T const& x,bool output){
    if constexpr(std::is_same_v<T,double>) {
        if(!std::isfinite(x))throw std::domain_error("nonfinite coordinate");
        if(output)return static_cast<std::uint64_t>(std::llround(std::ldexp(x,grid)));
        std::uint64_t z;std::memcpy(&z,&x,8);return z;
    } else return text_hash(x.str());
}
template<class R> std::uint64_t path(R const& r,bool output) {
    words e;std::uint64_t last=0,h=r.size();bool first=true;
    for(auto const& p:r) {
        auto k=mix(coord(bg::get<0>(p),output),coord(bg::get<1>(p),output));
        if(output && !first)e.push_back(mix(std::min(k,last),std::max(k,last)));
        first=false;last=k;h=mix(h,k);
    }
    return output?mix(r.size(),r.size()==1?last:group(e)):h;
}
template<class P> std::uint64_t polygon_key(P const& p,bool out) {
    words holes;auto h=path(p.outer(),out);
    for(auto const& r:p.inners()){auto k=path(r,out);holes.push_back(k);if(!out)h=mix(h,k);}
    return mix(h,out?group(holes):holes.size());
}
std::uint64_t line_key(lines const& ls,bool out) {
    words cs;std::uint64_t h=ls.size();
    for(auto const& l:ls){auto k=path(l,out);cs.push_back(k);h=mix(h,k);}
    return out?group(cs):h;
}
std::string file(char const* name) {
    std::ifstream f(folder+"/"+name); std::string s;if(!std::getline(f,s))throw std::runtime_error("missing input");return s;
}
template<class G> G read(std::string const& text){G g;bg::read_wkt(text,g);return g;}
auto polygon(std::string const& s) {
#ifdef FILTERED
    return bg::prepare_filtered_decimal_polygon(read<poly<bg::exact_coordinate<std::string>>>(s));
#elif defined(RATIONAL)
    // Exact decimal token to rational fraction, never via floating point.
    static std::regex const decimal("(-?)([0-9]+)\\.([0-9]+)");
    std::string dst;std::size_t pos=0;
    for(std::sregex_iterator i(s.begin(),s.end(),decimal),end;i!=end;++i) {
        auto m=*i;dst+=s.substr(pos,m.position()-pos);auto n=m[2].str()+m[3].str();
        auto j=n.find_first_not_of('0');
        dst+=m[1].str()+(j==std::string::npos?"0":n.substr(j))+"/1"+std::string(m[3].length(),'0');pos=m.position()+m.length();
    }
    return read<poly<C>>(dst+s.substr(pos));
#else
    return read<poly<C>>(s);
#endif
}
template<class G> bool valid(char const* key,G const& g) {
    std::string why;auto ok=bg::is_valid(g,why);emit(key,ok);emit(std::string(key)+"_reason",why);return ok;
}
template<class F> BOOST_NOINLINE std::size_t invoke(F& f){return f();}
template<class F,class Obs> void bench(F call,Obs observe) {
    sink=invoke(call);
    if(!observe())return;
    if(phase=="check"){emit("status","ok");return;}
    auto c=std::clock();auto start=std::chrono::steady_clock::now();std::uint64_t n=0;double wall=0;
    do{sink=invoke(call);++n;wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();}while(wall<.15);
    double cpu=double(std::clock()-c)/CLOCKS_PER_SEC;
    emit("calls",n);emit("seconds",wall);emit("cpu_seconds",cpu);emit("ns_per_call",wall*1e9/n);emit("status","ok");
}
int main(int argc,char**argv) {
    try {
        if(argc!=6)return 2;
        folder=argv[1];std::string kind=argv[2];op=argv[3];phase=argv[4];grid=std::stoi(argv[5]);
        std::cout<<std::setprecision(17);strategy s;
        auto b=polygon(file("b.wkt"));
        if(phase=="check")valid("b_source_valid",b);
        bg::correct(b);
        auto setup=[&](auto const& a,std::uint64_t key){
            emit("a_points",bg::num_points(a));emit("b_points",bg::num_points(b));
            emit("input_control",mix(key,polygon_key(b,false)));
            if(phase!="check")return true;
            bool va=valid("a_valid",a),vb=valid("b_valid",b);
            if(!va||!vb){emit("status","invalid-input");return false;}return true;
        };
        if(kind=="polygon"){
            auto a=polygon(file("a.wkt"));if(phase=="check")valid("a_source_valid",a);bg::correct(a);
            if(!setup(a,polygon_key(a,false)))return 0;
            polygons out;std::string matrix;
            bench([&]()->std::size_t {
                out.clear();
                if(op=="intersection")bg::intersection(a,b,out,s);
                else if(op=="union")bg::union_(a,b,out,s);
                else if(op=="difference")bg::difference(a,b,out,s);
                else if(op=="sym_difference")bg::sym_difference(a,b,out,s);
                else if(op=="relate_aa"){matrix=bg::relation(a,b,s).str();return matrix[0];}
                else throw std::invalid_argument(op);
                return out.size();
            },[&](){
                if(op=="relate_aa"){emit("output_control",matrix);return true;}
                words keys;for(auto const& p:out)keys.push_back(polygon_key(p,true));
                emit("output_control",group(keys));emit("components",out.size());emit("points",bg::num_points(out));
                if(phase=="check" && !valid("output_valid",out)){emit("status","invalid-output");return false;}return true;
            });
        }
#ifndef RATIONAL
        else if(kind=="track"){
            auto a=read<lines>(file("a.wkt"));if(phase=="check")valid("a_source_valid",a);
            if(!setup(a,line_key(a,false)))return 0;
            lines out;std::string matrix;
            bench([&]()->std::size_t {
                out.clear();if(op=="clip"){bg::intersection(a,b,out,s);return out.size();}
                if(op!="relate_la")throw std::invalid_argument(op);
                matrix=bg::relation(a,b,s).str();return matrix[0];
            },[&](){
                if(op=="relate_la"){emit("output_control",matrix);return true;}
                emit("output_control",line_key(out,true));emit("components",out.size());emit("points",bg::num_points(out));
                if(phase=="check" && !valid("output_valid",out)){emit("status","invalid-output");return false;}return true;
            });
        }
#endif
        else throw std::invalid_argument(kind);
        return 0;
    } catch(std::exception const&e){std::cerr<<e.what()<<'\n';return 1;}
}
