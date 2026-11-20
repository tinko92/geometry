#include <boost/geometry.hpp>
#include <boost/property_tree/json_parser.hpp>
#ifdef CERTIFIED
#include <boost/geometry/algorithms/buffer_linestring_certified.hpp>
#include <boost/geometry/algorithms/buffer_polygon_certified.hpp>
#endif
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
namespace bg=boost::geometry;
using pt=bg::model::d2::point_xy<double>;
using line=bg::model::linestring<pt>;
using poly=bg::model::polygon<pt>;
#ifdef CERTIFIED
using outpoly=bg::model::polygon<bg::model::d2::point_xy<boost::multiprecision::cpp_rational>>;
#else
using outpoly=poly;
#endif
using multi=bg::model::multi_polygon<outpoly>;
using record=boost::property_tree::ptree;
using timer=std::chrono::steady_clock;
enum class result { success, cannot_certify, invalid_input };
char const* name(result r) { return r==result::success?"success":r==result::cannot_certify?"cannot_certify":"invalid_input"; }

template<class G,class Out> result operation(G const& g,Out& out,double r,double tol) {
#ifdef CERTIFIED
    bg::certified_buffer_status st;
    if constexpr(std::is_same_v<G,line>) st=bg::buffer_linestring_certified(g,out,r,tol);
    else st=bg::buffer_polygon_certified(g,out,std::abs(r),tol,r<0?bg::certified_polygon_buffer_mode::erosion:bg::certified_polygon_buffer_mode::dilation);
    if(st==bg::certified_buffer_status::success)return result::success;
    if(st==bg::certified_buffer_status::cannot_certify)return result::cannot_certify;
    return result::invalid_input;
#else
    (void)tol;
    bg::buffer(g,out,bg::strategy::buffer::distance_symmetric<double>(r),bg::strategy::buffer::side_straight{},
        bg::strategy::buffer::join_round(16),bg::strategy::buffer::end_round(16),bg::strategy::buffer::point_circle(16));
    return result::success;
#endif
}
template<class G> void bench(G g,double r,double tol,unsigned n,record& info) {
#ifdef CERTIFIED
    using Out=std::conditional_t<std::is_same_v<G,line>,outpoly,multi>;
#else
    using Out=multi;
#endif
    if constexpr(std::is_same_v<G,poly>) bg::correct(g);
    if(!n) {
        std::string why;
        info.put("input_valid",bg::is_valid(g,why)); info.put("input_reason",why);
        info.put("input_points",bg::num_points(g));
        Out out;
        std::cerr<<"entered_buffer\n";
        auto const t=timer::now();
        auto const st=operation(g,out,r,tol);
        auto const duration=timer::now()-t;
        std::cerr<<"returned_"<<name(st)<<"\n";
        info.put("first_ns",std::chrono::duration<double,std::nano>(duration).count());
        info.put("status",name(st));
        if(st==result::success) {
            info.put("valid",bg::is_valid(out,why));info.put("reason",why);
            info.put("points",bg::num_points(out));info.put("holes",bg::num_interior_rings(out));
            if constexpr(std::is_same_v<Out,multi>) info.put("components",out.size());
            else info.put("components",bg::is_empty(out)?0:1);
            info.put("area_m2",static_cast<double>(bg::area(out)));
            std::ostringstream w;w<<std::setprecision(17)<<bg::wkt(out);info.put("wkt",w.str());
        }
        return;
    }
    std::uint64_t checksum=0;
    unsigned counts[3]={0,0,0};
    auto const c=std::clock(); auto const t=timer::now();
    for(unsigned i=0;i<n;++i) {
        Out out;
        auto const st=operation(g,out,r,tol);++counts[static_cast<unsigned>(st)];
        checksum+=bg::num_points(out);
    }
    auto const dt=timer::now()-t; auto const dc=std::clock()-c;
    info.put("wall_ns",std::chrono::duration<double,std::nano>(dt).count());
    info.put("cpu_ns",double(dc)*1e9/CLOCKS_PER_SEC); info.put("calls",n);
    info.put("success",counts[0]);info.put("cannot_certify",counts[1]);info.put("invalid_input",counts[2]);
    info.put("checksum",checksum);
}
int main(int argc,char** argv) {
    if(argc!=6)return 2;
    try {
        std::ifstream f(argv[2]);std::string wkt; if(!std::getline(f,wkt))return 3;
        auto run=[&](auto g){bg::read_wkt(wkt,g); record obs;
            bench(g,std::stod(argv[3]),std::stod(argv[4]),std::stoul(argv[5]),obs);
            boost::property_tree::write_json(std::cout,obs,false);};
        if(std::string(argv[1])=="road")run(line{}); else run(poly{});
    } catch(std::exception const& e){std::cerr<<e.what()<<'\n'; return 1;}
}
