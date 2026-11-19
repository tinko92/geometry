// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
// Always build with the frozen before headers, never a baseline's validator.
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
namespace bg = boost::geometry;
using Q = boost::multiprecision::cpp_rational;
using P = bg::model::d2::point_xy<Q>;
using L = bg::model::linestring<P>;
using ML = bg::model::multi_linestring<L>;
using Poly = bg::model::polygon<P>;
using MP = bg::model::multi_polygon<Poly>;
using Points = bg::model::multi_point<P>;
using Tree = boost::property_tree::ptree;
bg::strategies::relate::cartesian<Q> validation;
template <typename G> G read(std::string const& s) { G g; bg::read_wkt(s,g); return g; }
MP polygons(std::string const& s)
{ if(s.find("MULTIPOLYGON")==0) return read<MP>(s); MP m; m.push_back(read<Poly>(s)); return m; }
ML lines(std::string const& s)
{ if(s.find("MULTILINESTRING")==0) return read<ML>(s); ML m; m.push_back(read<L>(s)); return m; }
template <typename G> void summarize(G const& g,Tree& t)
{
    std::string reason;
    t.put("rational_valid",bg::is_valid(g,reason,validation)); t.put("rational_reason",reason);
    t.put("empty",bg::is_empty(g));
}
void summarize_areal(MP const& g,Tree& t)
{
    summarize(g,t); t.put("components",g.size()); t.put("holes",bg::num_interior_rings(g));
    Q area=bg::area(g,validation); t.put("area",area); t.put("area_double",area.convert_to<double>());
}
template <typename A,typename B> bool equal_set(A const& a,B const& b)
{
    if(bg::is_empty(a)||bg::is_empty(b)) return bg::is_empty(a)&&bg::is_empty(b);
    return bg::relate(a,b,bg::de9im::mask("T*F**FFF*"),validation);
}
Q x(P const& p) { return bg::get<0>(p); }
Q y(P const& p) { return bg::get<1>(p); }
P sub(P const& a,P const& b) { return P(x(a)-x(b),y(a)-y(b)); }
Q cross(P const& a,P const& b) { return x(a)*y(b)-y(a)*x(b); }
bool same(P const& a,P const& b) { return x(a)==x(b)&&y(a)==y(b); }
bool on(P const& p,P const& a,P const& b)
{
    return cross(sub(p,a),sub(b,a))==0 && x(p)>=(std::min)(x(a),x(b))
        && x(p)<=(std::max)(x(a),x(b)) && y(p)>=(std::min)(y(a),y(b)) && y(p)<=(std::max)(y(a),y(b));
}
// Independent rational segment solve; not BG intersection/turn generation.
std::vector<P> intersections(P const& a,P const& b,P const& c,P const& d)
{
    P r=sub(b,a),s=sub(d,c); Q den=cross(r,s); std::vector<P> out;
    if(den==0)
    {
        for(auto p:{a,b,c,d}) if(on(p,a,b)&&on(p,c,d)
            &&std::none_of(out.begin(),out.end(),[&](P const& q){return same(p,q);})) out.push_back(p);
    }
    else
    {
        Q t=cross(sub(c,a),s)/den,u=cross(sub(c,a),r)/den;
        if(t>=0&&t<=1&&u>=0&&u<=1) out.emplace_back(x(a)+t*x(r),y(a)+t*y(r));
    }
    return out;
}
template <typename Ring> int inside_ring(P const& p,Ring const& ring)
{
    bool inside=false;
    for(std::size_t i=1;i<ring.size();++i)
    {
        P const& a=ring[i-1]; P const& b=ring[i];
        if(on(p,a,b)) return 0;
        if((y(a)>y(p))!=(y(b)>y(p)))
        { Q ix=x(a)+(y(p)-y(a))*(x(b)-x(a))/(y(b)-y(a)); if(ix>x(p)) inside=!inside; }
    }
    return inside?1:-1;
}
bool inside(P const& p,MP const& polygons)
{
    for(auto const& poly:polygons)
    {
        int outer=inside_ring(p,poly.outer()); if(outer==0) return true; if(outer<0) continue;
        bool hole=false;
        for(auto const& ring:poly.inners()) { int in=inside_ring(p,ring); if(in==0) return true; if(in>0) hole=true; }
        if(!hole) return true;
    }
    return false;
}
ML clip_oracle(ML const& input,MP const& polygons)
{
    ML out;
    for(auto const& line:input) for(std::size_t i=1;i<line.size();++i)
    {
        P a=line[i-1],b=line[i],delta=sub(b,a);
        if(same(a,b)) continue;
        std::vector<Q> cuts{Q(0),Q(1)};
        auto ring=[&](auto const& r)
        {
            for(std::size_t j=1;j<r.size();++j) for(auto p:intersections(a,b,r[j-1],r[j]))
                cuts.push_back(x(delta)!=0?(x(p)-x(a))/x(delta):(y(p)-y(a))/y(delta));
        };
        for(auto const& p:polygons) { ring(p.outer()); for(auto const& h:p.inners()) ring(h); }
        std::sort(cuts.begin(),cuts.end()); cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        auto point=[&](Q t){return P(x(a)+t*x(delta),y(a)+t*y(delta));};
        for(std::size_t j=1;j<cuts.size();++j)
            if(inside(point((cuts[j-1]+cuts[j])/2),polygons)) out.push_back(L{point(cuts[j-1]),point(cuts[j])});
    }
    return out;
}
void evaluate(Tree& r)
{
    auto& geoms=r.get_child("geometries");
    for(auto& entry:geoms)
    {
        std::string s=entry.second.get<std::string>("wkt");
        if(s.find("MULTIPOLYGON")==0||s.find("POLYGON")==0) summarize_areal(polygons(s),entry.second);
        else if(s.find("MULTILINESTRING")==0||s.find("LINESTRING")==0)
        { auto g=lines(s); summarize(g,entry.second); entry.second.put("components",g.size()); }
        else if(s.find("MULTIPOINT")==0)
        { auto g=read<Points>(s); summarize(g,entry.second); entry.second.put("components",g.size()); }
        else if(s.find("BOX")==0) summarize(read<bg::model::box<P>>(s),entry.second);
    }
    auto wkt=[&](char const* key){return geoms.get<std::string>(std::string(key)+".wkt");};
    auto valid=[&](char const* key){return geoms.get<bool>(std::string(key)+".rational_valid",false);};
    unsigned issue=r.get<unsigned>("issue");
    if(geoms.count("result")&&geoms.count("reverse")&&valid("result")&&valid("reverse"))
    {
        if(wkt("result").find("MULTIPOLYGON")==0) r.put("checks.reverse_set",equal_set(polygons(wkt("result")),polygons(wkt("reverse"))));
        else r.put("checks.reverse_set",equal_set(lines(wkt("result")),lines(wkt("reverse"))));
    }
    if(issue==528)
    {
        auto a=polygons(wkt("a")); auto b=lines(wkt("b")); std::vector<P> expected;
        for(auto const& poly:a) for(std::size_t i=1;i<poly.outer().size();++i)
            for(auto p:intersections(poly.outer()[i-1],poly.outer()[i],b[0][0],b[0][1]))
                if(std::none_of(expected.begin(),expected.end(),[&](P const& q){return same(p,q);})) expected.push_back(p);
        auto actual=read<Points>(wkt("result")); bool match=actual.size()==expected.size();
        for(auto p:actual) match=match&&std::any_of(expected.begin(),expected.end(),[&](P const& q){return same(p,q);});
        r.put("checks.boundary_oracle",match); r.put("reference.points",expected.size());
    }
    if(issue==440)
    {
        // a and b are lines, not outputs. Compare all eight direction/order calls.
        bool consistent=true,reference=true; Points previous;
        for(auto const& entry:geoms) if(entry.first.find("result")==0)
        {
            auto actual=read<Points>(entry.second.get<std::string>("wkt"));
            if(!previous.empty()) consistent=consistent&&actual.size()==previous.size()&&same(actual[0],previous[0]);
            reference=reference&&actual.size()==1&&x(actual[0])==4700000&&y(actual[0])==Q(-2665034)-Q(1351)*10043/24634;
            previous=actual;
        }
        r.put("checks.direction_consistent",consistent); r.put("checks.rational_point",reference);
    }
    if(issue==510||issue==511||issue==1253||issue==656||issue==858)
    {
        auto expected=clip_oracle(lines(wkt("a")),polygons(wkt("b")));
        auto actual=lines(wkt("result"));
        r.put("reference.linear_components",expected.size());
        r.put("checks.linear_oracle",equal_set(actual,expected));
        if(issue==858)
        {
            Q shared=Q(9)-Q(2)*y(lines(wkt("a"))[0][0])/9;
            r.put("reference.shared_x",shared);
            if(!actual.empty()&&!actual[0].empty())
            {
                Q left=x(actual[0].front()),right=x(actual[0].back());
                if(left>right) std::swap(left,right);
                r.put("result_shared_x",std::string(r.get<std::string>("case")).back()=='0'?right:left);
            }
        }
    }
    if(issue==1232)
    {
        auto actual=lines(wkt("result"));
        r.put("checks.head",bg::covered_by(P(5,-8),actual,validation));
        r.put("checks.tail",bg::covered_by(P(7,4),actual,validation));
    }
    if(issue==1044)
    {
        // The upper chain can run in either orientation after bg::correct.
        auto b=polygons(wkt("b"));
        Q bottom=y(b[0].outer()[0]);
        for(auto p:b[0].outer()) bottom=(std::min)(bottom,y(p));
        Q xmin(-5.8),xmax(-5.5),area=0;
        auto const& ring=b[0].outer();
        for(std::size_t i=1;i<ring.size();++i)
        {
            P a=ring[i-1],b=ring[i];
            if(y(a)<Q(5.1)||y(b)<Q(5.1)||x(a)==x(b)) continue;
            if(x(a)>x(b)) std::swap(a,b);
            Q lo=(std::max)(xmin,x(a)),hi=(std::min)(xmax,x(b));
            if(lo>=hi) continue;
            Q slope=(y(b)-y(a))/(x(b)-x(a));
            Q yl=y(a)+(lo-x(a))*slope,yh=y(a)+(hi-x(a))*slope;
            area+=(hi-lo)*((yl+yh)/2-bottom);
        }
        r.put("reference.top_chain_area",area); r.put("reference.top_chain_area_double",area.convert_to<double>());
        r.put("checks.top_chain_area",bg::area(polygons(wkt("result")),validation)==area);
    }
    if(issue==1236&&valid("result")&&valid("reference_container"))
        r.put("checks.rational_within",bg::within(polygons(wkt("result")),polygons(wkt("reference_container")),validation));
    if(issue==588&&valid("result"))
        r.put("checks.rational_within",bg::within(polygons(wkt("result")),polygons(wkt("input0")),validation));
    if(issue==900&&valid("result")&&valid("identity_reference"))
        r.put("checks.pipeline_identity",equal_set(polygons(wkt("result")),polygons(wkt("identity_reference"))));
}
int main()
{
    Tree r;
    try
    {
        boost::property_tree::read_json(std::cin,r);
        if(r.count("geometries")) evaluate(r);
        r.put("evaluation","ok");
    }
    catch(std::exception const& e) { r.put("evaluation","exception"); r.put("evaluation_exception",e.what()); }
    boost::property_tree::write_json(std::cout,r,false);
}
