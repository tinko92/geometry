// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
// COMPARISON_GROUP: 1 (areal), 2 (API), 3 (pipeline), 4 (downloaded API).
// Group 1 is split by COMPARISON_COORDINATE to bound compiler memory.
// COMPARISON_EXACT selects the frozen exact strategy; otherwise use defaults.
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/relate/cartesian.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/property_tree/json_parser.hpp>
#ifdef COMPARISON_EXACT
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#endif
#include "exact_additional_issue_cases.hpp"
#include <algorithm>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <tuple>

namespace bg = boost::geometry;
namespace cases = exact_additional_issue_cases;
using rational = boost::multiprecision::cpp_rational;
using record = boost::property_tree::ptree;
template <typename P, bool Clockwise = true> struct models
{
    using point = P;
    using polygon = bg::model::polygon<P, Clockwise>;
    using polygons = bg::model::multi_polygon<polygon>;
    using line = bg::model::linestring<P>;
    using lines = bg::model::multi_linestring<line>;
    using points = bg::model::multi_point<P>;
    using box = bg::model::box<P>;
};
template <typename C, bool Clockwise = true> using native = models<bg::model::d2::point_xy<C>, Clockwise>;
using rational_models = native<rational>;
#ifdef COMPARISON_EXACT
namespace ec = bg::strategies::detail::exact_cartesian;
template <typename C, bool Clockwise = true> using retained = models<ec::exact_point<C>, Clockwise>;
#ifndef COMPARISON_EXACT_STRATEGY
#if defined(COMPARISON_NEAREST_OUTPUT)
#ifndef BOOST_GEOMETRY_EXACT_TYPED_CONTEXT
#error COMPARISON_NEAREST_OUTPUT requires the redesigned exact strategy
#endif
#define COMPARISON_EXACT_STRATEGY(C) bg::strategies::relate::cartesian_exact<C, bg::strategies::relate::nearest_output>
#else
#define COMPARISON_EXACT_STRATEGY(C) bg::strategies::relate::cartesian_exact<C>
#endif
#endif

#ifdef BOOST_GEOMETRY_EXACT_TYPED_CONTEXT
template <typename C, typename Action>
void comparison_in_context(Action const& action, record& r)
{
    ec::point_cache_scope<C> operation;
    action(r);
}
#endif
#endif

template <typename C> struct api
{
#ifdef COMPARISON_EXACT
    COMPARISON_EXACT_STRATEGY(C) strategy;
#define COMPARISON_STRATEGY , strategy
#else
#define COMPARISON_STRATEGY
#endif
    template <typename A, typename B, typename R>
    void overlay(std::string const& op, A const& a, B const& b, R& r) const
    {
        if (op == "intersection") bg::intersection(a, b, r COMPARISON_STRATEGY);
        else if (op == "union") bg::union_(a, b, r COMPARISON_STRATEGY);
        else if (op == "difference") bg::difference(a, b, r COMPARISON_STRATEGY);
        else if (op == "sym_difference") bg::sym_difference(a, b, r COMPARISON_STRATEGY);
    }
    template <typename A, typename B, typename R>
    void intersection(A const& a, B const& b, R& r) const
    { bg::intersection(a, b, r COMPARISON_STRATEGY); }
    template <typename A, typename B, typename R>
    void difference(A const& a, B const& b, R& r) const
    { bg::difference(a, b, r COMPARISON_STRATEGY); }
    template <typename A, typename B> std::string relation(A const& a, B const& b) const
    { return bg::relation(a, b COMPARISON_STRATEGY).str(); }
    template <typename A, typename B> bool within(A const& a, B const& b) const
    { return bg::within(a, b COMPARISON_STRATEGY); }
    template <typename A, typename B> bool covered(A const& a, B const& b) const
    { return bg::covered_by(a, b COMPARISON_STRATEGY); }
    template <typename A, typename B> bool crosses(A const& a, B const& b) const
    { return bg::crosses(a, b COMPARISON_STRATEGY); }
    template <typename A, typename B> bool overlaps(A const& a, B const& b) const
    { return bg::overlaps(a, b COMPARISON_STRATEGY); }
    template <typename A, typename B> bool touches(A const& a, B const& b) const
    { return bg::touches(a, b COMPARISON_STRATEGY); }
};
#undef COMPARISON_STRATEGY

template <typename G> G read(char const* text)
{ G g; bg::read_wkt(text, g); return g; }
template <typename C> typename native<C>::polygons read_polygons(char const* text, bool correct = true)
{
    typename native<C>::polygons g;
    if (std::string(text).find("MULTIPOLYGON") == 0) bg::read_wkt(text, g);
    else g.push_back(read<typename native<C>::polygon>(text));
    // Orientation/closure normalization only, the same native coordinates on
    // all builds. Do not let an arithmetic-strategy choice alter the fixtures.
    if (correct) bg::correct(g);
    return g;
}

template <typename Tag> struct rational_geometry;
template <> struct rational_geometry<bg::polygon_tag> { using type = rational_models::polygon; };
template <> struct rational_geometry<bg::ring_tag> { using type = rational_models::polygon; };
template <> struct rational_geometry<bg::multi_polygon_tag> { using type = rational_models::polygons; };
template <> struct rational_geometry<bg::linestring_tag> { using type = rational_models::line; };
template <> struct rational_geometry<bg::multi_linestring_tag> { using type = rational_models::lines; };
template <> struct rational_geometry<bg::multi_point_tag> { using type = rational_models::points; };
template <> struct rational_geometry<bg::box_tag> { using type = rational_models::box; };

template <typename G> auto lift(G const& input) -> typename rational_geometry<typename bg::tag<G>::type>::type
{
    typename rational_geometry<typename bg::tag<G>::type>::type result;
    // cpp_rational's floating constructor preserves the binary value, not its
    // shortest decimal rendering. The execution script tests this invariant.
    bg::convert(input, result);
    return result;
}
template <typename G> void geometry(record& r, std::string const& key, G const& g)
{
    auto q = lift(g);
    std::ostringstream wkt;
    wkt << bg::wkt(q);
    r.put("geometries." + key + ".wkt", wkt.str());
    r.put("geometries." + key + ".points", bg::num_points(g));
}
template <typename G> void output(record& r, std::string const& key, G const& g)
{
    geometry(r, key, g);
    std::string reason;
    // This is the selected build's public validator, separate from the frozen
    // rational validator used by exact_issue_comparison_evaluate.cpp.
    r.put("geometries." + key + ".public_valid", bg::is_valid(g, reason));
    r.put("geometries." + key + ".public_reason", reason);
}
void emit(record const& r)
{ boost::property_tree::write_json(std::cout, r, false); std::cout.flush(); }

#if COMPARISON_GROUP == 1
template <typename C, typename Out, typename A, typename B>
void areal(exact_issue_cases::wkt_pair const& item, A a, B b, record& r)
{
    if(item.correct_inputs) { bg::correct(a); bg::correct(b); }
    geometry(r, "a", a); geometry(r, "b", b);
    api<C> call;
    std::string op(item.operation);
    if (op == "covered_by" || op == "within")
    {
        r.put("predicate", op == "within" ? call.within(a,b) : call.covered(a,b));
        r.put("matrix", call.relation(a,b));
        return;
    }
    typename Out::polygons result;
    call.overlay(op, a, b, result);
    output(r, "result", result);
    if (op != "difference")
    {
        result.clear(); call.overlay(op, b, a, result);
        output(r, "reverse", result);
    }
    if (item.issue == 704)
    {
        r.put("overlaps", call.overlaps(a,b));
        r.put("matrix", call.relation(a,b));
    }
}
template <typename C, typename Out>
void areal_shapes(exact_issue_cases::wkt_pair const& item, record& r)
{
    using N=native<C>;
    bool am=std::string(item.first).find("MULTIPOLYGON")==0;
    bool bm=std::string(item.second).find("MULTIPOLYGON")==0;
    if(am&&bm) areal<C,Out>(item,read<typename N::polygons>(item.first),read<typename N::polygons>(item.second),r);
    else if(am) areal<C,Out>(item,read<typename N::polygons>(item.first),read<typename N::polygon>(item.second),r);
    else if(bm) areal<C,Out>(item,read<typename N::polygon>(item.first),read<typename N::polygons>(item.second),r);
    else areal<C,Out>(item,read<typename N::polygon>(item.first),read<typename N::polygon>(item.second),r);
}
template <typename C> void areal_dispatch(exact_issue_cases::wkt_pair const& item, record& r, bool keep)
{
#ifdef COMPARISON_EXACT
    if (keep) return areal_shapes<C, retained<C>>(item, r);
#endif
    areal_shapes<C, native<C>>(item, r);
}
#endif

#if COMPARISON_GROUP == 2
template <typename C, typename Out, typename A, typename B>
void linear_clip(cases::typed_pair const& item, A a, B b, record& r)
{
    api<C> call;
    bg::correct(b);
    geometry(r, "a", a); geometry(r, "b", b);
    typename Out::lines result;
    bool const polygon_first = item.data.issue==510 || item.data.issue==511
        || item.data.issue==1253 || item.data.issue==858;
    if (polygon_first) call.intersection(b,a,result); else call.intersection(a,b,result);
    output(r, "result", result);
    bg::reverse(a); result.clear();
    if (polygon_first) call.intersection(b,a,result); else call.intersection(a,b,result);
    output(r, "reverse", result);
}
template <typename Out>
void double_api(cases::typed_pair const& item, record& r)
{
    using N = native<double>;
    api<double> call;
    using kind = cases::geometry_type;
    if (item.data.issue == 1235)
    {
        auto a = read<N::polygon>(item.data.first);
        auto b = read<N::polygons>(item.data.second);
        bg::correct(a); bg::correct(b);
        geometry(r,"a",a); geometry(r,"b",b);
        std::tuple<typename Out::points, typename Out::lines, typename Out::polygons> out;
        call.intersection(a,b,out);
        geometry(r,"isolated",std::get<0>(out));
        output(r,"lines",std::get<1>(out)); output(r,"result",std::get<2>(out));
    }
    else if (item.data.issue == 1406)
    {
        auto a = read<N::line>(item.data.first);
        auto b = read<N::polygon>(item.data.second); bg::correct(b);
        geometry(r,"a",a); geometry(r,"b",b);
        for (unsigned reverse=0; reverse<2; ++reverse)
        {
            std::tuple<typename Out::points, typename Out::lines> out;
            call.intersection(b,a,out);
            std::string suffix = reverse ? "_reverse" : "";
            geometry(r,"isolated"+suffix,std::get<0>(out));
            output(r,"lines"+suffix,std::get<1>(out));
            bg::reverse(a);
        }
    }
    else if (item.data.issue == 1232)
    {
        auto a = read<N::line>(item.data.first);
        auto b = read<N::lines>(item.data.second);
        geometry(r,"a",a); geometry(r,"b",b);
        typename Out::lines result;
        call.difference(a,b,result); output(r,"result",result);
        bg::reverse(a); result.clear();
        call.difference(a,b,result); output(r,"reverse",result);
    }
    else if (item.first_type == kind::multi_linestring)
        linear_clip<double,Out>(item,read<N::lines>(item.data.first),read<N::polygons>(item.data.second),r);
    else linear_clip<double,Out>(item,read<N::line>(item.data.first),read<N::polygon>(item.data.second),r);
}
template <typename Out> void integer_api(cases::typed_pair const& item, record& r)
{
    using N = native<std::int32_t>;
    api<std::int32_t> call;
    auto a = read<N::line>(item.data.first);
    if (item.data.issue == 440)
    {
        auto b = read<N::line>(item.data.second);
        geometry(r,"a",a); geometry(r,"b",b);
        for (unsigned ra=0; ra<2; ++ra) for (unsigned rb=0; rb<2; ++rb)
        {
            bg::model::segment<N::point> sa(a[ra],a[1-ra]), sb(b[rb],b[1-rb]);
            for (bool swap : {false,true})
            {
                typename Out::points result;
                call.intersection(swap?sb:sa,swap?sa:sb,result);
                geometry(r,"result"+std::to_string(ra)+std::to_string(rb)+std::to_string(swap),result);
            }
        }
    }
    else if (item.data.issue == 1448)
        linear_clip<std::int32_t,Out>(item,a,read<N::box>(item.data.second),r);
    else
    {
        auto polygon=read<N::polygon>(item.data.second); bg::correct(polygon);
        auto ring=polygon.outer();
        linear_clip<std::int32_t,Out>(item,a,ring,r);
    }
}
void crosses(cases::typed_pair const& item, record& r)
{
    using N = native<float>;
    auto a = read<N::line>(item.data.first);
    auto b = read<N::polygon>(item.data.second); bg::correct(b);
    geometry(r,"a",a); geometry(r,"b",b);
    api<float> call;
    r.put("predicate",call.crosses(a,b)); r.put("matrix",call.relation(a,b));
    bg::reverse(a);
    r.put("reverse_predicate",call.crosses(a,b)); r.put("reverse_matrix",call.relation(a,b));
}
double binary64(std::uint64_t u) { double d; std::memcpy(&d,&u,sizeof d); return d; }
template <typename Out> void boundary528(record& r)
{
    using N = native<double>;
    N::polygon a; N::line b;
    for (auto p : cases::issue_528_polygon()) a.outer().emplace_back(binary64(p.x),binary64(p.y));
    for (auto p : cases::issue_528_line()) b.emplace_back(binary64(p.x),binary64(p.y));
    bg::correct(a); geometry(r,"a",a); geometry(r,"b",b);
    typename Out::points result;
    api<double>().intersection(a,b,result); geometry(r,"result",result);
}
#endif

#if COMPARISON_GROUP == 3
template <typename Out> void pipeline(cases::pipeline_case const& item, record& r)
{
    using M = typename Out::polygons;
    api<double> call;
    std::vector<M> input;
    for (auto wkt : item.inputs)
    {
        auto n = read_polygons<double>(wkt);
        geometry(r,"input"+std::to_string(input.size()),n);
        M g; bg::convert(n,g); input.push_back(std::move(g));
    }
    unsigned steps=0;
    auto overlay = [&](std::string op, M const& a, M const& b)
    {
        M out; call.overlay(op,a,b,out);
        output(r,"step"+std::to_string(++steps),out);
        return out;
    };
    auto diff = [&](M const& a,M const& b) { return overlay("difference",a,b); };
    auto sym = [&](M const& a,M const& b) { return overlay("sym_difference",a,b); };
    M result;
    if (item.issue==1419) result=diff(diff(input[0],input[1]),input[2]);
    else if (item.issue==1236)
    {
        auto a=diff(input[0],input[1]); result=diff(a,input[3]); auto c=diff(input[2],input[3]);
        r.put("predicate",call.within(result,c));
        geometry(r,"reference_container",c);
    }
    else if (item.issue==900)
    {
        auto a=sym(input[0],input[1]); auto b=sym(a,input[2]); auto c=sym(input[3],input[2]);
        result=sym(b,c); auto expected=sym(a,input[3]); geometry(r,"identity_reference",expected);
    }
    else if (item.issue==588 || item.issue==1290)
    {
        result=diff(input[0],input[1]);
        r.put("predicate",call.within(result,input[0]));
        r.put("matrix",call.relation(result,input[0]));
        if (item.issue==588) r.put("frozen_predicate",call.within(input[2],input[0]));
        else
        {
            r.put("subtractor_matrix",call.relation(result,input[1]));
            r.put("overlaps",call.overlaps(result,input[1]));
            r.put("touches",call.touches(result,input[1]));
        }
    }
    else if (item.issue==888)
    {
        std::reverse(input.begin(),input.end());
        while (input.size()>1)
        {
            std::size_t k=(input.size()+1)/2;
            for (std::size_t i=0;i<input.size()/2;++i) input[i]=sym(input[i+k],input[i]);
            input.resize(k);
        }
        result=input[0];
    }
    else if (item.issue==1174)
    {
        auto n=read_polygons<double>(item.inputs[0]);
        native<double>::box bounds; bg::envelope(n,bounds);
        double x0=bg::get<bg::min_corner,0>(bounds), x1=bg::get<bg::max_corner,0>(bounds);
        double y0=bg::get<bg::min_corner,1>(bounds), y1=bg::get<bg::max_corner,1>(bounds);
        double cx=(x0+x1)/2, cy=(y0+y1)/2;
        for (bool right : {false,true}) for (bool upper : {false,true})
        {
            native<double>::box box({right?cx:x0,upper?cy:y0},{right?x1:cx,upper?y1:cy});
            typename Out::box qbox; bg::convert(box,qbox);
            M piece; call.intersection(input[0].front(),qbox,piece);
            output(r,"step"+std::to_string(++steps),piece);
            result.insert(result.end(),piece.begin(),piece.end());
        }
    }
    output(r,"result",result); r.put("steps",steps);
    if(item.issue==1174) r.put("aggregate_note","Synthetic concatenation, not an output multipolygon contract: quadrant pieces share edges. Validate each step separately.");
}
#endif

#if COMPARISON_GROUP == 4
template <typename Multi> void rotate_starts(Multi& m)
{
    auto rotate=[](auto& ring)
    {
        if(ring.size()<4) return;
        ring.pop_back(); std::rotate(ring.begin(),ring.begin()+1,ring.end()); ring.push_back(ring.front());
    };
    for(auto& p:m) { rotate(p.outer()); for(auto& ring:p.inners()) rotate(ring); }
}
template <typename Out> void api1044(bool control,bool shifted,record& r)
{
    auto a=read_polygons<double>(cases::data::issue_1044_0);
    auto b=read_polygons<double>(control?cases::data::issue_1044_control_1:cases::data::issue_1044_1);
    if(shifted) bg::for_each_point(b,[](native<double>::point& p){bg::set<1>(p,bg::get<1>(p)+1e-13);});
    geometry(r,"a",a); geometry(r,"b",b);
    api<double> call; typename Out::polygons result;
    call.intersection(a,b,result); output(r,"result",result);
    result.clear(); call.intersection(b,a,result); output(r,"reverse",result);
}
template <typename Out> void api669(unsigned branch,record& r)
{
    using N=native<double,false>;
    auto a=read<typename N::polygons>(cases::data::issue_669_0);
    auto b=read<typename N::polygons>(cases::data::issue_669_1);
    bg::correct(a); bg::correct(b);
    geometry(r,"a",a); geometry(r,"b",b);
    typename Out::polygons source;
    api<double> call;
    if(branch==0) bg::convert(a,source);
    else
    {
        typename N::box box; bg::envelope(b,box);
        bg::set<bg::min_corner,0>(box,bg::get<bg::min_corner,0>(box)-1);
        bg::set<bg::min_corner,1>(box,bg::get<bg::min_corner,1>(box)-1);
        bg::set<bg::max_corner,0>(box,bg::get<bg::max_corner,0>(box)+1);
        bg::set<bg::max_corner,1>(box,bg::get<bg::max_corner,1>(box)+1);
        if(branch==1) call.intersection(box,a,source);
        else
        {
            typename N::polygons legacy;
            bg::intersection(box,a,legacy);
            bg::convert(legacy,source);
        }
        output(r,"crop",source);
    }
    geometry(r,"source",source);
    for(unsigned rotate=0;rotate<2;++rotate)
    {
        if(rotate) rotate_starts(source);
        typename Out::polygons result;
        call.difference(source,b,result);
        output(r,rotate?"rotated":"result",result);
    }
}
#endif

int main(int argc,char** argv)
{
#if defined(COMPARISON_EXACT) && !defined(BOOST_GEOMETRY_EXACT_TYPED_CONTEXT)
    ec::point_cache_scope operation;
#endif
    if (argc==2 && std::string(argv[1])=="--binary-check")
    {
        rational q(0.1);
        if (q != rational(3602879701896397LL)/36028797018963968LL) return 1;
        std::cout << q << '\n'; return 0;
    }
    bool keep=argc>2 && std::string(argv[2])=="retained";
    std::string wanted=argc>1?argv[1]:"--list";
    auto run=[&](unsigned issue,char const* name,char const* coord,char const* op,auto action)
    {
        if (wanted=="--list") { std::cout<<name<<'\t'<<issue<<'\t'<<coord<<'\n'; return; }
        if (wanted!=name) return;
        record r; r.put("issue",issue); r.put("case",name); r.put("coordinate",coord); r.put("operation",op);
        r.put("mode",keep?"exact_retained":"native_output");
        r.put("status","ok");
        try
        {
#if defined(COMPARISON_EXACT) && defined(BOOST_GEOMETRY_EXACT_TYPED_CONTEXT)
            std::string const coordinate(coord);
            if (coordinate == "float") comparison_in_context<float>(action, r);
            else if (coordinate == "int32_t") comparison_in_context<std::int32_t>(action, r);
            else if (coordinate == "int64_t") comparison_in_context<std::int64_t>(action, r);
            else if (coordinate == "double") comparison_in_context<double>(action, r);
            else throw std::invalid_argument("Unknown exact comparison coordinate type");
#else
            action(r);
#endif
        }
#if defined(COMPARISON_EXACT) && defined(BOOST_GEOMETRY_EXACT_TYPED_CONTEXT)
        catch(boost::numeric::bad_numeric_cast const& e)
        {
            r.put("status","exception"); r.put("exception",e.what());
            r.put("exception_kind","numeric_conversion");
        }
#endif
        catch(std::exception const& e) { r.put("status","exception"); r.put("exception",e.what()); }
        emit(r);
    };
#if COMPARISON_GROUP == 1
#ifndef COMPARISON_COORDINATE
#define COMPARISON_COORDINATE double
#define COMPARISON_COORDINATE_NAME "double"
#endif
    auto selected=[&](exact_issue_cases::wkt_pair const& item)
    {
        if (std::string(item.operation)=="within" && item.issue==1238) return;
        if (std::string(item.coordinate_type)!=COMPARISON_COORDINATE_NAME) return;
        run(item.issue,item.name,item.coordinate_type,item.operation,[&](record& r)
        {
            areal_dispatch<COMPARISON_COORDINATE>(item,r,keep);
        });
    };
    for(auto const& i:exact_issue_cases::cases()) selected(i);
    for(auto const& i:cases::pairs())
        if ((i.first_type==cases::geometry_type::polygon || i.first_type==cases::geometry_type::multi_polygon)
            && std::string(i.data.operation)!="intersection_tuple") selected(i.data);
#elif COMPARISON_GROUP == 2
    for (auto const& i:cases::pairs())
    {
        if (i.first_type==cases::geometry_type::polygon || i.first_type==cases::geometry_type::multi_polygon)
            if(i.data.issue!=1235) continue;
        run(i.data.issue,i.data.name,i.data.coordinate_type,i.data.operation,[&](record& r)
        {
            if(i.data.issue==831) return crosses(i,r);
            bool integral=std::string(i.data.coordinate_type)=="int32_t";
#ifdef COMPARISON_EXACT
            if(keep)
            {
                if(integral) integer_api<retained<std::int32_t>>(i,r);
                else double_api<retained<double>>(i,r);
                return;
            }
#endif
            if(integral) integer_api<native<std::int32_t>>(i,r);
            else double_api<native<double>>(i,r);
        });
    }
    run(528,"issue_528_binary64","double","intersection_points",[&](record& r)
    {
#ifdef COMPARISON_EXACT
        if(keep) return boundary528<retained<double>>(r);
#endif
        boundary528<native<double>>(r);
    });
    for(auto const& i:exact_issue_cases::cases()) if(i.issue==1238)
        run(i.issue,i.name,i.coordinate_type,i.operation,[&](record& r)
        {
            auto a=read<native<double>::lines>(i.first), b=read<native<double>::lines>(i.second);
            geometry(r,"a",a); geometry(r,"b",b);
            api<double> call; r.put("predicate",call.within(a,b)); r.put("matrix",call.relation(a,b));
        });
#elif COMPARISON_GROUP == 3
    for(auto const& i:cases::pipelines()) run(i.issue,i.name,i.coordinate_type,"pipeline",[&](record& r)
    {
#ifdef COMPARISON_EXACT
        if(keep) return pipeline<retained<double>>(i,r);
#endif
        pipeline<native<double>>(i,r);
    });
#elif COMPARISON_GROUP == 4
    for(bool control:{false,true}) for(bool shifted:{false,true})
    {
        std::string name=std::string("issue_1044_api_")+(control?"control":"bug")+(shifted?"_shifted":"");
        run(1044,name.c_str(),"double","intersection",[&](record& r)
        {
#ifdef COMPARISON_EXACT
            if(keep) return api1044<retained<double>>(control,shifted,r);
#endif
            api1044<native<double>>(control,shifted,r);
        });
    }
    for(unsigned branch=0;branch<3;++branch)
    {
        std::string name=std::string("issue_669_api_")+(branch==0?"full":branch==1?"strategy_crop_feedback":"source_default_crop");
        run(669,name.c_str(),"double","difference",[&](record& r)
        {
#ifdef COMPARISON_EXACT
            if(keep) return api669<retained<double,false>>(branch,r);
#endif
            api669<native<double,false>>(branch,r);
        });
    }
#endif
}
