// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Distributed under the Boost Software License, Version 1.0.
#define BOOST_GEOMETRY_EXACT_FILTER_STATISTICS
#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/algorithms/disjoint.hpp>
#include <boost/geometry/algorithms/intersection.hpp>
#include <boost/geometry/algorithms/assign.hpp>
#include <boost/geometry/strategies/cartesian/exact/static_box.hpp>
#include <boost/geometry/strategies/cartesian/box_in_box.hpp>
#include <boost/geometry/strategies/cartesian/point_in_box.hpp>
#include <boost/geometry/strategies/cartesian/disjoint_box_box.hpp>
#include <boost/geometry/strategies/detail/box_topology.hpp>
#include <boost/geometry/strategies/relate/cartesian_exact.hpp>
#include <array>
#include <cfenv>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <new>
#include <random>
#include <vector>

namespace bg = boost::geometry;
namespace ec = bg::strategies::detail::exact_cartesian;
using q = ec::exact_number;
template <typename C> using point = bg::model::d2::point_xy<C>;
template <typename C> using box = bg::model::box<point<C>>;
template <typename C> using exact_box = bg::model::box<ec::exact_point<C>>;

namespace
{
bool count_allocations = false;
std::size_t allocations = 0;
}
void* operator new(std::size_t size)
{
    if (count_allocations) ++allocations;
    if (void* result = std::malloc(size ? size : 1)) return result;
    throw std::bad_alloc();
}
void operator delete(void* value) noexcept { std::free(value); }
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete[](void* value) noexcept { ::operator delete(value); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* value, std::size_t) noexcept { ::operator delete(value); }
void operator delete[](void* value, std::size_t) noexcept { ::operator delete(value); }
#endif

template <typename F> void no_allocations(F function)
{
    function(); // Grow indexed point caches before counting numerical work.
    allocations = 0;
    count_allocations = true;
    function();
    count_allocations = false;
    BOOST_TEST_EQ(allocations, 0u);
}

template <typename C, typename Tag>
struct opaque_view
{
    std::array<C,4> const& values;
    unsigned& reads;
    opaque_view(std::array<C,4> const& v,unsigned& r):values(v),reads(r){}
    opaque_view(opaque_view const&)=delete;
    opaque_view& operator=(opaque_view const&)=delete;
};
namespace boost { namespace geometry { namespace traits {
template <typename C,typename Tag> struct tag<opaque_view<C,Tag>> { using type=Tag; };
template <typename C,typename Tag> struct point_type<opaque_view<C,Tag>> { using type=::point<C>; };
template <typename C,typename Tag,std::size_t I,std::size_t D>
struct indexed_access<opaque_view<C,Tag>,I,D>
{
    static C get(opaque_view<C,Tag> const& view) { ++view.reads;return view.values[2*I+D]; }
};
}}}

template <typename C> q value(C const& x) { return q(x); }
template <typename P> point<q> oracle_point(P const& p)
{ return {value(bg::get<0>(p)),value(bg::get<1>(p))}; }
template <typename B> box<q> oracle_box(B const& b)
{
    return {{value(bg::get<bg::min_corner,0>(b)),value(bg::get<bg::min_corner,1>(b))},
            {value(bg::get<bg::max_corner,0>(b)),value(bg::get<bg::max_corner,1>(b))}};
}
template <typename B> void same_bounds(B const& actual,box<q> const& expected)
{
    auto const a=oracle_box(actual);
    BOOST_TEST((bg::get<bg::min_corner,0>(a)==bg::get<bg::min_corner,0>(expected)));
    BOOST_TEST((bg::get<bg::min_corner,1>(a)==bg::get<bg::min_corner,1>(expected)));
    BOOST_TEST((bg::get<bg::max_corner,0>(a)==bg::get<bg::max_corner,0>(expected)));
    BOOST_TEST((bg::get<bg::max_corner,1>(a)==bg::get<bg::max_corner,1>(expected)));
}
void encloses(ec::search_box const& actual,box<q> const& expected)
{
    double const lx=bg::get<bg::min_corner,0>(actual),ly=bg::get<bg::min_corner,1>(actual);
    double const ux=bg::get<bg::max_corner,0>(actual),uy=bg::get<bg::max_corner,1>(actual);
    BOOST_TEST((lx==-(std::numeric_limits<double>::infinity)() || q(lx)<=bg::get<bg::min_corner,0>(expected)));
    BOOST_TEST((ly==-(std::numeric_limits<double>::infinity)() || q(ly)<=bg::get<bg::min_corner,1>(expected)));
    BOOST_TEST((ux==(std::numeric_limits<double>::infinity)() || q(ux)>=bg::get<bg::max_corner,0>(expected)));
    BOOST_TEST((uy==(std::numeric_limits<double>::infinity)() || q(uy)>=bg::get<bg::max_corner,1>(expected)));
}
template <typename F> void no_q(F function)
{
    ec::filtering_statistics()={};
    function();
    BOOST_TEST_EQ(ec::filtering_statistics().native_rationals,0u);
    BOOST_TEST_EQ(ec::filtering_statistics().lazy_exact_values,0u);
    BOOST_TEST_EQ(ec::filtering_statistics().point_coordinate_values,0u);
    BOOST_TEST_EQ(ec::filtering_statistics().homogeneous_values,0u);
}
template <typename A,typename B> void box_predicates(A const& a,B const& b)
{
    auto const x=oracle_box(a),y=oracle_box(b);
    using pred=ec::static_box_predicates<double>;
    bool covered,within,disjoint,touches,overlaps;
    no_q([&]{covered=pred::covered_by(a,b);within=pred::within(a,b);disjoint=pred::disjoint(a,b);
        touches=pred::touches(a,b);overlaps=pred::overlaps(a,b);});
    BOOST_TEST_EQ(covered,bg::strategy::covered_by::cartesian_box_box::apply(x,y));
    BOOST_TEST_EQ(within,bg::strategy::within::cartesian_box_box::apply(x,y));
    BOOST_TEST_EQ(disjoint,bg::strategy::disjoint::cartesian_box_box::apply(x,y));
    BOOST_TEST_EQ(touches,bg::strategies::detail::legacy_box_touches::apply(x,y));
    BOOST_TEST_EQ(overlaps,bg::strategies::detail::legacy_box_overlaps::apply(x,y));
}
template <typename P,typename B> void point_predicates(P const& p,B const& b)
{
    auto const x=oracle_point(p);auto const y=oracle_box(b);
    bool covered,within;
    no_q([&]{covered=ec::static_box_predicates<double>::point_covered_by(p,b);
        within=ec::static_box_predicates<double>::point_covered_by(p,b,true);});
    BOOST_TEST_EQ(covered,bg::strategy::covered_by::cartesian_point_box::apply(x,y));
    BOOST_TEST_EQ(within,bg::strategy::within::cartesian_point_box::apply(x,y));
}
template <typename C> void native_bounds()
{
    using p=point<C>;
    std::array<p,4> const points{{{-8,3},{4,-5},{12,7},{2,0}}};
    box<q> const expected({-8,-5},{12,7});
    exact_box<C> semantic;ec::search_box search;
    no_q([&]{ec::static_envelope_range<C>::apply(points,semantic);ec::static_envelope_range<C>::apply(points,search);});
    same_bounds(semantic,expected);encloses(search,expected);
    point_predicates(points[3],semantic);point_predicates(points[0],semantic);
    exact_box<C> expanded;
    ec::static_envelope_point<C>::apply(points[3],expanded);
    for(auto const& p:points)ec::static_expand_point<C>::apply(expanded,p);
    same_bounds(expanded,expected);
    bg::model::segment<p> const segment(points[0],points[1]);
    ec::static_envelope_segment<C>::apply(segment,expanded);
    same_bounds(expanded,box<q>({-8,-5},{4,3}));
    bg::model::referring_segment<p const> const referring(points[2],points[3]);
    ec::static_expand_segment<C>::apply(expanded,referring);same_bounds(expanded,expected);
    typename ec::static_envelope_boxes<C>::template state<ec::search_box> state;
    ec::static_envelope_boxes<C>::apply(state,semantic);
    ec::static_envelope_boxes<C>::apply(state,expanded);
    ec::static_envelope_boxes<C>::result(state,search);encloses(search,expected);
    std::vector<p> empty;ec::static_envelope_range<C>::apply(empty,search);
    BOOST_TEST((bg::get<bg::min_corner,0>(search)>bg::get<bg::max_corner,0>(search)));
    ec::static_expand_point<C>::apply(search,points[0]);encloses(search,box<q>({-8,3},{-8,3}));
}
template <typename C> void opaque_snapshots(std::array<C,4> input)
{
    using kernel_coordinate = std::conditional_t<ec::fixed_native_supported<C>::value, C, double>;
    ec::point_cache_scope<kernel_coordinate> operation;
    box<q> const expected({q(input[0]),q(input[1])},{q(input[2]),q(input[3])});
    exact_box<kernel_coordinate> retained,segment,expanded;ec::search_box search;
    unsigned box_reads=0,segment_reads=0;
    auto coordinates=input;
    std::array<C,4> segment_coordinates{{input[2],input[3],input[0],input[1]}};
    opaque_view<C,bg::box_tag> const box_view(coordinates,box_reads);
    opaque_view<C,bg::segment_tag> const segment_view(segment_coordinates,segment_reads);
    no_q([&]{ec::static_envelope_box<kernel_coordinate>::apply(box_view,retained);ec::static_envelope_box<kernel_coordinate>::apply(box_view,search);});
    BOOST_TEST(box_reads > 0);
    {
        ec::static_envelope_segment<kernel_coordinate>::apply(segment_view,segment);
        ec::static_envelope_point<kernel_coordinate>::apply(point<C>(0,0),expanded);
        ec::static_expand_segment<kernel_coordinate>::apply(expanded,segment_view);
        BOOST_TEST(segment_reads > 0);
    }
    same_bounds(retained,expected);same_bounds(segment,expected);encloses(search,expected);
    same_bounds(expanded,box<q>({0,0},{q(input[2]),q(input[3])}));
}
ec::exact_point<double> third(double x,double y)
{
    using p=ec::exact_point<double>;
    return ec::construct_intersection(p(x,y),p(x+1,y+1),p(x,y+1),p(x+1,y-1));
}
ec::exact_point<double> implicit_point(double x,double y)
{
    using p=ec::exact_point<double>;
    return ec::construct_intersection(p(x-1,y),p(x+1,y),p(x,y-1),p(x,y+1));
}

void segment_adapter_snapshots()
{
    using p=ec::exact_point<double>;
    bg::model::segment<p> segment(third(0,0),third(3,3));
    no_q([&]
    {
        bg::detail::segment_as_subrange<decltype(segment)> snapshot(segment);
        segment.first=p(100,100);
        BOOST_TEST(ec::static_equal_points(snapshot.at(0),third(0,0)));
        BOOST_TEST(ec::static_equal_points(snapshot.at(1),segment.second));
    });
    p a=third(0,0),b=third(3,3);
    bg::model::referring_segment<p> const reference(a,b);
    no_q([&]
    {
        bg::detail::segment_as_subrange<decltype(reference)> snapshot(reference);
        a=p(100,100);
        BOOST_TEST(ec::static_equal_points(snapshot.at(0),third(0,0)));
        BOOST_TEST(ec::static_equal_points(snapshot.at(1),b));
    });
    std::array<int,4> values{{1,2,3,4}};
    unsigned reads=0;
    opaque_view<int,bg::segment_tag> const adapted(values,reads);
    bg::detail::segment_as_subrange<decltype(adapted)> snapshot(adapted);
    BOOST_TEST_EQ(reads,4u);
    values.fill(100);
    BOOST_TEST_EQ(bg::get<0>(snapshot.at(0)),1);
    BOOST_TEST_EQ(bg::get<1>(snapshot.at(0)),2);
    BOOST_TEST_EQ(bg::get<0>(snapshot.at(1)),3);
    BOOST_TEST_EQ(bg::get<1>(snapshot.at(1)),4);
    BOOST_TEST_EQ(reads,4u);
}

template <typename Segment> void segment_box_routing(Segment const& segment)
{
    using qp=point<q>;
    bg::model::segment<qp> const oracle(oracle_point(ec::static_box_detail::corner<0>(segment)),
                                      oracle_point(ec::static_box_detail::corner<1>(segment)));
    box<int> const bounds({0,0},{1,1});
    bool const expected=bg::disjoint(oracle,oracle_box(bounds),bg::strategies::relate::cartesian<q>());
    bg::strategies::relate::cartesian_exact<double> const strategy;
    no_q([&]{no_allocations([&]
    {
        BOOST_TEST_EQ(bg::disjoint(segment,bounds,strategy),expected);
        BOOST_TEST_EQ(bg::disjoint(bounds,segment,strategy),expected);
    });});
}

template <typename A,typename B> void segment_segment_routing(A const& a,B const& b)
{
    using qp=point<q>;
    auto lift=[](auto const& segment)
    {
        return bg::model::segment<qp>(oracle_point(ec::static_box_detail::corner<0>(segment)),
                                     oracle_point(ec::static_box_detail::corner<1>(segment)));
    };
    auto const qa=lift(a),qb=lift(b);
    bg::strategies::relate::cartesian<q> const validation;
    bool const expected=bg::disjoint(qa,qb,validation);
    std::vector<qp> oracle;
    bg::intersection(qa,qb,oracle,validation);
    std::vector<point<double>> actual;
    actual.reserve(2);
    bg::strategies::relate::cartesian_exact<double> const strategy;
    no_q([&]{no_allocations([&]
    {
        BOOST_TEST_EQ(bg::disjoint(a,b,strategy),expected);
        actual.clear();
        bg::detail::intersection::intersection_segment_segment_point<point<double>>
            ::apply(a,b,std::back_inserter(actual),strategy);
    });});
    BOOST_TEST_EQ(actual.size(),oracle.size());
    for(std::size_t i=0;i<(std::min)(actual.size(),oracle.size());++i)
    {
        BOOST_TEST_EQ(bg::get<0>(actual[i]),ec::point_detail::output_rounder<double>::apply(bg::get<0>(oracle[i])));
        BOOST_TEST_EQ(bg::get<1>(actual[i]),ec::point_detail::output_rounder<double>::apply(bg::get<1>(oracle[i])));
    }
    // The public wrapper owns a working result vector; only its math is no-Q.
    no_q([&]{actual.clear();bg::intersection(a,b,actual,strategy);});
    BOOST_TEST_EQ(actual.size(),oracle.size());
    for(std::size_t i=0;i<(std::min)(actual.size(),oracle.size());++i)
    {
        BOOST_TEST_EQ(bg::get<0>(actual[i]),ec::point_detail::output_rounder<double>::apply(bg::get<0>(oracle[i])));
        BOOST_TEST_EQ(bg::get<1>(actual[i]),ec::point_detail::output_rounder<double>::apply(bg::get<1>(oracle[i])));
    }
}

void implicit_segment_routing()
{
    using p=ec::exact_point<double>;
    using segment=bg::model::segment<p>;
    std::array<segment,7> const segments{{
        {third(-1,-1),third(1,1)}, {third(2,0),third(3,0)},
        {implicit_point(0,0),third(-1,-1)}, {implicit_point(0,-1),implicit_point(0,2)},
        {third(0,0),third(0,0)}, {implicit_point(1,1),implicit_point(1,1)},
        {implicit_point(2,2),implicit_point(2,2)}
    }};
    for(auto const& s:segments)
    {
        segment_box_routing(s);
        segment_box_routing(bg::model::referring_segment<p const>(s.second,s.first));
    }
    double const tiny=std::ldexp(1.0,-150);
    p const near=ec::construct_intersection(p(0,0),p(1,0),p(tiny,-1),p(tiny,1));
    segment const horizontal(implicit_point(0,0),implicit_point(1,0));
    std::array<segment,6> const others{{
        {third(0,-1),third(0,1)}, {third(2,0),third(3,0)},
        {implicit_point(1,0),implicit_point(2,1)}, {near,near},
        {implicit_point(1,0),near}, {implicit_point(0,1),implicit_point(1,1)}
    }};
    for(auto const& s:others)
    {
        segment_segment_routing(horizontal,s);
        segment_segment_routing(bg::model::referring_segment<p const>(horizontal.second,horizontal.first),
                                bg::model::referring_segment<p const>(s.second,s.first));
    }
}
void retained_implicit()
{
    exact_box<double> retained;ec::search_box search;
    {
        auto a=third(0,0),b=third(3,3);
        no_q([&]{ec::static_envelope_segment<double>::apply(a,b,retained);
            ec::static_envelope_segment<double>::apply(a,b,search);});
        bg::set<0>(a,q(100));bg::set<1>(b,q(-100));
    }
    q const t=q(1)/3;box<q> const expected({t,t},{3+t,3+t});
    same_bounds(retained,expected);encloses(search,expected);
    point_predicates(third(1,1),retained);point_predicates(point<float>(1,1),retained);
    box_predicates(retained,box<std::int64_t>({0,0},{2,2}));
    ec::static_expand_box<double>::apply(retained,box<int>({-1,-2},{0,0}));
    same_bounds(retained,box<q>({-1,-2},{3+t,3+t}));
}
void integer_extremes()
{
    using p=point<std::int64_t>;
    auto const lo=(std::numeric_limits<std::int64_t>::min)(),hi=(std::numeric_limits<std::int64_t>::max)();
    std::array<p,2> const input{{{lo,lo},{hi,hi}}};exact_box<std::int64_t> exact;ec::search_box search;
    no_q([&]{ec::static_envelope_range<std::int64_t>::apply(input,exact);
        ec::static_envelope_range<std::int64_t>::apply(input,search);});
    box<q> const expected({q(lo),q(lo)},{q(hi),q(hi)});same_bounds(exact,expected);encloses(search,expected);
    point_predicates(p(hi,hi),exact);point_predicates(p(lo+1,hi-1),exact);
    std::int64_t const n=9007199254740993LL;
    opaque_snapshots<std::int64_t>({{n,n,n+2,n+4}});
    box_predicates(box<std::int64_t>({n,n},{n+1,n+1}),box<std::int64_t>({n+1,n},{n+2,n+1}));
    auto const top=(std::numeric_limits<std::uint64_t>::max)();
    std::array<point<std::uint64_t>,2> const unsigned_input{{{0,0},{top,top}}};
    exact_box<std::uint64_t> unsigned_exact;
    ec::static_envelope_range<std::uint64_t>::apply(unsigned_input,unsigned_exact);
    ec::static_envelope_range<std::uint64_t>::apply(unsigned_input,search);
    same_bounds(unsigned_exact,box<q>({0,0},{q(top),q(top)}));encloses(search,box<q>({0,0},{q(top),q(top)}));
}
void output_rounding()
{
    using p=ec::exact_point<double>;
    q const h=q(1)/2;
    exact_box<double> const source(p::literal(h,3*h),p::literal(5*h,7*h));box<int> rounded;
    ec::static_envelope_box<double>::apply(source,rounded);same_bounds(rounded,box<q>({0,2},{2,4}));
    box<double> native;
    ec::static_envelope_box<double>::apply(exact_box<double>(third(0,0),third(3,3)),native);
    BOOST_TEST((bg::get<bg::min_corner,0>(native)==ec::point_detail::output_rounder<double>::apply(q(1)/3)));
    for(std::int64_t outside:{-32769LL,32768LL})
    {
        box<std::int64_t> input({outside,0},{outside,1});box<std::int16_t> narrow;
        BOOST_TEST_THROWS(ec::static_envelope_box<std::int64_t>::apply(input,narrow),boost::numeric::bad_numeric_cast);
    }
}
void spatial_grid()
{
    std::vector<box<int>> boxes;
    for(int x=-1;x<=1;++x)for(int xx=x;xx<=1;++xx)
    for(int y=-1;y<=1;++y)for(int yy=y;yy<=1;++yy)boxes.emplace_back(point<int>(x,y),point<int>(xx,yy));
    for(auto const& a:boxes)
    {
        for(auto const& b:boxes)box_predicates(a,b);
        for(double x:{-1.0,0.0,1.0})for(double y:{-1.0,0.0,1.0})point_predicates(point<double>(x,y),a);
    }
}
int main()
{
    ec::point_cache_scope operation;
    static_assert(sizeof(ec::search_box)==4*sizeof(double),"Search envelopes remain four native scalars");
    segment_adapter_snapshots();
    int const original=std::fegetround();
    for(int mode:{FE_TONEAREST})
    {
        if(std::fesetround(mode)!=0)continue;
        native_bounds<double>();native_bounds<std::int64_t>();retained_implicit();integer_extremes();output_rounding();
        opaque_snapshots<double>({{1,2,7,8}});
        opaque_snapshots<q>({{q(1)/3,q(2)/3,q(7)/3,q(8)/3}});
        implicit_segment_routing();
    }
    BOOST_TEST_EQ(std::fesetround(original),0);spatial_grid();
    return boost::report_errors();
}
