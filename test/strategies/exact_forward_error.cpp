// Boost.Geometry
// Copyright (c) 2026 Tinko Bartels.
// Licensed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt

#include <boost/core/lightweight_test.hpp>
#include <boost/geometry/strategies/cartesian/exact/forward_error.hpp>

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <new>
#include <random>
#include <vector>

namespace ec = boost::geometry::strategies::detail::exact_cartesian;
namespace fe = ec::forward_error;
using q = ec::exact_number;
using point = fe::point_approximation;
using scalar = fe::scalar_approximation;
using exact_point = std::array<q, 2>;

namespace
{
thread_local bool counting = false;
thread_local std::size_t allocations = 0;
}
void* operator new(std::size_t size)
{
    if (counting) ++allocations;
    if (void* p = std::malloc(size == 0 ? 1 : size)) return p;
    throw std::bad_alloc();
}
void operator delete(void* p) noexcept { std::free(p); }
void* operator new[](std::size_t size) { return ::operator new(size); }
void operator delete[](void* p) noexcept { ::operator delete(p); }
#if defined(__cpp_sized_deallocation)
void operator delete(void* p, std::size_t) noexcept { std::free(p); }
void operator delete[](void* p, std::size_t) noexcept { std::free(p); }
#endif

template <typename F> bool allocation_free(F f)
{
    allocations = 0; counting = true;
    bool const result = f();
    counting = false;
    BOOST_TEST_EQ(allocations, 0u);
    return result;
}

q absolute(q x) { return x < 0 ? q(-x) : x; }
void bound(scalar const& a, q const& exact)
{
    BOOST_TEST(a.valid());
    if (a.valid()) BOOST_TEST(absolute(ec::exact_value(a.value) - exact) <= ec::exact_value(a.error));
}
void bound(point const& a, exact_point const& exact)
{
    BOOST_TEST(a.valid());
    for (unsigned d = 0; d < 2; ++d) bound({a.value[d], a.error[d]}, exact[d]);
}
q cross(exact_point const& a, exact_point const& b) { return a[0]*b[1]-a[1]*b[0]; }
q dot(exact_point const& a, exact_point const& b) { return a[0]*b[0]+a[1]*b[1]; }
exact_point difference(exact_point const& a, exact_point const& b) { return {{a[0]-b[0], a[1]-b[1]}}; }
bool intersection(std::array<exact_point, 4> const& p, exact_point& result)
{
    auto const u = difference(p[1],p[0]), v = difference(p[3],p[2]), w = difference(p[2],p[0]);
    q const denominator = cross(u,v);
    if (denominator == 0) return false;
    q const t = cross(w,v)/denominator;
    result = {{p[0][0]+t*u[0], p[0][1]+t*u[1]}};
    return true;
}
q parameter(std::array<exact_point, 4> const& p, bool second)
{
    auto const u=difference(p[1],p[0]),v=difference(p[3],p[2]),w=difference(p[2],p[0]);
    return cross(w,second?u:v)/cross(u,v);
}
point sentinel() { return {{{17,19}},{{23,29}}}; }
void unchanged(point const& p) { BOOST_TEST(p.value == sentinel().value); BOOST_TEST(p.error == sentinel().error); }

template <typename C> void conversion(C value)
{
    double a = 17, e = 19;
    BOOST_TEST(allocation_free([&]{return fe::try_native_coordinate(value,a,e);}));
    bound({a,e},ec::exact_value(value));
}

void predicates(std::array<point,4> const& p,std::array<exact_point,4> const& x);

template <typename C> bool native_case(std::array<C,8> const& values)
{
    std::array<exact_point,4> exact;
    std::array<point,4> points;
    std::array<double,8> rounded, errors;
    for (unsigned i=0;i<4;++i)
    {
        BOOST_TEST(allocation_free([&]{return fe::try_native_point(values[2*i],values[2*i+1],points[i]);}));
        exact[i]={{ec::exact_value(values[2*i]),ec::exact_value(values[2*i+1])}};
        bound(points[i],exact[i]);
        for(unsigned d=0;d<2;++d) {rounded[2*i+d]=points[i].value[d];errors[2*i+d]=points[i].error[d];}
    }
    point fast=sentinel(), generic=sentinel(), already=sentinel();
    bool const accepted=allocation_free([&]{return fe::native_intersection(values,fast);});
    bool const perturbed=allocation_free([&]{return fe::try_intersection(points,generic);});
    bool const converted=allocation_free([&]{return fe::try_intersection(rounded,errors,already);});
    exact_point expected;
    bool const unique=intersection(exact,expected);
    if(accepted) {BOOST_TEST(unique);if(unique)bound(fast,expected);} else unchanged(fast);
    if(perturbed) {BOOST_TEST(unique);if(unique)bound(generic,expected);} else unchanged(generic);
    if(converted) {BOOST_TEST(unique);if(unique)bound(already,expected);} else unchanged(already);
    for(bool second:{false,true})
    {
        scalar direct{17,19},generic_parameter{17,19};
        bool const native_ok=allocation_free([&]{return fe::try_segment_parameter(values,second,direct);});
        bool const generic_ok=allocation_free([&]{return fe::try_segment_parameter(points,second,generic_parameter);});
        BOOST_TEST(!accepted||native_ok);
        if(native_ok) {BOOST_TEST(unique);if(unique)bound(direct,parameter(exact,second));}
        else {BOOST_TEST_EQ(direct.value,17);BOOST_TEST_EQ(direct.error,19);}
        if(generic_ok) {BOOST_TEST(unique);if(unique)bound(generic_parameter,parameter(exact,second));}
        else {BOOST_TEST_EQ(generic_parameter.value,17);BOOST_TEST_EQ(generic_parameter.error,19);}
    }
    predicates(points,exact);
    return accepted;
}

template <typename F> bool certified(F f, q const& expected)
{
    int sign=17;
    bool const certain=allocation_free([&]{return f(sign);});
    BOOST_TEST_EQ(sign,certain?expected.sign():17);
    return certain;
}

void predicates(std::array<point,4> const& p,std::array<exact_point,4> const& x)
{
    scalar s;
    BOOST_TEST(allocation_free([&]{return fe::try_determinant(p[0],p[1],s);}));bound(s,cross(x[0],x[1]));
    certified([&](int& r){return fe::try_determinant_sign(p[0],p[1],r);},cross(x[0],x[1]));
    BOOST_TEST(allocation_free([&]{return fe::try_dot(p[0],p[1],s);}));bound(s,dot(x[0],x[1]));
    auto const ab=difference(x[1],x[0]),ac=difference(x[2],x[0]),bc=difference(x[2],x[1]);
    auto const cd=difference(x[3],x[2]);
    BOOST_TEST(allocation_free([&]{return fe::try_signed_determinant(p[0],p[1],p[2],s);}));bound(s,cross(ab,ac));
    certified([&](int& r){return fe::try_orientation(p[0],p[1],p[2],r);},cross(ab,ac));
    certified([&](int& r){return fe::try_direction(p[0],p[1],p[2],r);},dot(ab,bc));
    BOOST_TEST(allocation_free([&]{return fe::try_squared_distance(p[0],p[1],s);}));bound(s,dot(ab,ab));
    certified([&](int& r){return fe::try_squared_distance_compare(p[0],p[1],p[2],p[3],r);},dot(ab,ab)-dot(cd,cd));
    certified([&](int& r){return fe::try_squared_distance_compare(p[0],p[1],p[2],r);},dot(ab,ab)-dot(ac,ac));
    for(unsigned d=0;d<2;++d)
    {
        certified([&](int& r){return fe::try_compare_coordinates(p[0],p[1],d,r);},x[0][d]-x[1][d]);
        scalar const target{p[1].value[d],p[1].error[d]};
        certified([&](int& r){return fe::try_compare_coordinate(p[0],d,target,r);},x[0][d]-x[1][d]);
        certified([&](int& r){return fe::try_compare_scalars(p[0].value[d],p[0].error[d],target.value,target.error,r);},x[0][d]-x[1][d]);
    }
    BOOST_TEST(allocation_free([&]{return fe::try_box_measure(p[0],p[1],s);}));bound(s,ab[0]*ab[1]);
}

void native_final_coefficients()
{
    q const u=ec::exact_value(std::ldexp(1.0,-53));
    q const one_minus_u=1-u;
    q const reverse_three=1/(one_minus_u*one_minus_u*one_minus_u)-1;
    BOOST_TEST(reverse_three<4*u);
    BOOST_TEST(8*u*one_minus_u>reverse_three+u);
    BOOST_TEST(16*u*one_minus_u*one_minus_u>reverse_three+2*u+u*u);

    point const a{{{0,0}},{{0,0}}},b{{{3,4}},{{0,0}}},p{{{1,2}},{{0,0}}};
    scalar result;
    BOOST_TEST(fe::try_signed_determinant(a,b,p,result));
    BOOST_TEST_EQ(result.value,2);
    BOOST_TEST_EQ(result.error,8*fe::detail::u*10);
    BOOST_TEST(fe::try_squared_distance(a,b,result));
    BOOST_TEST_EQ(result.value,25);
    BOOST_TEST_EQ(result.error,8*fe::detail::u*25);
    BOOST_TEST(fe::try_box_measure(a,b,result));
    BOOST_TEST_EQ(result.value,12);
    BOOST_TEST_EQ(result.error,8*fe::detail::u*12);

    auto perturbed=p;
    perturbed.error[0]=std::ldexp(1.0,-20);
    auto const generic=fe::detail::determinant(fe::detail::difference(b,a),fe::detail::difference(perturbed,a));
    BOOST_TEST(fe::try_signed_determinant(a,b,perturbed,result));
    BOOST_TEST_EQ(result.value,generic.value);
    BOOST_TEST_EQ(result.error,generic.error);
    BOOST_TEST(result.error>8*fe::detail::u*10);
    bound(result,q(2)-4*ec::exact_value(perturbed.error[0]));
    bound(result,q(2)+4*ec::exact_value(perturbed.error[0]));

    int sign=17;
    double const maximum=(std::numeric_limits<double>::max)();
    BOOST_TEST(fe::try_compare_scalars({maximum,0},{-maximum,0},sign));BOOST_TEST_EQ(sign,1);
    BOOST_TEST(fe::try_compare_scalars({-0.0,0},{0.0,0},sign));BOOST_TEST_EQ(sign,0);
    sign=17;
    BOOST_TEST(!fe::try_compare_scalars({1,4*fe::detail::u},{1,0},sign));BOOST_TEST_EQ(sign,17);

    // Normal-range cancellation exercises reverse-error bounds at many scales.
    for(int exponent:{-200,-100,0,100,200})
    {
        double const scale=std::ldexp(1.0,exponent);
        double const large=std::ldexp(1.0,52);
        native_case(std::array<double,8>{{0,0,large*scale,(large-1)*scale,
            (large-1)*scale,(large-2)*scale,scale,2*scale}});
        native_case(std::array<double,8>{{0,0,scale,scale,0,2*scale,
            scale,std::nextafter(scale,2*scale)}});
    }
    // Nearly equal distances and a dot product near its sign boundary.
    native_case(std::array<double,8>{{0,0,1,0,0,0,1,std::ldexp(1.0,-26)}});
    native_case(std::array<double,8>{{0,0,1,1,0,std::nextafter(2.0,3.0),1,0}});
}

void native_inputs()
{
    for(double v:{0.0,-0.0,1.0,-1.0,0.1,std::ldexp(1.0,100),std::ldexp(1.0,-100)})conversion(v);
    for(float v:{0.0f,-0.0f,0.1f,1e20f,1e-20f})conversion(v);
    for(std::int64_t v:{std::int64_t(0),std::int64_t(1),std::int64_t(9007199254740991LL),
        std::int64_t(9007199254740993LL),std::int64_t(-9007199254740993LL),
        (std::numeric_limits<std::int64_t>::min)(),(std::numeric_limits<std::int64_t>::max)()})conversion(v);
    conversion((std::numeric_limits<std::uint64_t>::max)());
    BOOST_TEST(native_case(std::array<double,8>{{0,0,3,3,0,3,3,0}}));
    BOOST_TEST(native_case(std::array<float,8>{{0.1f,0.2f,3.1f,3.2f,0.1f,3.2f,3.1f,0.2f}}));
    BOOST_TEST(native_case(std::array<std::int32_t,8>{{0,0,3,3,0,3,3,0}}));
    BOOST_TEST(native_case(std::array<std::int64_t,8>{{0,0,3,3,0,3,3,0}}));
    auto const large=std::int64_t(1)<<60;
    BOOST_TEST(native_case(std::array<std::int64_t,8>{{1,1,large-1,large-3,1,large-5,large-7,1}}));
    auto const maximum=(std::numeric_limits<std::uint64_t>::max)();
    BOOST_TEST(native_case(std::array<std::uint64_t,8>{{0,0,maximum,maximum,0,maximum,maximum,0}}));
    // Conversion uncertainty near a huge common origin must not disappear.
    native_case(std::array<std::int64_t,8>{{large+1,large+3,large+1001,large+2999,
        large+1,large+3001,large+2003,large+3}});
    BOOST_TEST(!native_case(std::array<double,8>{{0,0,1,1,0,1,1,2}}));
    BOOST_TEST(!native_case(std::array<double,8>{{0,0,1,1,0,1,1,std::nextafter(2.0,3.0)}}));
}

void perturbed_inputs()
{
    std::array<point,4> p{{{{{0,0}},{{1e-6,2e-6}}},{{{3,3}},{{2e-6,1e-6}}},
                           {{{0,3}},{{1e-6,1e-6}}},{{{3,0}},{{2e-6,2e-6}}}}};
    point approximation;
    BOOST_TEST(allocation_free([&]{return fe::try_intersection(p,approximation);}));
    std::array<scalar,2> parameters;
    for(unsigned second=0;second<2;++second)
        BOOST_TEST(allocation_free([&]{return fe::try_segment_parameter(p,second!=0,parameters[second]);}));
    // All 256 conversion/perturbation corners, not just nominal centers.
    for(unsigned mask=0;mask<256;++mask)
    {
        std::array<exact_point,4> x;
        for(unsigned i=0;i<4;++i)for(unsigned d=0;d<2;++d)
            x[i][d]=ec::exact_value(p[i].value[d])+((mask&(1u<<(2*i+d)))?1:-1)*ec::exact_value(p[i].error[d]);
        exact_point exact;
        BOOST_TEST(intersection(x,exact));bound(approximation,exact);predicates(p,x);
        for(unsigned second=0;second<2;++second)bound(parameters[second],parameter(x,second!=0));
    }
}

void randomized()
{
    std::mt19937_64 rng(624826);
    unsigned accepted=0;
    for(unsigned i=0;i<1024;++i)
    {
        std::array<double,8> values;
        for(auto& v:values)v=(double(rng()%2000001)-1000000)/1024;
        accepted+=native_case(values);
        std::array<point,4> p;
        std::array<exact_point,4> x;
        for(unsigned j=0;j<4;++j)for(unsigned d=0;d<2;++d)
        {
            double const v=values[2*j+d];
            double const e=std::ldexp(1.0,-40)*(1+std::abs(v));
            p[j].value[d]=v;p[j].error[d]=e;
            x[j][d]=ec::exact_value(v)+((rng()&1)?q(e)/2:-q(e)/2);
        }
        predicates(p,x);
        point result=sentinel();exact_point exact;
        if(allocation_free([&]{return fe::try_intersection(p,result);}))
        {BOOST_TEST(intersection(x,exact));bound(result,exact);}else unchanged(result);
    }
    BOOST_TEST(accepted>1000);
    std::cout<<"native forward intersections certified "<<accepted<<"/1024\n";
    // Near cancellation over many scales, including barely nonparallel lines.
    for(int exponent:{-100,-20,0,20,100})
    {
        double const scale=std::ldexp(1.0,exponent);
        for(int offset=-4;offset<=4;++offset)
        {
            double const a=std::ldexp(1.0,40)+offset;
            native_case(std::array<double,8>{{0,0,a*scale,(a+1)*scale,
                scale,0,(a+1)*scale,(a+2)*scale}});
        }
    }
}

q ring_reference(std::vector<exact_point> const& p,bool close)
{
    q sum=0;
    for(std::size_t i=1;i<p.size();++i)sum+=(p[i-1][0]+p[i][0])*(p[i-1][1]-p[i][1]);
    if(close&&!p.empty())sum+=(p.back()[0]+p.front()[0])*(p.back()[1]-p.front()[1]);
    return sum/2;
}
void ring_tests()
{
    for(bool closed:{false,true})for(bool close:{false,true})for(bool reverse:{false,true})
    {
        std::vector<std::array<std::int64_t,2>> native{{-3,-2},{-3,4},{5,4},{5,-2}};
        if(closed)native.push_back(native.front());
        if(reverse)std::reverse(native.begin(),native.end());
        std::vector<point> p;std::vector<exact_point> x;
        for(auto const& v:native)
        {point a;BOOST_TEST(fe::try_native_point(v[0],v[1],a));p.push_back(a);x.push_back({{q(v[0]),q(v[1])}});}
        scalar result;
        BOOST_TEST(allocation_free([&]{return fe::try_ring_measure(p.begin(),p.end(),result,close);}));
        bound(result,ring_reference(x,close));
        auto reader=[](auto const& v,point& a){return fe::try_native_point(v[0],v[1],a);};
        BOOST_TEST(allocation_free([&]{return fe::try_ring_measure(native.begin(),native.end(),reader,result,close);}));
        bound(result,ring_reference(x,close));
    }
    std::vector<point> p;std::vector<exact_point> x;
    std::mt19937_64 rng(9123);
    for(unsigned i=0;i<1024;++i)
    {
        std::int64_t const v=(std::int64_t(1)<<53)+std::int64_t(rng()%2000001)-1000000;
        std::int64_t const w=(std::int64_t(1)<<53)+std::int64_t(rng()%2000001)-1000000;
        point a;BOOST_TEST(fe::try_native_point(v,w,a));p.push_back(a);x.push_back({{q(v),q(w)}});
    }
    scalar result;
    BOOST_TEST(allocation_free([&]{return fe::try_ring_measure(p.begin(),p.end(),result);}));
    bound(result,ring_reference(x,true));
    std::vector<point> empty;
    BOOST_TEST(fe::try_ring_measure(empty.begin(),empty.end(),result));BOOST_TEST_EQ(result.value,0);BOOST_TEST_EQ(result.error,0);
}

void invalid_and_modes()
{
    std::array<double,8> const native{{0,0,3,3,0,3,3,0}};
    std::array<point,4> p;
    for(unsigned i=0;i<4;++i)BOOST_TEST(fe::try_native_point(native[2*i],native[2*i+1],p[i]));
    for(double invalid:{(std::numeric_limits<double>::infinity)(),(std::numeric_limits<double>::quiet_NaN)()})
    {
        point out=sentinel();BOOST_TEST(!fe::try_native_point(invalid,0,out));unchanged(out);
        auto values=native;values[0]=invalid;
        BOOST_TEST(!fe::try_intersection(values,out));unchanged(out);
        scalar result{17,19};
        BOOST_TEST(!fe::try_segment_parameter(values,false,result));
        BOOST_TEST_EQ(result.value,17);BOOST_TEST_EQ(result.error,19);
    }
    auto invalid=p;invalid[0].error[0]=-1;point out=sentinel();
    BOOST_TEST(!fe::try_intersection(invalid,out));unchanged(out);
    scalar invalid_parameter{17,19};
    BOOST_TEST(!fe::try_segment_parameter(invalid,true,invalid_parameter));
    BOOST_TEST_EQ(invalid_parameter.value,17);BOOST_TEST_EQ(invalid_parameter.error,19);
    scalar invalid_determinant{17,19};
    BOOST_TEST(!fe::try_signed_determinant(invalid[0],invalid[1],invalid[2],invalid_determinant));
    BOOST_TEST_EQ(invalid_determinant.value,17);BOOST_TEST_EQ(invalid_determinant.error,19);
    BOOST_TEST(!fe::try_native_point(1.0L,2.0L,out));unchanged(out);
}

int main()
{
    int const original=std::fegetround();
    if(std::fesetround(FE_TONEAREST)!=0)return 1;
    native_final_coefficients();native_inputs();perturbed_inputs();randomized();ring_tests();invalid_and_modes();
    std::fesetround(original);
    return boost::report_errors();
}
