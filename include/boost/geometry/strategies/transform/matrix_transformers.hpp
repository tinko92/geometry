// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2007-2015 Barend Gehrels, Amsterdam, the Netherlands.
// Copyright (c) 2008-2015 Bruno Lalande, Paris, France.
// Copyright (c) 2009-2015 Mateusz Loskot, London, UK.

// This file was modified by Oracle on 2015.
// Modifications copyright (c) 2015 Oracle and/or its affiliates.

// Contributed and/or modified by Menelaos Karavelas, on behalf of Oracle

// Parts of Boost.Geometry are redesigned from Geodan's Geographic Library
// (geolib/GGL), copyright (c) 1995-2010 Geodan, Amsterdam, the Netherlands.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_STRATEGIES_TRANSFORM_MATRIX_TRANSFORMERS_HPP
#define BOOST_GEOMETRY_STRATEGIES_TRANSFORM_MATRIX_TRANSFORMERS_HPP


#include <cstddef>

#include <boost/qvm/mat.hpp>
#include <boost/qvm/vec.hpp>
#include <boost/qvm/mat_access.hpp>
#include <boost/qvm/vec_access.hpp>
#include <boost/qvm/mat_operations.hpp>
#include <boost/qvm/vec_mat_operations.hpp>

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/util/math.hpp>
#include <boost/geometry/util/promote_floating_point.hpp>
#include <boost/geometry/util/select_coordinate_type.hpp>
#include <boost/geometry/util/select_most_precise.hpp>


namespace boost { namespace geometry
{

namespace strategy { namespace transform
{

namespace detail { namespace matrix_transformer {

template
<
    typename CalculationType,
    typename P2,
    std::size_t Dimension2,
    std::size_t Dimension2Count
>
struct set_point_from_vec {
    typedef CalculationType ct;
    static inline void apply(P2& p2, qvm::vec<ct,Dimension2 + 1>& p2temp)
    {
        typedef typename geometry::coordinate_type<P2>::type ct2;
        set<Dimension2Count>(p2, boost::numeric_cast<ct2>(qvm::A<Dimension2Count>(p2temp)));
        set_point_from_vec<CalculationType, P2, Dimension2, Dimension2Count-1>::apply(p2,p2temp);
    }
};

template
<
    typename CalculationType,
    typename P2,
    std::size_t Dimension2
>
struct set_point_from_vec<CalculationType, P2, Dimension2, -1> {
    typedef CalculationType ct;
    static inline void apply(P2& p2, qvm::vec<ct,Dimension2+1>& p2temp) {}
};

template
<
    typename CalculationType,
    typename P1,
    std::size_t Dimension1,
    std::size_t Dimension1Count
>
struct set_vec_from_point {
    typedef CalculationType ct;
    static inline void apply(const P1& p1, qvm::vec<ct,Dimension1+1>& p1temp)
    {
        qvm::A<Dimension1Count>(p1temp) = get<Dimension1Count>(p1);
        set_vec_from_point<CalculationType, P1, Dimension1, Dimension1Count-1>::apply(p1,p1temp);
    }
};

template
<
    typename CalculationType,
    typename P1,
    std::size_t Dimension1
>
struct set_vec_from_point<CalculationType, P1, Dimension1, -1> {
    typedef CalculationType ct;
    static inline void apply(const P1& p1, qvm::vec<ct,Dimension1+1>& p1temp) {}
};


template
<
    typename CalculationType,
    std::size_t Dimension1,
    std::size_t Dimension2
>
class matrix_transformer
{
protected :
    typedef CalculationType ct;
    typedef boost::qvm::mat<ct, Dimension2+1, Dimension1+1> matrix_type;
    matrix_type m_matrix;
public :

    matrix_type const& matrix() const { return m_matrix; }
    template <typename P1, typename P2>
    inline bool apply(P1 const& p1, P2& p2) const
    {
        assert_dimension<P1,Dimension1>();
        assert_dimension<P2,Dimension2>();
        qvm::vec<ct,Dimension1+1> p1temp;
        qvm::vec<ct,Dimension2+1> p2temp;
        set_vec_from_point<ct,P1,Dimension1, Dimension1-1>::apply(p1, p1temp);
        qvm::A<Dimension1>(p1temp) = 1;
        p2temp = m_matrix*p1temp;
        set_point_from_vec<ct,P2,Dimension2, Dimension2-1>::apply(p2, p2temp);
        return true;
    }

};

}} // namespace detail::matrix_transform

/*!
\brief Affine transformation strategy in Cartesian system.
\details The strategy serves as a generic definition of affine transformation matrix
         and procedure of application it to given point.
\see http://en.wikipedia.org/wiki/Affine_transformation
     and http://www.devmaster.net/wiki/Transformation_matrices
\ingroup strategies
\tparam Dimension1 number of dimensions to transform from
\tparam Dimension2 number of dimensions to transform to
 */
template
<
    typename CalculationType,
    std::size_t Dimension1,
    std::size_t Dimension2
>
class matrix_transformer : public detail::matrix_transformer::matrix_transformer<CalculationType,Dimension1,Dimension2>
{
protected:
    typedef CalculationType ct;
    typedef boost::qvm::mat<ct, Dimension2+1, Dimension1+1> matrix_type;
public:

    template<typename Matrix>
    inline matrix_transformer(Matrix const& matrix)
    {
        qvm::assign(this->m_matrix,matrix);
    }
    inline matrix_transformer() {}
};


template <typename CalculationType>
class matrix_transformer<CalculationType, 2, 2> : public detail::matrix_transformer::matrix_transformer<CalculationType, 2, 2>
{
protected :
    typedef CalculationType ct;
    typedef boost::qvm::mat<ct, 3, 3> matrix_type;
public :

    template<typename Matrix>
    inline matrix_transformer(Matrix const& matrix)
    {
        qvm::assign(this->m_matrix,matrix);
    }
    inline matrix_transformer() {}
    inline matrix_transformer(
                ct const& m_0_0, ct const& m_0_1, ct const& m_0_2,
                ct const& m_1_0, ct const& m_1_1, ct const& m_1_2,
                ct const& m_2_0, ct const& m_2_1, ct const& m_2_2)
    {
        qvm::A<0,0>(this->m_matrix) = m_0_0;   qvm::A<0,1>(this->m_matrix) = m_0_1;   qvm::A<0,2>(this->m_matrix) = m_0_2;
        qvm::A<1,0>(this->m_matrix) = m_1_0;   qvm::A<1,1>(this->m_matrix) = m_1_1;   qvm::A<1,2>(this->m_matrix) = m_1_2;
        qvm::A<2,0>(this->m_matrix) = m_2_0;   qvm::A<2,1>(this->m_matrix) = m_2_1;   qvm::A<2,2>(this->m_matrix) = m_2_2;
    }
};


// It IS possible to go from 3 to 2 coordinates
template <typename CalculationType>
class matrix_transformer<CalculationType, 3, 2> : public detail::matrix_transformer::matrix_transformer<CalculationType, 3, 2>
{
    typedef CalculationType ct;
    typedef boost::qvm::mat<ct, 3, 4> matrix_type;
public :

    template<typename Matrix>
    inline matrix_transformer(Matrix const& matrix)
    {
        qvm::assign(this->m_matrix,matrix);
    }
    inline matrix_transformer() {}
    inline matrix_transformer(
                ct const& m_0_0, ct const& m_0_1, ct const& m_0_2,
                ct const& m_1_0, ct const& m_1_1, ct const& m_1_2,
                ct const& m_2_0, ct const& m_2_1, ct const& m_2_2)
    {
        qvm::A<0,0>(this->m_matrix) = m_0_0;   qvm::A<0,1>(this->m_matrix) = m_0_1;   qvm::A<0,2>(this->m_matrix) = 0;   qvm::A<0,3>(this->m_matrix) = m_0_2;
        qvm::A<1,0>(this->m_matrix) = m_1_0;   qvm::A<1,1>(this->m_matrix) = m_1_1;   qvm::A<1,2>(this->m_matrix) = 0;   qvm::A<1,3>(this->m_matrix) = m_1_2;
        qvm::A<2,0>(this->m_matrix) = m_2_0;   qvm::A<2,1>(this->m_matrix) = m_2_1;   qvm::A<2,2>(this->m_matrix) = 0;   qvm::A<2,3>(this->m_matrix) = m_2_2;
    }
};


template <typename CalculationType>
class matrix_transformer<CalculationType, 3, 3> : public detail::matrix_transformer::matrix_transformer<CalculationType, 3, 3>
{
protected :
    typedef CalculationType ct;
    typedef boost::qvm::mat<ct, 4, 4> matrix_type;
public :

    template<typename Matrix>
    inline matrix_transformer(Matrix const& matrix)
    {
        qvm::assign(this->m_matrix,matrix);
    }
    inline matrix_transformer() {}
    inline matrix_transformer(
                ct const& m_0_0, ct const& m_0_1, ct const& m_0_2, ct const& m_0_3,
                ct const& m_1_0, ct const& m_1_1, ct const& m_1_2, ct const& m_1_3,
                ct const& m_2_0, ct const& m_2_1, ct const& m_2_2, ct const& m_2_3,
                ct const& m_3_0, ct const& m_3_1, ct const& m_3_2, ct const& m_3_3
                )
    {
        qvm::A<0,0>(this->m_matrix) = m_0_0; qvm::A<0,1>(this->m_matrix) = m_0_1; qvm::A<0,2>(this->m_matrix) = m_0_2; qvm::A<0,3>(this->m_matrix) = m_0_3;
        qvm::A<1,0>(this->m_matrix) = m_1_0; qvm::A<1,1>(this->m_matrix) = m_1_1; qvm::A<1,2>(this->m_matrix) = m_1_2; qvm::A<1,3>(this->m_matrix) = m_1_3;
        qvm::A<2,0>(this->m_matrix) = m_2_0; qvm::A<2,1>(this->m_matrix) = m_2_1; qvm::A<2,2>(this->m_matrix) = m_2_2; qvm::A<2,3>(this->m_matrix) = m_2_3;
        qvm::A<3,0>(this->m_matrix) = m_3_0; qvm::A<3,1>(this->m_matrix) = m_3_1; qvm::A<3,2>(this->m_matrix) = m_3_2; qvm::A<3,3>(this->m_matrix) = m_3_3;
    }
};


/*!
\brief Strategy of translate transformation in Cartesian system.
\details Translate moves a geometry a fixed distance in 2 or 3 dimensions.
\see http://en.wikipedia.org/wiki/Translation_%28geometry%29
\ingroup strategies
\tparam Dimension1 number of dimensions to transform from
\tparam Dimension2 number of dimensions to transform to
 */
template
<
    typename CalculationType,
    std::size_t Dimension1,
    std::size_t Dimension2
>
class translate_transformer
{
};


template<typename CalculationType>
class translate_transformer<CalculationType, 2, 2> : public matrix_transformer<CalculationType, 2, 2>
{
public :
    // To have translate transformers compatible for 2/3 dimensions, the
    // constructor takes an optional third argument doing nothing.
    inline translate_transformer(CalculationType const& translate_x,
                CalculationType const& translate_y,
                CalculationType const& = 0)
        : matrix_transformer<CalculationType, 2, 2>(
                1, 0, translate_x,
                0, 1, translate_y,
                0, 0, 1)
    {}
};


template <typename CalculationType>
class translate_transformer<CalculationType, 3, 3> : public matrix_transformer<CalculationType, 3, 3>
{
public :
    inline translate_transformer(CalculationType const& translate_x,
                CalculationType const& translate_y,
                CalculationType const& translate_z)
        : matrix_transformer<CalculationType, 3, 3>(
                1, 0, 0, translate_x,
                0, 1, 0, translate_y,
                0, 0, 1, translate_z,
                0, 0, 0, 1)
    {}

};


/*!
\brief Strategy of scale transformation in Cartesian system.
\details Scale scales a geometry up or down in all its dimensions.
\see http://en.wikipedia.org/wiki/Scaling_%28geometry%29
\ingroup strategies
\tparam Dimension1 number of dimensions to transform from
\tparam Dimension2 number of dimensions to transform to
*/
template
<
    typename CalculationType,
    std::size_t Dimension1,
    std::size_t Dimension2
>
class scale_transformer
{
};

template
<
    typename CalculationType,
    std::size_t Dimension1
>
class scale_transformer<CalculationType, Dimension1, Dimension1> : public matrix_transformer<CalculationType, Dimension1, Dimension1>
{
public:
    inline scale_transformer(CalculationType const& scale)
    {
        boost::qvm::set_identity(this->m_matrix);
        this->m_matrix*=scale;
        qvm::A<Dimension1,Dimension1>(this->m_matrix) = 1;
    }
};

template <typename CalculationType>
class scale_transformer<CalculationType, 2, 2> : public matrix_transformer<CalculationType, 2, 2>
{

public :
    inline scale_transformer(CalculationType const& scale_x,
                CalculationType const& scale_y,
                CalculationType const& = 0)
        : matrix_transformer<CalculationType, 2, 2>(
                scale_x, 0,       0,
                0,       scale_y, 0,
                0,       0,       1)
    {}


    inline scale_transformer(CalculationType const& scale)
        : matrix_transformer<CalculationType, 2, 2>(
                scale, 0,     0,
                0,     scale, 0,
                0,     0,     1)
    {}
};


template <typename CalculationType>
class scale_transformer<CalculationType, 3, 3> : public matrix_transformer<CalculationType, 3, 3>
{
public :
    inline scale_transformer(CalculationType const& scale_x,
                CalculationType const& scale_y,
                CalculationType const& scale_z)
        : matrix_transformer<CalculationType, 3, 3>(
                scale_x, 0,       0,       0,
                0,       scale_y, 0,       0,
                0,       0,       scale_z, 0,
                0,       0,       0,       1)
    {}


    inline scale_transformer(CalculationType const& scale)
        : matrix_transformer<CalculationType, 3, 3>(
                scale, 0,     0,     0,
                0,     scale, 0,     0,
                0,     0,     scale, 0,
                0,     0,     0,     1)
    {}
};


#ifndef DOXYGEN_NO_DETAIL
namespace detail
{


template <typename DegreeOrRadian>
struct as_radian
{};


template <>
struct as_radian<radian>
{
    template <typename T>
    static inline T get(T const& value)
    {
        return value;
    }
};

template <>
struct as_radian<degree>
{
    template <typename T>
    static inline T get(T const& value)
    {
        typedef typename promote_floating_point<T>::type promoted_type;
        return value * math::d2r<promoted_type>();
    }

};


template
<
    typename CalculationType,
    std::size_t Dimension1,
    std::size_t Dimension2
>
class rad_rotate_transformer
    : public matrix_transformer<CalculationType, Dimension1, Dimension2>
{
public :
    inline rad_rotate_transformer(CalculationType const& angle)
        : matrix_transformer<CalculationType, Dimension1, Dimension2>(
                 cos(angle), sin(angle), 0,
                -sin(angle), cos(angle), 0,
                 0,          0,          1)
    {}
};


} // namespace detail
#endif // DOXYGEN_NO_DETAIL


/*!
\brief Strategy for rotate transformation in Cartesian coordinate system.
\details Rotate rotates a geometry of specified angle about a fixed point (e.g. origin).
\see http://en.wikipedia.org/wiki/Rotation_%28mathematics%29
\ingroup strategies
\tparam DegreeOrRadian degree/or/radian, type of rotation angle specification
\note A single angle is needed to specify a rotation in 2D.
      Not yet in 3D, the 3D version requires special things to allow
      for rotation around X, Y, Z or arbitrary axis.
\todo The 3D version will not compile.
 */
template
<
    typename DegreeOrRadian,
    typename CalculationType,
    std::size_t Dimension1,
    std::size_t Dimension2
>
class rotate_transformer : public detail::rad_rotate_transformer<CalculationType, Dimension1, Dimension2>
{

public :
    inline rotate_transformer(CalculationType const& angle)
        : detail::rad_rotate_transformer
            <
                CalculationType, Dimension1, Dimension2
            >(detail::as_radian<DegreeOrRadian>::get(angle))
    {}
};


}} // namespace strategy::transform


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_STRATEGIES_TRANSFORM_MATRIX_TRANSFORMERS_HPP
