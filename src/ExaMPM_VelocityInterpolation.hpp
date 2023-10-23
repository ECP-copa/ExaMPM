/****************************************************************************
 * Copyright (c) 2018-2020 by the ExaMPM authors                            *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the ExaMPM library. ExaMPM is distributed under a   *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef EXAMPM_VELOCITYINTERPOLATION_HPP
#define EXAMPM_VELOCITYINTERPOLATION_HPP

#include <ExaMPM_DenseLinearAlgebra.hpp>
#include <ExaMPM_Types.hpp>

#include <Cabana_Grid.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

#include <cmath>
#include <type_traits>

namespace ExaMPM
{
namespace APIC
{
//---------------------------------------------------------------------------//
// Inertial tensor scale factor.
template <class SplineDataType>
KOKKOS_INLINE_FUNCTION typename SplineDataType::scalar_type inertialScaling(
    const SplineDataType& sd,
    typename std::enable_if<( 2 == SplineDataType::order ), void*>::type = 0 )
{
    return 4.0 / ( sd.dx[0] * sd.dx[0] );
}

template <class SplineDataType>
KOKKOS_INLINE_FUNCTION typename SplineDataType::scalar_type inertialScaling(
    const SplineDataType& sd,
    typename std::enable_if<( 3 == SplineDataType::order ), void*>::type = 0 )
{
    return 3.0 / ( sd.dx[0] * sd.dx[0] );
}

//---------------------------------------------------------------------------//
// Interpolate particle momentum to the nodes. (Second and Third order
// splines)
template <class SplineDataType, class MomentumView>
KOKKOS_INLINE_FUNCTION void
p2g( const typename MomentumView::original_value_type m_p,
     const typename MomentumView::original_value_type u_p[3],
     const typename MomentumView::original_value_type B_p[3][3],
     const SplineDataType& sd, const MomentumView& node_momentum,
     typename std::enable_if<
         ( Cabana::Grid::isNode<typename SplineDataType::entity_type>::value &&
           ( SplineDataType::order == 2 || SplineDataType::order == 3 ) ),
         void*>::type = 0 )
{
    static_assert( Cabana::Grid::P2G::is_scatter_view<MomentumView>::value,
                   "P2G requires a Kokkos::ScatterView" );
    auto momentum_access = node_momentum.access();

    using value_type = typename MomentumView::original_value_type;

    // Scaling factor from inertial tensor with quadratic shape
    // functions.
    value_type D_p_inv = inertialScaling( sd );

    // Project momentum.
    value_type distance[3];
    value_type B_p_d[3];
    value_type wm_ip;
    for ( int i = 0; i < SplineDataType::num_knot; ++i )
        for ( int j = 0; j < SplineDataType::num_knot; ++j )
            for ( int k = 0; k < SplineDataType::num_knot; ++k )
            {
                // Physical distance to entity.
                distance[Dim::I] = sd.d[Dim::I][i];
                distance[Dim::J] = sd.d[Dim::J][j];
                distance[Dim::K] = sd.d[Dim::K][k];

                // Compute the action of B_p on the distance.
                DenseLinearAlgebra::matVecMultiply( B_p, distance, B_p_d );

                // Weight times mass.
                wm_ip =
                    sd.w[Dim::I][i] * sd.w[Dim::J][j] * sd.w[Dim::K][k] * m_p;

                // Interpolate particle momentum to the entity.
                for ( int d = 0; d < 3; ++d )
                    momentum_access( sd.s[Dim::I][i], sd.s[Dim::J][j],
                                     sd.s[Dim::K][k], d ) +=
                        wm_ip * ( u_p[d] + D_p_inv * B_p_d[d] );
            }
}

//---------------------------------------------------------------------------//
// Interpolate particle momentum to the nodes. (First order splines)
template <class SplineDataType, class MomentumView>
KOKKOS_INLINE_FUNCTION void
p2g( const typename MomentumView::original_value_type m_p,
     const typename MomentumView::original_value_type u_p[3],
     const typename MomentumView::original_value_type B_p[3][3],
     const SplineDataType& sd, const MomentumView& node_momentum,
     typename std::enable_if<
         ( Cabana::Grid::isNode<typename SplineDataType::entity_type>::value &&
           ( SplineDataType::order == 1 ) ),
         void*>::type = 0 )
{
    static_assert( Cabana::Grid::P2G::is_scatter_view<MomentumView>::value,
                   "P2G requires a Kokkos::ScatterView" );
    auto momentum_access = node_momentum.access();

    using value_type = typename MomentumView::original_value_type;

    // Project momentum.
    value_type B_g_d[3];
    value_type wm_ip;
    value_type gm_ip[3];
    for ( int i = 0; i < SplineDataType::num_knot; ++i )
        for ( int j = 0; j < SplineDataType::num_knot; ++j )
            for ( int k = 0; k < SplineDataType::num_knot; ++k )
            {
                // Weight times mass.
                wm_ip =
                    sd.w[Dim::I][i] * sd.w[Dim::J][j] * sd.w[Dim::K][k] * m_p;

                // Weight gradient times mass.
                gm_ip[0] =
                    sd.g[Dim::I][i] * sd.w[Dim::J][j] * sd.w[Dim::K][k] * m_p;
                gm_ip[1] =
                    sd.w[Dim::I][i] * sd.g[Dim::J][j] * sd.w[Dim::K][k] * m_p;
                gm_ip[2] =
                    sd.w[Dim::I][i] * sd.w[Dim::J][j] * sd.g[Dim::K][k] * m_p;

                // Compute the action of B_p on the gradient.
                DenseLinearAlgebra::matVecMultiply( B_p, gm_ip, B_g_d );

                // Interpolate particle momentum to the entity.
                for ( int d = 0; d < 3; ++d )
                    momentum_access( sd.s[Dim::I][i], sd.s[Dim::J][j],
                                     sd.s[Dim::K][k], d ) +=
                        wm_ip * u_p[d] + B_g_d[d];
            }
}

//---------------------------------------------------------------------------//
// Interpolate grid node velocity to the particle.
template <class SplineDataType, class VelocityView>
KOKKOS_INLINE_FUNCTION void
g2p( const VelocityView& node_velocity, const SplineDataType& sd,
     typename VelocityView::value_type u_p[3],
     typename VelocityView::value_type B_p[3][3],
     typename std::enable_if<
         Cabana::Grid::isNode<typename SplineDataType::entity_type>::value,
         void*>::type = 0 )
{
    using value_type = typename VelocityView::value_type;

    for ( int d = 0; d < 3; ++d )
        u_p[d] = 0.0;

    for ( int d0 = 0; d0 < 3; ++d0 )
        for ( int d1 = 0; d1 < 3; ++d1 )
            B_p[d0][d1] = 0.0;

    value_type distance[3];
    value_type w_ip;

    for ( int i = 0; i < SplineDataType::num_knot; ++i )
        for ( int j = 0; j < SplineDataType::num_knot; ++j )
            for ( int k = 0; k < SplineDataType::num_knot; ++k )
            {
                // Projection weight.
                w_ip = sd.w[Dim::I][i] * sd.w[Dim::J][j] * sd.w[Dim::K][k];

                // Update velocity.
                for ( int d = 0; d < 3; ++d )
                    u_p[d] +=
                        w_ip * node_velocity( sd.s[Dim::I][i], sd.s[Dim::J][j],
                                              sd.s[Dim::K][k], d );

                // Physical distance to entity.
                distance[Dim::I] = sd.d[Dim::I][i];
                distance[Dim::J] = sd.d[Dim::J][j];
                distance[Dim::K] = sd.d[Dim::K][k];

                // Update affine matrix.
                for ( int d0 = 0; d0 < 3; ++d0 )
                    for ( int d1 = 0; d1 < 3; ++d1 )
                        B_p[d0][d1] +=
                            w_ip *
                            node_velocity( sd.s[Dim::I][i], sd.s[Dim::J][j],
                                           sd.s[Dim::K][k], d0 ) *
                            distance[d1];
            }
}

//---------------------------------------------------------------------------//

} // end namespace APIC
} // end namespace ExaMPM

#endif // end EXAMPM_VELOCITYINTERPOLATION_HPP
