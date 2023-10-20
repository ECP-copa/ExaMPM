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

#ifndef EXAMPM_TIMEINTEGRATOR_HPP
#define EXAMPM_TIMEINTEGRATOR_HPP

#include <ExaMPM_BoundaryConditions.hpp>
#include <ExaMPM_ProblemManager.hpp>
#include <ExaMPM_VelocityInterpolation.hpp>

#include <Cabana_Grid.hpp>

#include <Kokkos_Core.hpp>

#include <cmath>

namespace ExaMPM
{
namespace TimeIntegrator
{
//---------------------------------------------------------------------------//
// Particle-to-grid.
template <class ProblemManagerType, class ExecutionSpace>
void p2g( const ExecutionSpace& exec_space, const ProblemManagerType& pm )
{
    // Get the particle data we need.
    auto m_p = pm.get( Location::Particle(), Field::Mass() );
    auto u_p = pm.get( Location::Particle(), Field::Velocity() );
    auto B_p = pm.get( Location::Particle(), Field::Affine() );
    auto x_p = pm.get( Location::Particle(), Field::Position() );
    auto v_p = pm.get( Location::Particle(), Field::Volume() );
    auto j_p = pm.get( Location::Particle(), Field::J() );

    // Get the views we need.
    auto m_i = pm.get( Location::Node(), Field::Mass() );
    auto mu_i = pm.get( Location::Node(), Field::Momentum() );
    auto f_i = pm.get( Location::Node(), Field::Force() );

    // Reset write views.
    Kokkos::deep_copy( m_i, 0.0 );
    Kokkos::deep_copy( mu_i, 0.0 );
    Kokkos::deep_copy( f_i, 0.0 );

    // Create the scatter views we need.
    auto m_i_sv = Kokkos::Experimental::create_scatter_view( m_i );
    auto mu_i_sv = Kokkos::Experimental::create_scatter_view( mu_i );
    auto f_i_sv = Kokkos::Experimental::create_scatter_view( f_i );

    // Get the fluid properties.
    double bulk_mod = pm.bulkModulus();
    double gamma = pm.gamma();

    // Build the local mesh.
    auto local_mesh = Cabana::Grid::createLocalMesh<ExecutionSpace>(
        *( pm.mesh()->localGrid() ) );

    // Loop over particles.
    Kokkos::parallel_for(
        "p2g",
        Kokkos::RangePolicy<ExecutionSpace>( exec_space, 0, pm.numParticle() ),
        KOKKOS_LAMBDA( const int p ) {
            // Get the particle position.
            double x[3] = { x_p( p, 0 ), x_p( p, 1 ), x_p( p, 2 ) };

            // Setup interpolation to the nodes.
            Cabana::Grid::SplineData<double, 2, 3, Cabana::Grid::Node> sd;
            Cabana::Grid::evaluateSpline( local_mesh, x, sd );

            // Compute the pressure on the particle with an equation of
            // state.
            double pressure = -bulk_mod * ( pow( j_p( p ), -gamma ) - 1.0 );

            // Project the pressure gradient to the grid.
            Cabana::Grid::P2G::gradient( -v_p( p ) * j_p( p ) * pressure, sd,
                                         f_i_sv );

            // Extract the particle velocity
            double vel_p[3] = { u_p( p, 0 ), u_p( p, 1 ), u_p( p, 2 ) };

            // Extract the affine particle matrix.
            double aff_p[3][3];
            for ( int d0 = 0; d0 < 3; ++d0 )
                for ( int d1 = 0; d1 < 3; ++d1 )
                    aff_p[d0][d1] = B_p( p, d0, d1 );

            // Project momentum to the grid.
            APIC::p2g( m_p( p ), vel_p, aff_p, sd, mu_i_sv );

            // Project mass to the grid.
            Cabana::Grid::P2G::value( m_p( p ), sd, m_i_sv );
        } );

    // Complete local scatter.
    Kokkos::Experimental::contribute( m_i, m_i_sv );
    Kokkos::Experimental::contribute( mu_i, mu_i_sv );
    Kokkos::Experimental::contribute( f_i, f_i_sv );

    // Complete global scatter.
    pm.scatter( Location::Node() );
}

//---------------------------------------------------------------------------//
// Field solve.
template <class ProblemManagerType, class ExecutionSpace>
void fieldSolve( const ExecutionSpace& exec_space, const ProblemManagerType& pm,
                 const double delta_t, const double gravity,
                 const BoundaryCondition& bc )
{
    // Get the views we need.
    auto m_i = pm.get( Location::Node(), Field::Mass() );
    auto mu_i = pm.get( Location::Node(), Field::Momentum() );
    auto f_i = pm.get( Location::Node(), Field::Force() );
    auto u_i = pm.get( Location::Node(), Field::Velocity() );

    // Velocity increment from gravity.
    double delta_g = delta_t * gravity;

    // Node mass epsilon. Masses smaller than this will be ignored.
    double mass_epsilon = 1.0e-12;

    // Compute the velocity.
    auto l2g = Cabana::Grid::IndexConversion::createL2G(
        *( pm.mesh()->localGrid() ), Cabana::Grid::Node() );
    auto local_nodes = pm.mesh()->localGrid()->indexSpace(
        Cabana::Grid::Ghost(), Cabana::Grid::Node(), Cabana::Grid::Local() );
    Kokkos::parallel_for(
        Cabana::Grid::createExecutionPolicy( local_nodes, exec_space ),
        KOKKOS_LAMBDA( const int li, const int lj, const int lk ) {
            int gi, gj, gk;
            l2g( li, lj, lk, gi, gj, gk );

            // Only compute velocity if a node has mass
            u_i( li, lj, lk, 0 ) = ( m_i( li, lj, lk, 0 ) > mass_epsilon )
                                       ? ( mu_i( li, lj, lk, 0 ) +
                                           delta_t * f_i( li, lj, lk, 0 ) ) /
                                             m_i( li, lj, lk, 0 )
                                       : 0.0;
            u_i( li, lj, lk, 1 ) = ( m_i( li, lj, lk, 0 ) > mass_epsilon )
                                       ? ( mu_i( li, lj, lk, 1 ) +
                                           delta_t * f_i( li, lj, lk, 1 ) ) /
                                             m_i( li, lj, lk, 0 )
                                       : 0.0;
            u_i( li, lj, lk, 2 ) = ( m_i( li, lj, lk, 0 ) > mass_epsilon )
                                       ? ( mu_i( li, lj, lk, 2 ) +
                                           delta_t * f_i( li, lj, lk, 2 ) ) /
                                                 m_i( li, lj, lk, 0 ) -
                                             delta_g
                                       : 0.0;

            // Apply the boundary condition.
            bc( gi, gj, gk, u_i( li, lj, lk, 0 ), u_i( li, lj, lk, 1 ),
                u_i( li, lj, lk, 2 ) );
        } );
}

//---------------------------------------------------------------------------//
// Grid-to-particle.
template <class ProblemManagerType, class ExecutionSpace>
void g2p( const ExecutionSpace& exec_space, const ProblemManagerType& pm,
          const double delta_t )
{
    // Get the particle data we need.
    auto m_p = pm.get( Location::Particle(), Field::Mass() );
    auto B_p = pm.get( Location::Particle(), Field::Affine() );
    auto u_p = pm.get( Location::Particle(), Field::Velocity() );
    auto x_p = pm.get( Location::Particle(), Field::Position() );
    auto j_p = pm.get( Location::Particle(), Field::J() );

    // Get the views we need.
    auto u_i = pm.get( Location::Node(), Field::Velocity() );
    auto r_c = pm.get( Location::Cell(), Field::Density() );
    auto k_c = pm.get( Location::Cell(), Field::Mark() );

    // Reset write views.
    Kokkos::deep_copy( r_c, 0.0 );
    Kokkos::deep_copy( k_c, 0.0 );

    // Create the scatter views we need.
    auto r_c_sv = Kokkos::Experimental::create_scatter_view( r_c );
    auto k_c_sv = Kokkos::Experimental::create_scatter_view( k_c );

    // Build the local mesh.
    auto local_mesh = Cabana::Grid::createLocalMesh<ExecutionSpace>(
        *( pm.mesh()->localGrid() ) );
    auto cell_size =
        pm.mesh()->localGrid()->globalGrid().globalMesh().cellSize( 0 );
    auto cell_volume = cell_size * cell_size * cell_size;

    // Gather the data we need.
    pm.gather( Location::Node() );

    // Loop over particles.
    Kokkos::parallel_for(
        "g2p",
        Kokkos::RangePolicy<ExecutionSpace>( exec_space, 0, pm.numParticle() ),
        KOKKOS_LAMBDA( const int p ) {
            // Get the particle position.
            double x[3] = { x_p( p, 0 ), x_p( p, 1 ), x_p( p, 2 ) };

            // Setup interpolation from the nodes.
            Cabana::Grid::SplineData<double, 2, 3, Cabana::Grid::Node> sd_i;
            Cabana::Grid::evaluateSpline( local_mesh, x, sd_i );

            // Update particle velocity.
            double vel_p[3];
            double aff_p[3][3];
            APIC::g2p( u_i, sd_i, vel_p, aff_p );
            for ( int d = 0; d < 3; ++d )
                u_p( p, d ) = vel_p[d];
            for ( int d0 = 0; d0 < 3; ++d0 )
                for ( int d1 = 0; d1 < 3; ++d1 )
                    B_p( p, d0, d1 ) = aff_p[d0][d1];

            // Compute the velocity divergence (this is the trace of the
            // velocity gradient).
            double div_u;
            Cabana::Grid::G2P::divergence( u_i, sd_i, div_u );

            // Update the deformation gradient determinant.
            j_p( p ) *= exp( delta_t * div_u );

            // Move the particle
            for ( int d = 0; d < 3; ++d )
            {
                x[d] += delta_t * vel_p[d];
                x_p( p, d ) = x[d];
            }

            // Project density to cell.
            Cabana::Grid::SplineData<double, 1, 3, Cabana::Grid::Cell> sd_c1;
            Cabana::Grid::evaluateSpline( local_mesh, x, sd_c1 );
            Cabana::Grid::P2G::value( m_p( p ) / cell_volume, sd_c1, r_c_sv );

            // Mark cells. Indicates whether or not cells have particles.
            Cabana::Grid::SplineData<double, 0, 3, Cabana::Grid::Cell> sd_c0;
            Cabana::Grid::evaluateSpline( local_mesh, x, sd_c0 );
            Cabana::Grid::P2G::value( 1.0, sd_c0, k_c_sv );
        } );

    // Complete local scatter.
    Kokkos::Experimental::contribute( r_c, r_c_sv );
    Kokkos::Experimental::contribute( k_c, k_c_sv );

    // Complete global scatter.
    pm.scatter( Location::Cell() );
}

//---------------------------------------------------------------------------//
// Correct particle positions.
template <class ProblemManagerType, class ExecutionSpace>
void correctParticlePositions( const ExecutionSpace& exec_space,
                               const ProblemManagerType& pm,
                               const double delta_t,
                               const BoundaryCondition& bc )
{
    // Get the particle data we need.
    auto x_p = pm.get( Location::Particle(), Field::Position() );

    // Get the views we need.
    auto r_c = pm.get( Location::Cell(), Field::Density() );
    auto k_c = pm.get( Location::Cell(), Field::Mark() );
    auto x_i = pm.get( Location::Node(), Field::PositionCorrection() );

    // Reset write views.
    Kokkos::deep_copy( x_i, 0.0 );

    // Create the scatter views we need.
    auto x_i_sv = Kokkos::Experimental::create_scatter_view( x_i );

    // Get the fluid properties.
    double kappa = pm.kappa();
    double density = pm.density();

    // Build the local mesh.
    auto local_mesh = Cabana::Grid::createLocalMesh<ExecutionSpace>(
        *( pm.mesh()->localGrid() ) );

    // Compute nodal correction.
    auto local_cells = pm.mesh()->localGrid()->indexSpace(
        Cabana::Grid::Own(), Cabana::Grid::Cell(), Cabana::Grid::Local() );
    Kokkos::parallel_for(
        "compute_position_correction",
        Cabana::Grid::createExecutionPolicy( local_cells, exec_space ),
        KOKKOS_LAMBDA( const int i, const int j, const int k ) {
            // Get the cell center.
            int idx[3] = { i, j, k };
            double x[3];
            local_mesh.coordinates( Cabana::Grid::Cell(), idx, x );

            // Setup interpolation from cell center to nodes.
            Cabana::Grid::SplineData<double, 1, 3, Cabana::Grid::Node> sd_i;
            Cabana::Grid::evaluateSpline( local_mesh, x, sd_i );

            // Clamp the density outside the fluid.
            double rho = ( k_c( i, j, k, 0 ) > 0.0 )
                             ? r_c( i, j, k, 0 )
                             : fmax( r_c( i, j, k, 0 ), density );

            // Compute correction.
            double correction =
                -delta_t * delta_t * kappa * ( 1 - rho / density ) / density;
            Cabana::Grid::P2G::gradient( correction, sd_i, x_i_sv );
        } );

    // Complete local scatter.
    Kokkos::Experimental::contribute( x_i, x_i_sv );

    // Complete global scatter.
    pm.scatter( Location::Node(), Field::PositionCorrection() );

    // Gather the position correction.
    pm.gather( Location::Node(), Field::PositionCorrection() );

    // Apply boundary condition to position correction.
    // Compute the velocity.
    auto l2g = Cabana::Grid::IndexConversion::createL2G(
        *( pm.mesh()->localGrid() ), Cabana::Grid::Node() );
    auto local_nodes = pm.mesh()->localGrid()->indexSpace(
        Cabana::Grid::Ghost(), Cabana::Grid::Node(), Cabana::Grid::Local() );
    Kokkos::parallel_for(
        Cabana::Grid::createExecutionPolicy( local_nodes, exec_space ),
        KOKKOS_LAMBDA( const int li, const int lj, const int lk ) {
            int gi, gj, gk;
            l2g( li, lj, lk, gi, gj, gk );
            bc( gi, gj, gk, x_i( li, lj, lk, 0 ), x_i( li, lj, lk, 1 ),
                x_i( li, lj, lk, 2 ) );
        } );

    // Update particle positions.
    Kokkos::parallel_for(
        "correct_particles",
        Kokkos::RangePolicy<ExecutionSpace>( exec_space, 0, pm.numParticle() ),
        KOKKOS_LAMBDA( const int p ) {
            // Get the particle position.
            double x[3] = { x_p( p, 0 ), x_p( p, 1 ), x_p( p, 2 ) };

            // Setup interpolation from the nodes.
            Cabana::Grid::SplineData<double, 2, 3, Cabana::Grid::Node> sd_i;
            Cabana::Grid::evaluateSpline( local_mesh, x, sd_i );

            // Correct the particle position.
            double delta_x[3];
            Cabana::Grid::G2P::value( x_i, sd_i, delta_x );
            for ( int d = 0; d < 3; ++d )
                x_p( p, d ) += delta_x[d];
        } );
}

//---------------------------------------------------------------------------//
// Take a time step.
template <class ProblemManagerType, class ExecutionSpace>
void step( const ExecutionSpace& exec_space, const ProblemManagerType& pm,
           const double delta_t, const double gravity,
           const BoundaryCondition& bc )
{
    p2g( exec_space, pm );
    fieldSolve( exec_space, pm, delta_t, gravity, bc );
    g2p( exec_space, pm, delta_t );
    correctParticlePositions( exec_space, pm, delta_t, bc );
}

//---------------------------------------------------------------------------//

} // end namespace TimeIntegrator
} // end namespace ExaMPM

#endif // EXAMPM_TIMEINTEGRATOR_HPP
