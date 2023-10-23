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

#ifndef EXAMPM_TIMESTEPCONTROL_HPP
#define EXAMPM_TIMESTEPCONTROL_HPP

#include <ExaMPM_ProblemManager.hpp>

#include <Kokkos_Core.hpp>

#include <cmath>

namespace ExaMPM
{

//---------------------------------------------------------------------------//
// Find minimum timestep across MPI ranks.
double updateGlobalTimestep( MPI_Comm comm, const double current_dt,
                             const double local_min_dt )
{
    double global_min_dt = 0.0;
    MPI_Allreduce( &local_min_dt, &global_min_dt, 1, MPI_DOUBLE, MPI_MIN,
                   comm );

    if ( global_min_dt < current_dt )
        return global_min_dt;
    else
        return current_dt;
}

//---------------------------------------------------------------------------//
// Momentum CFL condition
template <class ExecutionSpace, class ProblemManagerType>
double momentumCFL( MPI_Comm comm, ExecutionSpace, const ProblemManagerType& pm,
                    const double current_dt, const double cfl )
{
    using Kokkos::abs;
    using Kokkos::sqrt;

    // Get the particle data we need.
    auto m_p = pm.get( Location::Particle(), Field::Mass() );
    auto u_p = pm.get( Location::Particle(), Field::Velocity() );
    auto x_p = pm.get( Location::Particle(), Field::Position() );
    auto v_p = pm.get( Location::Particle(), Field::Volume() );
    auto B = pm.bulkModulus();
    auto cell_size =
        pm.mesh()->localGrid()->globalGrid().globalMesh().cellSize( 0 );

    double min_dt;
    std::size_t num_p = pm.numParticle();
    Kokkos::parallel_reduce(
        "momentum_cfl_timestep",
        Kokkos::RangePolicy<ExecutionSpace>( 0, num_p ),
        KOKKOS_LAMBDA( const int p, double& local_min ) {
            // local particle velocity
            auto u = u_p( p, 0 );
            auto v = u_p( p, 1 );
            auto w = u_p( p, 2 );

            // local particle density
            double rho = m_p( p ) / v_p( p );

            // local wave speed
            auto c = sqrt( B / rho );

            // local dt
            double local_dt =
                cfl * cell_size /
                ( abs( u ) + abs( v ) + abs( w ) + c * sqrt( 3.0 ) );

            if ( local_dt < local_min )
                local_min = local_dt;
        },
        Kokkos::Min<double>( min_dt ) );

    return updateGlobalTimestep( comm, current_dt, min_dt );
}

//---------------------------------------------------------------------------//
// Max velocity condition
template <class ExecutionSpace, class ProblemManagerType>
double maxVelocity( MPI_Comm comm, ExecutionSpace, const ProblemManagerType& pm,
                    const double current_dt, const double cfl )
{
    using Kokkos::sqrt;

    // Get the particle data we need.
    auto u_p = pm.get( Location::Particle(), Field::Velocity() );
    auto cell_size =
        pm.mesh()->localGrid()->globalGrid().globalMesh().cellSize( 0 );

    double min_dt;
    std::size_t num_p = pm.numParticle();
    Kokkos::parallel_reduce(
        "max_velocity_timestep",
        Kokkos::RangePolicy<ExecutionSpace>( 0, num_p ),
        KOKKOS_LAMBDA( const int p, double& local_min ) {
            // local particle velocity
            auto u = u_p( p, 0 );
            auto v = u_p( p, 1 );
            auto w = u_p( p, 2 );

            // local dt
            double local_dt =
                cfl * cell_size / ( 2.0 * sqrt( u * u + v * v + w * w ) );

            if ( local_dt < local_min )
                local_min = local_dt;
        },
        Kokkos::Min<double>( min_dt ) );

    return updateGlobalTimestep( comm, current_dt, min_dt );
}

// Set new restricted timestep based on CFL and maximum velocity critera
// Note this will return the user-specified timestep if it is below the
// restictions.
template <class ExecutionSpace, class ProblemManagerType>
double timeStepControl( MPI_Comm comm, const ExecutionSpace& exec_space,
                        const ProblemManagerType& pm, const double current_dt,
                        const double cfl = 0.5 )
{
    double momentum_cfl_dt =
        momentumCFL( comm, exec_space, pm, current_dt, cfl );
    double max_velocity_dt =
        maxVelocity( comm, exec_space, pm, current_dt, cfl );

    return std::min( { momentum_cfl_dt, max_velocity_dt } );
}

} // end namespace ExaMPM

#endif // EXAMPM_TIMESTEPCONTROL_HPP
