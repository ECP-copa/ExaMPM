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

#ifndef EXAMPM_PARTICLECOMMUNICATION_HPP
#define EXAMPM_PARTICLECOMMUNICATION_HPP

#include <Cajita.hpp>

#include <ExaMPM_Types.hpp>
#include <ExaMPM_Mesh.hpp>

#include <Cabana_Core.hpp>

#include <Kokkos_Core.hpp>

#include <mpi.h>

#include <vector>
#include <numeric>
#include <type_traits>
#include <algorithm>

namespace ExaMPM
{
namespace ParticleCommunication
{
//---------------------------------------------------------------------------//
// Check for the global number of particles that must be communicated.
//---------------------------------------------------------------------------//
template<class LocalGridType,
         class CoordSliceType>
int communicationCount( const LocalGridType& local_grid,
                        const CoordSliceType& coords,
                        const int minimum_halo_width )
{
    using execution_space = typename CoordSliceType::execution_space;

    // Locate the particles in the local mesh and count how many have left the
    // halo region.
    auto local_mesh = Cajita::createLocalMesh<Kokkos::HostSpace>( local_grid );
    auto dx = local_grid.globalGrid().globalMesh().cellSize(0);
    const Kokkos::Array<double,3> local_low =
        { local_mesh.lowCorner(Cajita::Ghost(),Dim::I) + minimum_halo_width*dx,
          local_mesh.lowCorner(Cajita::Ghost(),Dim::J) + minimum_halo_width*dx,
          local_mesh.lowCorner(Cajita::Ghost(),Dim::K) + minimum_halo_width*dx };
    const Kokkos::Array<double,3> local_high =
        { local_mesh.highCorner(Cajita::Ghost(),Dim::I) - minimum_halo_width*dx,
          local_mesh.highCorner(Cajita::Ghost(),Dim::J) - minimum_halo_width*dx,
          local_mesh.highCorner(Cajita::Ghost(),Dim::K) - minimum_halo_width*dx };
    int comm_count = 0;
    Kokkos::parallel_reduce(
        "redistribute_count",
        Kokkos::RangePolicy<execution_space>(0,coords.size()),
        KOKKOS_LAMBDA( const int p, int& result ){
            if ( coords(p,Dim::I) < local_low[Dim::I] ||
                 coords(p,Dim::I) > local_high[Dim::I] ||
                 coords(p,Dim::J) < local_low[Dim::J] ||
                 coords(p,Dim::J) > local_high[Dim::J] ||
                 coords(p,Dim::K) < local_low[Dim::K] ||
                 coords(p,Dim::K) > local_high[Dim::K] )
                result += 1;
        },
        comm_count );

    MPI_Allreduce( MPI_IN_PLACE, &comm_count, 1, MPI_INT, MPI_SUM,
                   local_grid.globalGrid().comm() );

    return comm_count;
}

//---------------------------------------------------------------------------//
// Compute particle destinations and shift periodic coordinates.
//---------------------------------------------------------------------------//
template<class LocalGridType,
         class CoordSliceType,
         class NeighborRankView,
         class DestinationRankView>
void prepareCommunication(
    const LocalGridType& local_grid,
    const NeighborRankView& neighbor_ranks,
    DestinationRankView& destinations,
    CoordSliceType& coords )
{
    using execution_space = typename CoordSliceType::execution_space;

    // Locate the particles in the global grid and get their destination
    // rank. The particle halo should be constructed such that particles will
    // only move to a location in the 26 neighbor halo or stay on this
    // rank. If the particle crosses a periodic boundary update it's
    // coordinates to represent the shift.
    auto local_mesh = Cajita::createLocalMesh<Kokkos::HostSpace>( local_grid );
    const auto& global_grid = local_grid.globalGrid();
    const auto& global_mesh = global_grid.globalMesh();
    const Kokkos::Array<double,3> local_low =
        { local_mesh.lowCorner(Cajita::Own(),Dim::I),
          local_mesh.lowCorner(Cajita::Own(),Dim::J),
          local_mesh.lowCorner(Cajita::Own(),Dim::K) };
    const Kokkos::Array<double,3> local_high =
        { local_mesh.highCorner(Cajita::Own(),Dim::I),
          local_mesh.highCorner(Cajita::Own(),Dim::J),
          local_mesh.highCorner(Cajita::Own(),Dim::K) };
    const Kokkos::Array<bool,3> period =
        { global_grid.isPeriodic(Dim::I),
          global_grid.isPeriodic(Dim::J),
          global_grid.isPeriodic(Dim::K) };
    const Kokkos::Array<double,3> global_low =
        { global_mesh.lowCorner(Dim::I),
          global_mesh.lowCorner(Dim::J),
          global_mesh.lowCorner(Dim::K) };
    const Kokkos::Array<double,3> global_high =
        { global_mesh.highCorner(Dim::I),
          global_mesh.highCorner(Dim::J),
          global_mesh.highCorner(Dim::K) };
    const Kokkos::Array<double,3> global_span =
        { global_mesh.extent(Dim::I),
          global_mesh.extent(Dim::J),
          global_mesh.extent(Dim::K) };
    Kokkos::parallel_for(
        "redistribute_locate_shift",
        Kokkos::RangePolicy<execution_space>(0,coords.size()),
        KOKKOS_LAMBDA( const int p ){
            // Compute the logical index of the neighbor we are sending to.
            int nid[3] = {1,1,1};
            for ( int d = 0; d < 3; ++d )
            {
                if ( coords(p,d) < local_low[d] ) nid[d] = 0;
                else if ( coords(p,d) > local_high[d] ) nid[d] = 2;
            }

            // Compute the destination MPI rank.
            destinations( p ) =
                neighbor_ranks( nid[Dim::I] + 3*(nid[Dim::J] + 3*nid[Dim::K]) );

            // Shift periodic coordinates if needed.
            for ( int d = 0; d < 3; ++d )
            {
                if ( period[d] )
                {
                    if ( coords(p,d) > global_high[d] )
                        coords(p,d) -= global_span[d];
                    else if ( coords(p,d) < global_low[d] )
                        coords(p,d) += global_span[d];
                }
            }
        });
}

//---------------------------------------------------------------------------//
// Particle redistribution
//---------------------------------------------------------------------------//
/*!
  \brief Redistribute particles to new owning ranks based on their location.

  \param local_grid The local_grid in which the particles are currently located.

  \param particles The particles to redistribute.

  \param Member index in the AoSoA of the particle coordinates.

  \param force_communication If true communication will always occur even if
  particles have not exited the halo.
 */
template<class LocalGridType, class ParticleContainer, std::size_t CoordIndex>
void redistribute( const LocalGridType& local_grid,
                   const int minimum_halo_width,
                   ParticleContainer& particles,
                   std::integral_constant<std::size_t,CoordIndex>,
                   const bool force_communication = false )
{
    using device_type = typename ParticleContainer::device_type;

    // Get the coordinates.
    auto coords = Cabana::slice<CoordIndex>( particles );

    // If we are not forcing communication check to see if we need to
    // communicate.
    if ( !force_communication )
    {
        // Check to see if we need to communicate.
        auto comm_count =
            communicationCount( local_grid, coords, minimum_halo_width );

        // If we have no particle communication to do then exit.
        if ( 0 == comm_count )
            return;
    }

    // Of the 27 potential local grids figure out which are in our topology. Some
    // of the ranks in this list may be invalid. We will update this list
    // after we compute destination ranks so it is unique and valid.
    std::vector<int> topology( 27, -1 );
    int nr = 0;
    for ( int k = -1; k < 2; ++k )
        for ( int j = -1; j < 2; ++j )
            for ( int i = -1; i < 2; ++i, ++nr )
                topology[nr] = local_grid.neighborRank(i,j,k);

    // Locate the particles in the global grid and get their destination
    // rank and shift periodic coordinates if necessary.
    Kokkos::View<int*,Kokkos::HostSpace,Kokkos::MemoryUnmanaged>
        neighbor_ranks( topology.data(), topology.size() );
    auto nr_mirror = Kokkos::create_mirror_view_and_copy(
        device_type(), neighbor_ranks );
    Kokkos::View<int*,device_type> destinations(
        Kokkos::ViewAllocateWithoutInitializing("destinations"),
        particles.size() );
    prepareCommunication( local_grid, nr_mirror, destinations, coords );

    // Make the topology a list of unique and valid ranks.
    auto remove_end = std::remove( topology.begin(), topology.end(), -1 );
    std::sort( topology.begin(), remove_end );
    auto unique_end = std::unique( topology.begin(), remove_end );
    topology.resize( std::distance(topology.begin(),unique_end) );

    // Create the Cabana distributor.
    Cabana::Distributor<device_type> distributor( local_grid.globalGrid().comm(),
                                                  destinations,
                                                  topology );

    // Redistribute the particles.
    Cabana::migrate( distributor, particles );
}

//---------------------------------------------------------------------------//

} // end namespace ParticleCommunication
} // end namespace ExaMPM

#endif // end EXAMPM_PARTICLECOMMUNICATION_HPP
