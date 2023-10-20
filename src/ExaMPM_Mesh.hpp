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

#ifndef EXAMPM_MESH_HPP
#define EXAMPM_MESH_HPP

#include <Cabana_Grid.hpp>

#include <Kokkos_Core.hpp>

#include <memory>

#include <mpi.h>

#include <limits>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*!
  \class Mesh
  \brief Logically and spatially uniform Cartesian mesh.
*/
template <class MemorySpace>
class Mesh
{
  public:
    using memory_space = MemorySpace;

    // Construct a mesh.
    Mesh( const Kokkos::Array<double, 6>& global_bounding_box,
          const std::array<int, 3>& global_num_cell,
          const std::array<bool, 3>& periodic,
          const Cabana::Grid::BlockPartitioner<3>& partitioner,
          const int halo_cell_width, const int minimum_halo_cell_width,
          MPI_Comm comm )
    {
        // Make a copy of the global number of cells so we can modify it.
        std::array<int, 3> num_cell = global_num_cell;

        // Compute the cell size.
        double cell_size =
            ( global_bounding_box[3] - global_bounding_box[0] ) / num_cell[0];

        // Because the mesh is uniform check that the domain is evenly
        // divisible by the cell size in each dimension within round-off
        // error. This will let us do cheaper math for particle location.
        for ( int d = 0; d < 3; ++d )
        {
            double extent = num_cell[d] * cell_size;
            if ( std::abs( extent - ( global_bounding_box[d + 3] -
                                      global_bounding_box[d] ) ) >
                 double( 10.0 ) * std::numeric_limits<double>::epsilon() )
                throw std::logic_error(
                    "Extent not evenly divisible by uniform cell size" );
        }

        // Create global mesh bounds.
        std::array<double, 3> global_low_corner = { global_bounding_box[0],
                                                    global_bounding_box[1],
                                                    global_bounding_box[2] };
        std::array<double, 3> global_high_corner = { global_bounding_box[3],
                                                     global_bounding_box[4],
                                                     global_bounding_box[5] };
        for ( int d = 0; d < 3; ++d )
        {
            _min_domain_global_node_index[d] = 0;
            _max_domain_global_node_index[d] = num_cell[d];
        }

        // For dimensions that are not periodic we pad by the minimum halo
        // cell width to allow for projections outside of the domain.
        for ( int d = 0; d < 3; ++d )
        {
            if ( !periodic[d] )
            {
                global_low_corner[d] -= cell_size * minimum_halo_cell_width;
                global_high_corner[d] += cell_size * minimum_halo_cell_width;
                num_cell[d] += 2 * minimum_halo_cell_width;
                _min_domain_global_node_index[d] += minimum_halo_cell_width;
                _max_domain_global_node_index[d] += minimum_halo_cell_width;
            }
        }

        // Create the global mesh.
        auto global_mesh = Cabana::Grid::createUniformGlobalMesh(
            global_low_corner, global_high_corner, num_cell );

        // Build the global grid.
        auto global_grid = Cabana::Grid::createGlobalGrid(
            comm, global_mesh, periodic, partitioner );

        // Build the local grid.
        int halo_width = std::max( minimum_halo_cell_width, halo_cell_width );
        _local_grid = Cabana::Grid::createLocalGrid( global_grid, halo_width );
    }

    // Get the local grid.
    const std::shared_ptr<
        Cabana::Grid::LocalGrid<Cabana::Grid::UniformMesh<double>>>&
    localGrid() const
    {
        return _local_grid;
    }

    // Get the cell size.
    double cellSize() const
    {
        return _local_grid->globalGrid().globalMesh().cellSize( 0 );
    }

    // Get the minimum node index in the domain.
    Kokkos::Array<int, 3> minDomainGlobalNodeIndex() const
    {
        return _min_domain_global_node_index;
    }

    // Get the maximum node index in the domain.
    Kokkos::Array<int, 3> maxDomainGlobalNodeIndex() const
    {
        return _max_domain_global_node_index;
    }

  public:
    std::shared_ptr<Cabana::Grid::LocalGrid<Cabana::Grid::UniformMesh<double>>>
        _local_grid;

    Kokkos::Array<int, 3> _min_domain_global_node_index;
    Kokkos::Array<int, 3> _max_domain_global_node_index;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_MESH_HPP
