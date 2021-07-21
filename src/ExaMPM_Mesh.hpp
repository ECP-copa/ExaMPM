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

#include <Cajita.hpp>

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
template<class MemorySpace>
class Mesh
{
  public:

    using memory_space = MemorySpace;

    // Construct a mesh.
    Mesh( const Kokkos::Array<double,6>& global_bounding_box,
          const std::array<int,3>& global_num_cell,
          const std::array<bool,3>& periodic,
          const Cajita::ManualBlockPartitioner<3>& partitioner, //todo(sschulz): should work with Cajita::BlockPartitioner<3>!
          const int halo_cell_width,
          const int minimum_halo_cell_width,
          MPI_Comm comm )
    : _global_bounding_box( global_bounding_box )
    , _global_num_cell( global_num_cell )
    , _periodic( periodic )
    , _partitioner( partitioner )
    , _halo_cell_width( halo_cell_width )
    , _minimum_halo_cell_width( minimum_halo_cell_width )
    , _comm( comm )
    {
        // Make a copy of the global number of cells so we can modify it.
        std::array<int,3> num_cell = global_num_cell;

        // Compute the cell size.
        double cell_size =
            (global_bounding_box[3] - global_bounding_box[0]) /
            num_cell[0];

        // Because the mesh is uniform check that the domain is evenly
        // divisible by the cell size in each dimension within round-off
        // error. This will let us do cheaper math for particle location.
        for ( int d = 0; d < 3; ++d )
        {
            double extent = num_cell[d] * cell_size;
            if ( std::abs(
                     extent - (global_bounding_box[d+3]-
                               global_bounding_box[d]) ) >
                 double( 10.0 ) * std::numeric_limits<double>::epsilon() )
                throw std::logic_error(
                    "Extent not evenly divisible by uniform cell size" );
        }

        // Create global mesh bounds.
        std::array<double,3> global_low_corner = { global_bounding_box[0],
                                                   global_bounding_box[1],
                                                   global_bounding_box[2] };
        std::array<double,3> global_high_corner = { global_bounding_box[3],
                                                    global_bounding_box[4],
                                                    global_bounding_box[5] };
        for ( int d = 0; d < 3; ++d )
        {
            _min_domain_global_node_index[d] = 0;
            _max_domain_global_node_index[d] = num_cell[d] + 1;
        }

        // For dimensions that are not periodic we pad by the minimum halo
        // cell width to allow for projections outside of the domain.
        for ( int d = 0; d < 3; ++d )
        {
            if ( !periodic[d] )
            {
                global_low_corner[d] -= cell_size*minimum_halo_cell_width;
                global_high_corner[d] += cell_size*minimum_halo_cell_width;
                num_cell[d] += 2*minimum_halo_cell_width;
                _min_domain_global_node_index[d] += minimum_halo_cell_width;
                _max_domain_global_node_index[d] -= minimum_halo_cell_width;
            }
        }

        // Create the global mesh.
        auto global_mesh = Cajita::createUniformGlobalMesh(
            global_low_corner, global_high_corner, num_cell );

        // Build the global grid.
        _global_grid = Cajita::createGlobalGrid(
            comm, global_mesh, periodic, partitioner );

        // Build the local grid.
        int halo_width = std::max( minimum_halo_cell_width, halo_cell_width );
        _local_grid = Cajita::createLocalGrid( _global_grid, halo_width );
    }

    void updateGeometry(
            std::array<int,3>& local_offset,
            std::array<int,3>& local_num_cell )
    {
        // Create global mesh bounds.
        std::array<double,3> global_low_corner = { _global_bounding_box[0],
                                                   _global_bounding_box[1],
                                                   _global_bounding_box[2] };
        std::array<double,3> global_high_corner = { _global_bounding_box[3],
                                                    _global_bounding_box[4],
                                                    _global_bounding_box[5] };
        // For dimensions that are not periodic we pad by the minimum halo
        // cell width to allow for projections outside of the domain.
        double cell_size = _local_grid->globalGrid().globalMesh().cellSize(0);
        for ( int d = 0; d < 3; ++d )
        {
            if ( !_periodic[d] )
            {
                global_low_corner[d] -= cell_size*_minimum_halo_cell_width;
                global_high_corner[d] += cell_size*_minimum_halo_cell_width;
            }
        }

        // Create the global mesh.
        auto global_mesh = Cajita::createUniformGlobalMesh(
            global_low_corner, global_high_corner, _global_num_cell );

        // Build the global grid.
        _global_grid = Cajita::createGlobalGrid(
            _comm, global_mesh, _periodic, _partitioner, local_offset, local_num_cell );

        // Build the local grid.
        int halo_width = std::max( _minimum_halo_cell_width, _halo_cell_width );
        _local_grid = Cajita::createLocalGrid( _global_grid, halo_width );
        // todo(sschulz): Can we just reassing without causing a memory leak?
    }

    
    // Get the local grid.
    const std::shared_ptr<Cajita::LocalGrid<Cajita::UniformMesh<double>>>& localGrid()
    {
        return _local_grid;
    }

    // Get the global grid.
    const std::shared_ptr<Cajita::GlobalGrid<Cajita::UniformMesh<double>>>& globalGrid()
    {
        return _global_grid;
    }

    // Get the cell size.
    double cellSize() const
    {
        return _local_grid->globalGrid().globalMesh().cellSize(0);
    }

    // Get the minimum node index in the domain.
    Kokkos::Array<int,3> minDomainGlobalNodeIndex() const
    {
        return _min_domain_global_node_index;
    }

    // Get the maximum node index in the domain.
    Kokkos::Array<int,3> maxDomainGlobalNodeIndex() const
    {
        return _max_domain_global_node_index;
    }

  public:

    Kokkos::Array<double,6> _global_bounding_box;
    std::array<int,3> _global_num_cell;
    std::array<bool,3> _periodic;
    Cajita::ManualBlockPartitioner<3> _partitioner;
    int _halo_cell_width;
    int _minimum_halo_cell_width;
    MPI_Comm _comm;

    std::shared_ptr<
      Cajita::LocalGrid<Cajita::UniformMesh<double>>> _local_grid;
    std::shared_ptr<
      Cajita::GlobalGrid<Cajita::UniformMesh<double>>> _global_grid;

    Kokkos::Array<int,3> _min_domain_global_node_index;
    Kokkos::Array<int,3> _max_domain_global_node_index;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_MESH_HPP
