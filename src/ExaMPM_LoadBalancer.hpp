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

#ifndef EXAMPM_LOADBALANCER_HPP
#define EXAMPM_LOADBALANCER_HPP

#include <ExaMPM_Mesh.hpp>
#include <ExaMPM_ProblemManager.hpp>
#include <ExaMPM_VTKDomainWriter.hpp>

#include <ALL.hpp>

#include <string>

#include <mpi.h>

namespace ExaMPM
{
template <class MemorySpace>
class LoadBalancer
{
  public:
    using mesh_type = Mesh<MemorySpace>;
    // todo(sschulz): Allow for arbitrary dimension
    LoadBalancer( MPI_Comm comm, std::shared_ptr<mesh_type>& mesh )
        : _comm( comm )
        , _mesh( mesh )
    {
        MPI_Comm_rank( comm, &_rank );
        _liball =
            std::make_shared<ALL::ALL<double, double>>( ALL::TENSOR, 3, 0 );
        // todo(sschulz): Check if the code still crashes.
        // For some reason only(!) the following line causes the code to crash
        // auto global_grid = _mesh->localGrid()->globalGrid();
        std::vector<int> block_id( 3, 0 );
        for ( std::size_t i = 0; i < 3; ++i )
            block_id.at( i ) = _mesh->globalGrid().dimBlockId( i );
        std::vector<int> blocks_per_dim( 3, 0 );
        for ( std::size_t i = 0; i < 3; ++i )
            blocks_per_dim.at( i ) = _mesh->globalGrid().dimNumBlock( i );
        _liball->setProcGridParams( block_id, blocks_per_dim );
        std::vector<double> min_domain_size( 3, 0 );
        for ( std::size_t i = 0; i < 3; ++i )
            min_domain_size.at( i ) = 3. * _mesh->cellSize();
        _liball->setMinDomainSize( min_domain_size );
        _liball->setCommunicator( _comm );
        _liball->setProcTag( _rank );
        _liball->setup();
        std::vector<ALL::Point<double>> lb_vertices( 2,
                                                     ALL::Point<double>( 3 ) );
        for ( std::size_t d = 0; d < 3; ++d )
            lb_vertices.at( 0 )[d] =
                _mesh->localGrid()->globalGrid().globalOffset( d ) *
                _mesh->cellSize();
        for ( std::size_t d = 0; d < 3; ++d )
            lb_vertices.at( 1 )[d] =
                lb_vertices.at( 0 )[d] +
                _mesh->localGrid()->globalGrid().ownedNumCell( d ) *
                    _mesh->cellSize();
        _liball->setVertices( lb_vertices );
    }

    // This will update the domain decomposition and also update the mesh
    void balance( std::shared_ptr<ProblemManager<MemorySpace>> pm )
    {
        _liball->setWork( pm->numParticle() );
        _liball->balance();
        std::vector<ALL::Point<double>> updated_vertices =
            _liball->getVertices();
        std::array<int, 3> cell_index_lo, cell_index_hi;
        for ( std::size_t d = 0; d < 3; ++d )
            cell_index_lo[d] =
                std::rint( updated_vertices.at( 0 )[d] / _mesh->cellSize() );
        for ( std::size_t d = 0; d < 3; ++d )
            cell_index_hi[d] =
                std::rint( updated_vertices.at( 1 )[d] / _mesh->cellSize() );
        std::array<int, 3> num_cell;
        for ( std::size_t d = 0; d < 3; ++d )
            num_cell[d] = cell_index_hi[d] - cell_index_lo[d];
        _mesh->globalGrid().setNumCellAndOffset( num_cell, cell_index_lo );
        _liball->setVertices( updated_vertices );
        pm->updateMesh( _mesh );
    }

    // Output the actual and internal load balancing grid to VTK files
    void output( const int t ) const
    {
        std::string vtk_actual_domain_basename( "domain_act" );
        std::string vtk_lb_domain_basename( "domain_lb" );
        std::vector<ALL::Point<double>> updated_vertices =
            _liball->getVertices();
        // todo(sschulz): The official VTK routine seems to create mangled files
        // on my system.
        // _liball->printVTKoutlines( t );
        std::array<double, 6> vertices;
        for ( std::size_t d = 0; d < 3; ++d )
            vertices[d] =
                static_cast<double>( _mesh->globalGrid().globalOffset( d ) ) *
                _mesh->cellSize();
        for ( std::size_t d = 3; d < 6; ++d )
            vertices[d] = vertices[d - 3] +
                          static_cast<double>(
                              _mesh->globalGrid().ownedNumCell( d - 3 ) ) *
                              _mesh->cellSize();
        VTKDomainWriter::writeDomain( _comm, t, vertices,
                                      vtk_actual_domain_basename );
        for ( std::size_t d = 0; d < 3; ++d )
            vertices[d] = updated_vertices.at( 0 )[d];
        for ( std::size_t d = 3; d < 6; ++d )
            vertices[d] = updated_vertices.at( 1 )[d - 3];
        VTKDomainWriter::writeDomain( _comm, t, vertices,
                                      vtk_lb_domain_basename );
    }

  private:
    MPI_Comm _comm;
    std::shared_ptr<mesh_type> _mesh;
    std::shared_ptr<ALL::ALL<double, double>> _liball;
    int _rank;
};
} // end namespace ExaMPM
#endif // EXAMPM_LOADBALANCER_HPP
