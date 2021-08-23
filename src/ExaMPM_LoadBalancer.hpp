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

#include <Cajita_LoadBalancer.hpp>

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
        _lb = std::make_shared<Cajita::LoadBalancer<UniformMesh<double>>>( _mesh->localGrid()->globalGrid(), 3.*_mesh->cellSize() );
    }

    // This will update the domain decomposition and also update the mesh
    void balance( std::shared_ptr<ProblemManager<MemorySpace>> pm )
    {
        double work = pm->numParticle();
        // todo(sschulz): How to save the partitioner, without requiring a shared ptr everywhere?
        auto partitioner = 
        auto global_grid = _lb->createBalancedGlobalGrid( _mesh->globalMesh(), partitioner, work )
        _mesh->newGlobalGrid(global_grid);
        pm->updateMesh( _mesh );
    }

    // Output the actual and internal load balancing grid to VTK files
    void output( const int t ) const
    {
        std::string vtk_actual_domain_basename( "domain_act" );
        std::string vtk_lb_domain_basename( "domain_lb" );
        VTKDomainWriter::writeDomain( _comm, t, _lb->getVertices(),
                                      vtk_actual_domain_basename );
        VTKDomainWriter::writeDomain( _comm, t, _lb->getInternalVertices(),
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
