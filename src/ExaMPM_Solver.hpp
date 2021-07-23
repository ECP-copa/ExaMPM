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

#ifndef EXAMPM_SOLVER_HPP
#define EXAMPM_SOLVER_HPP

#include <ExaMPM_ProblemManager.hpp>
#include <ExaMPM_TimeIntegrator.hpp>
#include <ExaMPM_Mesh.hpp>
#include <ExaMPM_BoundaryConditions.hpp>
#include <ExaMPM_SiloParticleWriter.hpp>
#include <ExaMPM_VTKDomainWriter.hpp>

#include <Kokkos_Core.hpp>

#include <ALL.hpp>

#include <memory>
#include <string>

#include <mpi.h>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
class SolverBase
{
  public:
    virtual ~SolverBase() = default;
    virtual void solve( const double t_final, const int write_freq ) = 0;
};

//---------------------------------------------------------------------------//
template<class MemorySpace, class ExecutionSpace>
class Solver : public SolverBase
{
  public:

    template<class InitFunc>
    Solver( MPI_Comm comm,
            const Kokkos::Array<double,6>& global_bounding_box,
            const std::array<int,3>& global_num_cell,
            const std::array<bool,3>& periodic,
            const Cajita::ManualBlockPartitioner<3>& partitioner, //todo(sschulz): should work with Cajita::BlockPartitioner<3>!
            const int halo_cell_width,
            const InitFunc& create_functor,
            const int particles_per_cell,
            const double bulk_modulus,
            const double density,
            const double gamma,
            const double kappa,
            const double delta_t,
            const double gravity,
            const BoundaryCondition& bc )
    : _dt( delta_t )
    , _gravity( gravity )
    , _bc( bc )
    , _halo_min( 3 )
    , _comm( comm )
    {
        _mesh = std::make_shared<Mesh<MemorySpace>>(
            global_bounding_box, global_num_cell, periodic, partitioner,
            halo_cell_width, _halo_min, comm );

        _bc.min = _mesh->minDomainGlobalNodeIndex();
        _bc.max = _mesh->maxDomainGlobalNodeIndex();

        _pm = std::make_shared<ProblemManager<MemorySpace>>(
            ExecutionSpace(), _mesh, create_functor, particles_per_cell,
            bulk_modulus, density, gamma, kappa );

        MPI_Comm_rank( comm, &_rank );

        _liball = std::make_shared<ALL::ALL<double, double>>(ALL::TENSOR, 3, 0);
        // For some reason only(!) the following line causes the code to crash
        //auto global_grid = _mesh->localGrid()->globalGrid();
        std::vector<int> block_id(3,0);
        for(std::size_t i=0; i<3; ++i)
            block_id.at(i) = _mesh->mutGlobalGrid().dimBlockId(i);
        std::vector<int> blocks_per_dim(3,0);
        for(std::size_t i=0; i<3; ++i)
            blocks_per_dim.at(i) = _mesh->mutGlobalGrid().dimNumBlock(i);
        _liball->setProcGridParams(block_id, blocks_per_dim);
        std::vector<double> min_domain_size(3,0);
        for(std::size_t i=0; i<3; ++i)
            min_domain_size.at(i) = 3.*_mesh->cellSize();
        _liball->setMinDomainSize(min_domain_size);
        _liball->setCommunicator(MPI_COMM_WORLD);
        _liball->setProcTag(_rank);
        _liball->setup();
    }

    void solve( const double t_final, const int write_freq ) override
    {
        std::string vtk_actual_domain_basename("domain");
        SiloParticleWriter::writeTimeStep(
            _mesh->localGrid()->globalGrid(),
            0,
            0.0,
            _pm->get( Location::Particle(), Field::Position() ),
            _pm->get( Location::Particle(), Field::Velocity() ),
            _pm->get( Location::Particle(), Field::J() ) );

        std::vector<ALL::Point<double>> lb_vertices(2, ALL::Point<double>(3));
        for(std::size_t d=0; d<3; ++d)
            lb_vertices.at(0)[d] = _mesh->localGrid()->globalGrid().globalOffset(d) * _mesh->cellSize();
        for(std::size_t d=0; d<3; ++d)
            lb_vertices.at(1)[d] = lb_vertices.at(0)[d] + _mesh->localGrid()->globalGrid().ownedNumCell(d) * _mesh->cellSize();
        _liball->setVertices(lb_vertices);
        printf(">> %d, %g %g %g | %g %g %g\n", _rank,
                lb_vertices.at(0)[0],lb_vertices.at(0)[1],lb_vertices.at(0)[2],
                lb_vertices.at(1)[0],lb_vertices.at(1)[1],lb_vertices.at(1)[2]);

        int num_step = t_final / _dt;
        double delta_t = t_final / num_step;
        double time = 0.0;
        for ( int t = 0; t < num_step; ++t )
        {
            if ( 0 == _rank && 0 == t % write_freq )
                printf( "Step %d / %d\n", t+1, num_step );

            TimeIntegrator::step( ExecutionSpace(), *_pm, delta_t, _gravity, _bc );

            _liball->setWork(_pm->numParticle());
            _liball->balance();
            std::vector<ALL::Point<double>> updated_vertices = _liball->getVertices();
            for(std::size_t d=0; d<3; ++d)
                _mesh->mutGlobalGrid().setGlobalOffset(d,
                        std::rint( updated_vertices.at(0)[d] / _mesh->cellSize() ));
            for(std::size_t d=0; d<3; ++d)
                _mesh->mutGlobalGrid().setOwnedNumCell(d,
                        std::rint( (updated_vertices.at(1)[d] - updated_vertices.at(0)[d])/_mesh->cellSize() ));
            _liball->setVertices(updated_vertices);
            printf(">> %d, balanced vertices: %g %g %g, %g %g %g\n", _rank,
                updated_vertices.at(0)[0],updated_vertices.at(0)[1],updated_vertices.at(0)[2],
                updated_vertices.at(1)[0],updated_vertices.at(1)[1],updated_vertices.at(1)[2]);

            _pm->updateMesh(_mesh);


            _pm->communicateParticles( _halo_min );

            if ( 0 == t % write_freq )
            {
                SiloParticleWriter::writeTimeStep(
                    _mesh->localGrid()->globalGrid(),
                    t+1,
                    time,
                    _pm->get( Location::Particle(), Field::Position() ),
                    _pm->get( Location::Particle(), Field::Velocity() ),
                    _pm->get( Location::Particle(), Field::J() ) );
                _liball->printVTKoutlines(t);
                std::array<double, 6> vertices;
                for(std::size_t d=0; d<3; ++d)
                    vertices[d] = _mesh->mutGlobalGrid().globalOffset(d) * _mesh->cellSize();
                for(std::size_t d=3; d<6; ++d)
                    vertices[d] = vertices[d-3] + _mesh->mutGlobalGrid().ownedNumCell(d) * _mesh->cellSize();
                printf(">> %d, %g %g %g | %g %g %g\n", _rank,
                        vertices[0], vertices[1], vertices[2],
                        vertices[3], vertices[4], vertices[5]);
                VTKDomainWriter::writeDomain(_comm, t, vertices);
            }

           time += delta_t;
        }
        printf("Finished\n");
    }

  private:

    double _dt;
    double _gravity;
    BoundaryCondition _bc;
    int _halo_min;
    MPI_Comm _comm;
    std::shared_ptr<Mesh<MemorySpace>> _mesh;
    std::shared_ptr<ProblemManager<MemorySpace>> _pm;
    int _rank;
    std::shared_ptr<ALL::ALL<double, double>> _liball;
};

//---------------------------------------------------------------------------//
// Creation method.
template<class InitFunc>
std::shared_ptr<SolverBase>
createSolver( const std::string& device,
              MPI_Comm comm,
              const Kokkos::Array<double,6>& global_bounding_box,
              const std::array<int,3>& global_num_cell,
              const std::array<bool,3>& periodic,
              const Cajita::ManualBlockPartitioner<3>& partitioner,
              const int halo_cell_width,
              const InitFunc& create_functor,
              const int particles_per_cell,
              const double bulk_modulus,
              const double density,
              const double gamma,
              const double kappa,
              const double delta_t,
              const double gravity,
              const BoundaryCondition& bc )
{
    if ( 0 == device.compare("serial") )
    {
#ifdef KOKKOS_ENABLE_SERIAL
        return std::make_shared<ExaMPM::Solver<Kokkos::HostSpace,Kokkos::Serial>>(
            comm,
            global_bounding_box,
            global_num_cell,
            periodic,
            partitioner,
            halo_cell_width,
            create_functor,
            particles_per_cell,
            bulk_modulus,
            density,
            gamma,
            kappa,
            delta_t,
            gravity,
            bc );
#else
        throw std::runtime_error( "Serial Backend Not Enabled" );
#endif
    }
    else if ( 0 == device.compare("openmp") )
    {
#ifdef KOKKOS_ENABLE_OPENMP
        return std::make_shared<ExaMPM::Solver<Kokkos::HostSpace,Kokkos::OpenMP>>(
            comm,
            global_bounding_box,
            global_num_cell,
            periodic,
            partitioner,
            halo_cell_width,
            create_functor,
            particles_per_cell,
            bulk_modulus,
            density,
            gamma,
            kappa,
            delta_t,
            gravity,
            bc );
#else
        throw std::runtime_error( "OpenMP Backend Not Enabled" );
#endif
    }
    else if ( 0 == device.compare("cuda") )
    {
#ifdef KOKKOS_ENABLE_CUDA
        return std::make_shared<ExaMPM::Solver<Kokkos::CudaSpace,Kokkos::Cuda>>(
            comm,
            global_bounding_box,
            global_num_cell,
            periodic,
            partitioner,
            halo_cell_width,
            create_functor,
            particles_per_cell,
            bulk_modulus,
            density,
            gamma,
            kappa,
            delta_t,
            gravity,
            bc );
#else
        throw std::runtime_error( "CUDA Backend Not Enabled" );
#endif
    }
    else if ( 0 == device.compare("hip") )
    {
#ifdef KOKKOS_ENABLE_HIP
        return std::make_shared<ExaMPM::Solver<Kokkos::Experimental::HIPSpace,Kokkos::Experimental::HIP>>(
            comm,
            global_bounding_box,
            global_num_cell,
            periodic,
            partitioner,
            halo_cell_width,
            create_functor,
            particles_per_cell,
            bulk_modulus,
            density,
            gamma,
            kappa,
            delta_t,
            gravity,
            bc );
#else
        throw std::runtime_error( "HIP Backend Not Enabled" );
#endif
    }
    else
    {
        throw std::runtime_error( "invalid backend" );
        return nullptr;
    }
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_SOLVER_HPP
