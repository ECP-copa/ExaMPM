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

#include <ExaMPM_BoundaryConditions.hpp>
#include <ExaMPM_Mesh.hpp>
#include <ExaMPM_ProblemManager.hpp>
#include <ExaMPM_SiloParticleWriter.hpp>
#include <ExaMPM_TimeIntegrator.hpp>
#include <ExaMPM_VTKDomainWriter.hpp>

#ifdef ExaMPM_ENABLE_LB
#include <Cajita_LoadBalancer.hpp>
#endif

#include <Kokkos_Core.hpp>

#include <memory>

#include <mpi.h>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
class SolverBase
{
  public:
    virtual ~SolverBase() = default;
    virtual void solve( const double t_final, const int write_freq ) = 0;
    virtual void solve( const double t_final, const int write_freq,
                        const int lb_freq ) = 0;
};

//---------------------------------------------------------------------------//
template <class MemorySpace, class ExecutionSpace>
class Solver : public SolverBase
{
  public:
    template <class InitFunc>
    Solver( MPI_Comm comm, const Kokkos::Array<double, 6>& global_bounding_box,
            const std::array<int, 3>& global_num_cell,
            const std::array<bool, 3>& periodic,
            const std::shared_ptr<Cajita::ManualPartitioner>& partitioner,
            const int halo_cell_width, const InitFunc& create_functor,
            const int particles_per_cell, const double bulk_modulus,
            const double density, const double gamma, const double kappa,
            const double delta_t, const double gravity,
            const BoundaryCondition& bc )
        : _dt( delta_t )
        , _gravity( gravity )
        , _bc( bc )
        , _halo_min( 3 )
        , _comm( comm )
        , _partitioner( partitioner )
    {
        _mesh = std::make_shared<Mesh<MemorySpace>>(
            global_bounding_box, global_num_cell, periodic, *partitioner,
            halo_cell_width, _halo_min, comm );

        _bc.min = _mesh->minDomainGlobalNodeIndex();
        _bc.max = _mesh->maxDomainGlobalNodeIndex();

        _pm = std::make_shared<ProblemManager<MemorySpace>>(
            ExecutionSpace(), _mesh, create_functor, particles_per_cell,
            bulk_modulus, density, gamma, kappa );

        MPI_Comm_rank( comm, &_rank );

#ifdef ExaMPM_ENABLE_LB
        _lb = std::make_shared<
            Cajita::Experimental::LoadBalancer<Cajita::UniformMesh<double>>>(
            _comm, _mesh->globalGrid(), 3. * _mesh->cellSize() );
#endif
    }

    void solve( const double t_final, const int write_freq ) override
    {
        solve( t_final, write_freq, 10 );
    }

    void solve( const double t_final, const int write_freq,
                const int lb_freq ) override
    {
        double output_time = 0;
        double integrate_time = 0;
        double lb_time = 0;
        double comm_part_time = 0;
        Kokkos::Timer timer, output_timer, integrate_timer, lb_timer,
            comm_part_timer, step_timer;

        output_timer.reset();
        SiloParticleWriter::writeTimeStep(
            _mesh->localGrid()->globalGrid(), 0, 0.0,
            _pm->get( Location::Particle(), Field::Position() ),
            _pm->get( Location::Particle(), Field::Velocity() ),
            _pm->get( Location::Particle(), Field::J() ) );

        std::string vtk_actual_domain_basename( "domain_act" );
        std::string vtk_lb_domain_basename( "domain_lb" );
        std::array<double, 6> vertices;

#ifdef ExaMPM_ENABLE_LB
        vertices = _lb->getVertices();
        VTKDomainWriter::writeDomain( _comm, 0, vertices,
                                      vtk_actual_domain_basename );
        vertices = _lb->getInternalVertices();
        VTKDomainWriter::writeDomain( _comm, 0, vertices,
                                      vtk_lb_domain_basename );
#endif
        output_time += output_timer.seconds();

        int total_num_particle;
        int num_particle = _pm->numParticle();
        MPI_Reduce( &num_particle, &total_num_particle, 1, MPI_INT, MPI_SUM, 0,
                    _comm );
        int comm_size;
        MPI_Comm_size( _comm, &comm_size );
        step_timer.reset();

        int num_step = t_final / _dt;
        double delta_t = t_final / num_step;
        double time = 0.0;
        for ( int t = 0; t < num_step; ++t )
        {
            if ( 0 == _rank && 0 == t % write_freq )
                printf( "Step %d / %d\n", t + 1, num_step );

            integrate_timer.reset();
            TimeIntegrator::step( ExecutionSpace(), *_pm, delta_t, _gravity,
                                  _bc );
            integrate_time += integrate_timer.seconds();

            lb_timer.reset();
#ifdef ExaMPM_ENABLE_LB
            if ( lb_freq > 0 && 0 == t % lb_freq )
            {
                double work = _pm->numParticle();
                auto global_grid = _lb->createBalancedGlobalGrid(
                    _mesh->globalMesh(), *_partitioner, work );
                _mesh->newGlobalGrid( global_grid );
                _pm->updateMesh( _mesh );
            }
#endif
            lb_time += lb_timer.seconds();

            comm_part_timer.reset();
            _pm->communicateParticles( _halo_min );
            comm_part_time += comm_part_timer.seconds();

            if ( 0 == t % write_freq )
            {
                output_timer.reset();
                SiloParticleWriter::writeTimeStep(
                    _mesh->localGrid()->globalGrid(), t + 1, time,
                    _pm->get( Location::Particle(), Field::Position() ),
                    _pm->get( Location::Particle(), Field::Velocity() ),
                    _pm->get( Location::Particle(), Field::J() ) );

#ifdef ExaMPM_ENABLE_LB
                vertices = _lb->getVertices();
                VTKDomainWriter::writeDomain( _comm, t, vertices,
                                              vtk_actual_domain_basename );
                vertices = _lb->getInternalVertices();
                VTKDomainWriter::writeDomain( _comm, t, vertices,
                                              vtk_lb_domain_basename );
                double imbalance = _lb->getImbalance();
                output_time += output_timer.seconds();
                output_timer.reset();
                if ( _rank == 0 )
                {
                    double step_time = step_timer.seconds();
                    std::cout
                        << std::fixed << std::setprecision( 5 ) << t << " "
                        << imbalance << " " << write_freq / step_time << " "
                        << write_freq / step_time * total_num_particle << " "
                        << write_freq / step_time * total_num_particle /
                               comm_size
                        << std::endl;
                }
                step_timer.reset();
                output_time += output_timer.seconds();
#endif
            }

            time += delta_t;
        }
        double total_time = timer.seconds();
        double steps_per_sec = 1. * num_step / total_time;
        if ( _rank == 0 )
        {
            std::cout
                << std::fixed << std::setprecision( 2 )
                << "\n#Procs Particles | Time T_Out T_lb T_Int T_Comm_part\n"
                << comm_size << " " << total_num_particle << " | " << total_time
                << " " << output_time << " " << lb_time << " " << integrate_time
                << " " << comm_part_time << " | PERFORMANCE\n"
                << comm_size << " " << total_num_particle << " | "
                << "1.0"
                << " " << output_time / total_time << " "
                << lb_time / total_time << " " << integrate_time / total_time
                << " " << comm_part_time / total_time << " | FRACTION\n\n"
                << "#Steps/s Particlesteps/s Particlesteps/(proc*s)\n"
                << std::scientific << steps_per_sec << " "
                << steps_per_sec * total_num_particle << " "
                << steps_per_sec * total_num_particle / comm_size << std::endl;
        }
    }

  private:
    double _dt;
    double _gravity;
    BoundaryCondition _bc;
    int _halo_min;
    std::shared_ptr<Mesh<MemorySpace>> _mesh;
    std::shared_ptr<ProblemManager<MemorySpace>> _pm;
    int _rank;
    MPI_Comm _comm;
    std::shared_ptr<Cajita::ManualPartitioner> _partitioner;
#ifdef ExaMPM_ENABLE_LB
    std::shared_ptr<
        Cajita::Experimental::LoadBalancer<Cajita::UniformMesh<double>>>
        _lb;
#endif
};

//---------------------------------------------------------------------------//
// Creation method.
template <class InitFunc>
std::shared_ptr<SolverBase>
createSolver( const std::string& device, MPI_Comm comm,
              const Kokkos::Array<double, 6>& global_bounding_box,
              const std::array<int, 3>& global_num_cell,
              const std::array<bool, 3>& periodic,
              const std::shared_ptr<Cajita::ManualPartitioner>& partitioner,
              const int halo_cell_width, const InitFunc& create_functor,
              const int particles_per_cell, const double bulk_modulus,
              const double density, const double gamma, const double kappa,
              const double delta_t, const double gravity,
              const BoundaryCondition& bc )
{
    if ( 0 == device.compare( "serial" ) )
    {
#ifdef KOKKOS_ENABLE_SERIAL
        return std::make_shared<
            ExaMPM::Solver<Kokkos::HostSpace, Kokkos::Serial>>(
            comm, global_bounding_box, global_num_cell, periodic, partitioner,
            halo_cell_width, create_functor, particles_per_cell, bulk_modulus,
            density, gamma, kappa, delta_t, gravity, bc );
#else
        throw std::runtime_error( "Serial Backend Not Enabled" );
#endif
    }
    else if ( 0 == device.compare( "openmp" ) )
    {
#ifdef KOKKOS_ENABLE_OPENMP
        return std::make_shared<
            ExaMPM::Solver<Kokkos::HostSpace, Kokkos::OpenMP>>(
            comm, global_bounding_box, global_num_cell, periodic, partitioner,
            halo_cell_width, create_functor, particles_per_cell, bulk_modulus,
            density, gamma, kappa, delta_t, gravity, bc );
#else
        throw std::runtime_error( "OpenMP Backend Not Enabled" );
#endif
    }
    else if ( 0 == device.compare( "cuda" ) )
    {
#ifdef KOKKOS_ENABLE_CUDA
        return std::make_shared<
            ExaMPM::Solver<Kokkos::CudaSpace, Kokkos::Cuda>>(
            comm, global_bounding_box, global_num_cell, periodic, partitioner,
            halo_cell_width, create_functor, particles_per_cell, bulk_modulus,
            density, gamma, kappa, delta_t, gravity, bc );
#else
        throw std::runtime_error( "CUDA Backend Not Enabled" );
#endif
    }
    else if ( 0 == device.compare( "hip" ) )
    {
#ifdef KOKKOS_ENABLE_HIP
        return std::make_shared<ExaMPM::Solver<Kokkos::Experimental::HIPSpace,
                                               Kokkos::Experimental::HIP>>(
            comm, global_bounding_box, global_num_cell, periodic, partitioner,
            halo_cell_width, create_functor, particles_per_cell, bulk_modulus,
            density, gamma, kappa, delta_t, gravity, bc );
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
