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
#include <ExaMPM_TimeIntegrator.hpp>
#include <ExaMPM_TimeStepControl.hpp>

#include <Cabana_Core.hpp>
#include <Kokkos_Core.hpp>

#include <memory>
#include <string>

#include <mpi.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/stat.h>

namespace ExaMPM
{
int nfork;
//---------------------------------------------------------------------------//
class SolverBase
{
  public:
    virtual ~SolverBase() = default;
    virtual void solve( const double t_final, const int write_freq ) = 0;
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
            const Cabana::Grid::BlockPartitioner<3>& partitioner,
            const int halo_cell_width, const InitFunc& create_functor,
            const int particles_per_cell, const double bulk_modulus,
            const double density, const double gamma, const double kappa,
            const double delta_t, const double gravity,
            const BoundaryCondition& bc )
        : _dt( delta_t )
        , _time( 0.0 )
        , _step( 0 )
        , _gravity( gravity )
        , _bc( bc )
        , _halo_min( 3 )
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
    }

    void solve( const double t_final, const int write_freq ) override
    {
        // Output initial state.
        outputParticles();

        while ( _time < t_final )
        {
            if ( 0 == _rank && 0 == _step % write_freq )
                printf( "Time %12.5e / %12.5e\n", _time, t_final );

            // Fixed timestep is guaranteed only when sufficently low dt
            // does not violate the CFL condition (otherwise user-set dt is
            // really a max_dt).
            _dt = timeStepControl( _mesh->localGrid()->globalGrid().comm(),
                                   ExecutionSpace(), *_pm, _dt );

            TimeIntegrator::step( ExecutionSpace(), *_pm, _dt, _gravity, _bc );

            _pm->communicateParticles( _halo_min );

            _time += _dt;
            _step++;

            // Output particles periodically.
            if ( 0 == ( _step ) % write_freq )
                outputParticles();
        }
        Cabana::Experimental::HDF5ParticleOutput::HDF5Config h5_config;
        const char* env_val = std::getenv("H5FUSE");
        if(env_val != NULL) {

          if(h5_config.subfiling) {

            MPI_Comm shmcomm;
            MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                   MPI_INFO_NULL, &shmcomm);

            int shmrank;
            MPI_Comm_rank(shmcomm, &shmrank);
            int status;
            if ( shmrank == 0) {
              //pid_t pid;
              for (int i = 0; i < nfork; i++) {
                waitpid(-1, NULL, 0);
              }
#if 0
              while ((pid = waitpid(-1, &status, WNOHANG)) > 0) {
                if (WIFEXITED(status)) {
                    int ret;

                    if ((ret = WEXITSTATUS(status)) != 0) {
                        printf("h5fuse process exited with error code %d\n", ret);
                        fflush(stdout);
                        MPI_Abort(MPI_COMM_WORLD, -1);
                    }
                }
                else {
                    printf("h5fuse process terminated abnormally\n");
                    fflush(stdout);
                    MPI_Abort(MPI_COMM_WORLD, -1);
                } 
              }
#endif
            }
            MPI_Comm_free(&shmcomm);
          }
       }

    }

    void outputParticles()
    {
        // Prefer HDF5 output over Silo. Only output if one is enabled.
#ifdef Cabana_ENABLE_HDF5
        Cabana::Experimental::HDF5ParticleOutput::HDF5Config h5_config;
        const char* env_val = std::getenv("H5FD_SUBFILING");
        if(env_val != NULL)
          h5_config.subfiling = true;

        // Sets the HDF5 alignment equal subfiling's stripe size
        env_val = std::getenv("H5FD_SUBFILING_STRIPE_SIZE");
        if(env_val != NULL) {
          h5_config.align = true;
          h5_config.threshold = 0;
          h5_config.alignment = std::atoi(env_val);
        } 

        Cabana::Experimental::HDF5ParticleOutput::writeTimeStep(
            h5_config, "particles", _mesh->localGrid()->globalGrid().comm(),
            _step, _time, _pm->numParticle(),
            _pm->get( Location::Particle(), Field::Position() ),
            _pm->get( Location::Particle(), Field::Velocity() ),
            _pm->get( Location::Particle(), Field::J() ) );
#if 1
        env_val = std::getenv("H5FUSE");
        if(env_val != NULL) {
      
        if(h5_config.subfiling) {

          MPI_Comm shmcomm;
          MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                 MPI_INFO_NULL, &shmcomm);

          int shmrank;
          MPI_Comm_rank(shmcomm, &shmrank);

          if ( shmrank == 0) {
             pid_t pid = 0;
             //pid_t tmppid;
             int   status;
             //std::cout << "DOING OUTPUT " << _rank << std::endl;
             //signal(SIGHUP, SIG_IGN);
             //not interested in its childs, prevent zombies 
             //signal(SIGCHLD, SIG_IGN);
             //signal(SIGQUIT, SIG_IGN);

             pid = fork();
             nfork++;
             if (pid == 0) {
              std::stringstream filename_hdf5;
              filename_hdf5 << "particles" << "_" << _step << ".h5";

              // Directory containing the subfiling configuration file
              std::stringstream config_dir;
              if(const char* env_value = std::getenv(H5FD_SUBFILING_CONFIG_FILE_PREFIX))
                config_dir << env_value;
              else
                config_dir << ".";
              // Find the name of the subfiling configuration file
              struct stat file_info;
              stat(filename_hdf5.str().c_str(), &file_info);

              char config_filename[PATH_MAX];
              snprintf(config_filename, PATH_MAX, "%s/" H5FD_SUBFILING_CONFIG_FILENAME_TEMPLATE, config_dir.str().c_str(),
                   filename_hdf5.str().c_str(), (uint64_t)file_info.st_ino);

              // Call the h5fuse utility
              char* args[] = { strdup("./h5fuse.sh"), strdup("-r"), strdup("-f"), config_filename, NULL};
              //char* args[] = { strdup("./h5fuse.sh"), strdup("-r"), strdup("-v"), strdup("-f"), config_filename, NULL};

              execvp(args[0], args);
              //exit(0);

#if 0
             //  char *tmp_filename;
               char *args[7];

            //tmp_filename = (char *)malloc(PATH_MAX);

            /* Generate name for configuration file */
            //HDsnprintf(tmp_filename, PATH_MAX, "%s/" H5FD_SUBFILING_CONFIG_FILENAME_TEMPLATE, config_dir,
            //           SUBF_FILENAME, file_inode);
            args[0] = strdup("ls");
            args[1] = strdup("-l");
            args[2] = NULL;
            execvp("ls", args);
#endif
        }
#if 0
        else {
            //tmppid = waitpid(pid, &status, WNOHANG);
            waitpid(pid, &status, 0);
            if (WIFEXITED(status)) {
                int ret;

                if ((ret = WEXITSTATUS(status)) != 0) {
                    printf("h5fuse process exited with error code %d\n", ret);
                    fflush(stdout);
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
            }
            else {
                printf("h5fuse process terminated abnormally\n");
                fflush(stdout);
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
        }
#endif
       }
}
}
#endif  
#if 0
            pid_t tmppid;
            int   status;

            pid_t pid = fork();
            if (pid == -1) {
                perror("fork");
                exit(EXIT_FAILURE);
            }
#endif
#if 0
            if (pid == 0) {
              // stub HDF5 filename
              std::stringstream filename_hdf5;
              filename_hdf5 << "particles" << "_" << _step << ".h5";

              // Directory containing the subfiling configuration file
              std::stringstream config_dir;
              if(const char* env_value = std::getenv(H5FD_SUBFILING_CONFIG_FILE_PREFIX))
                config_dir << env_value;
              else
                config_dir << ".";
              // Find the name of the subfiling configuration file
              struct stat file_info;
              stat(filename_hdf5.str().c_str(), &file_info);

              char config_filename[PATH_MAX];
              snprintf(config_filename, PATH_MAX, "%s/" H5FD_SUBFILING_CONFIG_FILENAME_TEMPLATE, config_dir.str().c_str(),
                   filename_hdf5.str().c_str(), (uint64_t)file_info.st_ino);

              // Call the h5fuse utility
              char* args[] = { strdup("./h5fuse.sh"), strdup("-r"), strdup("-v"), strdup("-f"), config_filename, NULL};

              execvp(args[0], args);
            }
            else {
              tmppid = waitpid(pid, &status, 0);

              if (WIFEXITED(status)) {
                int ret;

                if ((ret = WEXITSTATUS(status)) != 0) {
                  printf("h5fuse process exited with error code %d\n", ret);
                  fflush(stdout);
                  MPI_Abort( _mesh->localGrid()->globalGrid().comm(), -1);
                }
              }
              else {
                printf("h5fuse process terminated abnormally\n");
                fflush(stdout);
                MPI_Abort( _mesh->localGrid()->globalGrid().comm(), -1);
              }
            }
#endif
//          }
//        }
//#endif
#else
#ifdef Cabana_ENABLE_SILO
        Cabana::Grid::Experimental::SiloParticleOutput::writeTimeStep(
            "particles", _mesh->localGrid()->globalGrid(), _step, _time,
            _pm->get( Location::Particle(), Field::Position() ),
            _pm->get( Location::Particle(), Field::Velocity() ),
            _pm->get( Location::Particle(), Field::J() ) );
#else
        if ( _rank == 0 )
            std::cout << "No particle output enabled in Cabana. Add "
                         "Cabana_REQUIRE_HDF5=ON or Cabana_REQUIRE_SILO=ON to "
                         "the Cabana build if needed.";
#endif
#endif
    }

  private:
    double _dt;
    double _time;
    int _step;
    double _gravity;
    BoundaryCondition _bc;
    int _halo_min;
    std::shared_ptr<Mesh<MemorySpace>> _mesh;
    std::shared_ptr<ProblemManager<MemorySpace>> _pm;
    int _rank;
};

//---------------------------------------------------------------------------//
// Creation method.
template <class InitFunc>
std::shared_ptr<SolverBase>
createSolver( const std::string& device, MPI_Comm comm,
              const Kokkos::Array<double, 6>& global_bounding_box,
              const std::array<int, 3>& global_num_cell,
              const std::array<bool, 3>& periodic,
              const Cabana::Grid::BlockPartitioner<3>& partitioner,
              const int halo_cell_width, const InitFunc& create_functor,
              const int particles_per_cell, const double bulk_modulus,
              const double density, const double gamma, const double kappa,
              const double delta_t, const double gravity,
              const BoundaryCondition& bc )
{
    if ( 0 == device.compare( "serial" ) || 0 == device.compare( "Serial" ) ||
         0 == device.compare( "SERIAL" ) )
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
    else if ( 0 == device.compare( "openmp" ) ||
              0 == device.compare( "OpenMP" ) ||
              0 == device.compare( "OPENMP" ) )
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
    else if ( 0 == device.compare( "cuda" ) || 0 == device.compare( "Cuda" ) ||
              0 == device.compare( "CUDA" ) )
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
    else if ( 0 == device.compare( "hip" ) || 0 == device.compare( "Hip" ) ||
              0 == device.compare( "HIP" ) )
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
