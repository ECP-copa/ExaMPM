//---------------------------------------------------------------------------//
/*!
 * \file ProblemManager.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_PROBLEMMANAGER_HH
#define EXAMPM_PROBLEMMANAGER_HH

#include "Mesh.hh"
#include "Geometry.hh"
#include "Particle.hh"
#include "BoundaryCondition.hh"
#include "StressModel.hh"

#include <memory>
#include <vector>
#include <array>
#include <functional>
#include <string>
#include <chrono>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*!
 * \class ProblemManger
 * \brief Simulation manager.
 */
class ProblemManager
{
  public:

    // Constructor.
    ProblemManager( const int mesh_num_cells_x,
                    const int mesh_num_cells_y,
                    const int mesh_num_cells_z,
                    const double mesh_cell_width,
                    const bool has_gravity,
                    const int thread_count );

    // Set boundary conditions.
    void setBoundaryConditions(
        const std::array<std::shared_ptr<ExaMPM::BoundaryCondition>,6>& bc );

    // Set material models.
    void setMaterialModels(
        const std::vector<std::shared_ptr<StressModel> >& materials );

    // Initialize the problem with a given order over a set of geometries.
    void initialize( const std::vector<std::shared_ptr<Geometry> >& geometry,
                     const int order );

    // Solve the problem for a given number of time steps.
    void solve( const int num_time_steps,
                const double time_step_size,
                const std::string& output_file,
                const int write_frequency );

  private:

    // Locate the particles and compute grid values.
    void locateParticles();
//        std::vector<std::array<double,8> >& cell_basis,
//        std::vector<std::array<std::array<double,3>,8> >& cell_grads);

    // Update the grid state.
    void updateGrid( const double delta_t,
                     std::vector<double>& node_m,
                     std::vector<std::array<double,3> >& node_v,
                     std::vector<std::array<double,3> >& node_a );

    // Update the particles state.
    void updateParticles( const double delta_t,
                          const std::vector<std::array<double,3> >& node_v,
                          const std::vector<std::array<double,3> >& node_a );

    // Write a time step to file.
    void writeTimeStepToFile(
        const std::string& output_file, const int step	) const;

    // Write runtime results to file.
    void displayRuntime(
        const std::vector<double> step_times ) const;

  private:

    // Gravity boolean.
    bool d_has_gravity;

    // Mesh
    std::shared_ptr<Mesh> d_mesh;

    // Boundary conditions.
    std::array<std::shared_ptr<ExaMPM::BoundaryCondition>,6> d_bc;

    // Material models.
    std::vector<std::shared_ptr<StressModel> > d_materials;

    // Particles.
    std::vector<Particle> d_particles;

    // OpenMP Threads.
    int d_thread_count;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_PROBLEMMANAGER_HH

//---------------------------------------------------------------------------//
// end ProblemManager.hh
//---------------------------------------------------------------------------//
