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

    // Calculate the mass at the mesh nodes.
    void calculateNodalMass( std::vector<double>& node_m );

    // Calculate the nodal momentum.
    void calculateNodalMomentum( const std::vector<double>& node_m,
                                 std::vector<std::array<double,3> >& node_p );

    // Calculate the nodal velocity.
    void calculateNodalVelocity(
        const std::vector<std::array<double,3> >& node_p,
        const std::vector<std::array<double,3> >& node_imp,
        const std::vector<double>& node_m,
        std::vector<std::array<double,3> >& node_v );

    // Update the particle deformation gradient and velocity gradient.
    void updateParticleGradients(
        const std::vector<std::array<double,3> >& node_v,
        const double delta_t );

    // Update particle strain and stress tensors.
    void updateParticleStressStrain();

    // Calculate internal forces at mesh nodes.
    void calculateInternalNodalForces(
        std::vector<std::array<double,3> >& node_f_int	);

    // Calculate nodal impulse.
    void calculateNodalImpulse(
        const std::vector<std::array<double,3> >& node_f_int,
        const std::vector<double>& node_m,
        const double delta_t,
        std::vector<std::array<double,3> >& node_imp );

    // Update particle position and velocity.
    void updateParticlePositionAndVelocity(
        const std::vector<std::array<double,3> >& node_imp,
        const std::vector<std::array<double,3> >& node_p,
        const std::vector<double>& node_m,
        const double delta_t );

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
