//---------------------------------------------------------------------------//
/*!
 * \file ProblemManager.cc
 */
//---------------------------------------------------------------------------//

#include "ProblemManager.hh"
#include "TensorTools.hh"

#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>

#include <omp.h>
#include <chrono>
namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Constructor.
ProblemManager::ProblemManager( const int mesh_num_cells_x,
                                const int mesh_num_cells_y,
                                const int mesh_num_cells_z,
                                const double mesh_cell_width,
                                const bool has_gravity,
                                const int thread_count )
    : d_has_gravity( has_gravity )
    , d_thread_count( thread_count )
{
    // Create the mesh.
    d_mesh = std::make_shared<Mesh>( mesh_num_cells_x,
                                     mesh_num_cells_y,
                                     mesh_num_cells_z,
                                     mesh_cell_width,
                                     d_thread_count );
}

//---------------------------------------------------------------------------//
// Set boundary conditions.
void ProblemManager::setBoundaryConditions(
    const std::array<std::shared_ptr<ExaMPM::BoundaryCondition>,6>& bc )
{
    d_bc = bc;
}

//---------------------------------------------------------------------------//
// Set material models.
void ProblemManager::setMaterialModels(
    const std::vector<std::shared_ptr<StressModel> >& materials )
{
    d_materials = materials;
}

//---------------------------------------------------------------------------//
// Initialize the problem with a given order over a set of geometries.
void ProblemManager::initialize(
    const std::vector<std::shared_ptr<Geometry> >& geometry,
    const int order )
{
    // Create a vector of candidate particles.
    int ppcell = d_mesh->particlesPerCell( order );
    std::vector<ExaMPM::Particle> candidates( ppcell );

    // Loop through each cell and check its candidates against the geometries.
    int num_cells = d_mesh->totalNumCells();
    for ( int c = 0; c < num_cells; ++c )
    {
        // Create the candidates.
        d_mesh->initializeParticles( c, order, candidates );

        // Check each candidate
        for ( int p = 0; p < ppcell; ++p )
        {
            // Check each geometry
            for ( const auto& g : geometry )
            {
                // If the candidate is in the geometry add it.
                if ( g->particleInGeometry(candidates[p]) )
                {
                    g->initializeParticle(candidates[p]);
                    d_particles.push_back( candidates[p] );
                    break;
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Solve the problem for a given number of time steps.
void ProblemManager::solve( const int num_time_steps,
                            const double time_step_size,
                            const std::string& output_file,
                            const int write_frequency )
{
    // Allocate mesh fields.
    int num_nodes = d_mesh->totalNumNodes();

    // Nodal mass
    std::vector<double> node_m( num_nodes, 0.0 );

    // Nodal velocity
    std::vector<std::array<double,3> > node_v(
        num_nodes, std::array<double,3>() );

    // Nodal Acceleration
    std::vector<std::array<double,3> > node_a(
        num_nodes, std::array<double,3>() );

    // Nodal momentum
    std::vector<std::array<double,3> > node_p(
        num_nodes, std::array<double,3>() );

    // Nodal impulse
    std::vector<std::array<double,3> > node_imp(
        num_nodes, std::array<double,3>() );

    // Nodal internal force.
    std::vector<std::array<double,3> > node_f_int(
        num_nodes, std::array<double,3>() );

    // Time step parameters
    double time = 0.0;

    std::vector<double> step_times( num_time_steps, 0.0 );

    // Write the initial state.
    int write_step = 0;
    writeTimeStepToFile( output_file, write_step );

    // Time step
    for ( int step = 0; step < num_time_steps; ++step )
    {
        // Increment time.
        time += time_step_size;

        // Print time step info.
        if ( (step+1) % write_frequency == 0 )
        {
            std::cout << "Time Step " << step+1 << "/" << num_time_steps << ": "
                      <<  time <<  " (s)" << std::endl;
        }

        auto start = std::chrono::system_clock::now();

        // 1) Locate particles and evaluate basis functions.
        locateParticles();

        // 2) Update the grid state.
        updateGrid( time_step_size, node_m, node_v, node_a );

        // 3) Update the particles.
        updateParticles( time_step_size, node_v, node_a );

        auto stop = std::chrono::system_clock::now();
        std::chrono::duration<double> runtime = stop-start;
        // Write the time step to file.
        if ( (step+1) % write_frequency == 0 )
        {
            ++write_step;
            writeTimeStepToFile( output_file, write_step );
        }
        step_times[step] = runtime.count();
    }
    displayRuntime( step_times );

    // Write the end state.
    writeTimeStepToFile( output_file, write_step+1 );
}

//---------------------------------------------------------------------------//
// Locate the particles and compute grid values.
void ProblemManager::locateParticles()
{

    // Locate the particles
#pragma omp parallel for num_threads(d_thread_count)
    for ( std::size_t i = 0; i < d_particles.size(); ++i )
    {
        auto& p = d_particles.at(i);

        std::array<int,3> cell_id;
        std::array<double,3> ref_coords;

        // Locate the particle.
        d_mesh->locateParticle( p, cell_id );

        // Get the node ids local to the particle.
        d_mesh->cellNodeIds( cell_id, p.node_ids );

        // Map the particle to the reference frame of the cell.
        d_mesh->mapPhysicalToReferenceFrame( p, cell_id, ref_coords );

        // Evaluate the cell basis function at the particle location.
        d_mesh->shapeFunctionValue( ref_coords, p.basis_values );

        // Evaluate the cell basis function gradient at the particle location.
        d_mesh->shapeFunctionGradient( ref_coords, p.basis_gradients );
    }
}

//---------------------------------------------------------------------------//
// Update the grid state.
void ProblemManager::updateGrid( const double delta_t,
                                 std::vector<double>& node_m,
                                 std::vector<std::array<double,3> >& node_v,
                                 std::vector<std::array<double,3> >& node_a )
{
    int space_dim = d_mesh->spatialDimension();
    int nodes_per_cell = d_mesh->nodesPerCell();
    int num_nodes = d_mesh->totalNumNodes();
    std::vector<std::vector<double> > node_m_local(
        d_thread_count, std::vector<double>(num_nodes, 0.0) );
    std::vector<std::vector<std::array<double,3> > > node_v_local(
        d_thread_count, std::vector<std::array<double,3> > (
            num_nodes, std::array<double,3>() ) );
    std::vector<std::vector<std::array<double,3> > > node_a_local(
        d_thread_count, std::vector<std::array<double,3> > (
            num_nodes, std::array<double,3>() ) );

    // Reset grid data.
    std::fill( node_m.begin(), node_m.end(), 0.0 );
    for ( auto& vel : node_v )
        std::fill( vel.begin(), vel.end(), 0.0 );
    for ( auto& acc : node_a )
        std::fill( acc.begin(), acc.end(), 0.0 );

    // Initialize local storage to zero.
    for ( auto& thread : node_v_local )
        for ( auto& node_vel : thread )
            std::fill( node_vel.begin(), node_vel.end(), 0.0 );
    for ( auto& thread : node_a_local )
        for ( auto& node_acc : thread )
            std::fill( node_acc.begin(), node_acc.end(), 0.0 );

    // Project mass and momentum to grid.
#pragma omp parallel for num_threads(d_thread_count)
    for ( std::size_t i=0; i < d_particles.size(); ++i )
    {
        // Get the particle
        const auto& p = d_particles.at(i);

        // Get the thread id.
        int th = omp_get_thread_num();

    int node_id = 0;

        // Project to each node.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            // Get the node id.
            node_id = p.node_ids[n];

            // Deposit particle mass
            node_m_local[th][node_id] += p.basis_values[n] * p.m;

            // Deposit particle momentum.
            for ( int d = 0; d < space_dim; ++d )
                node_v_local[th][node_id][d] +=
                    p.m * p.v[d] * p.basis_values[n];

            // Add impulse from particle stress to acceleration.
            for ( int i = 0; i < space_dim; ++i )
                for ( int j = 0; j < space_dim; ++j )
                    node_a_local[th][node_id][i] -=
                        p.volume * p.basis_gradients[n][j] * p.stress[j][i];

            // Add gravity if needed to acceleration. Gravity is in the Z direction.
            if ( d_has_gravity )
                node_a_local[th][node_id][2] -= p.m * p.basis_values[n] * 9.81;
        }
    }

    // Combine thread results.
#pragma omp parallel for num_threads(d_thread_count)
    for ( int n = 0; n < num_nodes; ++n )
    {
        for ( int th = 0; th < d_thread_count; ++th )
        {
            node_m[n] += node_m_local[th][n];

            for ( int d = 0; d < space_dim; ++d )
                node_v[n][d] += node_v_local[th][n][d];

            for ( int d = 0; d < space_dim; ++d )
                node_a[n][d] += node_a_local[th][n][d];
        }
    }

    // Complete the boundary sum.
    for ( int b = 0; b < 6; ++b )
    {
        d_bc[b]->completeBoundarySum( d_mesh, b, node_m );
        d_bc[b]->completeBoundarySum( d_mesh, b, node_v );
        d_bc[b]->completeBoundarySum( d_mesh, b, node_a );
    }

    // Apply boundary conditions.
    for ( int b = 0; b < 6; ++b )
    {
        d_bc[b]->evaluateMomentumCondition( d_mesh, b, node_m, node_v );
        d_bc[b]->evaluateMomentumCondition( d_mesh, b, node_m, node_a );
    }

    // Compute grid velocity and acceleration from momentum.
#pragma omp parallel for num_threads(d_thread_count)
    for ( int n = 0; n < num_nodes; ++n )
    {
        if ( node_m[n] > 0.0 )
        {
            for ( int d = 0; d < space_dim; ++d )
                node_v[n][d] = (node_v[n][d] + delta_t * node_a[n][d]) / node_m[n];
            for ( int d = 0; d < space_dim; ++d )
                node_a[n][d] /= node_m[n];
        }
        else
        {
            for ( int d = 0; d < space_dim; ++d )
                node_v[n][d] = 0.0;
            for ( int d = 0; d < space_dim; ++d )
                node_a[n][d] = 0.0;
        }
    }
}

//---------------------------------------------------------------------------//
// Update the particles state.
void ProblemManager::updateParticles( const double delta_t,
                                      const std::vector<std::array<double,3> >& node_v,
                                      const std::vector<std::array<double,3> >& node_a )
{
    int space_dim = d_mesh->spatialDimension();
    int nodes_per_cell = d_mesh->nodesPerCell();

    // Update the particles.
#pragma omp parallel for num_threads(d_thread_count)
    for ( std::size_t i = 0; i < d_particles.size(); ++i )
    {
        // Get the particle.
        auto& p = d_particles.at(i);
        int node_id = 0;

        // Gradient work vectors.
        std::array<std::array<double,3>,3> delta_F;
        std::array<std::array<double,3>,3> work;

        // Reset the deformation gradient increment and velocity gradient.
        for ( int d = 0; d < space_dim; ++d )
        {
            std::fill( p.grad_v[d].begin(), p.grad_v[d].end(), 0.0 );
            std::fill( delta_F[d].begin(), delta_F[d].end(), 0.0 );
            std::fill( work[d].begin(), work[d].end(), 0.0 );
        }

        // Gather from each node.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            // Get the node id.
            node_id = p.node_ids[n];

            // Increment the position.
            for ( int d = 0; d < space_dim; ++d )
                p.r[d] += delta_t * node_v[node_id][d] * p.basis_values[n];

            // Increment the velocity. (FLIP Update)
            for ( int d = 0; d < space_dim; ++d )
                p.v[d] += delta_t * node_a[node_id][d] * p.basis_values[n];

            // Update velocity gradient.
            for ( int i = 0; i < space_dim; ++i )
                for ( int j = 0; j < space_dim; ++j )
                    p.grad_v[i][j] +=
                        p.basis_gradients[n][i] * node_v[node_id][j];
        }

        // Scale the velocity gradient.
        for ( int i = 0; i < space_dim; ++i )
            for ( int j = 0; j < space_dim; ++j )
                work[i][j] = p.grad_v[i][j] * delta_t;

        // Compute the deformation gradient increment.
        for ( int i = 0; i < space_dim; ++i )
            for ( int j = 0; j < space_dim; ++j )
                for ( int k = 0; k < space_dim; ++k )
                    delta_F[i][j] += work[i][k] * p.F[k][j];

        // Update the deformation gradient.
        for ( int i = 0; i < space_dim; ++i )
            for ( int j = 0; j < space_dim; ++j )
                p.F[i][j] += delta_F[i][j];

        // Add the scaled velocity gradient to the identity.
        for ( int i = 0; i < space_dim; ++i )
            work[i][i] += 1.0;

        // Update the particle volume.
        p.volume *= TensorTools::determinant( work );

        // Compute the particle stress.
        d_materials[p.matid]->calculateStress( p );
    }

    // Boundary conditions.
    for ( int b = 0; b < 6; ++b )
        d_bc[b]->evaluateBoundaryPosition( d_mesh, b, d_particles );
}

//---------------------------------------------------------------------------//
// Write a time step to file.
void ProblemManager::writeTimeStepToFile(
    const std::string& output_file, const int step ) const
{

    // Open the time step file.
    std::string filename = output_file + ".csv." + std::to_string(step);
    std::ofstream file( filename );

    // Write the data header.
    file << "x, y, z, velocity magnitude, J" << std::endl;

    // Write the particle data.
    double vmag = 0.0;
    for ( auto& p : d_particles )
    {
        vmag = std::sqrt( p.v[0]*p.v[0] + p.v[1]*p.v[1] + p.v[2]*p.v[2] );
        file << p.r[0] << ", "
             << p.r[1] << ", "
             << p.r[2] << ", "
             << vmag << ", "
             << TensorTools::determinant( p.F ) << std::endl;
    }

    // Close the time step file
    file.close();
}

//---------------------------------------------------------------------------//
// Write runtimes to file.
void ProblemManager::displayRuntime(
    const std::vector<double> step_times ) const
{

    double step_time = 0.0;

    // Compute average time for one time step.
#pragma omp parallel for reduction(+: step_time)
    for ( std::size_t i = 0; i < step_times.size(); ++i )
    {
        step_time += step_times[i];
    }
    step_time = step_time / step_times.size();

    // Print out runtime details.
    std::cout << "Average time per step: " << step_time
              << std::endl;
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end ProblemManager.cc
//---------------------------------------------------------------------------//
