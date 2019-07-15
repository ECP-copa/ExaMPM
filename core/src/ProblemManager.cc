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
        //  1) Locate particles and evaluate basis functions.
        locateParticles();

        //  2) Calculate the nodal masses.
        calculateNodalMass( node_m );
	
        //  3) Calculate nodal momentum.
        calculateNodalMomentum( node_m, node_p );

        //  4) Calculate internal forces.
        calculateInternalNodalForces( node_f_int );

        //  5) Calculate node impulse.
        calculateNodalImpulse(
            node_f_int, node_m, time_step_size, node_imp );
	
        //  6) Update the particle position and velocity.
        updateParticlePositionAndVelocity( node_v, node_m, time_step_size );

        //  7) Calculate nodal velocity.
        calculateNodalVelocity( node_p, node_imp, node_m, node_v );
	
        //  8) Update the particle gradients.
        updateParticleGradients( node_v, time_step_size );

        //  9) Update the particle stress and strain.
        updateParticleStressStrain();

        // 10) Update the grid node velocities.
        updateGridVelocity( node_imp, node_p, node_m, node_v );

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
// Calculate the mass at the mesh nodes.
void ProblemManager::calculateNodalMass( std::vector<double>& node_m )
{
    int nodes_per_cell = d_mesh->nodesPerCell();
    int node_id = 0;
    int num_nodes = d_mesh->totalNumNodes();
    std::vector<std::vector<double> > node_m_local(
        d_thread_count, std::vector<double>(num_nodes, 0.0) );

    // Reset the nodal mass.
    std::fill( node_m.begin(), node_m.end(), 0.0 );

    // Compute the nodal mass.
    #pragma omp parallel for num_threads(d_thread_count)
    for ( std::size_t i=0; i < d_particles.size(); ++i )
    {
        auto& p = d_particles.at(i);
        int th = omp_get_thread_num();

        // Assemble the nodal mass.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = p.node_ids[n];
            node_m_local[th][ node_id ] += p.basis_values[n] * p.m;
        }
    }

    // Combine thread results.
    #pragma omp parallel for num_threads(d_thread_count)
    for ( int n = 0; n < num_nodes; ++n )
        for ( int th = 0; th < d_thread_count; ++th )
            node_m[ n ] += node_m_local[th][n];

    // Boundary conditions.
    for ( int b = 0; b < 6; ++b )
        d_bc[b]->completeBoundarySum( d_mesh, b, node_m );
}

//---------------------------------------------------------------------------//
// Calculate the nodal momentum.
void ProblemManager::calculateNodalMomentum(
    const std::vector<double>& node_m,
    std::vector<std::array<double,3> >& node_p )
{
    int space_dim = d_mesh->spatialDimension();
    int nodes_per_cell = d_mesh->nodesPerCell();
    int node_id = 0;
    int num_nodes = d_mesh->totalNumNodes();

    std::vector<std::vector<std::array<double,3> > > node_p_local(
        d_thread_count, std::vector<std::array<double,3> > (
        num_nodes, std::array<double,3>() ) );

    // Initialize local storage to zero.
    for ( auto& thread : node_p_local )
        for ( auto& node_mom : thread )
            std::fill( node_mom.begin(), node_mom.end(), 0.0 );
    
    // Reset the momentum.
    for ( auto& mom : node_p )
        std::fill( mom.begin(), mom.end(), 0.0 );
  
    // Update the momentum.
    #pragma omp parallel for num_threads(d_thread_count)
    for ( std::size_t i = 0; i < d_particles.size(); ++i )
    {
        auto& p = d_particles.at(i);
        int th = omp_get_thread_num();

        std::array<double,8> s;
        std::array<double,3> coords;

        // Calculate momentum.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = p.node_ids[n];
            d_mesh->nodeCoordinates( node_id, coords );

            // Define modes.
            s[0] = 1;
            s[1] = coords[0] - p.r[0];
            s[2] = coords[1] - p.r[1];
            s[3] = coords[2] - p.r[2];
            s[4] = s[1] * s[2];
            s[5] = s[1] * s[3];
            s[6] = s[2] * s[3];
            s[7] = s[1] * s[2] * s[3];

            for ( std::size_t r = 0; r < s.size(); ++r )
                for ( int d = 0; d < space_dim; ++d )
                {
                    node_p_local[th][node_id][d] +=
                        p.m * p.basis_values[n] *
                        s[r] * p.c[r][d];
                }
        }
    }

    // Combine thread results.
    #pragma omp parallel for num_threads(d_thread_count)
    for ( int n = 0; n < num_nodes; ++n )
    {
        for ( int th = 0; th < d_thread_count; ++th )
        {
            for ( int d = 0; d < space_dim; ++d )
            {
                node_p[n][d] += node_p_local[th][n][d];
            }
        }
    }

    // Boundary conditions.
    for ( int b = 0; b < 6; ++b )
        d_bc[b]->completeBoundarySum( d_mesh, b, node_p );

    // Boundary conditions.
    for ( int b = 0; b < 6; ++b )
        d_bc[b]->evaluateMomentumCondition( d_mesh, b, node_m, node_p );
}

//---------------------------------------------------------------------------//
// Calculate the nodal velocity.
void ProblemManager::calculateNodalVelocity(
    const std::vector<std::array<double,3> >& node_p,
    const std::vector<std::array<double,3> >& node_imp,
    const std::vector<double>& node_m,
    std::vector<std::array<double,3> >& node_v )
{
    int space_dim = d_mesh->spatialDimension();
    int num_nodes = d_mesh->totalNumNodes();
    int nodes_per_cell = d_mesh->nodesPerCell();
    int node_id = 0;

    // Reset the velocity
    for ( auto& vel : node_v )
        std::fill( vel.begin(), vel.end(), 0.0 );

    // Update the momentum
    #pragma omp parallel for num_threads(d_thread_count)
    for ( std::size_t n = 0; n < num_nodes; ++n )
    {
	    if ( node_m[n] > 0.0 )
            for ( int d = 0; d < space_dim; ++d )
                node_v[n][d] = node_p[n][d] / node_m[n];
	    
	    // Otherwise no mass or momentum so no velocity.
	    else
		    for ( int d = 0; d < space_dim; ++d )
		        node_v[n][d] = 0.0;
    }

    /*
    // Boundary conditions.
    for ( int b = 0; b < 6; ++b )
        d_bc[b]->completeBoundarySum( d_mesh, b, node_v );
    */

    // Boundary conditions.
    for ( int b = 0; b < 6; ++b )
        d_bc[b]->evaluateMomentumCondition( d_mesh, b, node_m, node_v );
}

//---------------------------------------------------------------------------//
// Update the particle deformation gradient and velocity gradient.
void ProblemManager::updateParticleGradients(
    const std::vector<std::array<double,3> >& node_v,
    const double delta_t )
{
    int space_dim = d_mesh->spatialDimension();
    int nodes_per_cell = d_mesh->nodesPerCell();
    int node_id = 0;

    // Update the particles.
    #pragma omp parallel for num_threads(d_thread_count)
    for ( std::size_t index = 0; index < d_particles.size(); ++index )
    {
        auto& p = d_particles.at(index);
        std::array<std::array<double,3>,3> delta_F;
        std::array<std::array<double,3>,3> work;

        // Reset the deformation gradient increment and velocity gradient.
        for ( int d = 0; d < space_dim; ++d )
        {
            std::fill( p.grad_v[d].begin(), p.grad_v[d].end(), 0.0 );
            std::fill( delta_F[d].begin(), delta_F[d].end(), 0.0 );
            std::fill( work[d].begin(), work[d].end(), 0.0 );
        }

        // Compute the velocity gradient at the particle.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = p.node_ids[n];

            for ( int i = 0; i < space_dim; ++i )
            {
                for ( int j = 0; j < space_dim; ++j )
                    p.grad_v[i][j] +=
                        p.basis_gradients[n][i] * node_v[node_id][j];
            }
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
    }
}

//---------------------------------------------------------------------------//
// Update particle strain and stress tensors.
void ProblemManager::updateParticleStressStrain()
{
    // Compute the particle stress.
    #pragma omp parallel for num_threads(d_thread_count)
    for ( std::size_t i = 0; i < d_particles.size(); ++i )
    {
        auto& p = d_particles.at(i);
        d_materials[p.matid]->calculateStress( p );
    }
}

//---------------------------------------------------------------------------//
// Calculate internal forces at mesh nodes.
void ProblemManager::calculateInternalNodalForces(
    std::vector<std::array<double,3> >& node_f_int )
{
    int space_dim = d_mesh->spatialDimension();
    int nodes_per_cell = d_mesh->nodesPerCell();
    int node_id = 0;
    int num_nodes = d_mesh->totalNumNodes();

    std::vector<std::vector<std::array<double,3> > > node_f_local(
        d_thread_count, std::vector<std::array<double,3> > (
        num_nodes, std::array<double,3>() ) );

    // Initialize local storage to zero.
    for ( auto& thread : node_f_local )
        for ( auto& node_f : thread )
            std::fill( node_f.begin(), node_f.end(), 0.0 );

    // Reset the forces.
    for ( auto& f_int : node_f_int )
        std::fill( f_int.begin(), f_int.end(), 0.0 );

    // Compute forces.
    #pragma omp parallel for num_threads(d_thread_count)
    for ( std::size_t i = 0; i < d_particles.size(); ++i )
    {
        auto& p = d_particles.at(i);
        int th = omp_get_thread_num();

        assert( 0 <= p.matid && p.matid < d_materials.size() );

        // Project the particle stress gradients.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = p.node_ids[n];

            // Compute the force via the stress divergence.
            for ( int i = 0; i < space_dim; ++i )
                for ( int j = 0; j < space_dim; ++j )
                    node_f_local[th][node_id][i] -=
                        p.volume * p.basis_gradients[n][j] * p.stress[j][i];
        }
    }

    // Combine thread results.
    #pragma omp parallel for num_threads(d_thread_count)
    for ( int n = 0; n < num_nodes; ++n )
    {
	    for ( int th = 0; th < d_thread_count; ++th )
	    {
	        for ( int d = 0; d < space_dim; ++d )
	        {
		        node_f_int[n][d] += node_f_local[th][n][d];
	        }
	    }
    }

    // Boundary conditions.
    for ( int b = 0; b < 6; ++b )
        d_bc[b]->completeBoundarySum( d_mesh, b, node_f_int );
}

//---------------------------------------------------------------------------//
// Calculate nodal impulse.
void ProblemManager::calculateNodalImpulse(
    const std::vector<std::array<double,3> >& node_f_int,
    const std::vector<double>& node_m,
    const double delta_t,
    std::vector<std::array<double,3> >& node_imp )
{
    int space_dim = d_mesh->spatialDimension();
    int num_nodes = d_mesh->totalNumNodes();

    // Calculate impulse.
    #pragma omp parallel for num_threads(d_thread_count)
    for ( int n = 0; n < num_nodes; ++n )
        for ( int d = 0; d < space_dim; ++d )
            node_imp[n][d] = delta_t * node_f_int[n][d];

    // Add gravity if needed.
    if ( d_has_gravity )
        #pragma omp parallel for num_threads(d_thread_count)
        for ( int n = 0; n < num_nodes; ++n )
            node_imp[n][2] -= delta_t * node_m[n] * 9.81;
    
    // Boundary conditions.
    for ( int b = 0; b < 6; ++b )
        d_bc[b]->evaluateImpulseCondition( d_mesh, b, node_m, node_imp );
}

//---------------------------------------------------------------------------//
// Update grid node velocity.
void ProblemManager::updateGridVelocity(
    const std::vector<std::array<double,3> >& node_imp,
    const std::vector<std::array<double,3> >& node_p,
    const std::vector<double>& node_m ,
    std::vector<std::array<double,3> >& node_v )
{
    int num_nodes = d_mesh->totalNumNodes();
    int space_dim = d_mesh->spatialDimension();

    // Update the grid node velocities.
    #pragma omp parallel for num_threads(d_thread_count)
    for ( int n = 0; n < num_nodes; ++n )
        if ( node_m[n] > 0 )
            for ( int d = 0; d < space_dim; ++d )
                node_v[n][d] = ( node_p[n][d] + node_imp[n][d] ) / node_m[n];
}

//---------------------------------------------------------------------------//
// Update particle position and velocity.
void ProblemManager::updateParticlePositionAndVelocity(
    const std::vector<std::array<double,3> >& node_v,
    const std::vector<double>& node_m ,
    const double delta_t )
{
    int space_dim = d_mesh->spatialDimension();
    int nodes_per_cell = d_mesh->nodesPerCell();
    int node_id = 0;
    double width = d_mesh->cellWidth();

    // Update the particles.
    #pragma omp parallel for num_threads(d_thread_count)
    for ( std::size_t i = 0; i < d_particles.size(); ++i )
    {
        auto& p = d_particles.at(i);
        std::array<double,3> coords;

        for ( auto& mode : p.c )
            std::fill(mode.begin(),mode.end(),0.0);

        // Loop over adjacent nodes.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = p.node_ids[n];
            d_mesh->nodeCoordinates( node_id, coords );

            // Only add a contribution from an adjacent node if it has mass.
            if ( node_m[node_id] > 0.0 )
            {
                for ( int d = 0; d < space_dim; ++d )
                {
                    // Increment the velocity. (PolyPIC Update)
                    p.c[0][d] += p.basis_values[n] * node_v[node_id][d];
                    p.c[1][d] += p.basis_values[n] * ( coords[0] - p.r[0] ) *
                                 node_v[node_id][d] / ( pow(width, 2) / 4. );
                    p.c[2][d] += p.basis_values[n] * ( coords[1] - p.r[1] ) *
                                 node_v[node_id][d] / ( pow(width, 2) / 4. );
                    p.c[3][d] += p.basis_values[n] * ( coords[2] - p.r[2] ) *
                                 node_v[node_id][d] / ( pow(width, 2) / 4. );
                    p.c[4][d] += p.basis_values[n] * ( coords[0] - p.r[0] ) *
                                 ( coords[1] - p.r[1] ) * node_v[node_id][d] /
                                 ( pow(width, 4) / 16. );
                    p.c[5][d] += p.basis_values[n] * ( coords[0] - p.r[0] ) *
                                 ( coords[2] - p.r[2] ) * node_v[node_id][d] /
                                 ( pow(width, 4) / 16. );
                    p.c[6][d] += p.basis_values[n] * ( coords[1] - p.r[1] ) *
                                 ( coords[2] - p.r[2] ) * node_v[node_id][d] /
                                 ( pow(width, 4) / 16. );
                    p.c[7][d] += p.basis_values[n] * ( coords[0] - p.r[0] ) *
                                 ( coords[1] - p.r[1] ) * ( coords[2] - p.r[2] ) *
                                 node_v[node_id][d] / ( pow(width, 6 ) / 64. );
                }
            }
        }
        for ( int d = 0; d < space_dim; ++d )
        {
            // Increment the position.
            p.r[d] +=  delta_t * p.c[0][d];
        }
//        for ( int i = 0; i < 3; ++i )
//            for ( int j = 0; j < 3; ++j )
//                std::cout << p.grad_v[i][j] << "\t" << p.c[i+1][j] <<std::endl;
//        std::cout << std::endl;
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
    file << "x, y, z, velocity magnitude" << std::endl;

    // Write the particle data.
    double vmag = 0.0;
    for ( auto& p : d_particles )
    {
        vmag = std::sqrt( p.c[0][0]*p.c[0][0] + p.c[0][1]*p.c[0][1] + p.c[0][2]*p.c[0][2] );
        file << p.r[0] << ", "
             << p.r[1] << ", "
             << p.r[2] << ", "
             << vmag << std::endl;
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

