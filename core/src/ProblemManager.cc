//---------------------------------------------------------------------------//
/*!
 * \file ProblemManager.cc
 */
//---------------------------------------------------------------------------//

#include "ProblemManager.hh"
#include "FileIO.hh"
#include "TensorTools.hh"

#include <iostream>
#include <cassert>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Constructor.
ProblemManager::ProblemManager( const int mesh_num_cells_x,
                                const int mesh_num_cells_y,
                                const int mesh_num_cells_z,
                                const double mesh_cell_width )
{
    // Create the mesh.
    d_mesh = std::make_shared<Mesh>( mesh_num_cells_x,
                                     mesh_num_cells_y,
                                     mesh_num_cells_z,
                                     mesh_cell_width );

    // Create a zero specific body force.
    d_body_force = []( const std::array<double,3>& r, std::array<double,3>& f )
                   { std::fill( f.begin(), f.end(), 0.0 ); };
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
        const std::vector<MaterialModel>& materials )
{
    d_materials = materials;
}

//---------------------------------------------------------------------------//
// Set the specific body force.
void ProblemManager::setSpecificBodyForce(
    std::function<void(const std::array<double,3>& r,
                       std::array<double,3>& f)> force_field )
{
    d_body_force = force_field;
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

    // Nodal external force
    std::vector<std::array<double,3> > node_f_ext(
        num_nodes, std::array<double,3>() );

    // Setup an output writer.
    ExaMPM::FileIO file_io( output_file );

    // Time step parameters
    double time = 0.0;

    // Write the initial state.
    int write_step = 0;
    file_io.writeTimeStep( write_step, time, d_particles );

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

        // 1) Locate particles and evaluate basis functions.
        locateParticles();

        // 2) Calculate the nodal masses.
        calculateNodalMass( node_m );

        // 3) Calculate nodal momentum.
        calculateNodalMomentum( node_m, node_p );

        // 4) Calculate nodal velocity.
        calculateNodalVelocity( node_p, node_m, node_v );

        // 5) Update the particle gradients.
        updateParticleGradients( node_v, time_step_size );

        // 6) Update the particle stress and strain.
        updateParticleStressStrain();

        // 7) Calculate internal forces.
        calculateInternalNodalForces( node_f_int );

        // 8) Calculate external forces.
        calculateExternalNodalForces( node_f_ext );

        // 9) Calculate node impulse.
        calculateNodalImpulse(
            node_f_int, node_f_ext, node_m, time_step_size, node_imp );

        // 10) Update the particle position and velocity.
        updateParticlePositionAndVelocity(
            node_imp, node_p, node_m, time_step_size );

        // Write the time step to file.
        if ( (step+1) % write_frequency == 0 )
        {
            ++write_step;
            file_io.writeTimeStep( write_step, time, d_particles );
        }
    }

    // Write the end state.
    file_io.writeTimeStep( write_step+1, time, d_particles );
}

//---------------------------------------------------------------------------//
// Locate the particles and compute grid values.
void ProblemManager::locateParticles()
{
    std::array<int,3> cell_id;
    std::array<double,3> ref_coords;

    for ( auto& p : d_particles )
    {
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

    // Reset the nodal mass.
    std::fill( node_m.begin(), node_m.end(), 0.0 );

    // Compute the nodal mass.
    for ( auto& p : d_particles )
    {
        // Assemble the nodal mass.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = p.node_ids[n];
            node_m[ node_id ] += p.basis_values[n] * p.m;
        }
    }
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

    // Reset the momentum
    for ( auto& mom : node_p )
        std::fill( mom.begin(), mom.end(), 0.0 );

    // Update the momentum
    for ( auto& p : d_particles )
    {
        // Calculate momentum.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = p.node_ids[n];

            for ( int d = 0; d < space_dim; ++d )
            {
                node_p[node_id][d] +=
                    p.m * p.v[d] * p.basis_values[n];
            }
        }
    }

    // Boundary conditions.
    for ( int b = 0; b < 6; ++b )
        d_bc[b]->evaluateMomentumCondition( d_mesh, b, node_m, node_p );
}

//---------------------------------------------------------------------------//
// Calculate the nodal velocity.
void ProblemManager::calculateNodalVelocity(
    const std::vector<std::array<double,3> >& node_p,
    const std::vector<double>& node_m,
    std::vector<std::array<double,3> >& node_v )
{
    int space_dim = d_mesh->spatialDimension();
    int num_nodes = d_mesh->totalNumNodes();

    // Update the velocity
    for ( int n = 0; n < num_nodes; ++n )
    {
        // If we have nodal mass do the update.
        if ( node_m[n] > 0.0 )
        {
            for ( int d = 0; d < space_dim; ++d )
            {
                node_v[n][d] = node_p[n][d] / node_m[n];
            }
        }

        // Otherwise we have no mass or momentum so no velocity.
        else
        {
            for ( int d = 0; d < space_dim; ++d )
            {
                node_v[n][d] = 0.0;
            }
        }
    }
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
    std::array<std::array<double,3>,3> delta_F;
    std::array<std::array<double,3>,3> work;

    // Update the particles
    for ( auto& p : d_particles )
    {
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
    // Compute the particle strain.
    for ( auto& p : d_particles )
    {
        assert( 0 <= p.matid && p.matid < d_materials.size() );
        d_materials[p.matid].strain_model->calculateStrain( p );
    }

    // // Compute cell strain smoothing parameter.
    // int num_cells = d_mesh->totalNumCells();
    // std::vector<double> alpha( num_cells, 0.0 );
    // std::vector<double> cell_mass( num_cells, 0.0 );
    // double strain_trace;
    // for ( auto& p : d_particles )
    // {
    //     strain_trace = 0.0;
    //     for ( int i = 0; i < 3; ++i )
    //         strain_trace += p.strain[i][i];
    //     alpha[ p.cell_id ] += p.m * strain_trace;
    //     cell_mass[ p.cell_id ] += p.m;
    // }
    // for ( int c = 0; c < num_cells; ++c )
    // {
    //     if ( cell_mass[c] > 0.0 )
    //         alpha[c] /= cell_mass[c];
    // }

    // // Compute the smoothed particle strain.
    // for ( auto& p : d_particles )
    // {
    //     for ( int i = 0; i < 3; ++i )
    //         p.strain[i][i] =
    //             (2.0 * p.strain[i][i] + alpha[p.cell_id]) / 3.0;
    // }

    // Compute the particle stress.
    for ( auto& p : d_particles )
    {
        d_materials[p.matid].stress_model->calculateStress( p );
    }

    // // Compute cell stress smoothing parameter.
    // std::vector<double> beta( num_cells, 0.0 );
    // double stress_trace;
    // for ( auto& p : d_particles )
    // {
    //     stress_trace = 0.0;
    //     for ( int i = 0; i < 3; ++i )
    //         stress_trace += p.stress[i][i];
    //     beta[ p.cell_id ] += p.m * stress_trace;
    // }
    // for ( int c = 0; c < num_cells; ++c )
    // {
    //     if ( cell_mass[c] > 0.0 )
    //         beta[c] /= cell_mass[c];
    // }

    // // Compute the smoothed particle stress.
    // for ( auto& p : d_particles )
    // {
    //     for ( int i = 0; i < 3; ++i )
    //         p.stress[i][i] =
    //             (2.0 * p.stress[i][i] + beta[p.cell_id]) / 3.0;
    // }
}

//---------------------------------------------------------------------------//
// Calculate external forces at mesh nodes.
void ProblemManager::calculateExternalNodalForces(
    std::vector<std::array<double,3> >& node_f_ext )
{
    int space_dim = d_mesh->spatialDimension();
    int nodes_per_cell = d_mesh->nodesPerCell();
    std::array<double,3> local_acceleration;
    int node_id = 0;

    // Reset the forces.
    for ( auto& f_ext : node_f_ext )
        std::fill( f_ext.begin(), f_ext.end(), 0.0 );

    // Calculate forces.
    for ( auto& p : d_particles )
    {
        // Calculate force.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            d_body_force( p.r, local_acceleration );

            node_id = p.node_ids[n];
            for ( int d = 0; d < space_dim; ++d )
            {
                node_f_ext[node_id][d] +=
                    p.m * local_acceleration[d] * p.basis_values[n];
            }
        }
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

    // Reset the forces.
    for ( auto& f_int : node_f_int )
        std::fill( f_int.begin(), f_int.end(), 0.0 );

    // Compute forces.
    for ( auto& p : d_particles )
    {
        assert( 0 <= p.matid && p.matid < d_materials.size() );

        // Project the particle stress gradients.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = p.node_ids[n];

            // Compute the force via the stress divergence.
            for ( int i = 0; i < space_dim; ++i )
                for ( int j = 0; j < space_dim; ++j )
                    node_f_int[node_id][i] -=
                        p.m * p.basis_gradients[n][j] * p.stress[j][i];
        }
    }
}

//---------------------------------------------------------------------------//
// Calculate nodal impulse.
void ProblemManager::calculateNodalImpulse(
    const std::vector<std::array<double,3> >& node_f_int,
    const std::vector<std::array<double,3> >& node_f_ext,
    const std::vector<double>& node_m,
    const double delta_t,
    std::vector<std::array<double,3> >& node_imp )
{
    int space_dim = d_mesh->spatialDimension();
    int num_nodes = d_mesh->totalNumNodes();

    // Calculate impulse.
    for ( int n = 0; n < num_nodes; ++n )
        for ( int d = 0; d < space_dim; ++d )
            node_imp[n][d] = delta_t * (node_f_int[n][d] + node_f_ext[n][d]);

    // Boundary conditions.
    for ( int b = 0; b < 6; ++b )
        d_bc[b]->evaluateImpulseCondition( d_mesh, b, node_m, node_imp );
}

//---------------------------------------------------------------------------//
// Update particle position and velocity.
void ProblemManager::updateParticlePositionAndVelocity(
    const std::vector<std::array<double,3> >& node_imp,
    const std::vector<std::array<double,3> >& node_p,
    const std::vector<double>& node_m,
    const double delta_t )
{
    int space_dim = d_mesh->spatialDimension();
    int nodes_per_cell = d_mesh->nodesPerCell();
    int node_id = 0;

    // Update the particles.
    for ( auto& p : d_particles )
    {
        // Loop over adjacent nodes.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = p.node_ids[n];

            // Only add a contribution from an adjacent node if it has mass.
            if ( node_m[node_id] > 0.0 )
            {
                for ( int d = 0; d < space_dim; ++d )
                {
                    // Increment the position.
                    p.r[d] += delta_t *
                              (node_p[node_id][d] + node_imp[node_id][d]) *
                              p.basis_values[n] / node_m[node_id];

                    // Increment the velocity. (FLIP Update)
                    p.v[d] += node_imp[node_id][d] * p.basis_values[n] /
                              node_m[node_id];
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end ProblemManager.cc
//---------------------------------------------------------------------------//
