//---------------------------------------------------------------------------//
/*!
 * \file exampm.cc
 * \brief driver
 */
//---------------------------------------------------------------------------//

#include "Particle.hh"
#include "Mesh.hh"
#include "Mesh2d.hh"
#include "Square.hh"
#include "FileIO.hh"

#include <memory>
#include <vector>
#include <algorithm>
#include <functional>
#include <cassert>
#include <iostream>

//---------------------------------------------------------------------------//
// Calculate the nodal mass.
void calculateNodalMass( const std::shared_ptr<ExaMPM::Mesh>& mesh,
                         const std::vector<ExaMPM::Particle>& particles,
                         std::vector<double>& node_m )
{
    int num_nodes = node_m.size();
    int num_p = particles.size();
    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    std::vector<int> cell_id( space_dim );
    std::vector<int> cell_nodes( nodes_per_cell );
    std::vector<double> ref_coords( space_dim );
    std::vector<double> shape_values( nodes_per_cell );
    int node_id = 0;

    // Reset the nodal mass.
    std::fill( node_m.begin(), node_m.end(), 0.0 );

    // Compute the nodal mass.
    for ( int p = 0; p < num_p; ++p )
    {
        // Get the particle location.
        mesh->locateParticle( particles[p], cell_id );

        // Get the cell nodes.
        mesh->cellNodeIds( cell_id, cell_nodes );

        // Project the particle to the reference frame of the cell.
        mesh->mapPhysicalToReferenceFrame( particles[p],
                                           cell_id,
                                           ref_coords );

        // Evaluate the cell shape function at the particle location.
        mesh->shapeFunctionValue( ref_coords, shape_values );

        // Assemble the nodal mass.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = cell_nodes[n];
            node_m[node_id] += shape_values[n] * particles[p].m;
        }
    }
}

//---------------------------------------------------------------------------//
// Calculate external forces.
void calculateExternalForces(
    std::function<void(const std::vector<double>& r,std::vector<double>& v)> field,
    const std::shared_ptr<ExaMPM::Mesh>& mesh,
    const std::vector<ExaMPM::Particle>& particles,
    std::vector<std::vector<double> >& node_f_ext )
{
    int num_p = particles.size();
    int num_nodes = mesh->totalNumNodes();
    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    std::vector<int> cell_id( space_dim );
    std::vector<int> cell_nodes( nodes_per_cell );
    std::vector<double> ref_coords( space_dim );
    std::vector<double> shape_values( nodes_per_cell );
    std::vector<double> local_acceleration( space_dim, 0.0 );
    int node_id = 0;

    // Reset the forces.
    for ( int n = 0; n < num_nodes; ++n )
        std::fill( node_f_ext[n].begin(), node_f_ext[n].end(), 0.0 );

    // Calculate forces.
    for ( int p = 0; p < num_p; ++p )
    {
        // Get the particle location.
        mesh->locateParticle( particles[p], cell_id );

        // Get the cell nodes.
        mesh->cellNodeIds( cell_id, cell_nodes );

        // Project the particle to the reference frame of the cell.
        mesh->mapPhysicalToReferenceFrame( particles[p],
                                           cell_id,
                                           ref_coords );

        // Evaluate the cell shape function at the particle location.
        mesh->shapeFunctionValue( ref_coords, shape_values );

        // Calculate force.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            field( particles[p].r, local_acceleration );

            node_id = cell_nodes[n];
            for ( int d = 0; d < space_dim; ++d )
            {
                node_f_ext[node_id][d] +=
                    particles[p].m * local_acceleration[d] * shape_values[n];
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Calculate internal forces.
void calculateInternalForces( const std::shared_ptr<ExaMPM::Mesh>& mesh,
                              const std::vector<ExaMPM::Particle>& particles,
                              std::vector<std::vector<double> >& node_f_int )
{

    int num_p = particles.size();
    int num_nodes = mesh->totalNumNodes();
    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    std::vector<int> cell_id( space_dim );
    std::vector<int> cell_nodes( nodes_per_cell );
    std::vector<double> ref_coords( space_dim );
    std::vector<std::vector<double> > shape_gradients(
        nodes_per_cell, std::vector<double>(space_dim) );
    int node_id = 0;

    // Reset the forces.
    for ( int n = 0; n < num_nodes; ++n )
        std::fill( node_f_int[n].begin(), node_f_int[n].end(), 0.0 );

    // Compute forces.
    for ( int p = 0; p < num_p; ++p )
    {
        // Get the particle location.
        mesh->locateParticle( particles[p], cell_id );

        // Get the cell nodes.
        mesh->cellNodeIds( cell_id, cell_nodes );

        // Project the particle to the reference frame of the cell.
        mesh->mapPhysicalToReferenceFrame( particles[p],
                                           cell_id,
                                           ref_coords );

        // Evaluate the cell shape function gradient at the particle location.
        mesh->shapeFunctionGradient( ref_coords, shape_gradients );

        // Project the stress gradients.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = cell_nodes[n];

            for ( int d = 0; d < space_dim; ++d )
            {
                for ( int j = 0; j < space_dim; ++j )
                {
                    node_f_int[node_id][d] -=
                        particles[p].m * shape_gradients[n][j] *
                        particles[p].stress[d][j];
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Calculate nodal momentum.
void calculateNodalMomentum(
    const std::shared_ptr<ExaMPM::Mesh>& mesh,
    const std::vector<ExaMPM::Particle>& particles,
    const std::vector<std::vector<double> >& node_f_int,
    const std::vector<std::vector<double> >& node_f_ext,
    const double delta_t,
    std::vector<std::vector<double> >& node_p )
{
    int num_p = particles.size();
    int num_nodes = mesh->totalNumNodes();
    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    std::vector<int> cell_id( space_dim );
    std::vector<int> cell_nodes( nodes_per_cell );
    std::vector<double> ref_coords( space_dim );
    std::vector<double> shape_values( nodes_per_cell );
    int node_id = 0;

    // Reset the momentum
    for ( int n = 0; n < num_nodes; ++n )
        std::fill( node_p[n].begin(), node_p[n].end(), 0.0 );

    // Update the momentum
    for ( int p = 0; p < num_p; ++p )
    {
        // Get the particle location.
        mesh->locateParticle( particles[p], cell_id );

        // Get the cell nodes.
        mesh->cellNodeIds( cell_id, cell_nodes );

        // Project the particle to the reference frame of the cell.
        mesh->mapPhysicalToReferenceFrame( particles[p],
                                           cell_id,
                                           ref_coords );

        // Evaluate the cell shape function at the particle location.
        mesh->shapeFunctionValue( ref_coords, shape_values );

        // Calculate momentum.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = cell_nodes[n];

            for ( int d = 0; d < space_dim; ++d )
            {
                node_p[node_id][d] +=
                    particles[p].m * particles[p].v[d] * shape_values[n];
            }
        }
    }

    // Update with forces.
    for ( int n = 0; n < num_nodes; ++n )
    {
        for ( int d = 0; d < space_dim; ++d )
        {
            node_p[n][d] += delta_t * (node_f_int[n][d] + node_f_ext[n][d]);
        }
    }

    // Apply boundary conditions. No slip for now.
    std::vector<int> boundary_nodes;
    for ( int d = 0; d < space_dim; ++d )
    {
        std::vector<int> bid( space_dim, 0 );

        // low boundary.
        bid[d] = -1;
        mesh->getBoundaryNodes( bid, boundary_nodes );
        for ( auto n : boundary_nodes )
            for ( int d = 0; d < space_dim; ++d )
                node_p[n][d] = 0.0;

        // high boundary.
        bid[d] = 1;
        mesh->getBoundaryNodes( bid, boundary_nodes );
        for ( auto n : boundary_nodes )
            for ( int d = 0; d < space_dim; ++d )
                node_p[n][d] = 0.0;
    }
}

//---------------------------------------------------------------------------//
// Update particle velocity.
void updateParticleVelocity( const std::shared_ptr<ExaMPM::Mesh>& mesh,
                             const std::vector<std::vector<double> >& node_f_int,
                             const std::vector<std::vector<double> >& node_f_ext,
                             const std::vector<double>& node_m,
                             const double delta_t,
                             std::vector<ExaMPM::Particle>& particles )
{
    int num_p = particles.size();
    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    std::vector<int> cell_id( space_dim );
    std::vector<int> cell_nodes( nodes_per_cell );
    std::vector<double> ref_coords( space_dim );
    std::vector<double> shape_values( nodes_per_cell );
    int node_id = 0;

    // Update the velocity.
    for ( int p = 0; p < num_p; ++p )
    {
        // Get the particle location.
        mesh->locateParticle( particles[p], cell_id );

        // Get the cell nodes.
        mesh->cellNodeIds( cell_id, cell_nodes );

        // Project the particle to the reference frame of the cell.
        mesh->mapPhysicalToReferenceFrame( particles[p],
                                           cell_id,
                                           ref_coords );

        // Evaluate the cell shape function at the particle location.
        mesh->shapeFunctionValue( ref_coords, shape_values );

        // Increment the velocity
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = cell_nodes[n];

            assert( node_m[node_id] > 0.0 );

            for ( int d = 0; d < space_dim; ++d )
            {
                particles[p].v[d] +=
                    delta_t * (node_f_int[node_id][d] + node_f_ext[node_id][d]) *
                    shape_values[n] / node_m[node_id];
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Calculate nodal velocities.
void calculateNodalVelocities(
    const std::shared_ptr<ExaMPM::Mesh>& mesh,
    const std::vector<ExaMPM::Particle>& particles,
    const std::vector<double>& node_m,
    std::vector<std::vector<double> >& node_v )
{
    int num_p = particles.size();
    int num_nodes = mesh->totalNumNodes();
    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    std::vector<int> cell_id( space_dim );
    std::vector<int> cell_nodes( nodes_per_cell );
    std::vector<double> ref_coords( space_dim );
    std::vector<double> shape_values( nodes_per_cell );
    int node_id = 0;

    // Reset the velocities
    for ( int n = 0; n < num_nodes; ++n )
        std::fill( node_v[n].begin(), node_v[n].end(), 0.0 );

    // Update the velocity
    for ( int p = 0; p < num_p; ++p )
    {
        // Get the particle location.
        mesh->locateParticle( particles[p], cell_id );

        // Get the cell nodes.
        mesh->cellNodeIds( cell_id, cell_nodes );

        // Project the particle to the reference frame of the cell.
        mesh->mapPhysicalToReferenceFrame( particles[p],
                                           cell_id,
                                           ref_coords );

        // Evaluate the cell shape function at the particle location.
        mesh->shapeFunctionValue( ref_coords, shape_values );

        // Calculate velocity
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = cell_nodes[n];

            for ( int d = 0; d < space_dim; ++d )
            {
                node_v[node_id][d] +=
                    particles[p].m * particles[p].v[d] * shape_values[n] / node_m[node_id];
            }
        }
    }

    // Apply boundary conditions. No slip for now.
    std::vector<int> boundary_nodes;
    for ( int d = 0; d < space_dim; ++d )
    {
        std::vector<int> bid( space_dim, 0 );

        // low boundary.
        bid[d] = -1;
        mesh->getBoundaryNodes( bid, boundary_nodes );
        for ( auto n : boundary_nodes )
            for ( int d = 0; d < space_dim; ++d )
                node_v[n][d] = 0.0;

        // high boundary.
        bid[d] = 1;
        mesh->getBoundaryNodes( bid, boundary_nodes );
        for ( auto n : boundary_nodes )
            for ( int d = 0; d < space_dim; ++d )
                node_v[n][d] = 0.0;
    }
}

//---------------------------------------------------------------------------//
// Update particle stress and strain.
void updateParticleStressStrain(
    std::function<void(const std::vector<std::vector<double> >&,
                       std::vector<std::vector<double> >&)> material_model,
    const std::shared_ptr<ExaMPM::Mesh>& mesh,
    const std::vector<std::vector<double> >& node_v,
    const double delta_t,
    std::vector<ExaMPM::Particle>& particles )
{
    int num_p = particles.size();
    int num_nodes = mesh->totalNumNodes();
    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    std::vector<int> cell_id( space_dim );
    std::vector<int> cell_nodes( nodes_per_cell );
    std::vector<double> ref_coords( space_dim );
    std::vector<std::vector<double> > shape_gradients(
        nodes_per_cell, std::vector<double>(space_dim) );
    int node_id = 0;
    std::vector<std::vector<double> > strain_increment(
        space_dim, std::vector<double>(space_dim) );
    std::vector<std::vector<double> > stress_increment(
        space_dim, std::vector<double>(space_dim) );
    double strain_trace = 0.0;

    // Update the stress and strain.
    for ( int p = 0; p < num_p; ++p )
    {
        // Reset the strain incrememnt.
        for ( int d = 0; d < space_dim; ++d )
            std::fill( strain_increment[d].begin(), strain_increment[d].end(), 0.0 );

        // Get the particle location.
        mesh->locateParticle( particles[p], cell_id );

        // Get the cell nodes.
        mesh->cellNodeIds( cell_id, cell_nodes );

        // Project the particle to the reference frame of the cell.
        mesh->mapPhysicalToReferenceFrame( particles[p],
                                           cell_id,
                                           ref_coords );

        // Evaluate the cell shape function gradient at the particle location.
        mesh->shapeFunctionGradient( ref_coords, shape_gradients );

        // Calculate the local strain increment at the particle with small
        // deformation theory.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = cell_nodes[n];

            for ( int j = 0; j < space_dim; ++j )
            {
                for ( int i = 0; i < space_dim; ++i )
                {
                    strain_increment[i][j] +=
                        (shape_gradients[n][i] * node_v[node_id][j] +
                         shape_gradients[n][j] * node_v[node_id][i]) * delta_t / 2.0;
                }
            }
        }

        // Increment the particle strain.
        for ( int j = 0; j < space_dim; ++j )
        {
            for ( int i = 0; i < space_dim; ++i )
            {
                particles[p].strain[i][j] += strain_increment[i][j];
            }
        }

        // Compute stress increment.
        material_model( strain_increment, stress_increment );

        // Increment the particle stress.
        for ( int j = 0; j < space_dim; ++j )
        {
            for ( int i = 0; i < space_dim; ++i )
            {
                particles[p].stress[i][j] += stress_increment[i][j] / particles[p].rho;
            }
        }

        // Update the particle density.
        strain_trace = 0.0;
        for ( int d = 0; d < space_dim; ++d )
        {
            strain_trace += strain_increment[d][d];
        }
        particles[p].rho /= (1 + strain_trace);
    }
}

//---------------------------------------------------------------------------//
// Update particle position.
void updateParticlePosition( const std::shared_ptr<ExaMPM::Mesh>& mesh,
                             const std::vector<std::vector<double> >& node_p,
                             const std::vector<double>& node_m,
                             const double delta_t,
                             std::vector<ExaMPM::Particle>& particles )
{
    int num_p = particles.size();
    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    std::vector<int> cell_id( space_dim );
    std::vector<int> cell_nodes( nodes_per_cell );
    std::vector<double> ref_coords( space_dim );
    std::vector<double> shape_values( nodes_per_cell );
    int node_id = 0;

    // Update the position.
    for ( int p = 0; p < num_p; ++p )
    {
        // Get the particle location.
        mesh->locateParticle( particles[p], cell_id );

        // Get the cell nodes.
        mesh->cellNodeIds( cell_id, cell_nodes );

        // Project the particle to the reference frame of the cell.
        mesh->mapPhysicalToReferenceFrame( particles[p],
                                           cell_id,
                                           ref_coords );

        // Evaluate the cell shape function at the particle location.
        mesh->shapeFunctionValue( ref_coords, shape_values );

        // Increment the position.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = cell_nodes[n];

            assert( node_m[node_id] > 0.0 );

            for ( int d = 0; d < space_dim; ++d )
            {
                particles[p].r[d] +=
                    delta_t * node_p[node_id][d] *
                    shape_values[n] / node_m[node_id];
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Material model.
void materialModel( const std::vector<std::vector<double> >& strain_rate,
                    std::vector<std::vector<double> >& stress_rate )
{
    // youngs modulus
    double E = 0.05e9;

    // bulk modulus
    double B = 1.5e9;

    // poisson ratio
    double nu = 0.48;

    // 2d plane strain
    stress_rate[0][0] = B * (strain_rate[0][0] + nu*strain_rate[1][1]) /
                        ( 1 - nu*nu );

    stress_rate[1][1] = B * (nu*strain_rate[0][0] + strain_rate[1][1]) /
                        ( 1 - nu*nu );

    stress_rate[0][1] = B * strain_rate[0][1] / (1 + nu);

    stress_rate[1][0] = B * strain_rate[1][0] / (1 + nu);

    // stress_rate[0][0] = E*strain_rate[0][0];
    // stress_rate[0][1] = E*strain_rate[0][1];
    // stress_rate[1][0] = E*strain_rate[1][0];
    // stress_rate[1][1] = E*strain_rate[1][1];
}

//---------------------------------------------------------------------------//
int main( int argc, char *argv[] )
{
    // Set the spatial dimension.
    int space_dim = 2;

    // Create a mesh.
    int num_cells_x = 50;
    int num_cells_y = 30;
    double cell_width = 0.001;
    auto mesh = std::make_shared<ExaMPM::Mesh2d>(
        num_cells_x, num_cells_y, cell_width );

    // Get mesh data.
    int num_cells = mesh->totalNumCells();
    int num_nodes = mesh->totalNumNodes();

    // Create a geometry.
    std::vector<double> square_bnds = { 0.01, 0.03, 0.005, 0.025 };
    auto geom = std::make_shared<ExaMPM::Square>( square_bnds );

    // Assign properties to the geometry.
    int matid = 1;
    double density = 2000.0;
    double mass = 0.016;
    auto init_vf =
        [=](const std::vector<double>& r,std::vector<double>& v)
        { v[0] = 0.0; v[1] = 0.0; };
    geom->setMatId( matid );
    geom->setVelocityField( init_vf );
    geom->setDensity( density );
    geom->setMass( mass );

    // Initialize the particles in the geometry.
    std::vector<ExaMPM::Particle> particles;
    int order = 2;
    int ppcell = mesh->particlesPerCell( order );
    std::vector<ExaMPM::Particle> candidates( ppcell,
                                              ExaMPM::Particle(space_dim) );
    for ( int c = 0; c < num_cells; ++c )
    {
        // Create the particles.
        mesh->initializeParticles( c, order, candidates );

        // If they are in the geometry add them to the list.
        for ( int p = 0; p < ppcell; ++p )
        {
            if ( geom->particleInGeometry(candidates[p]) )
            {
                geom->initializeParticle(candidates[p]);
                particles.push_back( candidates[p] );
            }
        }
    }

    // Set the mass of the particles based on how many constructed the
    // geometry.
    int num_p = particles.size();
    double mass_p = geom->getMass() / num_p;
    for ( auto& p : particles ) p.m = mass_p;

    // Gravity acceleration function.
    auto gravity_field = [](const std::vector<double>& r,std::vector<double>& f)
                         { f.back() = -9.81; };

    // Nodal mass
    std::vector<double> node_m( num_nodes, 0.0 );

    // Nodal velocity
    std::vector<std::vector<double> > node_v(
        num_nodes, std::vector<double>(space_dim,0.0) );

    // Nodal momentum
    std::vector<std::vector<double> > node_p(
        num_nodes, std::vector<double>(space_dim,0.0) );

    // Nodal internal force.
    std::vector<std::vector<double> > node_f_int(
        num_nodes, std::vector<double>(space_dim,0.0) );

    // Nodal external force
    std::vector<std::vector<double> > node_f_ext(
        num_nodes, std::vector<double>(space_dim,0.0) );

    // Setup an output writer.
    ExaMPM::FileIO file_io( "particles.h5" );

    // Time step parameters
    double time = 0.0;

    // Write the initial state.
    int write_step = 0;
    file_io.writeTimeStep( write_step, time, particles );

    // Time step
    int num_step = 50000;
    double delta_t = 1.0e-6;
    int num_write = 200;
    int write_freq = num_step / num_write;
    for ( int step = 0; step < num_step; ++step )
    {
        time += delta_t;

        if ( (step+1) % write_freq == 0 )
        {
            std::cout << "Time Step " << step+1 <<  ": "
                      <<  time <<  " (s)" << std::endl;
        }

        // First calculate the nodal masses.
        calculateNodalMass( mesh, particles, node_m );

        // Calculate internal forces.
        calculateInternalForces( mesh, particles, node_f_int );

        // Calculate external forces.
        calculateExternalForces( gravity_field, mesh, particles, node_f_ext );

        // Calculate nodal momentum.
        calculateNodalMomentum( mesh, particles, node_f_int, node_f_ext, delta_t, node_p );

        // Update the particle velocity.
        updateParticleVelocity( mesh, node_f_int, node_f_ext, node_m, delta_t, particles );

        // Calculate nodal velocities.
        calculateNodalVelocities( mesh, particles, node_m, node_v );

        // Update the particle stress and strain.
        updateParticleStressStrain(
            materialModel, mesh, node_v, delta_t, particles );

        // Update the particle position.
        updateParticlePosition( mesh, node_p, node_m, delta_t, particles );

        // Write the time step.
        if ( (step+1) % write_freq == 0 )
        {
            ++write_step;
            file_io.writeTimeStep( write_step, time, particles );
        }
    }

    return 0;
}

//---------------------------------------------------------------------------//
// end exampm.cc
//---------------------------------------------------------------------------//
