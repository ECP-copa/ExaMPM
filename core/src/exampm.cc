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
#include <cmath>
#include <cassert>
#include <iostream>

//---------------------------------------------------------------------------//
// Locate the particles and compute grid values.
void locateParticles( const std::shared_ptr<ExaMPM::Mesh>& mesh,
                      std::vector<ExaMPM::Particle>& particles )
{
    int space_dim = mesh->spatialDimension();
    std::vector<int> cell_id( space_dim );
    std::vector<double> ref_coords( space_dim );

    for ( auto& p : particles )
    {
        // Locate the particle.
        mesh->locateParticle( p, cell_id );

        // Get the node ids local to the particle.
        mesh->cellNodeIds( cell_id, p.node_ids );

        // Map the particle to the reference frame of the cell.
        mesh->mapPhysicalToReferenceFrame( p, cell_id, ref_coords );

        // Evaluate the cell basis function at the particle location.
        mesh->shapeFunctionValue( ref_coords, p.basis_values );

        // Evaluate the cell basis function gradient at the particle location.
        mesh->shapeFunctionGradient( ref_coords, p.basis_gradients );
    }
}

//---------------------------------------------------------------------------//
// Calculate the nodal mass.
void calculateNodalMass( const std::shared_ptr<ExaMPM::Mesh>& mesh,
                         const std::vector<ExaMPM::Particle>& particles,
                         std::vector<double>& node_m )
{
    int nodes_per_cell = mesh->nodesPerCell();
    int node_id = 0;

    // Reset the nodal mass.
    std::fill( node_m.begin(), node_m.end(), 0.0 );

    // Compute the nodal mass.
    for ( auto& p : particles )
    {
        // Assemble the nodal mass.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = p.node_ids[n];
            node_m[node_id] += p.basis_values[n] * p.m;
        }
    }
}

//---------------------------------------------------------------------------//
// Calculate nodal momentum.
void calculateNodalMomentum(
    const std::shared_ptr<ExaMPM::Mesh>& mesh,
    const std::vector<ExaMPM::Particle>& particles,
    std::vector<std::vector<double> >& node_p )
{
    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    int node_id = 0;

    // Reset the momentum
    for ( auto& mom : node_p )
        std::fill( mom.begin(), mom.end(), 0.0 );

    // Update the momentum
    for ( auto& p : particles )
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
}

//---------------------------------------------------------------------------//
// Calculate external forces.
void calculateExternalForces(
    std::function<void(const std::vector<double>& r,std::vector<double>& v)> field,
    const std::shared_ptr<ExaMPM::Mesh>& mesh,
    const std::vector<ExaMPM::Particle>& particles,
    std::vector<std::vector<double> >& node_f_ext )
{
    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    std::vector<double> local_acceleration( space_dim, 0.0 );
    int node_id = 0;

    // Reset the forces.
    for ( auto& f_ext : node_f_ext )
        std::fill( f_ext.begin(), f_ext.end(), 0.0 );

    // Calculate forces.
    for ( auto& p : particles )
    {
        // Calculate force.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            field( p.r, local_acceleration );

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
// Calculate internal forces.
void calculateInternalForces(
    std::function<void(const ExaMPM::Particle&,
                       std::vector<std::vector<double> >&)> material_model,
    const std::shared_ptr<ExaMPM::Mesh>& mesh,
    const std::vector<ExaMPM::Particle>& particles,
    std::vector<std::vector<double> >& node_f_int )
{

    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    int node_id = 0;
    std::vector<std::vector<double> > stress(
        space_dim, std::vector<double>(space_dim) );

    // Reset the forces.
    for ( auto& f_int : node_f_int )
        std::fill( f_int.begin(), f_int.end(), 0.0 );

    // Compute forces.
    for ( auto& p : particles )
    {
        // Compute the first Piola-Kirchoff stress.
        material_model( p, stress );

        // Project the stress gradients.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = p.node_ids[n];

            for ( int d = 0; d < space_dim; ++d )
                for ( int j = 0; j < space_dim; ++j )
                    for ( int i = 0; i < space_dim; ++i )
                        node_f_int[node_id][d] -=
                            p.volume * stress[d][i] * p.F[j][i] * p.basis_gradients[n][j];
        }
    }
}

//---------------------------------------------------------------------------//
// Calculate nodal velocities.
void calculateNodalVelocities(
    const std::shared_ptr<ExaMPM::Mesh>& mesh,
    const std::vector<std::vector<double> >& node_p,
    const std::vector<double>& node_m,
    const std::vector<std::vector<double> >& node_f_int,
    const std::vector<std::vector<double> >& node_f_ext,
    const double delta_t,
    std::vector<std::vector<double> >& node_v )
{
    int space_dim = mesh->spatialDimension();
    int num_nodes = mesh->totalNumNodes();

    // Update the velocity
    for ( int n = 0; n < num_nodes; ++n )
    {
        // If we have nodal mass do the update.
        if ( node_m[n] > 0.0 )
        {
            for ( int d = 0; d < space_dim; ++d )
            {
                node_v[n][d] =
                    node_p[n][d] / node_m[n] +
                    delta_t * ( node_f_int[n][d] + node_f_ext[n][d] ) /
                    node_m[n];
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

    // Boundary conditions. Free slip for now.
    std::vector<int> boundary_nodes;
    for ( int d = 0; d < space_dim; ++d )
    {
        std::vector<int> bid( space_dim, 0 );

        // low boundary.
        bid[d] = -1;
        mesh->getBoundaryNodes( bid, boundary_nodes );
        for ( auto n : boundary_nodes )
            node_v[n][d] = 0.0;

        // high boundary.
        bid[d] = 1;
        mesh->getBoundaryNodes( bid, boundary_nodes );
        for ( auto n : boundary_nodes )
            node_v[n][d] = 0.0;
    }
}

//---------------------------------------------------------------------------//
// Update particle velocity. PIC update for now.
void updateParticles( const std::shared_ptr<ExaMPM::Mesh>& mesh,
                      const std::vector<std::vector<double> >& node_v,
                      const double delta_t,
                      std::vector<ExaMPM::Particle>& particles )
{
    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    int node_id = 0;
    std::vector<std::vector<double> > F_increment(
        space_dim, std::vector<double>(space_dim) );

    // Update the velocity.
    for ( auto& p : particles )
    {
        // Reset the velocity and deformation gradient increment.
        std::fill( p.v.begin(), p.v.end(), 0.0 );
        for ( int d = 0; d < space_dim; ++d )
            std::fill( F_increment[d].begin(), F_increment[d].end(), 0.0 );

        // Loop over adjacent nodes.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = p.node_ids[n];

            for ( int d = 0; d < space_dim; ++d )
            {
                // Increment the velocity
                p.v[d] += node_v[node_id][d] * p.basis_values[n];

                // Increment the position.
                p.r[d] += delta_t * node_v[node_id][d] * p.basis_values[n];

                // Increment the deformation gradient.
                for ( int j = 0; j < space_dim; ++j )
                    F_increment[d][j] += node_v[node_id][d] * p.basis_gradients[n][j];
            }
        }

        // Update the deformation gradient.
        for ( int i = 0; i < space_dim; ++i )
            for ( int j = 0; j < space_dim; ++j )
                p.F[i][j] *= ( i == j ? 1.0 : 0.0 ) + delta_t * F_increment[i][j];
    }
}

//---------------------------------------------------------------------------//
// Material model. Computes the first Piola-Kirchoff stress of a Neo-Hookean
// material.
void materialModel( const ExaMPM::Particle& particle,
                    std::vector<std::vector<double> >& stress )
{
    // youngs modulus
    double E = 1.0e8;

    // poisson ratio
    double nu = 0.4;

    // Calculate the jacobian of the deformation gradient.
    double J = particle.F[0][0] * particle.F[1][1] -
               particle.F[0][1] * particle.F[1][0];
    assert( J > 0.0 );

    // Calculate the transpse inverse of the deformation gradient.
    std::vector<std::vector<double> > F_T_inv = stress;
    F_T_inv[0][0] = particle.F[1][1] / J;
    F_T_inv[0][1] = -particle.F[1][0] / J;
    F_T_inv[1][0] = -particle.F[0][1] / J;
    F_T_inv[1][1] = particle.F[0][0] / J;

    // Compute the stress.
    double c = nu * std::log(J) - E;
    stress[0][0] = E * particle.F[0][0] + c * F_T_inv[0][0];
    stress[0][1] = E * particle.F[0][1] + c * F_T_inv[0][1];
    stress[1][0] = E * particle.F[1][0] + c * F_T_inv[1][0];
    stress[1][1] = E * particle.F[1][1] + c * F_T_inv[1][1];
}

//---------------------------------------------------------------------------//
int main( int argc, char *argv[] )
{
    // Set the spatial dimension.
    int space_dim = 2;

    // Create a mesh.
    int num_cells_x = 50;
    int num_cells_y = 50;
    double cell_width = 0.001;
    auto mesh = std::make_shared<ExaMPM::Mesh2d>(
        num_cells_x, num_cells_y, cell_width );

    // Get mesh data.
    int num_cells = mesh->totalNumCells();
    int num_nodes = mesh->totalNumNodes();
    int nodes_per_cell = mesh->nodesPerCell();

    // Create a geometry.
    std::vector<double> square_bnds = { 0.01, 0.03, 0.01, 0.03 };
    auto geom = std::make_shared<ExaMPM::Square>( square_bnds );

    // Assign properties to the geometry.
    int matid = 1;
    double density = 1.01e3;
    auto init_vf =
        [=](const std::vector<double>& r,std::vector<double>& v)
        { v[0] = 0.0; v[1] = -10.0; };
    geom->setMatId( matid );
    geom->setVelocityField( init_vf );
    geom->setDensity( density );

    // Initialize the particles in the geometry.
    std::vector<ExaMPM::Particle> particles;
    int order = 2;
    int ppcell = mesh->particlesPerCell( order );
    std::vector<ExaMPM::Particle> candidates(
        ppcell, ExaMPM::Particle(space_dim,nodes_per_cell) );
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
    int num_step = 10000;
    double delta_t = 1.0e-6;
    int num_write = 50;
    int write_freq = num_step / num_write;
    for ( int step = 0; step < num_step; ++step )
    {
        time += delta_t;

        if ( (step+1) % write_freq == 0 )
        {
            std::cout << "Time Step " << step+1 <<  ": "
                      <<  time <<  " (s)" << std::endl;
        }

        // 1) Locate particles and evaluate basis functions.
        locateParticles( mesh, particles );

        // 2) Calculate the nodal masses.
        calculateNodalMass( mesh, particles, node_m );

        // 3) Calculate nodal momentum.
        calculateNodalMomentum( mesh, particles, node_p );

        // 4) Calculate internal forces.
        calculateInternalForces( materialModel, mesh, particles, node_f_int );

        // 5) Calculate external forces.
        calculateExternalForces( gravity_field, mesh, particles, node_f_ext );

        // 6) Calculate nodal velocities.
        calculateNodalVelocities(
            mesh, node_p, node_m, node_f_int, node_f_ext, delta_t, node_v );

        // 7) Update the particle quantities.
        updateParticles( mesh, node_v, delta_t, particles );

        // Write the time step.
        if ( (step+1) % write_freq == 0 )
        {
            ++write_step;
            file_io.writeTimeStep( write_step, time, particles );
        }
    }

    file_io.writeTimeStep( write_step+1, time, particles );


    return 0;
}

//---------------------------------------------------------------------------//
// end exampm.cc
//---------------------------------------------------------------------------//
