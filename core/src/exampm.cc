//---------------------------------------------------------------------------//
/*!
 * \file exampm.cc
 * \brief driver
 */
//---------------------------------------------------------------------------//

#include "Particle.hh"
#include "Mesh.hh"
#include "Box.hh"
#include "FileIO.hh"

#include <memory>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <cmath>
#include <cassert>
#include <iostream>

//---------------------------------------------------------------------------//
// Calculate the determinant of a 3x3 matrix.
double determinant( const std::array<std::array<double,3>,3>& m )
{
    return
        m[0][0] * m[1][1] * m[2][2] +
        m[0][1] * m[1][2] * m[2][0] +
        m[0][2] * m[1][0] * m[2][1] -
        m[0][2] * m[1][1] * m[2][0] -
        m[0][1] * m[1][0] * m[2][2] -
        m[0][0] * m[1][2] * m[2][1];
}

//---------------------------------------------------------------------------//
// Initialize particles.
void initializeParticles(
    const std::vector<std::shared_ptr<ExaMPM::Geometry> >& geom,
    const std::shared_ptr<ExaMPM::Mesh>& mesh,
    const int order,
    std::vector<ExaMPM::Particle>& particles )
{
    int ppcell = mesh->particlesPerCell( order );
    std::vector<ExaMPM::Particle> candidates( ppcell );
    int num_cells = mesh->totalNumCells();
    for ( int c = 0; c < num_cells; ++c )
    {
        // Create the particles.
        mesh->initializeParticles( c, order, candidates );

        // If they are in the geometry add them to the list.
        for ( int p = 0; p < ppcell; ++p )
        {
            for ( const auto& g : geom )
            {
                if ( g->particleInGeometry(candidates[p]) )
                {
                    g->initializeParticle(candidates[p]);
                    particles.push_back( candidates[p] );
                    break;
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Locate the particles and compute grid values.
void locateParticles( const std::shared_ptr<ExaMPM::Mesh>& mesh,
                      std::vector<ExaMPM::Particle>& particles )
{
    int space_dim = mesh->spatialDimension();
    std::array<int,3> cell_id;
    std::array<double,3> ref_coords;

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

    // Boundary conditions. Free slip for now.
    std::vector<int> boundary_nodes;
    std::array<int,3> bid;
    for ( int d = 0; d < space_dim; ++d )
    {
        bid = {0,0,0};

        // low boundary.
        bid[d] = -1;
        mesh->getBoundaryNodes( bid, boundary_nodes );
        for ( auto n : boundary_nodes )
            node_p[n][d] = 0.0;

        // high boundary.
        bid[d] = 1;
        mesh->getBoundaryNodes( bid, boundary_nodes );
        for ( auto n : boundary_nodes )
            node_p[n][d] = 0.0;
    }
}

//---------------------------------------------------------------------------//
// Calculate nodal velocities.
void calculateNodalVelocity(
    const std::shared_ptr<ExaMPM::Mesh>& mesh,
    const std::vector<std::vector<double> >& node_p,
    const std::vector<double>& node_m,
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
// Update particle deformation gradient.
void updateDeformationGradient(
    const std::shared_ptr<ExaMPM::Mesh>& mesh,
    const std::vector<std::vector<double> >& node_v,
    const double delta_t,
    std::vector<ExaMPM::Particle>& particles )
{
    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    int node_id = 0;
    std::array<std::array<double,3>,3> V_grad;
    std::array<std::array<double,3>,3> delta_F;

    // Update the velocity.
    for ( auto& p : particles )
    {
        // Reset the deformation gradient increment and velocity gradient.
        for ( int d = 0; d < space_dim; ++d )
        {
            std::fill( V_grad[d].begin(), V_grad[d].end(), 0.0 );
            std::fill( delta_F[d].begin(), delta_F[d].end(), 0.0 );
        }

        // Loop over adjacent nodes.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = p.node_ids[n];

            for ( int i = 0; i < space_dim; ++i )
            {
                // Compute the velocity gradient.
                for ( int j = 0; j < space_dim; ++j )
                    V_grad[i][j] +=
                        p.basis_gradients[n][i] * node_v[node_id][j];
            }
        }

        // Scale the velocity gradient.
        for ( int i = 0; i < space_dim; ++i )
            for ( int j = 0; j < space_dim; ++j )
                V_grad[i][j] *= delta_t;

        // Compute the deformation gradient increment.
        for ( int i = 0; i < space_dim; ++i )
            for ( int j = 0; j < space_dim; ++j )
                for ( int k = 0; k < space_dim; ++k )
                    delta_F[i][j] += V_grad[i][k] * p.F[k][j];

        // Update the deformation gradient.
        for ( int i = 0; i < space_dim; ++i )
            for ( int j = 0; j < space_dim; ++j )
                p.F[i][j] += delta_F[i][j];

        // Add the velocity gradient to the identity.
        for ( int i = 0; i < space_dim; ++i )
            V_grad[i][i] += 1.0;

        // Update the particle volume.
        p.volume *= determinant( V_grad );
    }
}

//---------------------------------------------------------------------------//
// Calculate external forces.
void calculateExternalForces(
    std::function<void(const std::array<double,3>& r,std::array<double,3>& v)> field,
    const std::shared_ptr<ExaMPM::Mesh>& mesh,
    const std::vector<ExaMPM::Particle>& particles,
    std::vector<std::vector<double> >& node_f_ext )
{
    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    std::array<double,3> local_acceleration;
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
            local_acceleration = {0.0,0.0,0.0};
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
                       std::array<std::array<double,3>,3>&)> material_model,
    const std::shared_ptr<ExaMPM::Mesh>& mesh,
    const std::vector<ExaMPM::Particle>& particles,
    std::vector<std::vector<double> >& node_f_int )
{

    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    int node_id = 0;
    std::array<std::array<double,3>,3> stress;

    // Reset the forces.
    for ( auto& f_int : node_f_int )
        std::fill( f_int.begin(), f_int.end(), 0.0 );

    // Compute forces.
    for ( auto& p : particles )
    {
        // Compute the stress based on the particle state..
        material_model( p, stress );

        // Project the stress gradients.
        for ( int n = 0; n < nodes_per_cell; ++n )
        {
            node_id = p.node_ids[n];

            // Compute the force via the stress divergence.
            for ( int i = 0; i < space_dim; ++i )
                for ( int j = 0; j < space_dim; ++j )
                    node_f_int[node_id][i] -=
                        p.volume * p.basis_gradients[n][j] * stress[j][i];
        }
    }
}

//---------------------------------------------------------------------------//
// Calculate nodal impulse.
void calculateNodalImpulse(
    const std::shared_ptr<ExaMPM::Mesh>& mesh,
    const std::vector<std::vector<double> >& node_f_int,
    const std::vector<std::vector<double> >& node_f_ext,
    const double delta_t,
    std::vector<std::vector<double> >& node_imp )
{
    int space_dim = mesh->spatialDimension();
    int num_nodes = mesh->totalNumNodes();

    // Calculate impulse.
    for ( int n = 0; n < num_nodes; ++n )
        for ( int d = 0; d < space_dim; ++d )
            node_imp[n][d] = delta_t * (node_f_int[n][d] + node_f_ext[n][d]);

    // Boundary conditions. Free slip for now.
    std::vector<int> boundary_nodes;
    std::array<int,3> bid;
    for ( int d = 0; d < space_dim; ++d )
    {
        bid = {0,0,0};

        // low boundary.
        bid[d] = -1;
        mesh->getBoundaryNodes( bid, boundary_nodes );
        for ( auto n : boundary_nodes )
            node_imp[n][d] = 0.0;

        // high boundary.
        bid[d] = 1;
        mesh->getBoundaryNodes( bid, boundary_nodes );
        for ( auto n : boundary_nodes )
            node_imp[n][d] = 0.0;
    }
}

//---------------------------------------------------------------------------//
// Update particle velocity and position.
void updateParticles( const std::shared_ptr<ExaMPM::Mesh>& mesh,
                      const std::vector<std::vector<double> >& node_imp,
                      const std::vector<std::vector<double> >& node_p,
                      const std::vector<double>& node_m,
                      const double delta_t,
                      std::vector<ExaMPM::Particle>& particles )
{
    int space_dim = mesh->spatialDimension();
    int nodes_per_cell = mesh->nodesPerCell();
    int node_id = 0;

    // Update the velocity.
    for ( auto& p : particles )
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

                    // Increment the velocity.
                    p.v[d] += node_imp[node_id][d] * p.basis_values[n] /
                              node_m[node_id];
                }
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Material model. Computes the first Piola-Kirchoff stress of a Neo-Hookean
// material and then transforms it to the physical frame.
void materialModel( const ExaMPM::Particle& p,
                    std::array<std::array<double,3>,3>& stress )
{
    // youngs modulus
    double E = 1e4;

    // poisson ratio
    double nu = 0.3;

    // Calculate the determinant of the deformation gradient.
    double J = determinant( p.F );
    assert( J > 0.0 );

    // Reset the stress.
    for ( int d = 0; d < 3; ++d )
        stress[d] = { 0.0, 0.0, 0.0 };

    // Calculate F*F^T
    for ( int j = 0; j < 3; ++j )
        for ( int i = 0; i < 3; ++i )
            for ( int k = 0; k < 3; ++k )
                stress[i][j] += p.F[i][k] * p.F[j][k];

    // Subtract the identity.
    for ( int i = 0; i < 3; ++i )
        stress[i][i] -= 1.0;

    // Scale by youngs modulus.
    for ( int j = 0; j < 3; ++j )
        for ( int i = 0; i < 3; ++i )
            stress[i][j] *= E;

    // Add the scaled identity.
    double c = nu * std::log(J);
    for ( int i = 0; i < 3; ++i )
        stress[i][i] += c;

    // Scale by the determinant of the deformation gradient.
    for ( int j = 0; j < 3; ++j )
        for ( int i = 0; i < 3; ++i )
            stress[i][j] /= J;
}

//---------------------------------------------------------------------------//
int main( int argc, char *argv[] )
{
    // Set the spatial dimension.
    int space_dim = 3;

    // Create a mesh.
    int num_cells_x = 20;
    int num_cells_y = 20;
    int num_cells_z = 20;
    double cell_width = 0.05;
    auto mesh = std::make_shared<ExaMPM::Mesh>(
        num_cells_x, num_cells_y, num_cells_z, cell_width );

    // Get mesh data.
    int num_nodes = mesh->totalNumNodes();

    // Geometry
    std::vector<std::shared_ptr<ExaMPM::Geometry> > geom;

    // Create geometries.
    std::array<double,6> bnds = {0.1,0.3,0.5,0.7,0.5,0.7};
    geom.push_back( std::make_shared<ExaMPM::Box>(bnds) );
    bnds = {0.7,0.9,0.4,0.6,0.4,0.6};
    geom.push_back( std::make_shared<ExaMPM::Box>(bnds) );

    // Assign properties to the geometry.
    int matid = 1;
    double density = 1.0e3;
    auto init_vf1 =
        [=](const std::array<double,3>& r,std::array<double,3>& v)
        { v[0] = 0.1; v[1] = 0.0; v[2] = 0.0; };
    geom[0]->setMatId( matid );
    geom[0]->setVelocityField( init_vf1 );
    geom[0]->setDensity( density );
    auto init_vf2 =
        [=](const std::array<double,3>& r,std::array<double,3>& v)
        { v[0] = -0.1; v[1] = 0.0; v[2] = 0.0; };
    geom[1]->setMatId( matid );
    geom[1]->setVelocityField( init_vf2 );
    geom[1]->setDensity( density );

    // Initialize the particles in the geometry.
    std::vector<ExaMPM::Particle> particles;
    int order = 2;
    initializeParticles( geom, mesh, order, particles );

    // Gravity acceleration function.
    auto gravity_field = [](const std::array<double,3>& r,std::array<double,3>& f)
                         { f[0] = 0.0; f[1] = 0.0; f[2] = 0.0; };

    // Nodal mass
    std::vector<double> node_m( num_nodes, 0.0 );

    // Nodal velocity
    std::vector<std::vector<double> > node_v(
        num_nodes, std::vector<double>(space_dim,0.0) );

    // Nodal momentum
    std::vector<std::vector<double> > node_p(
        num_nodes, std::vector<double>(space_dim,0.0) );

    // Nodal impulse
    std::vector<std::vector<double> > node_imp(
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
    double delta_t = 1.0e-4;
    int num_write = 50;
    int write_freq = num_step / num_write;
    for ( int step = 0; step < num_step; ++step )
    {
        time += delta_t;

        if ( (step+1) % write_freq == 0 )
        {
            std::cout << "Time Step " << step+1 << "/" << num_step << ": "
                      <<  time <<  " (s)" << std::endl;
        }

        // 1) Locate particles and evaluate basis functions.
        locateParticles( mesh, particles );

        // 2) Calculate the nodal masses.
        calculateNodalMass( mesh, particles, node_m );

        // 3) Calculate nodal momentum.
        calculateNodalMomentum( mesh, particles, node_p );

        // 3) Calculate nodal velocity.
        calculateNodalVelocity( mesh, node_p, node_m, node_v );

        // 3) Update the particle deformation graadient.
        updateDeformationGradient( mesh, node_v, delta_t, particles );

        // 4) Calculate internal forces.
        calculateInternalForces( materialModel, mesh, particles, node_f_int );

        // 5) Calculate external forces.
        calculateExternalForces( gravity_field, mesh, particles, node_f_ext );

        // 6) Calculate node impulse.
        calculateNodalImpulse(
            mesh, node_f_int, node_f_ext, delta_t, node_imp );

        // 7) Update the particle quantities.
        updateParticles( mesh, node_imp, node_p, node_m, delta_t, particles );

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
