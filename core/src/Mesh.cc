//---------------------------------------------------------------------------//
/*!
 * \file Mesh.cc
 */
//---------------------------------------------------------------------------//

#include "Mesh.hh"
#include <cassert>
#include <cmath>

namespace ExaMPM
{

//---------------------------------------------------------------------------//
// Constructor.
Mesh::Mesh( const int num_cells_x,
            const int num_cells_y,
            const int num_cells_z,
            const double cell_width )
    : d_num_cells_x( num_cells_x )
    , d_num_cells_y( num_cells_y )
    , d_num_cells_z( num_cells_z )
    , d_cell_width( cell_width )
    , d_num_nodes_x( num_cells_x + 1 )
    , d_num_nodes_y( num_cells_y + 1 )
    , d_num_nodes_z( num_cells_z + 1 )
{
    // Create the boundary nodes.

    // -X boundary.
    {
        // Boundary id.
        int bid = 0;

        // X index.
        int i = 0;

        // Fill the boundary nodes.
        d_boundary_nodes[bid].resize( d_num_nodes_y*d_num_nodes_z );
        for ( int k = 0; k < d_num_nodes_z; ++k )
            for ( int j = 0; j < d_num_nodes_y; ++j )
                d_boundary_nodes[bid][k*d_num_nodes_y+j] = nodeId( i, j, k );
    }

    // +X boundary.
    {
        // Boundary id.
        int bid = 1;

        // X index
        int i = d_num_cells_x;

        // Fill the boundary nodes.
        d_boundary_nodes[bid].resize( d_num_nodes_y*d_num_nodes_z );
        for ( int k = 0; k < d_num_nodes_z; ++k )
            for ( int j = 0; j < d_num_nodes_y; ++j )
                d_boundary_nodes[bid][k*d_num_nodes_y+j] = nodeId( i, j, k );
    }


    // -Y boundary
    {
        // Boundary id.
        int bid = 2;

        // Y index.
        int j = 0;

        // Fill the boundary nodes.
        d_boundary_nodes[bid].resize( d_num_nodes_x*d_num_nodes_z );
        for ( int k = 0; k < d_num_nodes_z; ++k )
            for ( int i = 0; i < d_num_nodes_x; ++i )
                d_boundary_nodes[bid][k*d_num_nodes_x+i] = nodeId( i, j, k );

    }

    // +Y boundary
    {
        // Boundary id.
        int bid = 3;

        // Y index.
        int j = d_num_cells_y;

        // Fill the boundary nodes.
        d_boundary_nodes[bid].resize( d_num_nodes_x*d_num_nodes_z );
        for ( int k = 0; k < d_num_nodes_z; ++k )
            for ( int i = 0; i < d_num_nodes_x; ++i )
                d_boundary_nodes[bid][k*d_num_nodes_x+i] = nodeId( i, j, k );

    }

    // -Z boundary
    {
        // Boundary id.
        int bid = 4;

        // Z index.
        int k = 0;

        // Fill the boundary nodes.
        d_boundary_nodes[bid].resize( d_num_nodes_x*d_num_nodes_y );
        for ( int j = 0; j< d_num_nodes_y; ++j )
            for ( int i = 0; i < d_num_nodes_x; ++i )
                d_boundary_nodes[bid][j*d_num_nodes_x+i] = nodeId( i, j, k );
    }

    // +Z boundary
    {
        // Boundary id.
        int bid = 5;

        // Z index.
        int k = d_num_cells_z;

        // Fill the boundary nodes.
        d_boundary_nodes[bid].resize( d_num_nodes_x*d_num_nodes_y );
        for ( int j = 0; j< d_num_nodes_y; ++j )
            for ( int i = 0; i < d_num_nodes_x; ++i )
                d_boundary_nodes[bid][j*d_num_nodes_x+i] = nodeId( i, j, k );
    }
}

//---------------------------------------------------------------------------//
// Get the spatial dimension of the mesh.
int Mesh::spatialDimension() const
{
    return 3;
}

//---------------------------------------------------------------------------//
// Get the number of nodes per cell in the mesh.
int Mesh::nodesPerCell() const
{
    // Mesh consists of linear hex elements.
    return 8;
}

//---------------------------------------------------------------------------//
// Get the total number of cells in the mesh.
int Mesh::totalNumCells() const
{
    return d_num_cells_x * d_num_cells_y * d_num_cells_z;
}

//---------------------------------------------------------------------------//
// Get the total number of nodes in the mesh.
int Mesh::totalNumNodes() const
{
    return d_num_nodes_x * d_num_nodes_y * d_num_nodes_z;
}

//---------------------------------------------------------------------------//
// Given a node id get its coordinates.
void Mesh::nodeCoordinates( const int node_id,
                            std::array<double,3>& coords ) const
{
    assert( node_id < totalNumNodes() );

    // Get the ijk indices of the node.
    int nk = std::floor( node_id / (d_num_nodes_x*d_num_nodes_y) );
    int nj = std::floor(
        (node_id - d_num_nodes_x*d_num_nodes_y*nk) / d_num_nodes_x );
    int ni = node_id - d_num_nodes_x*d_num_nodes_y*nk -
                 d_num_nodes_x * nj;
    assert( node_id ==
            d_num_nodes_x*d_num_nodes_y*nk + d_num_nodes_x*nj + ni );

    // Get the coordinates.
    coords[0] = ni * d_cell_width;
    coords[1] = nj * d_cell_width;
    coords[2] = nk * d_cell_width;
}

//---------------------------------------------------------------------------//
// Given a boundary id get the ids of the nodes on that boundary.
const std::vector<int>&
Mesh::getBoundaryNodes( const int boundary_id ) const
{
    assert( 0 <= boundary_id && boundary_id < 6 );
    return d_boundary_nodes[boundary_id];
}

//---------------------------------------------------------------------------//
// Get the number of particles in a cell for a given order.
int Mesh::particlesPerCell( const int order ) const
{
    return order*order*order;
}

//---------------------------------------------------------------------------//
// Given a cardinal cell id intitalize a number of particles in that cell.
void Mesh::initializeParticles(
    const int cell_id,
    const int order,
    std::vector<Particle>& particles ) const
{
    assert( particles.size() == order*order*order );

    // Get the ijk cell indices.
    int ck = std::floor( cell_id / (d_num_cells_x*d_num_cells_y) );
    int cj = std::floor(
        (cell_id - d_num_cells_x*d_num_cells_y*ck) / d_num_cells_x );
    int ci = cell_id - d_num_cells_x*d_num_cells_y*ck -
             d_num_cells_x * cj;
    assert( cell_id ==
            d_num_cells_x*d_num_cells_y*ck + d_num_cells_x*cj + ci );

    // Calculate particle spacing. We want to space them evenly and away from
    // the edges of the cells.
    //
    // Order 1:
    //
    //        +----------------+
    //        |                |
    //        |                |
    //        |       o        |
    //        |                |
    //        |                |
    //        +----------------+
    //
    // Order 2:
    //
    //        +----------------+
    //        |                |
    //        |   o       o    |
    //        |                |
    //        |   o       o    |
    //        |                |
    //        +----------------+
    //
    // Order 3:
    //
    //        +----------------+
    //        | o     o     o  |
    //        |                |
    //        | o     o     o  |
    //        |                |
    //        | o     o     o  |
    //        +----------------+
    //
    // etc....


    double dp = d_cell_width / order;

    // Calculate particle volume.
    double cell_volume = d_cell_width * d_cell_width * d_cell_width;
    double p_vol = cell_volume / (order*order*order);

    // Initialize the particle locations.
    double x0 = ci*d_cell_width + dp / 2.0;
    double y0 = cj*d_cell_width + dp / 2.0;
    double z0 = ck*d_cell_width + dp / 2.0;
    int lid = 0;
    for ( int k = 0; k < order; ++k )
    {
        for ( int j = 0; j < order; ++j )
        {
            for ( int i = 0; i < order; ++i )
            {
                // Get the local particle id.
                lid = k*order*order + j*order + i;
                assert( lid < particles.size() );

                // Set the position.
                particles[lid].r[0] = x0 + i*dp;
                particles[lid].r[1] = y0 + j*dp;
                particles[lid].r[2] = z0 + k*dp;

                // Set the volume.
                particles[lid].volume = p_vol;
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Given a particle determine the cardinal index of the cell in which it
// is located.
void Mesh::locateParticle( const Particle& particle,
                           std::array<int,3>& cell_id ) const
{
    // Mesh spans (0.0,d_num_cells_x*cell_width) in x and
    // (0.0,d_num_cells_y*cell_width) in y
    cell_id[0] = std::floor( particle.r[0] / d_cell_width );
    cell_id[1] = std::floor( particle.r[1] / d_cell_width );
    cell_id[2] = std::floor( particle.r[2] / d_cell_width );

    assert( 0 <= cell_id[0] );
    assert( 0 <= cell_id[1] );
    assert( 0 <= cell_id[2] );
    assert( cell_id[0] < d_num_cells_x );
    assert( cell_id[1] < d_num_cells_y );
    assert( cell_id[2] < d_num_cells_z );
}

//---------------------------------------------------------------------------//
// Given a cell ids get the node ids of the cell.
void Mesh::cellNodeIds( const std::array<int,3>& cell_id,
                        std::array<int,8>& nodes ) const
{
    // Nodes are given in canonical order for a linear hex:
    //
    //   3 (i,j+1)---------(i+i,j+1) 2
    //        |                |
    //        |                |
    //        |                |
    //        |                |
    //    0 (i,j)------------(i+1,j) 1

    nodes[0] = nodeId( cell_id[0],     cell_id[1],   cell_id[2] );
    nodes[1] = nodeId( cell_id[0]+1,   cell_id[1],   cell_id[2] );
    nodes[2] = nodeId( cell_id[0]+1, cell_id[1]+1,   cell_id[2] );
    nodes[3] = nodeId( cell_id[0],   cell_id[1]+1,   cell_id[2] );

    nodes[4] = nodeId( cell_id[0],     cell_id[1], cell_id[2]+1 );
    nodes[5] = nodeId( cell_id[0]+1,   cell_id[1], cell_id[2]+1 );
    nodes[6] = nodeId( cell_id[0]+1, cell_id[1]+1, cell_id[2]+1 );
    nodes[7] = nodeId( cell_id[0],   cell_id[1]+1, cell_id[2]+1 );
}

//---------------------------------------------------------------------------//
// Map the coordinates of a particle from the physical frame to the
// reference frame of the cell in which it is located.
void Mesh::mapPhysicalToReferenceFrame(
    const Particle& particle,
    const std::array<int,3>& cell_id,
    std::array<double,3>& ref_coords ) const
{
    assert( cell_id[0] == std::floor( particle.r[0] / d_cell_width ) );
    assert( cell_id[1] == std::floor( particle.r[1] / d_cell_width ) );
    assert( cell_id[2] == std::floor( particle.r[2] / d_cell_width ) );

    // The reference cell spans -1 to 1 in all directions.
    //
    //   3 (-1,1)------------(1,1) 2
    //        |                |
    //        |                |
    //        |                |
    //        |                |
    //    0 (-1,-1)----------(1,-1) 1
    ref_coords[0] = (particle.r[0]/d_cell_width - cell_id[0])*2.0 - 1.0;
    ref_coords[1] = (particle.r[1]/d_cell_width - cell_id[1])*2.0 - 1.0;
    ref_coords[2] = (particle.r[2]/d_cell_width - cell_id[2])*2.0 - 1.0;

    assert( -1.0 <= ref_coords[0] && ref_coords[0] <= 1.0 );
    assert( -1.0 <= ref_coords[1] && ref_coords[1] <= 1.0 );
    assert( -1.0 <= ref_coords[2] && ref_coords[2] <= 1.0 );
}


//---------------------------------------------------------------------------//
// Given reference coordinates in a cell get the value of the shape
// function at those coordinates.
void Mesh::shapeFunctionValue( const std::array<double,3>& ref_coords,
                               std::array<double,8>& values ) const
{
    values[0] =
        (1.0 - ref_coords[0])*(1.0 - ref_coords[1])*(1.0 - ref_coords[2])/8.0;
    values[1] =
        (1.0 + ref_coords[0])*(1.0 - ref_coords[1])*(1.0 - ref_coords[2])/8.0;
    values[2] =
        (1.0 + ref_coords[0])*(1.0 + ref_coords[1])*(1.0 - ref_coords[2])/8.0;
    values[3] =
        (1.0 - ref_coords[0])*(1.0 + ref_coords[1])*(1.0 - ref_coords[2])/8.0;

    values[4] =
        (1.0 - ref_coords[0])*(1.0 - ref_coords[1])*(1.0 + ref_coords[2])/8.0;
    values[5] =
        (1.0 + ref_coords[0])*(1.0 - ref_coords[1])*(1.0 + ref_coords[2])/8.0;
    values[6] =
        (1.0 + ref_coords[0])*(1.0 + ref_coords[1])*(1.0 + ref_coords[2])/8.0;
    values[7] =
        (1.0 - ref_coords[0])*(1.0 + ref_coords[1])*(1.0 + ref_coords[2])/8.0;
}

//---------------------------------------------------------------------------//
// Given reference coordinates in a cell get the gradient of the shape
// function at those coordinates. Indexed as [Node][Dim].
void Mesh::shapeFunctionGradient(
    const std::array<double,3>& ref_coords,
    std::array<std::array<double,3>, 8>& gradients ) const
{
        gradients[0][0] = -(1.0 - ref_coords[1])*(1.0 - ref_coords[2])/8.0;
        gradients[0][1] = -(1.0 - ref_coords[0])*(1.0 - ref_coords[2])/8.0;
        gradients[0][2] = -(1.0 - ref_coords[0])*(1.0 - ref_coords[1])/8.0;

        gradients[1][0] =  (1.0 - ref_coords[1])*(1.0 - ref_coords[2])/8.0;
        gradients[1][1] = -(1.0 + ref_coords[0])*(1.0 - ref_coords[2])/8.0;
        gradients[1][2] = -(1.0 + ref_coords[0])*(1.0 - ref_coords[1])/8.0;

        gradients[2][0] =  (1.0 + ref_coords[1])*(1.0 - ref_coords[2])/8.0;
        gradients[2][1] =  (1.0 + ref_coords[0])*(1.0 - ref_coords[2])/8.0;
        gradients[2][2] = -(1.0 + ref_coords[0])*(1.0 + ref_coords[1])/8.0;

        gradients[3][0] = -(1.0 + ref_coords[1])*(1.0 - ref_coords[2])/8.0;
        gradients[3][1] =  (1.0 - ref_coords[0])*(1.0 - ref_coords[2])/8.0;
        gradients[3][2] = -(1.0 - ref_coords[0])*(1.0 + ref_coords[1])/8.0;

        gradients[4][0] = -(1.0 - ref_coords[1])*(1.0 + ref_coords[2])/8.0;
        gradients[4][1] = -(1.0 - ref_coords[0])*(1.0 + ref_coords[2])/8.0;
        gradients[4][2] =  (1.0 - ref_coords[0])*(1.0 - ref_coords[1])/8.0;

        gradients[5][0] =  (1.0 - ref_coords[1])*(1.0 + ref_coords[2])/8.0;
        gradients[5][1] = -(1.0 + ref_coords[0])*(1.0 + ref_coords[2])/8.0;
        gradients[5][2] =  (1.0 + ref_coords[0])*(1.0 - ref_coords[1])/8.0;

        gradients[6][0] =  (1.0 + ref_coords[1])*(1.0 + ref_coords[2])/8.0;
        gradients[6][1] =  (1.0 + ref_coords[0])*(1.0 + ref_coords[2])/8.0;
        gradients[6][2] =  (1.0 + ref_coords[0])*(1.0 + ref_coords[1])/8.0;

        gradients[7][0] = -(1.0 + ref_coords[1])*(1.0 + ref_coords[2])/8.0;
        gradients[7][1] =  (1.0 - ref_coords[0])*(1.0 + ref_coords[2])/8.0;
        gradients[7][2] =  (1.0 - ref_coords[0])*(1.0 + ref_coords[1])/8.0;
}

//---------------------------------------------------------------------------//
// Given ij indices of a node compute the node id.
int Mesh::nodeId( const int i, const int j, const int k ) const
{
    return d_num_nodes_x*d_num_nodes_y*k + d_num_nodes_x*j + i;
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end Mesh.cc
//---------------------------------------------------------------------------//
