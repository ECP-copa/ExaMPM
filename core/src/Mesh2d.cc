//---------------------------------------------------------------------------//
/*!
 * \file Mesh2d.cc
 */
//---------------------------------------------------------------------------//

#include "Mesh2d.hh"
#include <cassert>
#include <cmath>

namespace ExaMPM
{

//---------------------------------------------------------------------------//
// Constructor.
Mesh2d::Mesh2d( const int num_cells_x,
                const int num_cells_y,
                const double cell_width )
    : d_num_cells_x( num_cells_x )
    , d_num_cells_y( num_cells_y )
    , d_cell_width( cell_width )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Get the spatial dimension of the mesh.
int Mesh2d::spatialDimension() const
{
    return 2;
}

//---------------------------------------------------------------------------//
// Get the number of nodes per cell in the mesh.
int Mesh2d::nodesPerCell() const
{
    // Mesh consists of linear quadrilaterals.
    return 4;
}

//---------------------------------------------------------------------------//
// Get the total number of cells in the mesh.
int Mesh2d::totalNumCells() const
{
    return d_num_cells_x * d_num_cells_y;
}

//---------------------------------------------------------------------------//
// Get the total number of nodes in the mesh.
int Mesh2d::totalNumNodes() const
{
    return (d_num_cells_x+1) * (d_num_cells_y+1);
}

//---------------------------------------------------------------------------//
// Given a node id get its coordinates.
void Mesh2d::nodeCoordinates( const int node_id,
                              std::vector<double>& coords ) const
{
    assert( node_id < totalNumNodes() );
    assert( 2 == coords.size() );

    // Get the ij indices of the node.
    int node_j = std::floor( node_id / (d_num_cells_x+1) );
    int node_i = node_id - node_j * (d_num_cells_x+1);
    assert( node_id == node_j*(d_num_cells_x+1) + node_i );

    // Get the coordinates.
    coords[0] = node_i * d_cell_width;
    coords[1] = node_j * d_cell_width;
}

//---------------------------------------------------------------------------//
// Given a boundary id get the ids of the nodes on that boundary.
void Mesh2d::getBoundaryNodes( const std::vector<int>& boundary_id,
                               std::vector<int>& boundary_nodes ) const
{
    assert( 2 == boundary_id.size() );

    // X boundary.
    if ( boundary_id[0] )
    {
        // Make sure a Y boundary was also not input.
        assert( !boundary_id[1] );

        // Size the boundary nodes input vector.
        boundary_nodes.resize( d_num_cells_x + 1 );

        // Y index.
        int j = 0;

        // Low X boundary.
        if ( boundary_id[0] < 0 )
        {
            j = 0;
        }

        // High X boundary.
        else if ( boundary_id[0] > 0 )
        {
            j = d_num_cells_y;
        }

        // Fill the boundary nodes.
        for ( int i = 0; i < d_num_cells_x + 1; ++i )
            boundary_nodes[i] = nodeId( i, j );

    }

    // Y Boundary
    else if ( boundary_id[1] )
    {
        // Make sure an X boundary was also not input.
        assert( !boundary_id[0] );

        // Size the boundary nodes input vector.
        boundary_nodes.resize( d_num_cells_y + 1 );

        // X index.
        int i = 0;

        // Low Y boundary.
        if ( boundary_id[1] < 0 )
        {
            i = 0;
        }

        // High Y boundary.
        else if ( boundary_id[1] > 0 )
        {
            i = d_num_cells_x;
        }

        // Fill the boundary nodes.
        for ( int j = 0; j < d_num_cells_y + 1; ++j )
            boundary_nodes[j] = nodeId( i, j );
    }
}

//---------------------------------------------------------------------------//
// Given a particle determine the cardinal index of the cell in which it
// is located.
void Mesh2d::locateParticle( const Particle& particle,
                             std::vector<int>& cell_id ) const
{
    // Mesh spans (0.0,d_num_cells_x*cell_width) in x and
    // (0.0,d_num_cells_y*cell_width) in y
    cell_id[0] = std::floor( particle.r[0] / d_cell_width );
    cell_id[1] = std::floor( particle.r[1] / d_cell_width );

    assert( cell_id[0] < d_num_cells_x );
    assert( cell_id[1] < d_num_cells_y );
}

//---------------------------------------------------------------------------//
// Given a cell ids get the node ids of the cell.
void Mesh2d::cellNodeIds( const std::vector<int>& cell_id,
                          std::vector<int>& nodes ) const
{
    assert( 4 == nodes.size() );

    // Nodes are given in canonical order for a linear quadrilateral:
    //
    //   3 (i,j+1)---------(i+i,j+1) 2
    //        |                |
    //        |                |
    //        |                |
    //        |                |
    //    0 (i,j)------------(i+1,j) 1

    nodes[0] = nodeId( cell_id[0],   cell_id[1] );
    nodes[1] = nodeId( cell_id[0]+1, cell_id[1] );
    nodes[2] = nodeId( cell_id[0]+1, cell_id[1]+1 );
    nodes[3] = nodeId( cell_id[0],   cell_id[1]+1 );
}

//---------------------------------------------------------------------------//
// Map the coordinates of a particle from the physical frame to the
// reference frame of the cell in which it is located.
void Mesh2d::mapPhysicalToReferenceFrame(
    const Particle& particle,
    const std::vector<int>& cell_id,
    std::vector<double>& ref_coords ) const
{
    assert( 2 == ref_coords.size() );
    assert( cell_id[0] == std::floor( particle.r[0] / d_cell_width ) );
    assert( cell_id[1] == std::floor( particle.r[1] / d_cell_width ) );


    // The reference cell spans -1 to 1 in both directions.
    //
    //   3 (-1,1)------------(1,1) 2
    //        |                |
    //        |                |
    //        |                |
    //        |                |
    //    0 (-1,-1)----------(1,-1) 1
    ref_coords[0] = (particle.r[0]/d_cell_width - cell_id[0])*2.0 - 1.0;
    ref_coords[1] = (particle.r[1]/d_cell_width - cell_id[1])*2.0 - 1.0;

    assert( -1.0 <= ref_coords[0] && ref_coords[0] <= 1.0 );
    assert( -1.0 <= ref_coords[1] && ref_coords[1] <= 1.0 );
}


//---------------------------------------------------------------------------//
// Given reference coordinates in a cell get the value of the shape
// function at those coordinates.
void Mesh2d::shapeFunctionValue( const std::vector<double>& ref_coords,
                                 std::vector<double>& values ) const
{
    assert( 2 == ref_coords.size() );
    assert( 4 == values.size() );

    values[0] = (1.0 - ref_coords[0])*(1.0 - ref_coords[1])/4.0;
    values[1] = (1.0 + ref_coords[0])*(1.0 - ref_coords[1])/4.0;
    values[2] = (1.0 + ref_coords[0])*(1.0 + ref_coords[1])/4.0;
    values[3] = (1.0 - ref_coords[0])*(1.0 + ref_coords[1])/4.0;
}

//---------------------------------------------------------------------------//
// Given reference coordinates in a cell get the gradient of the shape
// function at those coordinates. Indexed as [Node][Dim].
void Mesh2d::shapeFunctionGradient(
    const std::vector<double>& ref_coords,
    std::vector<std::vector<double> >& gradients ) const
{
    assert( 2 == ref_coords.size() );
    assert( 4 == gradients.size() );
    assert( 2 == gradients[0].size() );
    assert( 2 == gradients[1].size() );
    assert( 2 == gradients[2].size() );
    assert( 2 == gradients[3].size() );

    gradients[0][0] = -(1.0 - ref_coords[1])/4.0;
    gradients[0][1] = -(1.0 - ref_coords[0])/4.0;

    gradients[1][0] =  (1.0 - ref_coords[1])/4.0;
    gradients[1][1] = -(1.0 + ref_coords[0])/4.0;

    gradients[2][0] =  (1.0 + ref_coords[1])/4.0;
    gradients[2][1] =  (1.0 + ref_coords[0])/4.0;

    gradients[3][0] = -(1.0 + ref_coords[1])/4.0;
    gradients[3][1] =  (1.0 - ref_coords[0])/4.0;
}

//---------------------------------------------------------------------------//
// Given ij indices of a node compute the node id.
int Mesh2d::nodeId( const int i, const int j ) const
{
    return j*(d_num_cells_x+1) + i;
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end Mesh2d.hh
//---------------------------------------------------------------------------//
