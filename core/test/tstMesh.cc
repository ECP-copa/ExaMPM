//---------------------------------------------------------------------------//
/*!
 * \file   test/tstMesh.cc
 */
//---------------------------------------------------------------------------//

#include "Mesh.hh"
#include "Mesh2d.hh"
#include "Particle.hh"

#include "gtest_main.hh"

//---------------------------------------------------------------------------//
class MeshTest : public ::testing::Test
{
  protected:
    void SetUp()
    { /* ... */ }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST_F(MeshTest, mesh_test_2d)
{
    // Create a 2d mesh.
    int num_cells_x = 145;
    int num_cells_y = 103;
    int num_nodes_x = num_cells_x + 1;
    int num_nodes_y = num_cells_y + 1;
    double cell_width = 0.4;
    std::shared_ptr<ExaMPM::Mesh> mesh =
        std::make_shared<ExaMPM::Mesh2d>(
            num_cells_x, num_cells_y, cell_width );

    // Check the mesh.
    EXPECT_EQ( 2, mesh->spatialDimension() );
    EXPECT_EQ( 4, mesh->nodesPerCell() );
    EXPECT_EQ( num_cells_x*num_cells_y, mesh->totalNumCells() );
    EXPECT_EQ( num_nodes_x*num_nodes_y, mesh->totalNumNodes() );

    // Create a particle and assign a location.
    ExaMPM::Particle p;
    p.r = { 32.3, 22.9 };

    // Check the particle location.
    std::vector<int> cell_id( 2 );
    mesh->locateParticle( p, cell_id );
    EXPECT_EQ( cell_id[0], 80 );
    EXPECT_EQ( cell_id[1], 57 );

    // Check getting the nodes of a cell.
    int node_0 = cell_id[1] * num_nodes_x + cell_id[0];
    int node_1 = cell_id[1] * num_nodes_x + cell_id[0] + 1;
    int node_2 = (cell_id[1]+1) * num_nodes_x + cell_id[0] + 1;
    int node_3 = (cell_id[1]+1) * num_nodes_x + cell_id[0];
    std::vector<int> cell_nodes( 4 );
    mesh->cellNodeIds( cell_id, cell_nodes );
    EXPECT_EQ( node_0, cell_nodes[0] );
    EXPECT_EQ( node_1, cell_nodes[1] );
    EXPECT_EQ( node_2, cell_nodes[2] );
    EXPECT_EQ( node_3, cell_nodes[3] );

    // Check the coordinates of a node.
    std::vector<double> node_coords(2);
    mesh->nodeCoordinates( cell_nodes[3], node_coords );
    EXPECT_FLOAT_EQ( node_coords[0], 32.0 );
    EXPECT_FLOAT_EQ( node_coords[1], 23.2 );

    // Check mapping to the reference frame.
    std::vector<double> ref_coords(2);
    mesh->mapPhysicalToReferenceFrame( p, cell_id, ref_coords );
    EXPECT_FLOAT_EQ( ref_coords[0], 0.5 );
    EXPECT_FLOAT_EQ( ref_coords[1], -0.5 );

    // Check the shape function value.
    std::vector<double> svals( 4 );
    mesh->shapeFunctionValue( ref_coords, svals );
    EXPECT_FLOAT_EQ( svals[0], 0.1875 );
    EXPECT_FLOAT_EQ( svals[1], 0.5625 );
    EXPECT_FLOAT_EQ( svals[2], 0.1875 );
    EXPECT_FLOAT_EQ( svals[3], 0.0625 );

    // Check the shape function gradient.
    std::vector<std::vector<double> > sgrads( 4, std::vector<double>(2) );
    mesh->shapeFunctionGradient( ref_coords, sgrads );

    EXPECT_FLOAT_EQ( sgrads[0][0], -0.375 );
    EXPECT_FLOAT_EQ( sgrads[0][1], -0.125 );

    EXPECT_FLOAT_EQ( sgrads[1][0], 0.375 );
    EXPECT_FLOAT_EQ( sgrads[1][1], -0.375 );

    EXPECT_FLOAT_EQ( sgrads[2][0], 0.125 );
    EXPECT_FLOAT_EQ( sgrads[2][1], 0.375 );

    EXPECT_FLOAT_EQ( sgrads[3][0], -0.125 );
    EXPECT_FLOAT_EQ( sgrads[3][1], 0.125 );

    // Check the boundaries.
    std::vector<int> boundary_id( 2 );
    std::vector<int> boundary_nodes;
    auto getNodeId = [=]( const int i, const int j )
                     { return j * num_nodes_x + i; };

    // lo x
    boundary_id[0] = -1;
    boundary_id[1] = 0;
    mesh->getBoundaryNodes( boundary_id, boundary_nodes );
    EXPECT_EQ( num_nodes_x, boundary_nodes.size() );
    for ( int i = 0; i < num_nodes_x; ++i )
        EXPECT_EQ( boundary_nodes[i], getNodeId(i,0) );

    // lo y
    boundary_id[0] = 0;
    boundary_id[1] = -1;
    mesh->getBoundaryNodes( boundary_id, boundary_nodes );
    EXPECT_EQ( num_nodes_y, boundary_nodes.size() );
    for ( int j = 0; j < num_nodes_y; ++j )
        EXPECT_EQ( boundary_nodes[j], getNodeId(0,j) );

    // hi x
    boundary_id[0] = 1;
    boundary_id[1] = 0;
    mesh->getBoundaryNodes( boundary_id, boundary_nodes );
    EXPECT_EQ( num_nodes_x, boundary_nodes.size() );
    for ( int i = 0; i < num_nodes_x; ++i )
        EXPECT_EQ( boundary_nodes[i], getNodeId(i,num_cells_y) );

    // hi y
    boundary_id[0] = 0;
    boundary_id[1] = 1;
    mesh->getBoundaryNodes( boundary_id, boundary_nodes );
    EXPECT_EQ( num_nodes_y, boundary_nodes.size() );
    for ( int j = 0; j < num_nodes_y; ++j )
        EXPECT_EQ( boundary_nodes[j], getNodeId(num_cells_x,j) );
}

//---------------------------------------------------------------------------//
// end of test/tstMesh.cc
//---------------------------------------------------------------------------//
