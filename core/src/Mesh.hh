//---------------------------------------------------------------------------//
/*!
 * \file Mesh.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_MESH_HH
#define EXAMPM_MESH_HH

#include "Particle.hh"

#include <array>
#include <vector>

namespace ExaMPM
{

//---------------------------------------------------------------------------//
/*!
 * \class Mesh
 *
 * \brief Uniform hexahedral mesh.
 */
class Mesh
{
  public:

    // Constructor.
    Mesh( const int num_cells_x,
          const int num_cells_y,
          const int num_cells_z,
          const double cell_width );

    // Get the spatial dimension of the mesh.
    int spatialDimension() const;

    // Get the number of nodes per cell in the mesh.
    int nodesPerCell() const;

    // Get the total number of cells in the mesh.
    int totalNumCells() const;

    // Get the total number of nodes in the mesh.
    int totalNumNodes() const;

    // Given a node id get its coordinates.
    void nodeCoordinates( const int node_id,
                          std::array<double,3>& coords ) const;

    // Given a boundary id get the ids of the nodes on that boundary.
    const std::vector<int>&
    getBoundaryNodes( const int boundary_id ) const;

    // Get the number of particles in a cell for a given order.
    int particlesPerCell( const int order ) const;

    // Given a cardinal cell id intitalize a number of particles in that cell.
    void initializeParticles(
        const int cell_id,
        const int order,
        std::vector<Particle>& particles ) const;

    // Given a particle determine the cardinal index of the cell in which it
    // is located.
    void locateParticle( const Particle& particle,
                         std::array<int,3>& cell_id ) const;

    // Given a cell ids get the node ids of the cell.
    void cellNodeIds( const std::array<int,3>& cell_id,
                      std::array<int,8>& nodes ) const;

    // Map the coordinates of a particle from the physical frame to the
    // reference frame of the cell in which it is located.
    void mapPhysicalToReferenceFrame(
        const Particle& particle,
        const std::array<int,3>& cell_id,
        std::array<double,3>& ref_coords ) const;

    // Given reference coordinates in a cell get the value of the shape
    // function at those coordinates.
    void shapeFunctionValue( const std::array<double,3>& ref_coords,
                             std::array<double,8>& values ) const;

    // Given reference coordinates in a cell get the gradient of the shape
    // function at those coordinates. Indexed as [Node][Dim].
    void shapeFunctionGradient(
        const std::array<double,3>& ref_coords,
        std::array<std::array<double,3>, 8>& gradients ) const;

  private:

    // Given ijk indices of a cell compute the cell id.
    int cellId( const int i, const int j, const int k ) const;

    // Given ijk indices of a node compute the node id.
    int nodeId( const int i, const int j, const int k ) const;

  private:

    // Number of cells in x direction.
    int d_num_cells_x;

    // Number of cells in y direction.
    int d_num_cells_y;

    // Number of cells in z direction.
    int d_num_cells_z;

    // Mesh cell width.
    double d_cell_width;

    // Number of nodes in the x direction.
    int d_num_nodes_x;

    // Number of nodes in the y direction.
    int d_num_nodes_y;

    // Number of nodes in the z direction.
    int d_num_nodes_z;

    // Boundary node ids.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    std::array<std::vector<int>,6> d_boundary_nodes;
};

} // end namespace ExaMPM

//---------------------------------------------------------------------------//

#endif // end EXAMPM_MESH_HH

//---------------------------------------------------------------------------//
// end Mesh.hh
//---------------------------------------------------------------------------//
