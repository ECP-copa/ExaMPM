//---------------------------------------------------------------------------//
/*!
 * \file Mesh2d.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_MESH2D_HH
#define EXAMPM_MESH2D_HH

#include "Particle.hh"
#include "Mesh.hh"

#include <vector>

namespace ExaMPM
{

//---------------------------------------------------------------------------//
/*!
 * \class Mesh2d
 *
 * \brief Uniform 2D quadrilateral mesh.
 */
class Mesh2d : public Mesh
{
  public:

    // Constructor.
    Mesh2d( const int num_cells_x,
            const int num_cells_y,
            const double cell_width );

    // Get the spatial dimension of the mesh.
    int spatialDimension() const override;

    // Get the number of nodes per cell in the mesh.
    int nodesPerCell() const override;

    // Get the total number of cells in the mesh.
    int totalNumCells() const override;

    // Get the total number of nodes in the mesh.
    int totalNumNodes() const override;

    // Given a node id get its coordinates.
    void nodeCoordinates( const int node_id,
                          std::vector<double>& coords ) const override;

    // Given a boundary id get the ids of the nodes on that boundary.
    void getBoundaryNodes( const std::vector<int>& boundary_id,
                           std::vector<int>& boundary_nodes ) const override;

    // Given a particle determine the cardinal index of the cell in which it
    // is located.
    void locateParticle( const Particle& particle,
                         std::vector<int>& cell_id ) const override;

    // Given a cell ids get the node ids of the cell.
    void cellNodeIds( const std::vector<int>& cell_id,
                      std::vector<int>& nodes ) const override;

    // Map the coordinates of a particle from the physical frame to the
    // reference frame of the cell in which it is located.
    void mapPhysicalToReferenceFrame(
        const Particle& particle,
        const std::vector<int>& cell_id,
        std::vector<double>& ref_coords ) const override;

    // Given reference coordinates in a cell get the value of the shape
    // function at those coordinates.
    void shapeFunctionValue( const std::vector<double>& ref_coords,
                             std::vector<double>& values ) const override;

    // Given reference coordinates in a cell get the gradient of the shape
    // function at those coordinates. Indexed as [Node][Dim].
    void shapeFunctionGradient(
        const std::vector<double>& ref_coords,
        std::vector<std::vector<double> >& gradients ) const override;

  private:

    // Given ij indices of a node compute the node id.
    int nodeId( const int i, const int j ) const;

  private:

    // Number of cells in x direction.
    int d_num_cells_x;

    // Number of cells in y direction.
    int d_num_cells_y;

    // Mesh cell width.
    double d_cell_width;
};

} // end namespace ExaMPM

//---------------------------------------------------------------------------//

#endif // end EXAMPM_MESH2D_HH

//---------------------------------------------------------------------------//
// end Mesh2d.hh
//---------------------------------------------------------------------------//
