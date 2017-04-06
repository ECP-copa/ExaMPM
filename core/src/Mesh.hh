//---------------------------------------------------------------------------//
/*!
 * \file Mesh.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_MESH_HH
#define EXAMPM_MESH_HH

#include "Particle.hh"

#include <vector>

namespace ExaMPM
{

//---------------------------------------------------------------------------//
/*!
 * \class Mesh
 *
 * \brief Interface for a logically structured grid.
 */
class Mesh
{
  public:

    // Destructor.
    virtual ~Mesh() = default;

    // Get the spatial dimension of the mesh.
    virtual int spatialDimension() const = 0;

    // Get the number of nodes per cell in the mesh.
    virtual int nodesPerCell() const = 0;

    // Get the total number of cells in the mesh.
    virtual int totalNumCells() const = 0;

    // Get the total number of nodes in the mesh.
    virtual int totalNumNodes() const = 0;

    // Given a node id get its coordinates.
    virtual void nodeCoordinates( const int node_id,
                                  std::vector<double>& coords ) const = 0;

    // Given a boundary id get the ids of the nodes on that boundary. Boundary
    // ids should be indicates as (+/- x, +/- y, ... )
    virtual void getBoundaryNodes( const std::vector<int>& boundary_id,
                                   std::vector<int>& boundary_nodes ) const = 0;

    // Given a particle determine the cardinal index of the cell in which it
    // is located.
    virtual void locateParticle( const Particle& particle,
                                 std::vector<int>& cell_id ) const = 0;

    // Given a cell ids get the node ids of the cell.
    virtual void cellNodeIds( const std::vector<int>& cell_id,
                              std::vector<int>& nodes ) const = 0;

    // Map the coordinates of a particle from the physical frame to the
    // reference frame of the cell in which it is located.
    virtual void mapPhysicalToReferenceFrame(
        const Particle& particle,
        const std::vector<int>& cell_id,
        std::vector<double>& ref_coords ) const = 0;

    // Given reference coordinates in a cell get the value of the shape
    // function at those coordinates.
    virtual void shapeFunctionValue( const std::vector<double>& ref_coords,
                                     std::vector<double>& values ) const = 0;

    // Given reference coordinates in a cell get the gradient of the shape
    // function at those coordinates. Indexed as [Node][Dim].
    virtual void shapeFunctionGradient(
        const std::vector<double>& ref_coords,
        std::vector<std::vector<double> >& gradients ) const = 0;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_MESH_HH

//---------------------------------------------------------------------------//
// end Mesh.hh
//---------------------------------------------------------------------------//
