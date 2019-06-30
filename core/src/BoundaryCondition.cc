//---------------------------------------------------------------------------//
/*!
 * \file BoundaryCondition.hh
 */
//---------------------------------------------------------------------------//

#include "BoundaryCondition.hh"

#include <cmath>
#include <cassert>

namespace ExaMPM
{

//---------------------------------------------------------------------------//
// No slip condition
//---------------------------------------------------------------------------//
// Evaluate the boundary condition for momentum.
// (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
void NoSlipBoundaryCondition::evaluateMomentumCondition(
    const std::shared_ptr<Mesh>& mesh,
    const int boundary,
    const std::vector<double>& mass,
    std::vector<std::array<double,3> >& momentum ) const
{
    // Get the boundary nodes.
    const auto& nodes = mesh->getBoundaryNodes( boundary );

    // Apply the boundary condition. No slip means the velocity and therefore
    // the momentum in all directions is zero.
    for ( auto n : nodes )
        std::fill( momentum[n].begin(), momentum[n].end(), 0.0 );
}

//---------------------------------------------------------------------------//
// Evaluate the boundary condition for impulse.
// (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
void NoSlipBoundaryCondition::evaluateImpulseCondition(
    const std::shared_ptr<Mesh>& mesh,
    const int boundary,
    const std::vector<double>& mass,
    std::vector<std::array<double,3> >& impulse ) const
{
    // Get the boundary nodes.
    const auto& nodes = mesh->getBoundaryNodes( boundary );

    // Apply the boundary condition. No slip means the velocity and therefore
    // the impulse in all directions is zero.
    for ( auto n : nodes )
        std::fill( impulse[n].begin(), impulse[n].end(), 0.0 );
}

//---------------------------------------------------------------------------//
// Evaluate the new particle positions on the boundary.
// Does not do anything unless periodic boundary conditions.
void NoSlipBoundaryCondition::evaluateBoundaryPosition(
    const std::shared_ptr<Mesh>& mesh, const int boundary,
    std::vector<Particle>& particles ) const {}

//---------------------------------------------------------------------------//
//
// Does not do anything unless periodic boundary conditions.
void NoSlipBoundaryCondition::completeBoundarySum(
    const std::shared_ptr<Mesh>& mesh, const int boundary,
    std::vector<double>& field ) const {}

//---------------------------------------------------------------------------//
// 
// Does not do anything unless periodic boundary conditions.
void NoSlipBoundaryCondition::completeBoundarySum(
    const std::shared_ptr<Mesh>& mesh, const int boundary,
    std::vector<std::array<double,3> >& field ) const {}

//---------------------------------------------------------------------------//
// Free slip condition
//---------------------------------------------------------------------------//
// Evaluate the boundary condition for momentum.
// (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
void FreeSlipBoundaryCondition::evaluateMomentumCondition(
    const std::shared_ptr<Mesh>& mesh,
    const int boundary,
    const std::vector<double>& mass,
    std::vector<std::array<double,3> >& momentum ) const
{
    // Determine if the boundary is on x, y, or z.
    int dim = std::floor( boundary / 2 );

    // Determine if it is positive or negative.
    int dir = (boundary % 2) ? 1 : -1;

    // Get the boundary nodes.
    const auto& nodes = mesh->getBoundaryNodes( boundary );

    // Apply the boundary condition. Free slip means the velocity component
    // normal to the grid boundary is zero and therefore so it the momentum.
    for ( auto n : nodes )
        momentum[n][dim] = 0.0;
}

//---------------------------------------------------------------------------//
// Evaluate the boundary condition for impulse.
// (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
void FreeSlipBoundaryCondition::evaluateImpulseCondition(
    const std::shared_ptr<Mesh>& mesh,
    const int boundary,
    const std::vector<double>& mass,
    std::vector<std::array<double,3> >& impulse ) const
{
    // Determine if the boundary is on x, y, or z.
    int dim = std::floor( boundary / 2 );

    // Determine if it is positive or negative.
    int dir = (boundary % 2) ? 1 : -1;

    // Get the boundary nodes.
    const auto& nodes = mesh->getBoundaryNodes( boundary );

    // Apply the boundary condition. Free slip means the velocity component
    // normal to the grid boundary is zero and therefore so is the impulse.
    for ( auto n : nodes )
        impulse[n][dim] = 0.0;
}

//---------------------------------------------------------------------------//
// Evaluate the new particle positions on the boundary.
// Does not do anything unless periodic boundary conditions.
void FreeSlipBoundaryCondition::evaluateBoundaryPosition(
    const std::shared_ptr<Mesh>& mesh, const int boundary,
    std::vector<Particle>& particles ) const {}

//---------------------------------------------------------------------------//
//
// Does not do anything unless periodic boundary conditions.
void FreeSlipBoundaryCondition::completeBoundarySum(
    const std::shared_ptr<Mesh>& mesh, const int boundary,
    std::vector<double>& field ) const {}

//---------------------------------------------------------------------------//
// 
// Does not do anything unless periodic boundary conditions.
void FreeSlipBoundaryCondition::completeBoundarySum(
    const std::shared_ptr<Mesh>& mesh, const int boundary,
    std::vector<std::array<double,3> >& field ) const {}

//---------------------------------------------------------------------------//
// Velocity condition
//---------------------------------------------------------------------------//
// Constructor.
VelocityBoundaryCondition::VelocityBoundaryCondition(
    const std::array<double,3>& vb )
    : d_vb( vb )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Evaluate the boundary condition for momentum.
// (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
void VelocityBoundaryCondition::evaluateMomentumCondition(
    const std::shared_ptr<Mesh>& mesh,
    const int boundary,
    const std::vector<double>& mass,
    std::vector<std::array<double,3> >& momentum ) const
{
    // Ensure that the normal velocity component is zero.
    assert( d_vb[std::floor(boundary/2)] == 0.0 );

    // Get the boundary nodes.
    const auto& nodes = mesh->getBoundaryNodes( boundary );

    // Apply the boundary condition by setting the prescribed velocity times
    // the mass to give the new momentum.
    for ( auto n : nodes )
        for ( int d = 0; d < 3; ++d )
            momentum[n][d] = d_vb[d] * mass[n];
}

//---------------------------------------------------------------------------//
// Evaluate the boundary condition for impulse.
// (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
void VelocityBoundaryCondition::evaluateImpulseCondition(
    const std::shared_ptr<Mesh>& mesh,
    const int boundary,
    const std::vector<double>& mass,
    std::vector<std::array<double,3> >& impulse ) const
{
    // Determine if the boundary is on x, y, or z.
    int dim = std::floor( boundary / 2 );

    // Get the boundary nodes.
    const auto& nodes = mesh->getBoundaryNodes( boundary );

    // Apply the boundary condition. We prescribed a velocity tangential to
    // the surface and force is constant with time so there is no impulse. The
    // mass must flow freely across the surface so this is effectively a free
    // slip condition.
    for ( auto n : nodes )
        impulse[n][dim] = 0.0;
}

//---------------------------------------------------------------------------//
// Evaluate the new particle positions on the boundary.
// Does not do anything unless periodic boundary conditions.
void VelocityBoundaryCondition::evaluateBoundaryPosition(
    const std::shared_ptr<Mesh>& mesh, const int boundary,
    std::vector<Particle>& particles ) const {}

//---------------------------------------------------------------------------//
//
// Does not do anything unless periodic boundary conditions.
void VelocityBoundaryCondition::completeBoundarySum(
    const std::shared_ptr<Mesh>& mesh, const int boundary,
    std::vector<double>& field ) const {}

//---------------------------------------------------------------------------//
// 
// Does not do anything unless periodic boundary conditions.
void VelocityBoundaryCondition::completeBoundarySum(
    const std::shared_ptr<Mesh>& mesh, const int boundary,
    std::vector<std::array<double,3> >& field ) const {}

//---------------------------------------------------------------------------//
// Periodic condition
//---------------------------------------------------------------------------//
// Evaluate the boundary condition for momentum.
// (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
void PeriodicBoundaryCondition::evaluateMomentumCondition(
    const std::shared_ptr<Mesh>& mesh,
    const int boundary,
    const std::vector<double>& mass,
    std::vector<std::array<double,3> >& momentum ) const
{
}

//---------------------------------------------------------------------------//
// Evaluate the boundary condition for impulse.
// (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
void PeriodicBoundaryCondition::evaluateImpulseCondition(
    const std::shared_ptr<Mesh>& mesh,
    const int boundary,
    const std::vector<double>& mass,
    std::vector<std::array<double,3> >& impulse ) const
{
}

//---------------------------------------------------------------------------//
// Evaluate the new particle positions on the boundary.
// (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
void PeriodicBoundaryCondition::evaluateBoundaryPosition(
    const std::shared_ptr<Mesh>& mesh, const int boundary,
    std::vector<Particle>& particles ) const
{
    // Determine if boundary is on x, y, or z.
    int dim = std::floor( boundary / 2 );

    // Get the length of the mesh in the current dimension.
    double length = mesh->getDimensionLength( dim );

    // If left boundary, move particles to right side.
    if ( boundary % 2 == 0 )
    {
        for ( auto& p : particles )
            if ( p.r[dim] < 0.0 )
		p.r[dim] += length;
    }
    // If right boundary, move particles to left side. 
    else
    {
	for ( auto& p : particles )
	    if ( p.r[dim] > length )
	        p.r[dim] -= length;
    }	
}

//---------------------------------------------------------------------------//
// Complete Boundary sum
// (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
void PeriodicBoundaryCondition::completeBoundarySum(
    const std::shared_ptr<Mesh>& mesh, const int boundary,
    std::vector<double>& field ) const
{
    // Update left and right boundaries to be their sum.
    if ( boundary % 2 == 0 )
    {
        // Get the boundary nodes.
        const auto& nodes_left = mesh->getBoundaryNodes( boundary );
        const auto& nodes_right = mesh->getBoundaryNodes( boundary+1 );

	assert( nodes_left.size() == nodes_right.size() );

        // Apply the boundary condition. Periodic means that the boundaries
	// are the sum of the left and right side. 
	for ( std::size_t i = 0; i < nodes_right.size(); ++i )
	{
            int right_id = nodes_right[i];
	    int left_id = nodes_left[i];
	    double val_left = field[left_id];
	    double val_right = field[right_id];
	    field[left_id] += val_right;
	    field[right_id] += val_left;
	}
    }
}

//---------------------------------------------------------------------------//
// Complete Boundary sum
// (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
void PeriodicBoundaryCondition::completeBoundarySum(
    const std::shared_ptr<Mesh>& mesh, const int boundary,
    std::vector<std::array<double,3> >& field ) const
{
    // Determine if the boundary is on x, y, or z.
    //int dim = std::floor( boundary / 2 );

    // Update left and right boundaries to be their sum.
    if ( boundary % 2 == 0 )
    {
        // Get the boundary nodes.
        const auto& nodes_left = mesh->getBoundaryNodes( boundary );
        const auto& nodes_right = mesh->getBoundaryNodes( boundary+1 );

	assert( nodes_left.size() == nodes_right.size() );

        // Apply the boundary condition. Periodic means that the boundaries
	// are the sum of the left and right side. 
	for ( std::size_t i = 0; i < nodes_right.size(); ++i )
	{
            int right = nodes_right[i];
	    int left = nodes_left[i];

            for ( int dim = 0; dim < 3; ++dim )
	    {
	        double val_left = field[left][dim];
	        double val_right = field[right][dim];
	        field[left][dim] += val_right;
	        field[right][dim] += val_left;
	    }
	}
    }
}
//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end BoundaryCondition.cc
//---------------------------------------------------------------------------//
