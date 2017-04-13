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

    // Get the boundary nodes.
    const auto& nodes = mesh->getBoundaryNodes( boundary );

    // Apply the boundary condition. Free slip means the velocity component
    // normal to the grid boundary is zero and therefore so is the impulse.
    for ( auto n : nodes )
        impulse[n][dim] = 0.0;
}

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

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end BoundaryCondition.cc
//---------------------------------------------------------------------------//
