//---------------------------------------------------------------------------//
/*!
 * \file Box.cc
 */
//---------------------------------------------------------------------------//

#include "Box.hh"

#include <cassert>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Constructor. Bounds give (-x,+x,-y,+y,-z,+z).
Box::Box( const std::array<double,6>& bounds )
    : d_bounds( bounds )
    , d_matid( -1 )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Set the initial material id of the geometry.
void Box::setMatId( const int matid )
{
    d_matid = matid;
}

//---------------------------------------------------------------------------//
// Set the initial velocity field of the geometry.
void Box::setVelocityField( VelocityField&& velocity_field )
{
    assert( velocity_field );
    d_velocity_field = velocity_field;
}

//---------------------------------------------------------------------------//
// Set the density.
void Box::setDensity( const double density )
{
    assert( density > 0.0 );
    d_density = density;
}

//---------------------------------------------------------------------------//
// Determine if a particle is in the geometry.
bool Box::particleInGeometry( const Particle& p ) const
{
    return
        ( d_bounds[0] <= p.r[0] && p.r[0] <= d_bounds[1] ) &&
        ( d_bounds[2] <= p.r[1] && p.r[1] <= d_bounds[3] ) &&
        ( d_bounds[4] <= p.r[2] && p.r[2] <= d_bounds[5] );
}

//---------------------------------------------------------------------------//
// Initialize the state of a particle in the geometry. The given particle will
// be in the geometry.
void Box::initializeParticle( Particle& p ) const
{
    assert( this->particleInGeometry(p) );

    // Assign the velocity.
    d_velocity_field( p.r, p.v );

    // Assign the material id.
    p.matid = d_matid;

    // Assign the mass
    p.m = d_density * p.volume;
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end Box.cc
//---------------------------------------------------------------------------//
