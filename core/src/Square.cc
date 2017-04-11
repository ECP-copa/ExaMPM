//---------------------------------------------------------------------------//
/*!
 * \file Square.cc
 */
//---------------------------------------------------------------------------//

#include "Square.hh"

#include <cassert>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Constructor. Bounds give (-x,+x,-y,+y).
Square::Square( const std::vector<double>& bounds )
    : d_bounds( bounds )
    , d_matid( -1 )
{
    assert( 4 == d_bounds.size() );
}

//---------------------------------------------------------------------------//
// Set the initial material id of the geometry.
void Square::setMatId( const int matid )
{
    d_matid = matid;
}

//---------------------------------------------------------------------------//
// Set the initial velocity field of the geometry.
void Square::setVelocityField( VelocityField&& velocity_field )
{
    assert( velocity_field );
    d_velocity_field = velocity_field;
}

//---------------------------------------------------------------------------//
// Set the density.
void Square::setDensity( const double density )
{
    assert( density > 0.0 );
    d_density = density;
}

//---------------------------------------------------------------------------//
// Determine if a particle is in the geometry.
bool Square::particleInGeometry( const Particle& p )
{
    return
        ( d_bounds[0] <= p.r[0] && p.r[0] <= d_bounds[1] ) &&
        ( d_bounds[2] <= p.r[1] && p.r[1] <= d_bounds[3] );
}

//---------------------------------------------------------------------------//
// Initialize the state of a particle in the geometry. The given particle will
// be in the geometry.
void Square::initializeParticle( Particle& p ) const
{
    assert( ( d_bounds[0] <= p.r[0] && p.r[0] <= d_bounds[1] ) &&
            ( d_bounds[2] <= p.r[1] && p.r[1] <= d_bounds[3] ) );

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
// end Square.cc
//---------------------------------------------------------------------------//
