//---------------------------------------------------------------------------//
/*!
 * \file Geometry.cc
 */
//---------------------------------------------------------------------------//

#include "Geometry.hh"

#include <cassert>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Constructor.
Geometry::Geometry()
{
    // Geometry is initially not in motion by default.
    d_velocity_field =
        []( const std::array<double,3>& r, std::array<std::array<double,3>,8>& c )
        { 
            for ( auto& mode : c )
                std::fill(mode.begin(),mode.end(),0.0); 
        };
}

//---------------------------------------------------------------------------//
// Set the initial material id of the geometry.
void Geometry::setMatId( const int matid )
{
    d_matid = matid;
}

//---------------------------------------------------------------------------//
// Set the color of the geometry.
void Geometry::setColor( const int color )
{
    d_color = color;
}

//---------------------------------------------------------------------------//
// Set the initial velocity field of the geometry.
void Geometry::setVelocityField( VelocityField&& velocity_field )
{
    assert( velocity_field );
    d_velocity_field = velocity_field;
}

//---------------------------------------------------------------------------//
// Set the density.
void Geometry::setDensity( const double density )
{
    assert( density > 0.0 );
    d_density = density;
}

//---------------------------------------------------------------------------//
// Initialize the state of a particle in the geometry. The given particle will
// be in the geometry.
void Geometry::initializeParticle( Particle& p ) const
{
    assert( this->particleInGeometry(p) );

    // Assign the velocity.
    d_velocity_field( p.r, p.c );

    // Assign the material id.
    p.matid = d_matid;

    // Assign the color
    p.color = d_color;

    // Assign the mass
    p.m = d_density * p.volume;
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end Geometry.cc
//---------------------------------------------------------------------------//
