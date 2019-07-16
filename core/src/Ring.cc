//---------------------------------------------------------------------------//
/*!
 * \file Ring.cc
 */
//---------------------------------------------------------------------------//

#include "Ring.hh"

#include <cassert>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Constructor. Bounds give (-x,+x,-y,+y,-z,+z).
Ring::Ring( const std::array<double,3>& center, const double radius )
    : d_center( center )
    , d_radius( radius )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Determine if a particle is in the geometry.
bool Ring::particleInGeometry( const Particle& p ) const
{
    std::array<double,2> ref_p = { p.r[0] - d_center[0],
                                   p.r[1] - d_center[1] };

    double dist = ref_p[0]*ref_p[0] + ref_p[1]*ref_p[1];
    double r1 = ( d_radius - 0.02 ) * ( d_radius - 0.02 );
    double r2 = d_radius * d_radius;
    return
        dist <= r2 && dist >= r1;
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end Ring.cc
//---------------------------------------------------------------------------//
