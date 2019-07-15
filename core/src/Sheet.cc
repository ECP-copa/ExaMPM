//---------------------------------------------------------------------------//
/*!
 * \file Sheet.cc
 */
//---------------------------------------------------------------------------//

#include "Sheet.hh"

#include <cassert>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Constructor. Bounds give (-x,+x,-y,+y,-z,+z).
Sheet::Sheet( const std::array<double,6>& bounds, 
	      const std::array<double,3>& center,
              const double radius )
    : d_bounds( bounds )
    , d_center( center )
    , d_radius( radius )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Determine if a particle is in the geometry.
bool Sheet::particleInGeometry( const Particle& p ) const
{
    std::array<double,3> ref_p = { p.r[0] - d_center[0],
                                   p.r[1] - d_center[1] };
    double dist = ref_p[0]*ref_p[0] + ref_p[1]*ref_p[1];
    return
        ( d_bounds[0] <= p.r[0] && p.r[0] <= d_bounds[1] ) &&
        ( d_bounds[2] <= p.r[1] && p.r[1] <= d_bounds[3] ) &&
        ( d_bounds[4] <= p.r[2] && p.r[2] <= d_bounds[5] ) &&
        ( dist >= d_radius*d_radius);
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end Sheet.cc
//---------------------------------------------------------------------------//
