//---------------------------------------------------------------------------//
/*!
 * \file Sheet.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_SHEET_HH
#define EXAMPM_SHEET_HH

#include "Geometry.hh"
#include "Particle.hh"

#include <array>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*!
 * \class Sheet
 */
class Sheet : public Geometry
{
  public:

    // Type aliases.
    using Base = Geometry;

    // Constructor. Bounds give (-x,+x,-y,+y,-z,+z).
    Sheet( const std::array<double,6>& bounds, 
           const std::array<double,3>& center,
	   const double radius );

    // Determine if a particle is in the geometry.
    bool particleInGeometry( const Particle& p ) const override;

  private:

    // Bounds.
    std::array<double,6> d_bounds;

    // Center.
    std::array<double,3> d_center;

    // Radius.
    double d_radius;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_SHEET_HH

//---------------------------------------------------------------------------//
// end Sheet.hh
//---------------------------------------------------------------------------//
