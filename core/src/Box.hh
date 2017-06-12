//---------------------------------------------------------------------------//
/*!
 * \file Box.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_BOX_HH
#define EXAMPM_BOX_HH

#include "Geometry.hh"
#include "Particle.hh"

#include <array>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*!
 * \class Box
 */
class Box : public Geometry
{
  public:

    // Type aliases.
    using Base = Geometry;

    // Constructor. Bounds give (-x,+x,-y,+y,-z,+z).
    Box( const std::array<double,6>& bounds );

    // Determine if a particle is in the geometry.
    bool particleInGeometry( const Particle& p ) const override;

  private:

    // Bounds.
    std::array<double,6> d_bounds;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_BOX_HH

//---------------------------------------------------------------------------//
// end Box.hh
//---------------------------------------------------------------------------//
