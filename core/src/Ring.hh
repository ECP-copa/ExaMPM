//---------------------------------------------------------------------------//
/*!
 * \file Ring.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_RING_HH
#define EXAMPM_RING_HH

#include "Geometry.hh"
#include "Particle.hh"

#include <array>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*!
 * \class Ring
 */
class Ring : public Geometry
{
  public:

    // Type aliases.
    using Base = Geometry;
    using VelocityField = typename Base::VelocityField;

    // Constructor.
    Ring( const std::array<double,3>& center, const double radius );

    // Determine if a particle is in the geometry.
    bool particleInGeometry( const Particle& p ) const override;

  private:

    // Center.
    std::array<double,3> d_center;

    // Radius.
    double d_radius;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_RING_HH

//---------------------------------------------------------------------------//
// end Ring.hh
//---------------------------------------------------------------------------//
