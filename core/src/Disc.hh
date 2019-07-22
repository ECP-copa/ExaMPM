//---------------------------------------------------------------------------//
/*!
 * \file Disc.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_DISC_HH
#define EXAMPM_DISC_HH

#include "Geometry.hh"
#include "Particle.hh"

#include <array>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*!
 * \class Disc
 */
class Disc : public Geometry
{
  public:

    // Type aliases.
    using Base = Geometry;
    using VelocityField = typename Base::VelocityField;

    // Constructor.
    Disc( const std::array<double,3>& center, const double radius );

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

#endif // end EXAMPM_DISC_HH

//---------------------------------------------------------------------------//
// end Disc.hh
//---------------------------------------------------------------------------//
