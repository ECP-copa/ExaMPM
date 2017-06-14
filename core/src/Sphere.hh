//---------------------------------------------------------------------------//
/*!
 * \file Sphere.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_SPHERE_HH
#define EXAMPM_SPHERE_HH

#include "Geometry.hh"
#include "Particle.hh"

#include <array>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*!
 * \class Sphere
 */
class Sphere : public Geometry
{
  public:

    // Type aliases.
    using Base = Geometry;
    using VelocityField = typename Base::VelocityField;

    // Constructor.
    Sphere( const std::array<double,3>& center, const double radius );

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

#endif // end EXAMPM_SPHERE_HH

//---------------------------------------------------------------------------//
// end Sphere.hh
//---------------------------------------------------------------------------//
