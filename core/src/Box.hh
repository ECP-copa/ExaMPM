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
    using VelocityField = typename Base::VelocityField;

    // Constructor. Bounds give (-x,+x,-y,+y,-z,+z).
    Box( const std::array<double,6>& bounds );

    // Set the initial material id of the geometry.
    void setMatId( const int matid ) override;

    // Set the initial velocity field of the geometry.
    void setVelocityField( VelocityField&& velocity_field ) override ;

    // Set the density.
    void setDensity( const double density ) override;

    // Determine if a particle is in the geometry.
    bool particleInGeometry( const Particle& p ) const override;

    // Initialize the state of a particle in the geometry. The given particle
    // will be in the geometry.
    void initializeParticle( Particle& p ) const override;

  private:

    // Bounds.
    std::array<double,6> d_bounds;

    // Material id.
    int d_matid;

    // Velocity field.
    VelocityField d_velocity_field;

    // Density.
    double d_density;

    // Total geometry mass.
    double d_total_mass;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_BOX_HH

//---------------------------------------------------------------------------//
// end Box.hh
//---------------------------------------------------------------------------//
