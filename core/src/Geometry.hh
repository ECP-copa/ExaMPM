//---------------------------------------------------------------------------//
/*!
 * \file Geometry.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_GEOMETRY_HH
#define EXAMPM_GEOMETRY_HH

#include "Particle.hh"

#include <array>
#include <functional>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*!
 * \class Geometry
 *
 * \brief Geometry interface for the creation of initial particle state.
 */
class Geometry
{
  public:

    // Velocity field function.
    using VelocityField = std::function<
      void(const std::array<double,3>& r,std::array<double,3>& v)>;

    // Destructor.
    virtual ~Geometry() = default;

    // Determine if a particle is in the geometry.
    virtual bool particleInGeometry( const Particle& p ) const = 0;

    // Set the initial material id of the geometry.
    void setMatId( const int matid );

    // Set the initial velocity field of the geometry.
    void setVelocityField( VelocityField&& velocity_field );

    // Set the density.
    void setDensity( const double density );

    // Initialize the state of a particle in the geometry. The given particle
    // will be in the geometry.
    void initializeParticle( Particle& p ) const;

  private:

    // Material id.
    int d_matid;

    // Velocity field.
    VelocityField d_velocity_field;

    // Density.
    double d_density;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_GEOMETRY_HH

//---------------------------------------------------------------------------//
// end Geometry.hh
//---------------------------------------------------------------------------//
