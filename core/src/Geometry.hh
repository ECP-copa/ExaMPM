//---------------------------------------------------------------------------//
/*!
 * \file Geometry.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_GEOMETRY_HH
#define EXAMPM_GEOMETRY_HH

#include "Particle.hh"

#include <vector>
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
      void(const std::vector<double>& r,std::vector<double>& v)>;

    // Destructor.
    virtual ~Geometry() = default;

    // Set the initial material id of the geometry.
    virtual void setMatId( const int matid ) = 0;

    // Set the initial velocity field of the geometry.
    virtual void setVelocityField( VelocityField&& velocity_field ) = 0;

    // Set the density.
    virtual void setDensity( const double density ) = 0;

    // Set the initial total mass of the geometry.
    virtual void setMass( const double total_mass ) = 0;

    // Set the initial total mass of the geometry.
    virtual double getMass() const = 0;

    // Determine if a particle is in the geometry. If it is, keep track of it
    // so we know how many total particles there are.
    virtual bool particleInGeometry( const Particle& p ) = 0;

    // Initialize the state of a particle in the geometry. This will be called
    // after we have counted how many particles are in the geometry with
    // particleInGeometry(). The given particle will be in the geometry.
    virtual void initializeParticle( Particle& p ) const = 0;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_GEOMETRY_HH

//---------------------------------------------------------------------------//
// end Geometry.hh
//---------------------------------------------------------------------------//
