//---------------------------------------------------------------------------//
/*!
 * \file NewtonianViscousStress.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_NEWTONIANVISCOUSSTRESS_HPP
#define EXAMPM_NEWTONIANVISCOUSSTRESS_HPP

#include "StressModel.hh"
#include "Particle.hh"

#include <array>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*!
 * \class NewtonianViscousStress
 * \brief Stress model for a viscous nearly-compressible newtonian fluid.
 *
 * Implements a stress model for viscous nearly-compressible newtonian fluids:
 *
 * sigma = (rho * k / rho_0) * theta * I + mu * ( grad V + grad V^T )
 *
 * where mu is the fluid dynamic viscosity, grad V the velocity gradient, rho
 * the fluid density, rho_0 the initial fluid density, k the fluid bulk
 * viscosity, theta the volumetric dilation, I the idenity, and sigma is the
 * resulting stress.
 */
class NewtonianViscousStress : public StressModel
{
  public:

    // Constructor.
    NewtonianViscousStress( const double dynamic_viscosity,
                            const double bulk_modulus,
                            const double initial_density );

    // Given a particle state calculate the particle stress.
    void calculateStress( ExaMPM::Particle& p ) const override;

  private:

    // Fluid dynamic viscosity scaled by the initial density.
    double d_viscosity;

    // Fluid bulk modulus scaled by the initial density.
    double d_bulk_modulus;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_NEWTONIANVISCOUSSTRESS_HPP

//---------------------------------------------------------------------------//
// end NewtonianViscousStress.hh
//---------------------------------------------------------------------------//
