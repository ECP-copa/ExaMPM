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
 *
 * \brief Stress model for a viscous nearly-incompressible viscous newtonian
 * fluid.
 *
 * Implements a stress model for viscous nearly-incompressible viscous
 * newtonian fluids based on Monaghan SPH model:
 *
 * sigma = -p * I + 2 * mu * d'
 *
 * where p is the fluid pressure, mu is the fluid dynamic viscosity, and d'
 * the deviatoric component of the rate of deformation tensor. Pressure is
 * calculated with an equation of state:
 *
 * p = K * ( (rho/rho_0)^gamma - 1 )
 *
 * where K is the bulk modulus of the fluid, rho the current fluid density,
 * rho_0 the initial fluid density, and gamma a constant.
 */
class NewtonianViscousStress : public StressModel
{
  public:

    // Constructor.
    NewtonianViscousStress( const double dynamic_viscosity,
                            const double bulk_modulus,
                            const double initial_density,
                            const double gamma );

    // Given a particle state calculate the particle stress.
    void calculateStress( ExaMPM::Particle& p ) const override;

  private:

    // Fluid dynamic viscosity
    double d_viscosity;

    // Fluid bulk modulus.
    double d_bulk_modulus;

    // Fluid initial density.
    double d_init_density;

    // Gamma constant.
    double d_gamma;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_NEWTONIANVISCOUSSTRESS_HPP

//---------------------------------------------------------------------------//
// end NewtonianViscousStress.hh
//---------------------------------------------------------------------------//
