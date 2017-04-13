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
 * \brief Stress model for a incompressible viscous newtonian fluid.
 *
 * Implements a stress model for incompressible viscous newtonian fluids:
 *
 * sigma = mu * ( grad V + grad V^T )
 *
 * where mu is the fluid viscosity, grad V the velocity gradient, and sigma is
 * the resulting stress.
 */
class NewtonianViscousStress : public StressModel
{
  public:

    // Constructor.
    NewtonianViscousStress( const double viscosity );

    // Given a particle state calculate the particle stress.
    void calculateStress(
        const ExaMPM::Particle& p,
        std::array<std::array<double,3>,3>& stress ) const override;

  private:

    // Fluid viscosity.
    double d_viscosity;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_NEWTONIANVISCOUSSTRESS_HPP

//---------------------------------------------------------------------------//
// end NewtonianViscousStress.hh
//---------------------------------------------------------------------------//
