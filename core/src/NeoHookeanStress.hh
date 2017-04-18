//---------------------------------------------------------------------------//
/*!
 * \file NeoHookeanStress.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_NEOHOOKEANSTRESS_HPP
#define EXAMPM_NEOHOOKEANSTRESS_HPP

#include "StressModel.hh"
#include "Particle.hh"

#include <array>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*!
 * \class NeoHookeanStress
 * \brief Stress model for a Neo-Hookean solid
 *
 * Implements an elastic stress model for large deformations:
 *
 * sigma = J^(-1) * mu * F * F^T + J^(-1) * (lambda * ln(J) - mu) * I
 *
 * where F is the deformation gradient, J is the determinant of the
 * deformation gradient, I is the identity matrix, lambda and mu are Lame's
 * parameters for the material, and sigma is the resulting stress.
 */
class NeoHookeanStress : public StressModel
{
  public:

    // Constructor.
    NeoHookeanStress( const double youngs_modulus,
                      const double poisson_ratio,
                      const double initial_density );

    // Given a particle state calculate the particle stress.
    void calculateStress( ExaMPM::Particle& p ) const override;

  private:

    // Bulk modulus scaled by the initial density.
    double d_K;

    // Shear modulus scaled by the initial density.
    double d_G;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_NEOHOOKEANSTRESS_HPP

//---------------------------------------------------------------------------//
// end NeoHookeanStress.hh
//---------------------------------------------------------------------------//
