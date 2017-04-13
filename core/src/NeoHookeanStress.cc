//---------------------------------------------------------------------------//
/*!
 * \file NeoHookeanStress.cc
 */
//---------------------------------------------------------------------------//

#include "NeoHookeanStress.hh"
#include "TensorTools.hh"

#include <array>
#include <cmath>
#include <cassert>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Constructor.
NeoHookeanStress::NeoHookeanStress( const double youngs_modulus,
                                  const double poisson_ratio )
{
    // First Lame parameter.
    d_lambda = youngs_modulus * poisson_ratio /
               ( (1+poisson_ratio)*(1-2*poisson_ratio) );

    // Second Lame parameter.
    d_mu = youngs_modulus / (2 *(1+poisson_ratio) );
}

//---------------------------------------------------------------------------//
// Given a particle state calculate the stress.
void NeoHookeanStress::calculateStress(
        const ExaMPM::Particle& p,
        std::array<std::array<double,3>,3>& stress ) const
{
    // Calculate the determinant of the deformation gradient.
    double J = TensorTools::determinant( p.F );
    assert( J > 0.0 );

    // Reset the stress.
    for ( int d = 0; d < 3; ++d )
        stress[d] = { 0.0, 0.0, 0.0 };

    // Calculate F*F^T
    for ( int j = 0; j < 3; ++j )
        for ( int i = 0; i < 3; ++i )
            for ( int k = 0; k < 3; ++k )
                stress[i][j] += p.F[i][k] * p.F[j][k];

    // Subtract the identity.
    for ( int i = 0; i < 3; ++i )
        stress[i][i] -= 1.0;

    // Scale by mu.
    for ( int j = 0; j < 3; ++j )
        for ( int i = 0; i < 3; ++i )
            stress[i][j] *= d_mu;

    // Add the scaled identity.
    double c = d_lambda * std::log(J);
    for ( int i = 0; i < 3; ++i )
        stress[i][i] += c;

    // Scale by the determinant of the deformation gradient.
    for ( int j = 0; j < 3; ++j )
        for ( int i = 0; i < 3; ++i )
            stress[i][j] /= J;
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end NeoHookeanStress.cc
//---------------------------------------------------------------------------//
