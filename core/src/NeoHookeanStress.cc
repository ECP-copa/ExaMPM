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
    // First parameter.
    d_lambda = youngs_modulus * poisson_ratio /
               ( (1 + poisson_ratio)*(1 - 2*poisson_ratio) );

    // Second parameter.
    d_mu = youngs_modulus / ( 2 * (1 + poisson_ratio) );
}

//---------------------------------------------------------------------------//
// Given a particle state calculate the stress.
void NeoHookeanStress::calculateStress( ExaMPM::Particle& p ) const
{
    // Reset the stress.
    for ( auto& s : p.stress )
        std::fill( s.begin(), s.end(), 0.0 );

    // Compute the determinant of F.
    double J = TensorTools::determinant( p.F );

    // Compute constants.
    double c1 = d_lambda * std::log(J) / J;
    double c2 = d_mu / J;

    // Compute F*F^T
    for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
            for ( int k = 0; k < 3; ++k )
                p.stress[i][j] += p.F[i][k] * p.F[j][k];

    // Scale F*F^T
    for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
            p.stress[i][j] *= c2;

    // Add the scaled identity.
    for ( int i = 0; i < 3; ++i )
        p.stress[i][i] += c1 - c2;
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end NeoHookeanStress.cc
//---------------------------------------------------------------------------//
