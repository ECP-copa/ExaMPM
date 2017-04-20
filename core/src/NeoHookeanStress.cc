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
#include <iostream>
namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Constructor.
NeoHookeanStress::NeoHookeanStress( const double youngs_modulus,
                                    const double poisson_ratio,
                                    const double initial_density )
{
    // Bulk modulus.
    d_K = youngs_modulus / (3 * (1 - 2*poisson_ratio) );
    d_K /= initial_density;

    // Shear modulus.
    d_G = youngs_modulus / (2 * (1+poisson_ratio) );
    d_G /= initial_density;
}

//---------------------------------------------------------------------------//
// Given a particle state calculate the stress.
void NeoHookeanStress::calculateStress( ExaMPM::Particle& p ) const
{
    // Calculate the volumetric dilation.
    double dilation = 0.0;
    for ( int i = 0; i < 3; ++i )
        dilation += p.strain[i][i];

    // Calculate the deviatoric stress.
    for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
            p.stress[i][j] = p.strain[i][j];
    for ( int i = 0; i < 3; ++i )
        p.stress[i][i] -= dilation / 3.0;
    for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
            p.stress[i][j] *= 2 * d_G;

    // Add the volumetric stress.
    for ( int i = 0; i < 3; ++i )
        p.stress[i][i] += d_K * dilation;
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end NeoHookeanStress.cc
//---------------------------------------------------------------------------//
