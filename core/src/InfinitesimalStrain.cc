//---------------------------------------------------------------------------//
/*!
 * \file InfinitesimalStrain.cc
 */
//---------------------------------------------------------------------------//

#include "InfinitesimalStrain.hh"

namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Given a particle state calculate the particle strain.
void InfinitesimalStrain::calculateStrain( ExaMPM::Particle& p ) const
{
    // Calculate F^T + F
    for ( int j = 0; j < 3; ++j )
        for ( int i = 0; i < 3; ++i )
            p.strain[i][j] = p.F[i][j] + p.F[j][i];

    // Scale by 1/2
    for ( int j = 0; j < 3; ++j )
        for ( int i = 0; i < 3; ++i )
            p.strain[i][j] *= 0.5;

    // Subtract the identity.
    for ( int i = 0; i < 3; ++i )
        p.strain[i][i] -= 1.0;
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end InfinitesimalStrain.cc
//---------------------------------------------------------------------------//
