//---------------------------------------------------------------------------//
/*!
 * \file LagrangianFiniteStrain.cc
 */
//---------------------------------------------------------------------------//

#include "LagrangianFiniteStrain.hh"

namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Given a particle state calculate the particle strain.
void LagrangianFiniteStrain::calculateStrain( ExaMPM::Particle& p ) const
{
    // Reset the particle strain.
    for ( int i = 0; i < 3; ++i )
        std::fill( p.strain[i].begin(), p.strain[i].end(), 0.0 );

    // Calculate F^T*F
    for ( int j = 0; j < 3; ++j )
        for ( int i = 0; i < 3; ++i )
            for ( int k = 0; k < 3; ++k )
                p.strain[i][j] += p.F[k][i] * p.F[k][j];

    // Subtract the identity.
    for ( int i = 0; i < 3; ++i )
        p.strain[i][i] -= 1.0;

    // Scale by 1/2
    for ( int j = 0; j < 3; ++j )
        for ( int i = 0; i < 3; ++i )
           p.strain[i][j] *= 0.5;
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end LagrangianFiniteStrain.cc
//---------------------------------------------------------------------------//
