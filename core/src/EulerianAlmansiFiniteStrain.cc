//---------------------------------------------------------------------------//
/*!
 * \file EulerianAlmansiFiniteStrain.cc
 */
//---------------------------------------------------------------------------//

#include "EulerianAlmansiFiniteStrain.hh"
#include "TensorTools.hh"

namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Given a particle state calculate the particle strain.
void EulerianAlmansiFiniteStrain::calculateStrain( ExaMPM::Particle& p ) const
{
    // Calculate F*F^T
    std::array<std::array<double,3>,3> fft;
    for ( int i = 0; i < 3; ++i )
        std::fill( fft[i].begin(), fft[i].end(), 0.0 );
    for ( int j = 0; j < 3; ++j )
        for ( int i = 0; i < 3; ++i )
            for ( int k = 0; k < 3; ++k )
                fft[i][j] += p.F[i][k] * p.F[j][k];

    // Invert F*F^T.
    TensorTools::inverse( fft, p.strain );

    // Subtract the identity.
    for ( int i = 0; i < 3; ++i )
        p.strain[i][i] -= 1.0;

    // Scale by -1/2
    for ( int j = 0; j < 3; ++j )
        for ( int i = 0; i < 3; ++i )
           p.strain[i][j] *= -0.5;
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end EulerianAlmansiFiniteStrain.cc
//---------------------------------------------------------------------------//
