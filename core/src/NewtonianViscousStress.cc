//---------------------------------------------------------------------------//
/*!
 * \file NewtonianViscousStress.cc
 */
//---------------------------------------------------------------------------//

#include "NewtonianViscousStress.hh"

namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Constructor.
NewtonianViscousStress::NewtonianViscousStress( const double dynamic_viscosity,
                                                const double bulk_modulus,
                                                const double initial_density )
    : d_viscosity( dynamic_viscosity )
    , d_bulk_modulus( bulk_modulus / initial_density )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Given a particle state calculate the stress.
void NewtonianViscousStress::calculateStress( ExaMPM::Particle& p ) const
{
    // Calculate the volumetric dilation. The dilation is the trace of the
    // particle strain. Strain = 0.5*( F^T * F  - I ) for large deformation
    // theory with lagrangian finite strain tensor.
    double dilation = 0.0;
    for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
            dilation += p.F[j][i] * p.F[j][i];
    dilation -= 3.0;
    dilation *= 0.5;

    // Calculate deviatoric stress.
    for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
            p.stress[i][j] = d_viscosity * ( p.grad_v[i][j] + p.grad_v[j][i] );

    // Add the volumetric stress.
    double density = p.m / p.volume;
    for ( int i = 0; i < 3; ++i )
        p.stress[i][i] += density * d_bulk_modulus * dilation;
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end NewtonianViscousStress.cc
//---------------------------------------------------------------------------//
