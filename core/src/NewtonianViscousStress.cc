//---------------------------------------------------------------------------//
/*!
 * \file NewtonianViscousStress.cc
 */
//---------------------------------------------------------------------------//

#include "NewtonianViscousStress.hh"

#include <cmath>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Constructor.
NewtonianViscousStress::NewtonianViscousStress( const double dynamic_viscosity,
                                                const double bulk_modulus,
                                                const double initial_density,
                                                const double gamma )
    : d_viscosity( dynamic_viscosity )
    , d_bulk_modulus( bulk_modulus )
    , d_init_density( initial_density )
    , d_gamma( gamma )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Given a particle state calculate the stress.
void NewtonianViscousStress::calculateStress( ExaMPM::Particle& p ) const
{
    // Calculate the rate of deformation tensor.
    for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
            p.stress[i][j] =
                0.5 * d_viscosity * ( p.grad_v[i][j] + p.grad_v[j][i] );

    // Compute the trace of the rate of deformation tensor.
    double tr_d = 0.0;
    for ( int i = 0; i < 3; ++i )
        tr_d += p.stress[i][i];

    // Remove 1/3 the trace to get the deviatoric component.
    for ( int i = 0; i < 3; ++i )
        p.stress[i][i] -= tr_d / 3.0;

    // Multiply by 2 * viscosity to get the viscous stress component.
    for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
            p.stress[i][j] *= 2.0 * d_viscosity;

    // Calculate pressure via equation of state.
    double density = p.m / p.volume;
    double pressure =
        d_bulk_modulus * ( std::pow(density/d_init_density,d_gamma) - 1.0 );

    // Add stress from the pressure.
    for ( int i = 0; i < 3; ++i )
        p.stress[i][i] -= pressure;
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end NewtonianViscousStress.cc
//---------------------------------------------------------------------------//
