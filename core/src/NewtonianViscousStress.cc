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
NewtonianViscousStress::NewtonianViscousStress( const double viscosity )
    : d_viscosity( viscosity )
{ /* ... */ }

//---------------------------------------------------------------------------//
// Given a particle state calculate the stress.
void NewtonianViscousStress::calculateStress(
        const ExaMPM::Particle& p,
        std::array<std::array<double,3>,3>& stress ) const
{
    // For now we assume a constant relative pressure of zero. This will
    // change when an energy equation is added.
    for ( int i = 0; i < 3; ++i )
        for ( int j = 0; j < 3; ++j )
            stress[i][j] = d_viscosity * ( p.grad_v[i][j] + p.grad_v[j][i] );
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end NewtonianViscousStress.cc
//---------------------------------------------------------------------------//
