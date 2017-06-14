//---------------------------------------------------------------------------//
/*!
 * \file StressModel.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_STRESSMODEL_HPP
#define EXAMPM_STRESSMODEL_HPP

#include "Particle.hh"

#include <array>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*
 * \class Stress model.
 * \brief Stress model for a given material.
 */
class StressModel
{
  public:

    // Destructor
    virtual ~StressModel() = default;

    // Given a particle state calculate the particle stress.
    virtual void calculateStress( ExaMPM::Particle& p ) const = 0;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_STRESSMODEL_HPP

//---------------------------------------------------------------------------//
// end StressModel.hh
//---------------------------------------------------------------------------//
