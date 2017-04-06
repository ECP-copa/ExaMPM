//---------------------------------------------------------------------------//
/*!
 * \file Particle.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_PARTICLE_HH
#define EXAMPM_PARTICLE_HH

#include <vector>

namespace ExaMPM
{

//---------------------------------------------------------------------------//
/*!
 * \class Particle
 */
class Particle
{
  public:

    // Physical location
    std::vector<double> r;

    // Velocity
    std::vector<double> v;

    // Mass
    double mass;

    // Material id
    int mat;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_PARTICLE_HH

//---------------------------------------------------------------------------//
// end Particle.hh
//---------------------------------------------------------------------------//
