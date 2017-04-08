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

    // Physical location.
    std::vector<double> r;

    // Velocity.
    std::vector<double> v;

    // Mass.
    double m;

    // Density
    double rho;

    // Material id.
    int matid;

    // Total strain (history dependent).
    std::vector<std::vector<double> > strain;

    // Total specific stress (history dependent).
    std::vector<std::vector<double> > stress;

  public:

    // Default constructor.
    Particle() = default;

    // Spatial dimension constructor.
    Particle( const int space_dim )
        : r( space_dim )
        , v( space_dim )
        , strain( space_dim, std::vector<double>(space_dim,0.0) )
        , stress( space_dim, std::vector<double>(space_dim,0.0) )
    { /* ... */ }
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_PARTICLE_HH

//---------------------------------------------------------------------------//
// end Particle.hh
//---------------------------------------------------------------------------//
