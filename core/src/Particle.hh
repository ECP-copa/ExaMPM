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

    //@{
    //! Particle State.

    //! Physical location.
    std::vector<double> r;

    //! Velocity.
    std::vector<double> v;

    //! Mass.
    double m;

    //! Density
    double rho;

    //! Volume.
    double volume;

    //! Material id.
    int matid;

    //! Total strain (history dependent).
    std::vector<std::vector<double> > strain;

    //! Total specific stress (history dependent).
    std::vector<std::vector<double> > stress;
    //@}

    //@{
    //! Grid state of particle.

    //! Adjacent node ids.
    std::vector<int> node_ids;

    //! Node basis functions.
    std::vector<double> basis_values;

    //! Node basis gradients.
    std::vector<std::vector<double> > basis_gradients;

    //@}

  public:

    // Default constructor.
    Particle() = default;

    // Spatial dimension constructor.
    Particle( const int space_dim, const int nodes_per_cell )
        : r( space_dim )
        , v( space_dim )
        , strain( space_dim, std::vector<double>(space_dim,0.0) )
        , stress( space_dim, std::vector<double>(space_dim,0.0) )
        , node_ids( nodes_per_cell )
        , basis_values( nodes_per_cell )
        , basis_gradients( nodes_per_cell, std::vector<double>(space_dim,0.0) )
    { /* ... */ }
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_PARTICLE_HH

//---------------------------------------------------------------------------//
// end Particle.hh
//---------------------------------------------------------------------------//
