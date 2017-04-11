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

    //! Volume.
    double volume;

    //! Material id.
    int matid;

    //! Deformation gradient.
    std::vector<std::vector<double> > F;
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
        , F( space_dim, std::vector<double>(space_dim) )
        , node_ids( nodes_per_cell )
        , basis_values( nodes_per_cell )
        , basis_gradients( nodes_per_cell, std::vector<double>(space_dim,0.0) )
    {
        // The initial deformation gradient is the identity.
        F[0][0] = 1.0;
        F[0][1] = 0.0;
        F[1][0] = 0.0;
        F[1][1] = 1.0;
    }
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_PARTICLE_HH

//---------------------------------------------------------------------------//
// end Particle.hh
//---------------------------------------------------------------------------//
