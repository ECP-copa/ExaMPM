//---------------------------------------------------------------------------//
/*!
 * \file Particle.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_PARTICLE_HH
#define EXAMPM_PARTICLE_HH

#include <array>

namespace ExaMPM
{

//---------------------------------------------------------------------------//
/*!
 * \class Particle
 * \brief 3d Particle
 */
class Particle
{
  public:

    //@{
    //! Particle State.

    //! Physical location.
    std::array<double,3> r;

    //! Velocity.
    std::array<double,3> v;

    //! Deformation gradient.
    std::array<std::array<double,3>,3> F;

    //! Velocity gradient.
    std::array<std::array<double,3>,3> grad_v;

    //! Strain.
    std::array<std::array<double,3>,3> strain;

    //! Stress.
    std::array<std::array<double,3>,3> stress;

    //! Mass.
    double m;

    //! Volume.
    double volume;

    //! Material id.
    int matid;
    //@}

    //@{
    //! Grid state of particle.

    //! Adjacent node ids.
    std::array<int,8> node_ids;

    //! Node basis functions.
    std::array<double,8> basis_values;

    //! Node basis gradients.
    std::array<std::array<double,3>,8> basis_gradients;

    //@}

    // Constructor.
    Particle()
    {
        // The initial deformation gradient is the identity.
        F[0][0] = 1.0;
        F[0][1] = 0.0;
        F[0][2] = 0.0;

        F[1][0] = 0.0;
        F[1][1] = 1.0;
        F[1][2] = 0.0;

        F[2][0] = 0.0;
        F[2][1] = 0.0;
        F[2][2] = 1.0;
    }
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_PARTICLE_HH

//---------------------------------------------------------------------------//
// end Particle.hh
//---------------------------------------------------------------------------//
