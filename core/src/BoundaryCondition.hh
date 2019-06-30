//---------------------------------------------------------------------------//
/*!
 * \file BoundaryCondition.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_BOUNDARYCONDITION_HH
#define EXAMPM_BOUNDARYCONDITION_HH

#include "Mesh.hh"

#include <memory>
#include <vector>
#include <array>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*!
 * \class BoundaryCondition
 * \brief Interface for rigid wall conditions on the mesh boundary.
 */
class BoundaryCondition
{
  public:

    // Destructor.
    virtual ~BoundaryCondition() = default;

    // Evaluate the boundary condition for momentum on the boundary nodes.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    virtual void evaluateMomentumCondition(
        const std::shared_ptr<Mesh>& mesh,
        const int boundary,
        const std::vector<double>& mass,
        std::vector<std::array<double,3> >& momentum ) const = 0;

    // Evaluate the boundary condition for impulse on the boundary nodes.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    virtual void evaluateImpulseCondition(
        const std::shared_ptr<Mesh>& mesh,
        const int boundary,
        const std::vector<double>& mass,
        std::vector<std::array<double,3> >& impulse ) const = 0;

    // Update position of particles based on boundary conditions.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    virtual void evaluateBoundaryPosition(
	const std::shared_ptr<Mesh>& mesh, const int boundary,
        std::vector<Particle>& particles ) const = 0;

    // Correct partial sums on boundaries for periodicity
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    virtual void completeBoundarySum(
	const std::shared_ptr<Mesh>& mesh, const int boundary,
	std::vector<double>& field ) const = 0;

    // Correct partial sums on boundaries for periodicity
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    virtual void completeBoundarySum(
	const std::shared_ptr<Mesh>& mesh, const int boundary,
	std::vector<std::array<double,3> >& field ) const = 0;

};

//---------------------------------------------------------------------------//
/*!
 * \class NoSlipBoundaryCondition
 *
 * V = 0 on boundary
 */
class NoSlipBoundaryCondition : public BoundaryCondition
{
  public:

    // Evaluate the boundary condition for momentum.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void evaluateMomentumCondition(
        const std::shared_ptr<Mesh>& mesh,
        const int boundary,
        const std::vector<double>& mass,
        std::vector<std::array<double,3> >& momentum ) const override;

    // Evaluate the boundary condition for impulse.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void evaluateImpulseCondition(
        const std::shared_ptr<Mesh>& mesh,
        const int boundary,
        const std::vector<double>& mass,
        std::vector<std::array<double,3> >& impulse ) const override;

    // Update position of particles based on boundary conditions.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void evaluateBoundaryPosition(
        const std::shared_ptr<Mesh>& mesh, const int boundary,
        std::vector<Particle>& particles ) const override;

    // Correct partial sums on boundaries for periodicity
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void completeBoundarySum(
	const std::shared_ptr<Mesh>& mesh, const int boundary,
	std::vector<double>& field ) const override;

    // Correct partial sums on boundaries for periodicity
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void completeBoundarySum(
	const std::shared_ptr<Mesh>& mesh, const int boundary,
	std::vector<std::array<double,3> >& field ) const override;

};

//---------------------------------------------------------------------------//
/*!
 * \class FreeSlipBoundaryCondition
 *
 * V dot n = 0 on boundary
 */
class FreeSlipBoundaryCondition : public BoundaryCondition
{
  public:

    // Evaluate the boundary condition for momentum.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void evaluateMomentumCondition(
        const std::shared_ptr<Mesh>& mesh,
        const int boundary,
        const std::vector<double>& mass,
        std::vector<std::array<double,3> >& momentum ) const override;

    // Evaluate the boundary condition for impulse.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void evaluateImpulseCondition(
        const std::shared_ptr<Mesh>& mesh,
        const int boundary,
        const std::vector<double>& mass,
        std::vector<std::array<double,3> >& impulse ) const override;

    // Update position of particles based on boundary conditions.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void evaluateBoundaryPosition(
        const std::shared_ptr<Mesh>& mesh, const int boundary,
        std::vector<Particle>& particles ) const override;

    // Correct partial sums on boundaries for periodicity
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void completeBoundarySum(
	const std::shared_ptr<Mesh>& mesh, const int boundary,
	std::vector<double>& field ) const override;

    // Correct partial sums on boundaries for periodicity
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void completeBoundarySum(
	const std::shared_ptr<Mesh>& mesh, const int boundary,
	std::vector<std::array<double,3> >& field ) const override;

};

//---------------------------------------------------------------------------//
/*!
 * \class VelocityBoundaryCondition
 *
 * V = V_b on boundary. The wall is rigid and so the velocity field must be
 * tangential.
 */
class VelocityBoundaryCondition : public BoundaryCondition
{
  public:

    // Constructor.
    VelocityBoundaryCondition( const std::array<double,3>& vb );

    // Evaluate the boundary condition for momentum.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void evaluateMomentumCondition(
        const std::shared_ptr<Mesh>& mesh,
        const int boundary,
        const std::vector<double>& mass,
        std::vector<std::array<double,3> >& momentum ) const override;

    // Evaluate the boundary condition for impulse.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void evaluateImpulseCondition(
        const std::shared_ptr<Mesh>& mesh,
        const int boundary,
        const std::vector<double>& mass,
        std::vector<std::array<double,3> >& impulse ) const override;

    // Update position of particles based on boundary conditions.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void evaluateBoundaryPosition(
        const std::shared_ptr<Mesh>& mesh, const int boundary,
        std::vector<Particle>& particles ) const override;

    // Correct partial sums on boundaries for periodicity
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void completeBoundarySum(
	const std::shared_ptr<Mesh>& mesh, const int boundary,
	std::vector<double>& field ) const override;

    // Correct partial sums on boundaries for periodicity
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void completeBoundarySum(
	const std::shared_ptr<Mesh>& mesh, const int boundary,
	std::vector<std::array<double,3> >& field ) const override;


  private:

    // Boundary velocity.
    std::array<double,3> d_vb;
};

//---------------------------------------------------------------------------//
/*!
 * \class PeriodicBoundaryCondition
 *
 * V dot n = 0 on boundary
 */
class PeriodicBoundaryCondition : public BoundaryCondition
{
  public:

    // Evaluate the boundary condition for momentum.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void evaluateMomentumCondition(
        const std::shared_ptr<Mesh>& mesh,
        const int boundary,
        const std::vector<double>& mass,
        std::vector<std::array<double,3> >& momentum ) const override;

    // Evaluate the boundary condition for impulse.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void evaluateImpulseCondition(
        const std::shared_ptr<Mesh>& mesh,
        const int boundary,
        const std::vector<double>& mass,
        std::vector<std::array<double,3> >& impulse ) const override;

    // Update position of particles based on boundary conditions.
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void evaluateBoundaryPosition(
	const std::shared_ptr<Mesh>& mesh, const int boundary,
        std::vector<Particle>& particles ) const override;

    // Correct partial sums on boundaries for periodicity
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void completeBoundarySum(
	const std::shared_ptr<Mesh>& mesh, const int boundary,
	std::vector<double>& field ) const override;

    // Correct partial sums on boundaries for periodicity
    // (-x = 0, +x = 1, -y = 2, +y = 3, -z = 4, +z = 5)
    void completeBoundarySum(
	const std::shared_ptr<Mesh>& mesh, const int boundary,
	std::vector<std::array<double,3> >& field ) const override;

};
//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_BOUNDARYCONDITION_HH

//---------------------------------------------------------------------------//
// end BoundaryCondition.hh
//---------------------------------------------------------------------------//
