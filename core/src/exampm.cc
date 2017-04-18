//---------------------------------------------------------------------------//
/*!
 * \file exampm.cc
 * \brief driver
 */
//---------------------------------------------------------------------------//

#include "Box.hh"
#include "NeoHookeanStress.hh"
#include "NewtonianViscousStress.hh"
#include "LagrangianFiniteStrain.hh"
#include "MaterialModel.hh"
#include "BoundaryCondition.hh"
#include "ProblemManager.hh"

#include <memory>
#include <vector>
#include <array>
#include <functional>

//---------------------------------------------------------------------------//
int main( int argc, char *argv[] )
{
    // Create the problem manager.
    int num_cells_x = 24;
    int num_cells_y = 2;
    int num_cells_z = 24;
    double cell_width = 0.25;
    ExaMPM::ProblemManager manager(
        num_cells_x, num_cells_y, num_cells_z, cell_width );

    // Create boundary conditions.
    std::array<std::shared_ptr<ExaMPM::BoundaryCondition>,6> bc;
    bc[0] = std::make_shared<ExaMPM::FreeSlipBoundaryCondition>();
    bc[1] = std::make_shared<ExaMPM::FreeSlipBoundaryCondition>();
    bc[2] = std::make_shared<ExaMPM::FreeSlipBoundaryCondition>();
    bc[3] = std::make_shared<ExaMPM::FreeSlipBoundaryCondition>();
    bc[4] = std::make_shared<ExaMPM::FreeSlipBoundaryCondition>();
    bc[5] = std::make_shared<ExaMPM::FreeSlipBoundaryCondition>();

    // Set boundary conditions with the manager.
    manager.setBoundaryConditions( bc );

    // Create materials.
    std::vector<ExaMPM::MaterialModel> materials( 1 );

    // Setup a strain model.
    materials[0].strain_model =
        std::make_shared<ExaMPM::LagrangianFiniteStrain>();

    // Setup a stress model.
    double viscosity = 1.0e-3;
    double density = 1.0e3;
    double bulk_modulus = 2.0e9;
    materials[0].stress_model =
        std::make_shared<ExaMPM::NewtonianViscousStress>(
            viscosity,bulk_modulus,density);

    // Set the materials with the manager.
    manager.setMaterialModels( materials );

    // Set gravity.
    auto gravity = [](const std::array<double,3>& r,
                      std::array<double,3>& f)
                   { f[0] = 0.0; f[1] = 0.0; f[2] = -9.81; };
    manager.setSpecificBodyForce( gravity );

    // Create geometries.
    std::vector<std::shared_ptr<ExaMPM::Geometry> > geom( 1 );

    std::array<double,6> bnds = {0.0,4.0,0.0,0.5,0.0,2.0};
    geom[0] = std::make_shared<ExaMPM::Box>(bnds);
    auto init_vf1 =
        [=](const std::array<double,3>& r,std::array<double,3>& v)
        { v[0] = 0.0; v[1] = 0.0; v[2] = 0.0; };
    geom[0]->setMatId( 0 );
    geom[0]->setVelocityField( init_vf1 );
    geom[0]->setDensity( density );

    // Initialize the manager.
    int order = 3;
    manager.initialize( geom, order );

    // Solve the problem.
    int num_steps = 13000;
    double delta_t = 1.0e-4;
    std::string output_file( "particles.h5" );
    int write_freq = 1000;
    manager.solve( num_steps, delta_t, output_file, write_freq );

    return 0;
}

//---------------------------------------------------------------------------//
// end exampm.cc
//---------------------------------------------------------------------------//
