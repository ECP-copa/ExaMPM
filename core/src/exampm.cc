//---------------------------------------------------------------------------//
/*!
 * \file exampm.cc
 * \brief driver
 */
//---------------------------------------------------------------------------//

#include "Box.hh"
#include "StressModel.hh"
#include "NeoHookeanStress.hh"
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
    int num_cells_x = 25;
    int num_cells_y = 25;
    int num_cells_z = 25;
    double cell_width = 0.04;
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

    // Setup a stress model.
    std::vector<std::shared_ptr<ExaMPM::StressModel> > materials( 1 );
    double youngs_modulus = 1.0e7;
    double poisson_ratio = 0.4;
    materials[0] = std::make_shared<ExaMPM::NeoHookeanStress>(
        youngs_modulus, poisson_ratio );

    // Set the materials with the manager.
    manager.setMaterialModels( materials );

    // Create geometries.
    std::vector<std::shared_ptr<ExaMPM::Geometry> > geom( 2 );

    double density = 1.0e3;

    std::array<double,6> bnds = {0.2,0.4,0.4,0.6,0.4,0.6};
    geom[0] = std::make_shared<ExaMPM::Box>(bnds);
    auto init_vf1 =
        [=](const std::array<double,3>& r,std::array<double,3>& v)
        { v[0] = 0.1; v[1] = 0.0; v[2] = 0.0; };
    geom[0]->setMatId( 0 );
    geom[0]->setVelocityField( init_vf1 );
    geom[0]->setDensity( density );

    bnds = {0.6,0.8,0.5,0.7,0.4,0.6};
    geom[1] = std::make_shared<ExaMPM::Box>(bnds);
    auto init_vf2 =
        [=](const std::array<double,3>& r,std::array<double,3>& v)
        { v[0] = -0.1; v[1] = 0.0; v[2] = 0.0; };
    geom[1]->setMatId( 0 );
    geom[1]->setVelocityField( init_vf2 );
    geom[1]->setDensity( density );

    // Initialize the manager.
    int order = 2;
    manager.initialize( geom, order );

    // Solve the problem.
    int num_steps = 10000;
    double delta_t = 1.0e-4;
    std::string output_file( "particles.h5" );
    int write_freq = 100;
    manager.solve( num_steps, delta_t, output_file, write_freq );

    return 0;
}

//---------------------------------------------------------------------------//
// end exampm.cc
//---------------------------------------------------------------------------//
