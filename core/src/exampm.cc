//---------------------------------------------------------------------------//
/*!
 * \file exampm.cc
 * \brief driver
 */
//---------------------------------------------------------------------------//

#include "Box.hh"
#include "Sphere.hh"
#include "NeoHookeanStress.hh"
#include "NewtonianViscousStress.hh"
#include "StressModel.hh"
#include "BoundaryCondition.hh"
#include "ProblemManager.hh"

#include <iostream>
#include <memory>
#include <vector>
#include <array>
#include <functional>
#include <cmath>

//---------------------------------------------------------------------------//
int main( int argc, char *argv[] )
{
    // Create the problem manager.
    int num_cells_x = 25;
    int num_cells_y = 25;
    int num_cells_z = 25;
    double cell_width = 0.004;
    ExaMPM::ProblemManager manager(
        num_cells_x, num_cells_y, num_cells_z, cell_width, true );

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
    std::vector<std::shared_ptr<ExaMPM::StressModel> > materials( 1 );

    // Setup a stress model.
    double youngs_modulus = 1.0e9;
    double poisson_ratio = 0.3;
    materials[0] = std::make_shared<ExaMPM::NeoHookeanStress>(
        youngs_modulus,poisson_ratio);

    // Set the materials with the manager.
    manager.setMaterialModels( materials );

    // Create geometries.
    std::vector<std::shared_ptr<ExaMPM::Geometry> > geom( 2 );

    double density = 1000.0;
    std::array<double,3> center = { 0.04, 0.04, 0.05 };
    double radius = 0.01;
    geom[0] = std::make_shared<ExaMPM::Sphere>(center,radius);
    auto init_vf1 =
        [](const std::array<double,3>& r,std::array<double,3>& v)
        { v[0] = 1.0; v[1] = 0.0; v[2] = 0.0; };
    geom[0]->setMatId( 0 );
    geom[0]->setColor( 0 );
    geom[0]->setVelocityField( init_vf1 );
    geom[0]->setDensity( density );

    std::array<double,6> bnds = {0.07,0.09,0.03,0.05,0.04,0.06};
    geom[1] = std::make_shared<ExaMPM::Box>(bnds);
    auto init_vf2 =
        [](const std::array<double,3>& r,std::array<double,3>& v)
        { v[0] = -0.1; v[1] = 0.0; v[2] = 0.0; };
    geom[1]->setMatId( 0 );
    geom[1]->setColor( 1 );
    geom[1]->setVelocityField( init_vf2 );
    geom[1]->setDensity( density );

    // Initialize the manager.
    int order = 3;
    manager.initialize( geom, order );

    // Calculate the time step paramters.
    double wave_speed = std::sqrt( youngs_modulus / density );
    double delta_t = cell_width / wave_speed;
    double final_time = 0.75;
    int num_steps = std::ceil( final_time / delta_t );
    std::cout << "Time step size: " << delta_t << std::endl;
    std::cout << "Num steps: " << num_steps << std::endl;

    // Solve the problem.
    std::string output_file( "particles.h5" );
    int num_write = 30;
    int write_freq = num_steps / num_write;
    manager.solve( num_steps, delta_t, output_file, write_freq );

    return 0;
}

//---------------------------------------------------------------------------//
// end exampm.cc
//---------------------------------------------------------------------------//
