//---------------------------------------------------------------------------//
/*!
 * \file test/tstTaylorGreen.cc
 */
//---------------------------------------------------------------------------//

#include "Box.hh"
#include "Sphere.hh"
#include "NewtonianViscousStress.hh"
#include "StressModel.hh"
#include "BoundaryCondition.hh"
#include "ProblemManager.hh"

#include <memory>
#include <vector>
#include <array>
#include <functional>
#include <cmath>

#include "gtest_main.hh"

//---------------------------------------------------------------------------//
class TaylorGreenTest : public ::testing::Test
{
  protected:
    void SetUp()
    { /* ... */ }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST_F(TaylorGreenTest, taylor_green_test)
{
    // Create the problem manager.
    int num_cells_x = 20;
    int num_cells_y = 20;
    int num_cells_z = 1;
    double cell_width = 0.314159;
    bool has_gravity = false;
    int thread_count = ( ExaMPM::test_argc == 2 ) ? std::atoi(ExaMPM::test_argv[1]) : 1;

    ExaMPM::ProblemManager manager(
        num_cells_x, num_cells_y, num_cells_z, cell_width, has_gravity, thread_count );

    // Create boundary conditions.
    std::array<std::shared_ptr<ExaMPM::BoundaryCondition>,6> bc;
    bc[0] = std::make_shared<ExaMPM::FreeSlipBoundaryCondition>();
    bc[1] = std::make_shared<ExaMPM::FreeSlipBoundaryCondition>();
    bc[2] = std::make_shared<ExaMPM::FreeSlipBoundaryCondition>();
    bc[3] = std::make_shared<ExaMPM::PeriodicBoundaryCondition>();
    bc[4] = std::make_shared<ExaMPM::PeriodicBoundaryCondition>();
    bc[5] = std::make_shared<ExaMPM::FreeSlipBoundaryCondition>();

    // Set boundary conditions with the manager.
    manager.setBoundaryConditions( bc );

    // Create materials.
    std::vector<std::shared_ptr<ExaMPM::StressModel> > materials( 1 );

    // Setup a stress model.
    double dynamic_viscosity = 3.5e-3;
    double density = 1100.0;
    double bulk_modulus = 2.0e9;
    double gamma = 7.0;
    materials[0] = std::make_shared<ExaMPM::NewtonianViscousStress>(
        dynamic_viscosity,bulk_modulus,density,gamma);

    // Set the materials with the manager.
    manager.setMaterialModels( materials );

    // Geometry
    std::vector<std::shared_ptr<ExaMPM::Geometry> > geom( 1 );

    // Create the initial fluid pool.
    std::array<double,6> bnds = {0.0,6.28318,0.0,6.28318,0.0,0.314159};
    geom[0] = std::make_shared<ExaMPM::Box>(bnds);
    geom[0]->setMatId( 0 );
    geom[0]->setColor( 1 );
    geom[0]->setDensity( density );

    // Give the droplet a downward velocity equal to a 300mm free fall
    auto init_vf1 =
        [](const std::array<double,3>& r,std::array<std::array<double,3>,8>& c)
        {
            for ( auto& mode : c )
                std::fill(mode.begin(),mode.end(),0.0);
            c[0][0] = cos(r[0]) * sin(r[1]); 
            c[0][1] = -1.0 * sin(r[0]) * cos(r[1]); 
        };
    geom[0]->setVelocityField( init_vf1 );

    // Initialize the manager.
    int order = 2;
    manager.initialize( geom, order );

    // Calculate the time step paramters.
    double delta_t = 1.0e-5;
    double final_time = 3.0;
    int num_steps = std::ceil( final_time / delta_t );
    std::cout << "Time step size: " << delta_t << std::endl;
    std::cout << "Num steps: " << num_steps << std::endl;

    // Solve the problem.
    std::string output_file( "taylor_green_particles" );
    int num_write = 50;
    int write_freq = num_steps / num_write;
    manager.solve( num_steps, delta_t, output_file, write_freq );
}

//---------------------------------------------------------------------------//
// end of test/tstTaylorGreen.cc
//---------------------------------------------------------------------------//
