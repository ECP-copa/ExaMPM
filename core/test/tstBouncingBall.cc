//---------------------------------------------------------------------------//
/*!
 * \file test/tstBouncingBall.cc
 */
//---------------------------------------------------------------------------//

#include "Sphere.hh"
#include "NeoHookeanStress.hh"
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
class BouncingBallTest : public ::testing::Test
{
  protected:
    void SetUp()
    { /* ... */ }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST_F(BouncingBallTest, bouncing_ball_test)
{
    // Create the problem manager.
    int num_cells_x = 20;
    int num_cells_y = 20;
    int num_cells_z = 20;
    double cell_width = 0.05;
    bool has_gravity = true;
    int thread_count = ( ExaMPM::test_argc == 2 ) ? std::atoi(ExaMPM::test_argv[1]) : 1;

    ExaMPM::ProblemManager manager(
        num_cells_x, num_cells_y, num_cells_z, cell_width, has_gravity, thread_count );

    // Create boundary conditions.
    std::array<std::shared_ptr<ExaMPM::BoundaryCondition>,6> bc;
    bc[0] = std::make_shared<ExaMPM::PeriodicBoundaryCondition>();
    bc[1] = std::make_shared<ExaMPM::PeriodicBoundaryCondition>();
    bc[2] = std::make_shared<ExaMPM::PeriodicBoundaryCondition>();
    bc[3] = std::make_shared<ExaMPM::PeriodicBoundaryCondition>();
    bc[4] = std::make_shared<ExaMPM::PeriodicBoundaryCondition>();
    bc[5] = std::make_shared<ExaMPM::PeriodicBoundaryCondition>();

    // Set boundary conditions with the manager.
    manager.setBoundaryConditions( bc );

    // Create materials.
    std::vector<std::shared_ptr<ExaMPM::StressModel> > materials( 1 );

    // Setup a stress model.
    double density = 1000.0;
    double youngs_modulus = 1.0e9;
    double poisson_ratio = 0.3;
    materials[0] = std::make_shared<ExaMPM::NeoHookeanStress>(
        youngs_modulus,poisson_ratio);

    // Set the materials with the manager.
    manager.setMaterialModels( materials );

    // Create the ball
    std::vector<std::shared_ptr<ExaMPM::Geometry> > geom( 1 );
    std::array<double,3> center = { 0.25, 0.5, 0.75 };
    double radius = 0.1;
    geom[0] = std::make_shared<ExaMPM::Sphere>(center,radius);
    geom[0]->setMatId( 0 );
    geom[0]->setColor( 0 );
    geom[0]->setDensity( density );

    // Set the ball moving to the right in X.
    auto init_vf1 =
        [](const std::array<double,3>& r,std::array<std::array<double,3>,8>& c)
        { 
            for ( auto& mode : c )
                std::fill(mode.begin(),mode.end(),0.0);
            c[0][0] = 2.0; c[0][1] = 0.0; c[0][2] = 0.0; 
        };
    geom[0]->setVelocityField( init_vf1 );

    // Initialize the manager.
    int order = 3;
    manager.initialize( geom, order );

    // Calculate the time step paramters.
    double delta_t = 1.0e-5;
    double final_time = 5.0;
    int num_steps = std::ceil( final_time / delta_t );
    std::cout << "Time step size: " << delta_t << std::endl;
    std::cout << "Num steps: " << num_steps << std::endl;

    // Solve the problem.
    std::string output_file( "bouncing_ball_particles" );
    int num_write = 100;
    int write_freq = num_steps / num_write;
    manager.solve( num_steps, delta_t, output_file, write_freq );
}

//---------------------------------------------------------------------------//
// end of test/tstBouncingBall.cc
//---------------------------------------------------------------------------//
