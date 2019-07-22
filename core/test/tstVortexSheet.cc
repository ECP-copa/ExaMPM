//---------------------------------------------------------------------------//
/*!
 * \file test/tstVortexSheet.cc
 */
//---------------------------------------------------------------------------//

#include "Sheet.hh"
#include "Ring.hh"
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
class VortexSheetTest : public ::testing::Test
{
  protected:
    void SetUp()
    { /* ... */ }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST_F(VortexSheetTest, vortex_sheet_test)
{
    // Create the problem manager.
    int num_cells_x = 20;
    int num_cells_y = 20;
    int num_cells_z = 1;
    double cell_width = 0.02;
    bool has_gravity = false;
    int thread_count = ( ExaMPM::test_argc == 2 ) ? std::atoi(ExaMPM::test_argv[1]) : 1;

    ExaMPM::ProblemManager manager(
        num_cells_x, num_cells_y, num_cells_z, cell_width, has_gravity, thread_count );

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
    double dynamic_viscosity = 0;//1.0e-3;
    double density = 1000.0;
    double bulk_modulus = 2.0e9;
    double gamma = 7.0;
    materials[0] = std::make_shared<ExaMPM::NewtonianViscousStress>(
        dynamic_viscosity,bulk_modulus,density,gamma);

    // Set the materials with the manager.
    manager.setMaterialModels( materials );

    // Geometry
    std::vector<std::shared_ptr<ExaMPM::Geometry> > geom( 2 );

    // Create the initial fluid pool.
    std::array<double,6> bnds = {0.0,0.4,0.0,0.4,0.0,0.02};
    std::array<double,3> center = { 0.2, 0.2, 0.01 };
    double radius = 0.1;
    geom[0] = std::make_shared<ExaMPM::Sheet>(bnds, center, radius);
    geom[0]->setMatId( 0 );
    geom[0]->setColor( 1 );
    geom[0]->setDensity( density );

    // Create the droplet.
    geom[1] = std::make_shared<ExaMPM::Ring>(center,radius);
    geom[1]->setMatId( 0 );
    geom[1]->setColor( 2 );
    geom[1]->setDensity( density );

    // Give the droplet a downward velocity equal to a 300mm free fall
    auto init_vf1 =
        [](const std::array<double,3>& r,std::array<std::array<double,3>,8>& c)
        {
            for ( auto& mode : c )
                std::fill(mode.begin(),mode.end(),0.0);
            std::array<double,3> center = { 0.2, 0.2, 0.01 };
            double x = r[0] - center[0];
            double y = r[1] - center[1];
            double rad = x*x + y*y;

            if ( y == 0 )
            {
                c[0][1] = ( x > 0 ) ? rad : -1.0*rad;
            }
            else
            {
                double den = sqrt(1.0 + pow( x/y , 2 ) );
                double sign = abs(y) / y;

                c[0][0] = -1.0 * rad * sign / den;
                c[0][1] = rad * x * sign / ( y * den );
            }
        };
    geom[1]->setVelocityField( init_vf1 );

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
    std::string output_file( "vortex_sheet_particles" );
    int num_write = 50;
    int write_freq = num_steps / num_write;
    manager.solve( num_steps, delta_t, output_file, write_freq );
}

//---------------------------------------------------------------------------//
// end of test/tstVortexSheet.cc
//---------------------------------------------------------------------------//
