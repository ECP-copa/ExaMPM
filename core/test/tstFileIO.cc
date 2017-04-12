//---------------------------------------------------------------------------//
/*!
 * \file   test/tstFileIO.cc
 */
//---------------------------------------------------------------------------//

#include "FileIO.hh"
#include "Geometry.hh"
#include "Box.hh"
#include "Particle.hh"

#include <vector>
#include <array>

#include "gtest_main.hh"

//---------------------------------------------------------------------------//
class FileIOTest : public ::testing::Test
{
  protected:
    void SetUp()
    { /* ... */ }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST_F(FileIOTest, file_io_test)
{
    // Create a box.
    std::array<double,6> bounds = {1.2, 2.3, 0.5, 1.1, 0.5, 1.1};
    std::shared_ptr<ExaMPM::Geometry> geometry =
        std::make_shared<ExaMPM::Box>( bounds );

    // Create some particles.
    std::vector<ExaMPM::Particle> particles( 2 );
    particles[0].r = { 2.1, 0.75, 33.2 };
    particles[0].volume = 2.0;
    particles[1].r = { 1.9, 1.0999, 1.30 };
    particles[1].volume = 2.0;

    // Assign values to the geometry.
    int matid = 3;
    double density = 3.43;
    auto vf = []( const std::array<double,3>& r,
                  std::array<double,3>& v)
              { std::copy( r.begin(), r.end(), v.begin() ); };

    geometry->setMatId( matid );
    geometry->setVelocityField( vf );
    geometry->setDensity( density );

    // Initialize particles
    geometry->initializeParticle( particles[0] );
    geometry->initializeParticle( particles[1] );

    // Write out results.
    ExaMPM::FileIO file_io( "test_file_io.h5" );
    file_io.writeTimeStep( 3, 3.332, particles );
}

//---------------------------------------------------------------------------//
// end of test/tstFileIO.cc
//---------------------------------------------------------------------------//
