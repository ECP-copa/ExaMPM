//---------------------------------------------------------------------------//
/*!
 * \file   test/tstFileIO.cc
 */
//---------------------------------------------------------------------------//

#include "FileIO.hh"
#include "Geometry.hh"
#include "Square.hh"
#include "Particle.hh"

#include <vector>

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
    // Create a square.
    std::vector<double> bounds = {1.2, 2.3, 0.5, 1.1};
    std::shared_ptr<ExaMPM::Geometry> geometry =
        std::make_shared<ExaMPM::Square>( bounds );

    // Create some particles.
    std::vector<ExaMPM::Particle> particles( 2, ExaMPM::Particle(2) );
    particles[0].r = { 2.1, 0.75 };
    particles[1].r = { 1.9, 1.0999 };

    // Assign values to the geometry.
    int matid = 3;
    double mass = 3.43;
    auto vf = []( const std::vector<double>& r,
                  std::vector<double>& v)
              { std::copy( r.begin(), r.end(), v.begin() ); };

    geometry->setMatId( matid );
    geometry->setVelocityField( vf );
    geometry->setMass( mass );

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
