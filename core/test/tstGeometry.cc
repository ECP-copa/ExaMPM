//---------------------------------------------------------------------------//
/*!
 * \file   test/tstGeometry.cc
 */
//---------------------------------------------------------------------------//

#include "Geometry.hh"
#include "Square.hh"
#include "Particle.hh"

#include <vector>

#include "gtest_main.hh"

//---------------------------------------------------------------------------//
class GeometryTest : public ::testing::Test
{
  protected:
    void SetUp()
    { /* ... */ }
};

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//
TEST_F(GeometryTest, square_test)
{
    // Create a square.
    std::vector<double> bounds = {1.2, 2.3, 0.5, 1.1};
    std::shared_ptr<ExaMPM::Geometry> geometry =
        std::make_shared<ExaMPM::Square>( bounds );

    // Create some particles.
    ExaMPM::Particle p1( 2 );
    p1.r = { 2.1, 0.75 };

    ExaMPM::Particle p2( 2 );
    p2.r = { 1.9, 1.0999 };

    ExaMPM::Particle p3( 2 );
    p3.r = { 32.3, 22.9 };

    // Check the geometric location.
    bool p1_inside = geometry->particleInGeometry(p1);
    EXPECT_TRUE( p1_inside );

    bool p2_inside = geometry->particleInGeometry(p2);
    EXPECT_TRUE( p2_inside );

    bool p3_inside = geometry->particleInGeometry(p3);
    EXPECT_FALSE( p3_inside );

    // Assign values to the geometry.
    int matid = 3;
    double mass = 3.43;
    auto vf = []( const std::vector<double>& r,
                  std::vector<double>& v)
              { std::copy( r.begin(), r.end(), v.begin() ); };

    geometry->setMatId( matid );
    geometry->setVelocityField( vf );
    geometry->setMass( mass );
    EXPECT_EQ( geometry->getMass(), mass );

    // Check the initialization.
    int np_in_geom = 2;

    geometry->initializeParticle( p1 );
    EXPECT_FLOAT_EQ( p1.v[0], p1.r[0] );
    EXPECT_FLOAT_EQ( p1.v[1], p1.r[1] );
    EXPECT_EQ( p1.matid, matid );

    geometry->initializeParticle( p2 );
    EXPECT_FLOAT_EQ( p2.v[0], p2.r[0] );
    EXPECT_FLOAT_EQ( p2.v[1], p2.r[1] );
    EXPECT_EQ( p2.matid, matid );
}

//---------------------------------------------------------------------------//
// end of test/tstGeometry.cc
//---------------------------------------------------------------------------//
