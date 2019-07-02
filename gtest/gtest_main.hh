//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   gtest/gtest_main.hh
 */
//---------------------------------------------------------------------------//

#ifndef gtest_gtest_main_hh
#define gtest_gtest_main_hh

#include "gtest.h"

namespace ExaMPM
{
//---------------------------------------------------------------------------//
static int test_argc;
static char** test_argv;

int gtest_main(int argc, char *argv[])
{
    // Initialize google test
    ::testing::InitGoogleTest(&argc, argv);

    // Capture command line args.
    test_argc = argc;
    test_argv = argv;

    // Run everything
    int failed = RUN_ALL_TESTS();

    // Print results
    if (argc)
        std::cout << "In " << argv[0] << ", ";
    std::cout << "overall test result: "
              << (failed ? "FAILED" : "PASSED")
              << std::endl;

    // Return 1 if any failure, 0 if all success
    return (failed ? 1 : 0);
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//


int main(int argc, char *argv[])
{
    ExaMPM::gtest_main(argc, argv);
}

//---------------------------------------------------------------------------//
#endif // gtest_gtest_main_hh

//---------------------------------------------------------------------------//
//              end of gtest/gtest_main.hh
//---------------------------------------------------------------------------//
