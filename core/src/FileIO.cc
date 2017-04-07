//---------------------------------------------------------------------------//
/*!
 * \file FileIO.cc
 */
//---------------------------------------------------------------------------//

#include "FileIO.hh"

#include <cassert>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
// Constructor.
FileIO::FileIO( const std::string& file_name )
{
    d_handle =
        H5Fcreate( file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
}

//---------------------------------------------------------------------------//
// Write particle data at a given time step to a file.
void FileIO::writeTimeStep( const int time_step,
                            const double time,
                            const std::vector<Particle>& particles )
{
    // Get the number of particles.
    int num_p = particles.size();

    // Check that we got some particles.
    assert( num_p > 0 );

    // Get the spatial dimension of the problem.
    int space_dim = particles[0].r.size();

    // Create the group name.
    std::string group_name = "TIME_STEP_" + std::to_string(time_step);

    // Make sure we didnt already write this time step.
    assert( !H5Lexists(d_handle, group_name.c_str(), H5P_DEFAULT) );

    // Create the group for the time step.
    hid_t group_handle = H5Gcreate( d_handle,
                                    group_name.c_str(),
                                    H5P_DEFAULT,
                                    H5P_DEFAULT,
                                    H5P_DEFAULT);

    // Write the time step number.
    writeScalar( group_handle, "TIME_STEP", time_step );

    // Write the current time.
    writeScalar( group_handle, "TIME", time );

    // Write the particle positions.
    std::vector<double> data( num_p );
    {
        for ( int d = 0; d < space_dim; ++d )
        {
            for ( int p = 0; p < num_p; ++p )
            {
                data[p] = particles[p].r[d];
            }

            std::string dset_name = "POS_" + std::to_string(d);
            writeVector( group_handle, dset_name, data );
        }
    }

    // Write the particle velocities.
    {
        for ( int d = 0; d < space_dim; ++d )
        {
            for ( int p = 0; p < num_p; ++p )
            {
                data[p] = particles[p].v[d];
            }

            std::string dset_name = "VEL_" + std::to_string(d);
            writeVector( group_handle, dset_name, data );
        }
    }

    // Write the particle material ids.
    {
        std::vector<int> matids( num_p );
        for ( int p = 0; p < num_p; ++p )
        {
            matids[p] = particles[p].matid;
        }
        std::string dset_name( "MATID" );
        writeVector( group_handle, dset_name, matids );
    }

    // Write the particle strains.
    {
        for ( int j = 0; j < space_dim; ++j )
        {
            for ( int i = 0; i < space_dim; ++i )
            {
                for ( int p = 0; p < num_p; ++p )
                {
                    data[p] = particles[p].strain[i][j];
                }

                std::string dset_name =
                    "STRAIN_" + std::to_string(i) + "_" + std::to_string(j);
                writeVector( group_handle, dset_name, data );
            }
        }
    }

    // Write the patricle stresses.
    {
        for ( int j = 0; j < space_dim; ++j )
        {
            for ( int i = 0; i < space_dim; ++i )
            {
                for ( int p = 0; p < num_p; ++p )
                {
                    data[p] = particles[p].stress[i][j];
                }

                std::string dset_name =
                    "STRESS_" + std::to_string(i) + "_" + std::to_string(j);
                writeVector( group_handle, dset_name, data );
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Write a scalar integer to file.
void FileIO::writeScalar( hid_t group_handle,
                          const std::string& dset_name,
                          const int data )
{
    hsize_t dims[1] = {1};
    H5LTmake_dataset_int( group_handle, dset_name.c_str(), 1, dims, &data);
}

//---------------------------------------------------------------------------//
// Write a scalar double to file.
void FileIO::writeScalar( hid_t group_handle,
                          const std::string& dset_name,
                          const double data )
{
    hsize_t dims[1] = {1};
    H5LTmake_dataset_double( group_handle, dset_name.c_str(), 1, dims, &data);
}

//---------------------------------------------------------------------------//
// Write an integer vector quantity to a file.
void FileIO::writeVector( hid_t group_handle,
                          const std::string& dset_name,
                          const std::vector<int>& data )
{
    hsize_t dims[1] = {data.size()};
    H5LTmake_dataset_int(
        group_handle, dset_name.c_str(), 1, dims, data.data() );
}

//---------------------------------------------------------------------------//
// Write a double vector quantity to a file.
void FileIO::writeVector( hid_t group_handle,
                          const std::string& dset_name,
                          const std::vector<double>& data )
{
    hsize_t dims[1] = {data.size()};
    H5LTmake_dataset_double(
        group_handle, dset_name.c_str(), 1, dims, data.data() );
}

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

//---------------------------------------------------------------------------//
// end FileIO.hh
//---------------------------------------------------------------------------//
