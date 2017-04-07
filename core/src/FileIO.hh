//---------------------------------------------------------------------------//
/*!
 * \file FileIO.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_FILEIO_HH
#define EXAMPM_FILEIO_HH

#include "Particle.hh"

#include <hdf5.h>
#include <hdf5_hl.h>

#include <string>
#include <vector>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*!
 * \brief File IO for reading and writing particle data.
 */
class FileIO
{
  public:

    // Constructor.
    FileIO( const std::string& file_name );

    // Write particle data at a given time step to a file.
    void writeTimeStep( const int time_step,
                        const double time,
                        const std::vector<Particle>& particles );

  private:

    // Write a scalar integer to file.
    void writeScalar( hid_t group_handle,
                      const std::string& dset_name,
                      const int data );

    // Write a scalar double to file.
    void writeScalar( hid_t group_handle,
                      const std::string& dset_name,
                      const double data );

    // Write an integer vector quantity to a file.
    void writeVector( hid_t group_handle,
                      const std::string& dset_name,
                      const std::vector<int>& data );

    // Write a double vector quantity to a file.
    void writeVector( hid_t group_handle,
                      const std::string& dset_name,
                      const std::vector<double>& data );

  private:

    // hdf5 file handle.
    hid_t d_handle;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_FILEIO_HH

//---------------------------------------------------------------------------//
// end FileIO.hh
//---------------------------------------------------------------------------//
