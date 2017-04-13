//---------------------------------------------------------------------------//
/*!
 * \file TensorTools.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_TENSORTOOLS_HH
#define EXAMPM_TENSORTOOLS_HH

#include <array>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*!
 * \brief Tools for tensor operations.
 */
class TensorTools
{
  public:

    //! Determinant of a 3x3 tensor
    static double determinant( const std::array<std::array<double,3>,3>& t )
    {
        return
            t[0][0] * t[1][1] * t[2][2] +
            t[0][1] * t[1][2] * t[2][0] +
            t[0][2] * t[1][0] * t[2][1] -
            t[0][2] * t[1][1] * t[2][0] -
            t[0][1] * t[1][0] * t[2][2] -
            t[0][0] * t[1][2] * t[2][1];
    }
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_TENSORTOOLS_HH

//---------------------------------------------------------------------------//
// end TensorTools.hh
//---------------------------------------------------------------------------//
