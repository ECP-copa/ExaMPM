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

    //! Inverse of a 3x3 tensor.
    static void inverse( const std::array<std::array<double,3>,3>& t,
                         std::array<std::array<double,3>,3>& t_inv )
    {
        double det_t = determinant( t );

        t_inv[0][0] = (t[1][1]*t[2][2] - t[1][2]*t[2][1]) / det_t;
        t_inv[0][1] = (t[0][2]*t[2][1] - t[0][1]*t[2][2]) / det_t;
        t_inv[0][2] = (t[0][1]*t[1][2] - t[0][2]*t[1][1]) / det_t;

        t_inv[1][0] = (t[1][2]*t[2][0] - t[1][0]*t[2][2]) / det_t;
        t_inv[1][1] = (t[0][0]*t[2][2] - t[0][2]*t[2][0]) / det_t;
        t_inv[1][2] = (t[0][2]*t[1][0] - t[0][0]*t[1][2]) / det_t;

        t_inv[2][0] = (t[1][0]*t[2][1] - t[1][1]*t[2][0]) / det_t;
        t_inv[2][1] = (t[0][1]*t[2][0] - t[0][0]*t[2][1]) / det_t;
        t_inv[2][2] = (t[0][0]*t[1][1] - t[0][1]*t[1][0]) / det_t;
    }
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_TENSORTOOLS_HH

//---------------------------------------------------------------------------//
// end TensorTools.hh
//---------------------------------------------------------------------------//
