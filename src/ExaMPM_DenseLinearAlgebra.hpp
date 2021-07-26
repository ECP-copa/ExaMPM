/****************************************************************************
 * Copyright (c) 2018-2020 by the ExaMPM authors                            *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of the ExaMPM library. ExaMPM is distributed under a   *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef EXAMPM_DENSELINEARALGEBRA_HPP
#define EXAMPM_DENSELINEARALGEBRA_HPP

#include <Kokkos_Core.hpp>

namespace ExaMPM
{
namespace DenseLinearAlgebra
{
//---------------------------------------------------------------------------//
// Compute the determinant of a 3x3 matrix.
template <class Real>
KOKKOS_INLINE_FUNCTION Real determinant( const Real m[3][3] )
{
    return m[0][0] * m[1][1] * m[2][2] + m[0][1] * m[1][2] * m[2][0] +
           m[0][2] * m[1][0] * m[2][1] - m[0][2] * m[1][1] * m[2][0] -
           m[0][1] * m[1][0] * m[2][2] - m[0][0] * m[1][2] * m[2][1];
}

//---------------------------------------------------------------------------//
// Compute the inverse of a 3x3 matrix with a precomputed determinant.
template <class Real>
KOKKOS_INLINE_FUNCTION void inverse( const Real m[3][3], Real det_m,
                                     Real m_inv[3][3] )
{
    Real det_m_inv = 1.0 / det_m;

    m_inv[0][0] = ( m[1][1] * m[2][2] - m[1][2] * m[2][1] ) * det_m_inv;
    m_inv[0][1] = ( m[0][2] * m[2][1] - m[0][1] * m[2][2] ) * det_m_inv;
    m_inv[0][2] = ( m[0][1] * m[1][2] - m[0][2] * m[1][1] ) * det_m_inv;

    m_inv[1][0] = ( m[1][2] * m[2][0] - m[1][0] * m[2][2] ) * det_m_inv;
    m_inv[1][1] = ( m[0][0] * m[2][2] - m[0][2] * m[2][0] ) * det_m_inv;
    m_inv[1][2] = ( m[0][2] * m[1][0] - m[0][0] * m[1][2] ) * det_m_inv;

    m_inv[2][0] = ( m[1][0] * m[2][1] - m[1][1] * m[2][0] ) * det_m_inv;
    m_inv[2][1] = ( m[0][1] * m[2][0] - m[0][0] * m[2][1] ) * det_m_inv;
    m_inv[2][2] = ( m[0][0] * m[1][1] - m[0][1] * m[1][0] ) * det_m_inv;
}

//---------------------------------------------------------------------------//
// Compute the inverse of a 3x3 matrix.
template <class Real>
KOKKOS_INLINE_FUNCTION void inverse( const Real m[3][3], Real m_inv[3][3] )
{
    Real det_m = determinant( m );
    inverse( m, det_m, m_inv );
}

//---------------------------------------------------------------------------//
// Matrix vector multiply. A*x = y
template <class Real>
KOKKOS_INLINE_FUNCTION void matVecMultiply( const Real a[3][3], const Real x[3],
                                            Real y[3] )
{
    for ( int i = 0; i < 3; ++i )
    {
        y[i] = 0.0;
        for ( int j = 0; j < 3; ++j )
            y[i] += a[i][j] * x[j];
    }
}

//---------------------------------------------------------------------------//
// Matrix Matrix  multiply. A*B = C
template <class Real>
KOKKOS_INLINE_FUNCTION void matMatMultiply( const Real a[3][3],
                                            const Real b[3][3], Real c[3][3] )
{
    for ( int i = 0; i < 3; ++i )
    {
        for ( int j = 0; j < 3; ++j )
        {
            c[i][j] = 0.0;
            for ( int k = 0; k < 3; ++k )
                c[i][j] += a[i][k] * b[k][j];
        }
    }
}

//---------------------------------------------------------------------------//
// Transpose Matrix A^{T}
template <class Real>
KOKKOS_INLINE_FUNCTION void transpose( const Real a[3][3],
                                       Real transpose_a[3][3] )
{
    for ( int i = 0; i < 3; ++i )
    {
        for ( int j = 0; j < 3; ++j )
            transpose_a[i][j] = a[j][i];
    }
}

//---------------------------------------------------------------------------//

} // end namespace DenseLinearAlgebra
} // end namespace ExaMPM

#endif // end EXAMPM_DENSELINEARALGEBRA_HPP
