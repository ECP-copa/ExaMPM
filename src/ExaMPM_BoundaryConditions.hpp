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

#ifndef EXAMPM_BOUNDARYCONDITIONS_HPP
#define EXAMPM_BOUNDARYCONDITIONS_HPP

#include <Kokkos_Core.hpp>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
struct BoundaryType
{
    enum Values
    {
        NONE = 0,
        NO_SLIP = 1,
        FREE_SLIP = 2
    };
};

//---------------------------------------------------------------------------//
struct BoundaryCondition
{
    KOKKOS_INLINE_FUNCTION
    void operator()( const int gi, const int gj, const int gk, double& ux,
                     double& uy, double& uz ) const
    {
        // Low x
        if ( gi <= min[0] )
        {
            if ( boundary[0] == BoundaryType::NO_SLIP )
            {
                ux = 0.0;
                uy = 0.0;
                uz = 0.0;
            }
            else if ( boundary[0] == BoundaryType::FREE_SLIP )
            {
                ux = 0.0;
            }
        }

        // Low y
        if ( gj <= min[1] )
        {
            if ( boundary[1] == BoundaryType::NO_SLIP )
            {
                ux = 0.0;
                uy = 0.0;
                uz = 0.0;
            }
            else if ( boundary[1] == BoundaryType::FREE_SLIP )
            {
                uy = 0.0;
            }
        }

        // Low z
        if ( gk <= min[2] )
        {
            if ( boundary[2] == BoundaryType::NO_SLIP )
            {
                ux = 0.0;
                uy = 0.0;
                uz = 0.0;
            }
            else if ( boundary[2] == BoundaryType::FREE_SLIP )
            {
                uz = 0.0;
            }
        }

        // High x
        if ( gi >= max[0] )
        {
            if ( boundary[3] == BoundaryType::NO_SLIP )
            {
                ux = 0.0;
                uy = 0.0;
                uz = 0.0;
            }
            else if ( boundary[3] == BoundaryType::FREE_SLIP )
            {
                ux = 0.0;
            }
        }

        // High y
        if ( gj >= max[1] )
        {
            if ( boundary[4] == BoundaryType::NO_SLIP )
            {
                ux = 0.0;
                uy = 0.0;
                uz = 0.0;
            }
            else if ( boundary[4] == BoundaryType::FREE_SLIP )
            {
                uy = 0.0;
            }
        }

        // High z
        if ( gk >= max[2] )
        {
            if ( boundary[5] == BoundaryType::NO_SLIP )
            {
                ux = 0.0;
                uy = 0.0;
                uz = 0.0;
            }
            else if ( boundary[5] == BoundaryType::FREE_SLIP )
            {
                uz = 0.0;
            }
        }
    }

    Kokkos::Array<int, 6> boundary;
    Kokkos::Array<int, 3> min;
    Kokkos::Array<int, 3> max;
};

//---------------------------------------------------------------------------//

} // end namespace ExaMPM

#endif // end EXAMPM_BOUNDARYCONDITIONS_HPP
