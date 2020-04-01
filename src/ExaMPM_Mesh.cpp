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

#include <ExaMPM_Mesh.hpp>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
template class Mesh<Kokkos::HostSpace>;

#ifdef KOKKOS_ENABLE_CUDA
template class Mesh<Kokkos::CudaSpace>;
#endif

//---------------------------------------------------------------------------//

} // end namespace ExaMPM
