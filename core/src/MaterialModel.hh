//---------------------------------------------------------------------------//
/*!
 * \file MaterialModel.hh
 */
//---------------------------------------------------------------------------//

#ifndef EXAMPM_MATERIALMODEL_HH
#define EXAMPM_MATERIALMODEL_HH

#include "StrainModel.hh"
#include "StressModel.hh"

#include <memory>

namespace ExaMPM
{
//---------------------------------------------------------------------------//
/*!
 * \class MaterialModel
 */
class MaterialModel
{
  public:

    // Strain model
    std::shared_ptr<StrainModel> strain_model;

    // Stress model
    std::shared_ptr<StressModel> stress_model;
};

//---------------------------------------------------------------------------//

};

#endif // end EXAMPM_MATERIALMODEL_HH

//---------------------------------------------------------------------------//
// end MaterialModel.hh
//---------------------------------------------------------------------------//
