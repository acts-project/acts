// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_MAGNETICFIELD_SHAREDBFIELD_H
#define ACTS_MAGNETICFIELD_SHAREDBFIELD_H

#include "ACTS/MagneticField/concept/AnyFieldLookup.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

/// @ingroup MagneticField
///
/// @brief allows to use a shared magnetic field
/// in several places and with multiple steppers
/// mainly targeted to save memory  
template <typename BField> class SharedBField {
  public:
    /// @brief the constructur with a shared pointer
    /// @tparam bField is the shared BField to be stored
    SharedBField(std::shared_ptr<BField> bField) :
      m_bField(bField)
     {}
  
    /// @brief retrieve magnetic field value
    ///
    /// @param [in] position global 3D position
    ///
    /// @return magnetic field vector at given position
    Vector3D
    getField(const Vector3D& position) const
    {
      return m_bField->getField(position);
    }

    /// @brief retrieve field cell for given position
    ///
    /// @param [in] position global 3D position
    /// @return field cell containing the given global position
    ///
    /// @pre The given @c position must lie within the range of the underlying
    ///      magnetic field map.
    concept::AnyFieldCell<>
    getFieldCell(const Vector3D& position) const
    {
      return m_bField->getFieldCell(position);
    }

    /// @brief retrieve magnetic field value & its gradient
    ///
    /// @param [in]  position   global 3D position
    /// @param [out] derivative gradient of magnetic field vector as (3x3) matrix
    /// @return magnetic field vector
    ///
    /// @note currently the derivative is not calculated
    /// @todo return derivative
    Vector3D
    getFieldGradient(const Vector3D& position,
                     ActsMatrixD<3, 3>& derivative) const
    {
      return m_bField->getField(position);
    }  
  private:
    std::shared_ptr<BField> m_bField;  
    
};

} // namespace Acts

#endif // ACTS_MAGNETICFIELD_SHAREDBFIELD_H
