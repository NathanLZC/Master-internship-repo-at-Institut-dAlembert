/*
 *  SPDX-License-Indentifier: AGPL-3.0-or-later
 *
 *  Copyright (©) 2016-2024 EPFL (École Polytechnique Fédérale de Lausanne),
 *  Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides),
 *  Laboratory (IJLRDA - Institut Jean Le Rond d'Alembert)
 *  Copyright (©) 2020-2024 Lucas Frérot
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as published
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#ifndef MAXWELL_VISCOELASTIC_HH
#define MAXWELL_VISCOELASTIC_HH
/* -------------------------------------------------------------------------- */
#include "polonsky_keer_rey.hh"
#include "westergaard.hh"
/* -------------------------------------------------------------------------- */
namespace tamaas {

class MaxwellViscoelastic : public PolonskyKeerRey {
public:
public:
  /// Constructor
  MaxwellViscoelastic(Model& model, const GridBase<Real>& surface,
                      Real tolerance, Real time_step,
                      std::vector<Real> shear_modulus_maxwell,
                      std::vector<Real> characteristic_time);
  ~MaxwellViscoelastic() override = default;

public:
  /// \cond DO_NOT_DOCUMENT
  using PolonskyKeerRey::solve;
  ///\endcond

  /// Compute the effective shear moduli for maxwell branches
  virtual Real
  computeGtilde(const Real& time_step,
                const std::vector<Real>& shear_modulus_maxwell,
                const std::vector<Real>& characteristic_time) const;
  /// Compute the partial displacement coefficient for each maxwell branch to
  /// update surface
  virtual std::vector<Real>
  computeGamma(const Real& time_step,
               const std::vector<Real>& shear_modulus_maxwell,
               const std::vector<Real>& characteristic_time) const;
  /// Main solve function with backward euler scheme
  Real solve(std::vector<Real> target) override;

protected:
  Real time_step_;
  std::vector<Real> shear_modulus_;
  std::vector<Real> characteristic_time_;
  // Two dimension array partial displacement M
  std::vector<GridBase<Real>> M;
  // Global displacement U
  GridBase<Real> U;
};

}  // namespace tamaas

#endif  // MAXWELL_VISCOELASTIC_HH