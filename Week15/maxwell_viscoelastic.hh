/*
 *  SPDX-License-Indentifier: AGPL-3.0-or-later
 *
 *  Copyright (©) 2016-2024 EPFL (École Polytechnique Fédérale de Lausanne),
 *  Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides),
 *  Laboratory (IJLRDA - Institut Jean Le Rond d'Alembert)
 *  Copyright (©) 2020-2024 Lucas Frérot
 *
 *  Part of codes in this script is insipred by Copilot and GPT-4.
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
#include "contact_solver.hh"
#include "grid_view.hh"
#include "meta_functional.hh"
#include "westergaard.hh"
/* -------------------------------------------------------------------------- */
namespace tamaas {

class MaxwellViscoelastic : public PolonskyKeerRey {
public:
  /// Types of algorithm (primal/dual) or constraint
  enum type { gap, pressure };

public:
  /// Constructor
  MaxwellViscoelastic(Model& model, const GridBase<Real>& surface,
                      Real tolerance, type variable_type, type constraint_type,
                      Real viscosity, Real time_step,
                      std::vector<Real> shear_modulus,
                      std::vector<Real> characteristic_time);
  ~MaxwellViscoelastic() override = default;

public:
  /// \cond DO_NOT_DOCUMENT
  using PolonskyKeerRey::solve;
  ///\endcond

  /// Solve
  Real solve(std::vector<Real> target) override;
  /// Mean on unsaturated constraint zone
  Real meanOnUnsaturated(const GridBase<Real>& field) const override;
  /// Compute squared norm
  Real computeSquaredNorm(const GridBase<Real>& field) const override;
  /// Update search direction
  void updateSearchDirection(Real factor) override;
  /// Compute critical step
  Real computeCriticalStep(Real target = 0) override;
  /// Update primal and check non-admissible state
  bool updatePrimal(Real step) override;
#endif  // MAXWELL_VISCOELASTIC_HH