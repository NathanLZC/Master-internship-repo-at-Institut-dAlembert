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

class MaxwellViscoelastic : public ContactSolver {
public:
  /// Types of algorithm (primal/dual) or constraint
  enum type { gap, pressure };

  pubilc :
      /// Constructor
      MaxwellViscoelastic(Model& model, const GridBase<Real>& surface,
                          Real tolerance, type variable_type,
                          type constraint_type, Real viscosity, Real time_step,
                          Real shear_modulus_elastic,
                          std::vector<Real> shear_modulus_maxwell,
                          std::vector<Real> characteristic_time);
  ~MaxwellViscoelastic() override = default;

public:
  /// \cond DO_NOT_DOCUMENT
  using PolonskyKeerRey::solve;
  ///\endcond

  /// Compute the viscosity \eta from shear moduli and characteristic
  /// time(relaxation time)
  virtual std::vector<Real>
  viscosityEta(const std::vector<Real>& shear_modulus_maxwell,
               const std::vector<Real>& characteristic_time) const;
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
  ///
