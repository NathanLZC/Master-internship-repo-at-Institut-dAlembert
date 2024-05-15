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
#include "maxwell_viscoelastic.hh"
#include "model_type.hh"
/* -------------------------------------------------------------------------- */

namespace tamaas {

MaxwellViscoelastic::MaxwellViscoelastic(Model& model,
                                         const GridBase<Real>& surface,
                                         Real tolerance, Real time_step,
                                         std::vector<Real> shear_modulus,
                                         std::vector<Real> characteristic_time)
    : PolonskyKeerRey(model, surface, tolerance, PolonskyKeerRey::pressure,
                      PolonskyKeerRey::pressure),
      time_step_(time_step), shear_modulus_(shear_modulus),
      characteristic_time_(characteristic_time) {
  // Check that the number of shear moduli and characteristic times is correct
  if (shear_modulus_.size() != characteristic_time_.size()) {
    throw std::invalid_argument(
        "The number of shear moduli and characteristic times must be equal.");
  }
}

/* ------------------------------------------------------------------------- */
Real MaxwellViscoelastic::computeGtilde(
    const Real& time_step, const std::vector<Real>& shear_modulus,
    const std::vector<Real>& characteristic_time) const {
  Real gtilde = 0;
  for (size_t i = 0; i < shear_modulus.size(); ++i) {
    gtilde += shear_modulus[i] * characteristic_time[i] /
              (characteristic_time[i] + time_step);
  }
  return gtilde;
}

/* ------------------------------------------------------------------------- */
std::vector<Real> MaxwellViscoelastic::computeGamma(
    const Real& time_step, const std::vector<Real>& shear_modulus,
    const std::vector<Real>& characteristic_time) const {
  std::vector<Real> gamma(shear_modulus.size());
  for (size_t i = 0; i < shear_modulus.size(); ++i) {
    gamma[i] = time_step / (characteristic_time[i] + time_step);
  }
  return gamma;
}

/* ------------------------------------------------------------------------- */
Real MaxwellViscoelastic::solve(std::vector<Real> target) {


  M_maxwell.fill(0);

  for (size_t k = 0; k < M.size(); ++k) {
    M_maxwell += gamma[k] * M[k];
  }

  H_new = alpha * Surface - beta * U + M_maxwell;
  std::vector<Real> target = {/* populate with needed target values based on H_new or other criteria */};

  // call PolonskyKeerRey::solve
  Real error = PolonskyKeerRey::solve(target);

  max_error = std::max(max_error, error);

  GridBase<Real> U_new = (1 / alpha) * (M_maxwell + beta * U - H_new);

  for (size_t k = 0; k < M.size(); ++k) {
    M[k] += gamma[k] * (U_new - U);
  }

  // Accumulate pressure or other metrics if needed
  Ac.push_back(
      calculateMetricBasedOnPressureOrOtherCriteria());  // Adjust as necessary
  U = U_new;

  return max_error;


}

}  // namespace tamaas
/* ------------------------------------------------------------------------- */