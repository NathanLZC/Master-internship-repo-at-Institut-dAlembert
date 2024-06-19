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
#include "logger.hh"
#include "loop.hh"
#include "model.hh"
#include "model_type.hh"
/* -------------------------------------------------------------------------- */

namespace tamaas {

MaxwellViscoelastic::MaxwellViscoelastic(Model& model,
                                         const GridBase<Real>& surface,
                                         Real tolerance, Real time_step,
                                         Real G_inf,
                                         std::vector<Real> shear_modulus,
                                         std::vector<Real> characteristic_time)
    : PolonskyKeerRey(model, surface, tolerance, PolonskyKeerRey::pressure,
                      PolonskyKeerRey::pressure),
      time_step_(time_step), G_inf(G_inf), shear_modulus_(shear_modulus),
      characteristic_time_(characteristic_time), M(characteristic_time.size()),
      U(model.getDisplacement().getNbPoints(),
        model.getDisplacement().getNbComponents()) {
  // Check that the number of shear moduli and characteristic times is correct
  if (shear_modulus_.size() != characteristic_time_.size()) {
    throw std::invalid_argument(
        "The number of shear moduli and characteristic times must be equal.");
  }

  for (UInt k = 0; k < M.size(); ++k) {
    M[k] = GridBase<Real>(model.getDisplacement().getNbPoints(),
                          model.getDisplacement().getNbComponents());
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
Real MaxwellViscoelastic::solve(std::vector<Real> target_viscoelastic) {

  auto& M_new = model.getDisplacement();

  auto gtilde = computeGtilde(time_step_, shear_modulus_, characteristic_time_);
  auto gamma = computeGamma(time_step_, shear_modulus_, characteristic_time_);

  Real alpha = gtilde + G_inf;
  Real beta = gtilde;

  GridBase<Real> M_maxwell(M_new.getNbPoints(), M_new.getNbComponents());
  GridBase<Real> H_new(M_new.getNbPoints(), M_new.getNbComponents());

  for (size_t k = 0; k < M.size(); ++k) {
    auto gamma_k = gamma[k];
    Loop::loop(
        [gamma_k](Real& M_maxwell, const Real& M) { M_maxwell += gamma_k * M; },
        M_maxwell, M[k]);
  }

  Loop::loop(
      [alpha, beta, gamma](Real& h_new, const Real& surface, const Real& u,
                           const Real& maxwell) {
        h_new = alpha * surface - beta * u + maxwell;
      },
      H_new, this->surface, U, M_maxwell);

  const Real target = target_viscoelastic.back();

  // call PolonskyKeerRey::solve
  Real error =
      PolonskyKeerRey::solve(target);  //???shall we return something as M_new?

  GridBase<Real> U_new(U.getNbPoints(), U.getNbComponents());

  Loop::loop(
      [alpha, beta](Real& U_new, const Real& M_new, const Real& M_maxwell,
                    const Real& U) {
        U_new = (1 / alpha) * (M_new - M_maxwell + beta * U);
      },
      U_new, M_new, M_maxwell, U);

  for (size_t k = 0; k < shear_modulus_.size(); ++k) {
    auto gamma_k = gamma[k];
    auto shear_modulus_k = shear_modulus_[k];
    Loop::loop(
        [gamma_k, shear_modulus_k](Real& M, const Real& U_new, const Real& U) {
          M = gamma_k * (M + shear_modulus_k * (U_new - U));
        },
        M[k], U_new, U);
  }

  U = U_new;

  return error;
}

}  // namespace tamaas
/* ------------------------------------------------------------------------- */
