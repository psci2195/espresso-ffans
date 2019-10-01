/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file brownian_inline.hpp */

#ifndef BROWNIAN_INLINE_HPP
#define BROWNIAN_INLINE_HPP

#include "thermostat.hpp"

#ifdef BROWNIAN_DYNAMICS
/** Propagate position: viscous drag driven by conservative forces.*/
/*********************************************************/
/** \name bd_drag */
/*********************************************************/
/**(Eq. (14.39) T. Schlick, https://doi.org/10.1007/978-1-4419-6351-2 (2010) for BD
 * and eq. (8a) Ermak & Buckholz, https://doi.org/10.1016/0021-9991(80)90084-4 (1980) for EB)
 * @param &p              Reference to the particle (Input)
 * @param dt              Time interval (Input)
 */
inline void bd_drag(Particle &p, double dt) {
  // The friction tensor Z from the Eq. (14.31) of Schlick2010:
  Thermostat::GammaType local_gamma;

  if (p.p.gamma >= Thermostat::GammaType{}) {
    local_gamma = p.p.gamma;
  } else {
    local_gamma = langevin_gamma;
  }

  for (int j = 0; j < 3; j++) {
#ifndef PARTICLE_ANISOTROPY
    double beta = local_gamma / p.p.mass;
#else
    double beta = local_gamma[j] / p.p.mass;
#endif // PARTICLE_ANISOTROPY
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
#ifdef ERMAK_BUCKHOLZ
      if (thermo_switch & THERMO_ERMAK_BUCKHOLZ) {
        // init p0 here.
        // It is used further by the velocity leap (8b) of Ermak1980.
        p.r.p0[j] = p.r.p[j];
      }
#endif // ERMAK_BUCKHOLZ
      // for BD:
      // Second (deterministic) term of the Eq. (14.39) of Schlick2010.
      // Only a conservative part of the force is used here
      // for EB:
      // same terms belong to the eq. (8a), Ermak1980
      if (!(thermo_switch & THERMO_LI)) {
        // for all BD-like thermostats except LI
        p.r.p[j] += p.f.f[j] * dt / (p.p.mass * beta);
      }
      // the remaining deterministic terms of the eq. (8a), Ermak1980.
      if (thermo_switch & THERMO_ERMAK_BUCKHOLZ) {
        double tmp_exp = (1. - exp(-beta * dt)) / beta;
        // velocity is taken from the previous step end.
        p.r.p[j] += tmp_exp * ((p.m.v[j]) - (p.f.f[j] / (p.p.mass * beta)));
      } else if ((thermo_switch & THERMO_EB_VELPOS) && (dt > 0.)) {
        p.r.p[j] += (1 / beta) * (p.m.v[j] + p.m.v0[j]
                    - 2.0 * p.f.f[j] / (p.p.mass * beta)) *
                    (1 - exp(-beta * dt)) / (1 + exp(-beta * dt));
      } else if ((thermo_switch & THERMO_LI) && (dt > 0.)) {
        double exp0;
        exp0 = exp(- beta * dt / 2.);
        p.r.p[j] += (1. - exp0 * exp0) * p.m.v[j] / (beta * exp0);
      }
    }
  }
}

/** Set the terminal velocity driven by the conservative forces drag.*/
/*********************************************************/
/** \name bd_drag_vel */
/*********************************************************/
/**(Eq. (14.34) T. Schlick, https://doi.org/10.1007/978-1-4419-6351-2 (2010) for BD
 * and eq. (8b) Ermak & Buckholz, https://doi.org/10.1016/0021-9991(80)90084-4 (1980) for EB)
 * @param &p              Reference to the particle (Input)
 * @param dt              Time interval (Input)
 */
inline void bd_drag_vel(Particle &p, double dt, bool start_flag = false, bool end_flag = false) {
  // The friction tensor Z from the eq. (14.31) of Schlick2010:
  Thermostat::GammaType local_gamma;

  if (p.p.gamma >= Thermostat::GammaType{}) {
    local_gamma = p.p.gamma;
  } else {
    local_gamma = langevin_gamma;
  }

  for (int j = 0; j < 3; j++) {
#ifndef PARTICLE_ANISOTROPY
    double beta = local_gamma / p.p.mass;
#else
    double beta = local_gamma[j] / p.p.mass;
#endif // PARTICLE_ANISOTROPY
#ifdef EXTERNAL_FORCES
    if (p.p.ext_flag & COORD_FIXED(j)) {
      p.m.v[j] = 0.0;
    } else
#endif
    {
      if (thermo_switch & THERMO_BROWNIAN) {
        // First (deterministic) term of the eq. (14.34) of Schlick2010 taking
        // into account eq. (14.35). Only conservative part of the force is used
        // here NOTE: velocity is assigned here and propagated by thermal part
        // further on top of it
        p.m.v[j] = p.f.f[j] / (p.p.mass * beta);
      } else if ((thermo_switch & THERMO_ERMAK_BUCKHOLZ) && (dt > 0.)) {
        // deterministic terms of the eq. (8b), Ermak1980
        double tmp_exp = (1. - exp(-beta * dt));
        double tmp_exp2 = (1. - exp(-2. * beta * dt));
        double C = 2 * beta * dt - 3. + 4. * exp(-beta * dt)
                   - exp(-2. * beta * dt);
        p.m.v[j] = (p.m.v[j] *
                  (2 * beta * dt * exp(-beta * dt) - tmp_exp2)
                  + beta * (p.r.p[j] - p.r.p0[j]) * pow(tmp_exp, 2)
                  + (p.f.f[j] / (p.p.mass * beta)) *
                  (beta * dt * tmp_exp2 - 2. * pow(tmp_exp, 2))) / C;
      } else if ((thermo_switch & THERMO_EB_VELPOS) && (dt > 0.)) {
        // init v0 here.
        // It is used further by the position leap (7b) of Ermak1980.
        p.m.v0[j] = p.m.v[j];
        double tmp_exp = (1. - exp(-beta * dt));
        p.m.v[j] = p.m.v[j] * exp(-beta * dt) + (p.f.f[j] / (p.p.mass * beta)) * tmp_exp;
      }  else if ((thermo_switch & THERMO_LI) && (dt > 0.)) {
        double exp0 = exp(- beta * dt);
        double wplus = (exp0 - 1. + beta * dt) / (beta * dt * (1. - exp0));
        double wminus = 1. - wplus;
        double S = 0.; // LI default
        // additional bit:
        if (thermo_switch & THERMO_VGB82) {
          S = 0.5 * (wplus - wminus);
        }
        if ((! start_flag) && (! end_flag)) {
          p.m.v[j] = exp(- beta * dt / 2.) * (exp(- beta * dt / 2.) *
                      p.m.v[j] + dt * p.f.f[j] / p.p.mass +
                      dt * S * (p.f.f[j] - p.f.f_saved[j]) / p.p.mass);
        } else if (start_flag) {
          p.m.v[j] = exp(- beta * dt / 2.) * (p.m.v[j] +
                      dt * wplus * p.f.f[j] / p.p.mass);
        } else if (end_flag) {
          p.m.v[j] = exp(- beta * dt / 2.) * p.m.v[j] +
                      dt * p.f.f[j] * wminus / p.p.mass +
                      dt * S * (p.f.f[j] - p.f.f_saved[j]) / p.p.mass;
        }
        p.f.f_saved[j] = p.f.f[j];
      } // else dt==0: is not needed, the original velocity is kept
    }
  }
}

/** Determine the velocities: random walk part.*/
/*********************************************************/
/** \name bd_random_walk_vel */
/*********************************************************/
/**(Eq. (10.2.16) N. Pottier, https://doi.org/10.1007/s10955-010-0114-6 (2010) for BD
 * and eq. (8b) Ermak & Buckholz,
 * https://doi.org/10.1016/0021-9991(80)90084-4 (1980) for EB)
 * @param &p              Reference to the particle (Input)
 * @param dt              Time interval (Input)
 */
inline void bd_random_walk_vel(Particle &p, double dt, bool start_flag = false, bool end_flag = false) {
  // skip the translation thermalizing for virtual sites unless enabled
  extern bool thermo_virtual;
  if (p.p.is_virtual && !thermo_virtual)
    return;
  // Just a square root of kT, see eq. (10.2.17) and comments in 2 paragraphs
  // afterwards, Pottier2010
  extern double brown_sigma_vel;
  // first, set defaults
  double brown_sigma_vel_temp;
  // The friction tensor Z from the eq. (14.31) of Schlick2010:
  Thermostat::GammaType local_gamma;

  if (p.p.gamma >= Thermostat::GammaType{}) {
    local_gamma = p.p.gamma;
  } else {
    local_gamma = langevin_gamma;
  }

  // Override defaults if per-particle values for T and gamma are given
#ifdef LANGEVIN_PER_PARTICLE
  auto const constexpr langevin_temp_coeff = 1.0;
  // Is a particle-specific temperature specified?
  if (p.p.T >= 0.) {
    brown_sigma_vel_temp = sqrt(langevin_temp_coeff * p.p.T);
  } else {
    brown_sigma_vel_temp = brown_sigma_vel;
  }
#endif /* LANGEVIN_PER_PARTICLE */

  //Utils::Vector3d noise = v_noise_g(p.p.identity, RNGSalt::BROWNIAN);
  Utils::Vector3d noise = {0.0, 0.0, 0.0};
  for (int j = 0; j < 3; j++) {
    noise[j] = gaussian_random();
  }
  for (int j = 0; j < 3; j++) {
#ifndef PARTICLE_ANISOTROPY
    double beta = local_gamma / p.p.mass;
#else
    double beta = local_gamma[j] / p.p.mass;
#endif // PARTICLE_ANISOTROPY
#ifdef EXTERNAL_FORCES
    if (!(p.p.ext_flag & COORD_FIXED(j)))
#endif
    {
      if (thermo_switch & THERMO_BROWNIAN) {
        // Random (heat) velocity is added here. It is already initialized in the
        // terminal drag part. See eq. (10.2.16) taking into account eq. (10.2.18)
        // and (10.2.29), Pottier2010. Note, that the Pottier2010 units system
        // (see Eq. (10.1.1) there) has been adapted towards the ESPResSo and the
        // referenced above Schlick2010 one, which is defined by the eq. (14.31)
        // of Schlick2010. A difference is the mass factor to the friction tensor.
        // The noise is Gaussian according to the convention at p. 237 (last
        // paragraph), Pottier2010.
        p.m.v[j] += brown_sigma_vel_temp * noise[j] / sqrt(p.p.mass);
      } else if ((thermo_switch & THERMO_ERMAK_BUCKHOLZ) && (dt > 0.)) {
        // the random terms of the (8b), Ermak1980.
        double tmp_exp = (1. - exp(-beta * dt));
        double tmp_exp2 = (1. - exp(-2. * beta * dt));
        double C = 2 * beta * dt - 3. + 4. * exp(-beta * dt)
                   - exp(-2. * beta * dt);
        p.m.v[j] += brown_sigma_vel_temp * sqrt((2. / (p.p.mass * C)) *
                    (beta * dt * tmp_exp2 - 2. * pow(tmp_exp, 2))) *
                    noise[j];
      } else if ((thermo_switch & THERMO_EB_VELPOS) && (dt > 0.)) {
        double tmp_exp2 = (1. - exp(-2. * beta * dt));
        p.m.v[j] += brown_sigma_vel_temp * noise[j] * sqrt(tmp_exp2 / p.p.mass);
      } else if ((thermo_switch & THERMO_LI) && (dt > 0.)) {
        double R, alpha, betacorr, a, b, c, pref, wplus, wminus, exp0;
        exp0 = exp(- beta * dt);
        pref = pow(brown_sigma_vel_temp, 2) / p.p.mass;
        wplus = (exp0 - 1. + beta * dt) / (beta * dt * (1. - exp0));
        wminus = 1. - wplus;
        a = pref * (2 * pow(wplus, 2) * beta * dt + wplus - wminus);
        b = pref * (2 * wplus * wminus * beta * dt + wminus - wplus);
        c = pref * (2 * pow(wminus, 2) * beta * dt + wplus - wminus);
        betacorr = b / (p.f.alpha_saved[j] + 1E-12);
        if ((! start_flag) && (! end_flag)) {
          alpha = sqrt(a + c - betacorr * betacorr + 1E-12);
          p.f.alpha_saved[j] = alpha; // for use in the next time step
          R = sqrt(exp0) * (betacorr * p.f.noise_saved[j] + alpha * noise[j]);
        } else if (start_flag) {
          alpha = sqrt(a + 1E-12);
          p.f.alpha_saved[j] = alpha; // for use in the next time step
          R = alpha * sqrt(exp0) * noise[j];
        } else if (end_flag) {
          alpha = sqrt(c - betacorr * betacorr + 1E-12);
          R = betacorr * p.f.noise_saved[j] + alpha * noise[j];
        }
        p.f.noise_saved[j] = noise[j]; // for use in the next time step
        p.m.v[j] += R;
      }
    }
  }
}
#endif // BROWNIAN_DYNAMICS

#endif // BROWNIAN_INLINE_HPP
